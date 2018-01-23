#include "Fitter.h"
#include "utils.h"

//RooFit includes
#include "RooDataSet.h"
#include "RooRealVar.h"
#include "RooWorkspace.h"
#include "RooFitResult.h"
#include "RooArgSet.h"
#include "RooAbsPdf.h"

//Root includes
#include "TTree.h"
#include "TFile.h"

//others
#include <memory>
#include <sstream>
#include <fstream>
#include <iostream>



const std::string Fitter::dataset_name = "dataset";
const std::string Fitter::snapshot_name = "results";


void Fitter::Fit(int numCPUs, bool enableMinos)
{
  scopeLog log("Fitter");
  if (!params_ok()) return;

  // Check if workspace already exists and delete it if force_refit is set
  if (!m_force_recreation && file_exists(m_ofname)) {

    TFile of(m_ofname.c_str(), "UPDATE");

    auto tmp_ws = dynamic_cast<RooWorkspace*>(of.Get(m_wsname.c_str()));
    if (tmp_ws) {
      if (m_force_refit) {
        std::cout << "Deleting existing workspace '" << m_wsname << "'." << std::endl;
        of.Delete((m_wsname + ";*").c_str());
        of.Write(0, TObject::kWriteDelete);
      }
      else {
        std::cout << "Workspace '" << m_wsname << "' already exists in file '" << m_ofname
          << "'.\nTo force a refit call SetWorkspace with true as second parameter." << std::endl;
        return;
      }
    }
  }

  // Open file
  TFile df(m_in_filename.c_str(), "READ");
  if (df.IsZombie()) {
    std::cout << "Fitter: Could not open file " << m_in_filename << "." << std::endl;
    return;
  }

  // Get dataset
  bool ok = false;
  auto ds = get_dataset(&df, &ok);
  if (!ok) return;

  // Setup workspace
  //ROOWORKSPACE HAS TO BE CREATED ON HEAP, warum auch immer, sonst gibts Problem mit Write();
  auto ws = new RooWorkspace(m_wsname.c_str());
  if (!setup_workspace(ws, &ds)) return;

  // Get model
  auto model = ws->pdf(m_modelname.c_str());
  if (!model) {
    std::cout << "Model '" << m_modelname << "' not found in workspace, check factory strings:\n";
    ws->Print();
    return;
  }

  //Fit
  std::unique_ptr<RooFitResult> fit_results;
  {
    scopeLog log("Fitting");
    stopwatch watch("Fitting");
    fit_results.reset(model->fitTo(ds, RooFit::NumCPU(numCPUs), RooFit::Minos(enableMinos), RooFit::Save(true)));
    fit_results->Print();
  }

  //Save results:
  ws->import(*fit_results);
  ws->saveSnapshot(snapshot_name.c_str(), RooArgSet(ws->allVars(), ws->allPdfs()));

  df.Close("R"); // necessary if in an out file are the same

  TFile f(m_ofname.c_str(), m_force_recreation ? "RECREATE" : "UPDATE");
  ws->Write(0, TObject::kWriteDelete);

  //TODO: What are the ProcessIDs used for? 
  f.Delete("ProcessID*;*"); // Seems that the ids are relicts

  //TODO:
  //write_logfile

  delete ws;
}

Fitter::~Fitter()
{
}

RooDataSet Fitter::get_dataset(TFile* f, bool *ok)
{

  scopeLog log("Creating dataset");
  stopwatch watch("Creating dataset");

  if (ok) *ok = false;

  // Setup RooArgSet
  RooArgSet arg_set;
  for (int i = 0, s = m_roo_vars.size(); i < s; ++i) {
    const auto & vn = m_roo_vars[i];
    Double_t min = m_roovar_borders[i].first;
    Double_t max = m_roovar_borders[i].second;

    if (min == max) { // Find out full range
      if (!m_in_treename.empty()) {
        TTree *tree = dynamic_cast<TTree *>(f->Get(m_in_treename.c_str()));
        if (tree) min = tree->GetMinimum(vn.c_str()), max = tree->GetMaximum(vn.c_str());
      }
      else if (!m_in_wsname.empty() && !m_in_dsname.empty()) {
        auto ws = dynamic_cast<RooWorkspace*>(f->Get(m_in_wsname.c_str()));
        RooRealVar* var = nullptr;
        if (ws && (var = ws->var(vn.c_str()))) min = var->getMin(), max = var->getMax();
      }
    }

    RooRealVar tempVar(vn.c_str(), vn.c_str(), min, max);
    arg_set.addClone(tempVar);
  }

  if (arg_set.getSize() == 0) {
    std::cout << "Problems creating RooRealVars, RooArgSet is empty." << std::endl;
    return RooDataSet();
  }

  // 1. Dataset from tree
  if (!m_in_treename.empty()) {

    // Get tree
    TTree *tree = dynamic_cast<TTree *>(f->Get(m_in_treename.c_str()));
    if (!tree) {
      std::cout << "Could not find tree '" << m_in_treename << "' in file '" << f->GetName() << "'." << std::endl;
      return RooDataSet();
    }

    // RooDataSet makes a FULL copy of the WHOLE tree into memory,
    // activating only the necessary branches makes this faster
    tree->SetBranchStatus("*", false);
    for (const auto &bn : m_roo_vars) {
      UInt_t check = 0;
      tree->SetBranchStatus(bn.c_str(), true, &check);
      if (!check) {
        std::cout << "Branch '" << bn << "' not found in tree '" << m_in_treename << "'." << std::endl;
        return RooDataSet();
      }
    }

    if (ok) *ok = true;
    return RooDataSet(dataset_name.c_str(), dataset_name.c_str(), tree, arg_set);
  }

  // 2. Dataset from existing dataset
  else if (!m_in_wsname.empty() && !m_in_dsname.empty()) {
    auto ws = dynamic_cast<RooWorkspace*>(f->Get(m_in_wsname.c_str()));
    if (!ws) {
      std::cout << "Could not find workspace '" << m_in_wsname << "' in file '" << f->GetName() << "'." << std::endl;
      return RooDataSet();
    }
    auto ds = dynamic_cast<RooDataSet*>(ws->data(m_in_dsname.c_str()));
    if (!ds) std::cout << "Could not find dataset '" << m_in_dsname << "' in workspace '" << m_in_wsname << "'." << std::endl;

    if (ok) *ok = true;
    return RooDataSet(dataset_name.c_str(), dataset_name.c_str(), ds, arg_set);
  }
  
  std::cout << "Could not create the RooDataset." << std::endl;
  return RooDataSet();
}
bool Fitter::setup_workspace(RooWorkspace* ws, RooDataSet *ds)
{
  scopeLog log("Setup workspace");
  stopwatch watch("Setup workspace");

  ws->import(*ds);

  // Create model
  for (auto f : m_factory_strings) {
    std::cout << "Running factory: " << f << std::endl;
    auto result = ws->factory(f.c_str());
    if (!result) {
      std::cout << "Problems creating factory " << f << std::endl;
      return false;
    }
  }
  return true;
}

bool Fitter::params_ok()
{
  bool isOk = true;
  if (m_wsname.empty()) isOk = false, std::cout << "No workspace name defined.\n";
  if (m_ofname.empty()) isOk = false, std::cout << "No outfile defined.\n";
  if (m_in_filename.empty()) isOk = false, std::cout << "You have to define a data input file.\n";
  if (m_in_treename.empty() && (m_in_wsname.empty() || m_in_dsname.empty())) isOk = false, std::cout << "You have to define an input tree or a workspace containing a dataset.\n";
  if (m_factory_strings.empty() || m_modelname.empty() || m_fitvarname.empty() || m_roovar_borders.empty() || m_roo_vars.empty()) isOk = false, std::cout << "Fit model not well defined.\n";

  std::cout << std::flush;
  return isOk;
}
