#include "FitAnalyser.h"

#include "TFile.h"
#include "RooHist.h"
#include "RooPlot.h"
#include "RooAddPdf.h"
#include "RooWorkspace.h"
#include "RooRealVar.h"
#include "RooAbsPdf.h"
#include "RooDataSet.h"
#include "RooFitResult.h"
#include "TCanvas.h"
#include "TAxis.h"
#include "TLine.h"
#include "TROOT.h"
#include "TLegend.h"
#include "TH1F.h"
#include "TPaveText.h"
#include "TFormula.h"
#include "RooStats/SPlot.h"

#include <sstream>
#include <memory>
#include <iomanip>
#include <regex>


FitAnalyser::FitAnalyser(const Fitter &f) :
  m_filename(f.m_ofname)
{
  SetWorkspaceName(f.m_wsname);
  SetModelName(f.m_modelname);
  SetFitVariableName(f.m_fitvarname);
}

FitAnalyser::FitAnalyser(csr filename, csr model_name, csr fit_variable_name, csr workspace_name, csr dataset_name, csr snapshot_name) :
  m_filename(filename)
{
  // Check if file is ok
  if (!get_file()) {
    std::cout << "EXITING PROGRAM: [FitAnalyser] could not open file. \n";
    exit(EXIT_FAILURE);
  }

  SetModelName(model_name);
  SetFitVariableName(fit_variable_name);
  SetWorkspaceName(workspace_name);
  if (!snapshot_name.empty()) SetSnapshotName(snapshot_name);
  if (!dataset_name.empty()) SetDatasetName(dataset_name);
}

void FitAnalyser::SetWorkspaceName(csr workspace_name)
{
  if (workspace_name.empty()) {
    // If there is only one workspace in the file take that one
    auto f = get_file();
    auto keys = f->GetListOfKeys();
    std::vector<RooWorkspace*> wspaces;
    for (auto k : *keys) {
      auto ws = dynamic_cast<RooWorkspace*>(f->Get(k->GetName()));
      if (ws) wspaces.emplace_back(ws);
    }
    if (wspaces.size() == 1) m_wsname = wspaces[0]->GetName();
    else if (wspaces.empty()) {
      std::cout << "No workspace found in file '" << m_filename << "'." << std::endl;
      return;
    }
    else if (wspaces.size() > 1) {
      std::cout << "Too many workspaces found in file '" << m_filename << "':\n";
      for (auto ws : wspaces) std::cout << '\t' << ws->GetName() << "; " << ws->GetTitle() << '\n';
      std::cout << "You have to provide a workspace name." << std::endl;
    }
  }
  else m_wsname = workspace_name;
}

RooWorkspace * FitAnalyser::GetWorkspace(bool lsnap)
{
  if (!m_workspace) {
    m_workspace = dynamic_cast<RooWorkspace*>(get_file()->Get(m_wsname.c_str()));
    if (!m_workspace) std::cout << "Workspace '" << m_wsname << "' not found in file '" << m_filename << "'." << std::endl;

  }

  if (m_workspace && lsnap && !m_snapshotname.empty()) m_workspace->loadSnapshot(m_snapshotname.c_str());

  return m_workspace;
}

RooDataSet * FitAnalyser::GetDataset(csr dataset_name)
{
  auto ws = GetWorkspace(true);
  if (!ws) return nullptr;

  auto dsname = dataset_name.empty() ? m_datasetname : dataset_name;
  auto ds = dynamic_cast<RooDataSet*>(ws->data(dsname.c_str()));
  if (!ds) std::cout << "No RooDataSet named '" << dsname << "' in RooWorkspace '" << m_wsname << "'." << std::endl;

  return ds;
}

RooRealVar * FitAnalyser::GetVariable(csr variable_name)
{
  auto ws = GetWorkspace(true);
  if (!ws) return nullptr;

  auto varname = variable_name.empty() ? m_fitvarname : variable_name;
  auto var = ws->var(varname.c_str());
  if (!var) std::cout << "No RooRealVar named '" << varname << "' in RooWorkspace '" << m_wsname << "'." << std::endl;

  return var;
}

RooFitResult * FitAnalyser::GetFitResult(csr fitresult_name)
{
  auto ws = GetWorkspace(true);
  if (!ws) return nullptr;

  if (!fitresult_name.empty()) {
    auto r = dynamic_cast<RooFitResult*>(ws->genobj(fitresult_name.c_str()));
    // if (!r) std::cout << "Could not find RooFitResult '" << fitresult_name << "' in RooWorkspace '" << m_wsname << "'." << std::endl;
    return r;
  }
  // else return fitresult if only one exists in workspace
  auto objs = ws->allGenericObjects();
  std::vector<RooFitResult*> allfitresults;
  for (auto o : objs) {
    auto r = dynamic_cast<RooFitResult*>(o);
    if (r) allfitresults.emplace_back(r);
  }
  if (allfitresults.size() == 1) return allfitresults[0];
  if (allfitresults.empty()) std::cout << "No RooFitResult found in RooWorkspace '" << m_wsname << "'." << std::endl;
  if (allfitresults.size() > 1) {
    std::cout << "Too many RooFitResults found in RooWorkspace '" << m_wsname << "':\n";
    for (auto res : allfitresults) std::cout << '\t' << res->GetName() << "; " << res->GetTitle() << '\n';
    std::cout << "You have to provide a RooFitResult name." << std::endl;
  }

  return nullptr;
}

RooAbsPdf * FitAnalyser::GetPdf(csr pdf_name)
{
  auto ws = GetWorkspace(true);
  if (!ws) return nullptr;

  return ws->pdf(pdf_name.empty() ? m_modelname.c_str() : pdf_name.c_str());
}

double FitAnalyser::GetVariableValue(csr variable_name, bool &ok)
{
  ok = false;

  auto var = GetVariable(variable_name);
  if (!var) return 0;

  return ok = true, var->getValV();
}

void FitAnalyser::PlotFitResult(csr output_file, const std::vector<double> & lines, bool addResultBox)
{
  if (!output_file.empty()) gROOT->SetBatch(kTRUE);
  auto canv = new TCanvas("fit_results", "Fit Results", 1200, 1000);

  auto var = GetVariable();
  auto ds = GetDataset();
  auto model = dynamic_cast<RooAddPdf*>(GetPdf());
  auto ws = GetWorkspace(true);
  if (!ds || !var || !model || !ws) return;

  auto plot = var->frame();
  auto pull = var->frame();
  plot->SetTitle(("Fit results from " + m_wsname + ";" + var->GetTitle() + ";").c_str());
  pull->SetTitle(";;Pull");
  plot->SetMinimum(0);

  // Legend - TODO: find a solution
  //https://root-forum.cern.ch/t/pdf-color-in-tlegend/26621

  // Print variables on the plot
  auto ibox = new TPaveText(0.7, 0.5, 0.99, 0.95, "BRNDC");
  ibox->SetTextSize(ibox->GetTextSize()*0.6);
  {
    auto chi2_ndof = get_chi2_ndof();
    std::stringstream ss;
    ss << std::fixed << "#chi^{2}/N_{dof}: "
      << std::setprecision(2) << chi2_ndof.first << "/"
      << std::setprecision(0) << chi2_ndof.second;
    ibox->AddText(ss.str().c_str());
  }
  ibox->AddText("");

  {
    auto vars = ws->allVars();
    auto it = vars.createIterator();
    TObject *var;
    while ((var = it->Next())) {
      auto v = dynamic_cast<RooRealVar*>(var);
      if (!v) continue;

      std::stringstream ss;
      ss << std::fixed << std::setprecision(4);
      ss << v->GetName() << ": " << v->getVal() << ", [ " << v->getMin() << ", " << v->getMax() << "]";

      ibox->AddText(ss.str().c_str());
    }
  }
  ibox->Print();



  // Pads for the two plots
  double pull_size = 0.15;
  double big_marg = 0.15;
  double small_marg = 0.05;
  TPad *plot_pad = new TPad("plot", "", 0, pull_size, 1, 1);
  TPad *pull_pad = new TPad("pull", "", 0, 0, 1, pull_size);
  plot_pad->SetMargin(big_marg, small_marg, 0.1, 0.1);
  pull_pad->SetMargin(big_marg, small_marg, small_marg, 0);
  plot_pad->SetGrid();
  pull_pad->SetGrid();
  pull_pad->SetTicks(1, 1);


  //auto fr = GetFitResult();
  //if (fr) model->plotOn(plot, RooFit::LineColor(kRed), RooFit::VisualizeError(*fr));

  // Plot model and dataset and get pull
  ds->plotOn(plot);
  model->plotOn(plot, RooFit::LineColor(kRed));
  auto pull_histo = plot->pullHist();

  auto y = pull->GetYaxis();
  y->SetLabelSize(big_marg*0.8);
  y->SetTitleOffset(0.3);
  y->SetTitleSize(0.2);
  y->SetNdivisions(6, 2, 0, true);
  y->SetLabelSize(0.2);
  y->SetTickSize(0.01);
  auto x = pull->GetXaxis();
  x->SetLabelSize(0);
  x->SetTickSize(0.15);
  pull->addPlotable(pull_histo, "P");
  pull->addObject(new TLine(var->getMin(), 0, var->getMax(), 0));


  // Plot single pdfs
  auto pdfs = model->pdfList();
  for (auto i = 0, s = pdfs.getSize(); i < s; ++i) {
    auto p = model->plotOn(plot, RooFit::Components(*pdfs.at(i)), RooFit::LineColor(kBlue + 2 - 2 * i), RooFit::LineWidth(2), RooFit::LineStyle(kDashed));
  }

  //plot dataset again, to be on top
  ds->plotOn(plot);

  for (double x : lines) plot->addObject(new TLine(x, plot->GetMinimum(), x, plot->GetMaximum()));

  if(addResultBox) plot->addObject(ibox);

  // Draw to canvas
  canv->cd();
  pull_pad->Draw();
  plot_pad->Draw();
  plot_pad->cd();
  plot->Draw();
  pull_pad->cd();
  pull->Draw();

  if (!output_file.empty()) canv->SaveAs(output_file.c_str());
}

double FitAnalyser::EvaluateFormula(csr formulastring, bool &ok)
{
  ok = false;
  auto ws = GetWorkspace(true);
  if (!ws) return 0;
  std::string formulaexpr(formulastring);
  auto params = ws->allVars();
  auto it = params.createIterator();

  while (auto param = dynamic_cast<RooRealVar*>(it->Next())) {
    std::string param_name(param->GetName());
    if (param_name != m_fitvarname) {
      std::regex reg("\\b" + param_name + "\\b");
      std::string param_val(std::to_string(param->getVal()));

      formulaexpr = std::regex_replace(formulaexpr, reg, param_val);
    }
  }

  // Create a TFormula and check if it's valid
  TFormula formula("temp_formula", formulaexpr.c_str());
  if (!formula.IsValid()) {
    std::cout << "Could not evaluate " << formulaexpr << " from workspace " << ws->GetName() << std::endl;
    return 0;
  }

  return ok = true, formula.Eval(0.0);
}

double FitAnalyser::GetXForQuantile(csr pdf_name, double q, bool &ok)
{
  ok = false;
  if (q > 1 || q < 0) {
    std::cout << "GetXForQuantile: The quantile has to be between 0 and 1, but is " << q << ".";
    return 0;
  }
  //TODO: see https://github.com/tmadlener/phys_utils/blob/unbinned-ml-fit/UnbinnedPolFW/interface/BkgWeightCalcChic.h

}

RooDataSet * FitAnalyser::AddSWeights(const std::vector<std::string>& yield_names)
{
  if (yield_names.empty()) return GetDataset();

  auto ws = GetWorkspace(true);
  auto ds = GetDataset();
  auto pdf = GetPdf();

  if (!ws || !ds || !pdf) return GetDataset();

  RooArgList yield_list;
  for (const auto &yld : yield_names) yield_list.add(*GetVariable(yld));
  yield_list.Print();

  //{
  //  // Set all variables, but yields and fitvar constant
  //  auto it = vars->createIterator();
  //  TObject *var;
  //  while ((var = it->Next())) {
  //    auto v = dynamic_cast<RooRealVar*>(var);
  //    if (!v) continue;
  //    std::string v_name(v->GetName());

  //    if (v_name == m_modelname || v_name == m_fitvarname) continue;
  //    bool yieldcheck_passed = true;
  //    for (const auto &y_name : yield_names) if (v_name == y_name) {
  //      yieldcheck_passed = false;
  //      break;
  //    }
  //    if (!yieldcheck_passed) continue;
  //    const_vars.push_back({ v,v->isConstant() });
  //    std::cout << "Setting variable " << v_name << " constant." << std::endl;
  //    v->setConstant();
  //  }
  //}

  RooStats::SPlot sData("sData", "sData", *ds, pdf, yield_list);

  get_file()->cd();
  ws->Write(0, TObject::kWriteDelete);



  return GetDataset();
  //
  // Next steps:
  // TreeMerger: RooDataset as argument and do the following:
  // TODO in calling program:
  // Set defaultstoragetype to ttree
  // Copy this dataset into a new one
  // get tree from dataset and merge it with the input tree
}

TFile * FitAnalyser::get_file()
{
  if (!m_file)  m_file.reset(new TFile(m_filename.c_str(), "update"));
  if (!m_file) std::cout << "FitAnalyser: Could not open file  '" << m_filename << "'." << std::endl;
  if (m_file->IsZombie()) {
    std::cout << "FitAnalyser: File  '" << m_filename << "' is a zombie." << std::endl;
    m_file.reset(nullptr);
  }

  return m_file.get();
}


std::pair<double, double> FitAnalyser::get_chi2_ndof()
{
  auto var = GetVariable();
  auto fr = GetFitResult();
  auto ds = GetDataset();
  auto model = GetPdf();
  if (!fr || !var || !ds || !model) return{ -1,-1 };

  auto frame = var->frame();
  ds->plotOn(frame);
  model->plotOn(frame);

  int fit_parameter_count = (fr->floatParsFinal()).getSize();
  int bin_count = frame->GetNbinsX();
  double reduced_chi2 = frame->chiSquare(fit_parameter_count);  //reduced chi-squared = chi2/ndof

  auto ndof = bin_count - fit_parameter_count;  //number of degrees of freedom
  auto chi2 = reduced_chi2*ndof;

  return{ chi2, ndof };


}