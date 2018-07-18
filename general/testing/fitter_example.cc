#include "Fitter.h"
#include "FitAnalyser.h"
#include "RooWorkspace.h"
#include "RooDataSet.h"
#include "TTree.h"
#include "RooStats/SPlot.h"

#include <iostream>


int main() {

  using FitFlag = Fitter::FitFlag;
  
  Fitter f;

  f.SetInputData("/afs/hephy.at/work/j/jnecker/data/bug_study/jpsi_2016.root", "data");
  f.SetOutfile("fitter_testoutput.root");
  f.SetWorkspaceName("test_workspace_pt_18_20",true);
  f.AddBinVariable("jpsiPt", 18, 20);
  f.SetModel({
//    "CBShape::jpsi_l(jpsiMass,mu[2.9,3.3],sigma[0.01,0.2],a_l[0.5,10],n_l[1,10])",
//    "CBShape::jpsi_r(jpsiMass,mu,sigma,a_r[-10,-0.5],n_r[1,10])",
//    "ExpGaussExp::jpsi(jpsiMass, mu[2.9,3.3],sigma[0.01,0.2],k_l[0.5,10],k_h[0.5,10])",
    "DoubleSidedCB::jpsi(jpsiMass, mu[3.05,3.15], sigma[0.01,0.05], a_l[0.5,7], n_l[1,10], a_h[0.5,7], n_h[1,10])",
//    "SUM::jpsi(frac_left[0,1]*jpsi_l,jpsi_r)",
    "Exponential::background(jpsiMass,l[-10,-0.001])",
    "SUM::model(N_signal[1e6,1e5,2e6]*jpsi,N_bkg[0,1e6]*background)"
  },
    "model", "jpsiMass", 2.95, 3.25);
  f.Fit(8, FitFlag::ExtendedFit | FitFlag::EnableMinos | FitFlag::SuppressOutput);
 
  // Another bin
/*
  f.SetWorkspaceName("test_workspace_pt_20_21");
  f.AddBinVariable("jpsiPt", 20,21);
  f.Fit(8, false,true);
*/


  FitAnalyser a(f);
  a.SetWorkspaceName("test_workspace_pt_18_20");

  bool ok = false;
  for (const auto &v : { "mu", "sigma", "blabla_test" }) {
    double val = a.GetVariableValue(v, ok);
    if (ok) std::cout << v << ": " << a.GetVariableValue(v, ok) << std::endl;
    else std::cout << "Could not read value for variable '" << v << "'." << std::endl;
  }
  a.PlotFitResult("fitresult_example.pdf");
  a.GetWorkspace()->Print();
  a.AddSWeights({ "N_signal", "N_bkg" });

  a.GetWorkspace()->Print(); // Now the weights should be in the dataset

  
  // Save final TTree with sWeights
  auto data = a.GetDataset();
  RooAbsData::setDefaultStorageType(RooAbsData::Tree);
  RooDataSet * tree_data = new RooDataSet("tree_with_weights", "Tree with sWeights", data, *data->get());
  auto tree = tree_data->tree();
  TFile treefile("fitter_exampletree_with_sweights.root", "recreate");
  tree->Write();

}
