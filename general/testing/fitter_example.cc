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

  f.SetInputData("/afs/hephy.at/work/j/jnecker/data/chib_results/preselected_data_chib2016_rereco.root", "data");
  f.SetOutfile("fitter_testoutput.root");
  f.SetWorkspaceName("test_workspace", true);
  f.AddBinVariable("dimuon_pt", 15, 16);
  f.SetModel({ "CBShape::ups1s_1(dimuon_mass, mu1s[9.5,9.3,9.6], sigma1s1[0.075,0.03,0.11],a_l[0.1,5], 3)",
      "CBShape::ups1s_2(dimuon_mass, mu1s, sigma1s2[0.075,0.03,0.11], a_r[-5,-0.1], 2)",
      "SUM::ups1s(tail_left1[0.5,0.,1.]*ups1s_1,ups1s_2)",
      "expr::mu2('mu1s*1.059508', mu1s)",
      "CBShape::ups2s(dimuon_mass, mu2, sigma2s[0.07,0.03,0.11], a2s[0.1,5], 3)",
      "expr::mu3('mu1s*1.094595', mu1s)",
      "CBShape::ups3s(dimuon_mass,  mu3, sigma3s[0.07,0.03,0.11], a3s[0.1,5], 3)",
     // "Exponential::background(dimuon_mass, l[-0.14,-1,1])",
      "Chebychev::background(dimuon_mass,{c1[-5,5]})",
      "SUM::model(N1[0,1e6]*ups1s, N2[0,1e6]*ups2s, N3[0,1e6]*ups3s, N_bkg[0,1e6]*background)"
  },
    "model", "dimuon_mass", 8.7, 11.1);
  f.SetBackground("background", { {8.7,9.1}, {10.65,11.1}, {9.7,9.8} });
  f.Fit(8, FitFlag::ExtendedFit | FitFlag::EnableMinos | FitFlag::SuppressOutput);

  // Another bin
/*
  f.SetWorkspaceName("test_workspace_pt_20_21");
  f.AddBinVariable("jpsiPt", 20,21);
  f.Fit(8, false,true);
*/


  FitAnalyser a(f);
  a.SetWorkspaceName("test_workspace");

  bool ok = false;
  for (const auto &v : { "mu1s", "sigma1s1", "blabla_test" }) {
    double val = a.GetVariableValue(v, ok);
    if (ok) std::cout << v << ": " << a.GetVariableValue(v, ok) << std::endl;
    else std::cout << "Could not read value for variable '" << v << "'." << std::endl;
  }
  a.PlotFitResult("fitresult_example.pdf");
  a.GetWorkspace()->Print();
  a.AddSWeights({ "N1","N2","N3", "N_bkg" });

  a.GetWorkspace()->Print(); // Now the weights should be in the dataset


  // Save final TTree with sWeights
  auto data = a.GetDataset();
  RooAbsData::setDefaultStorageType(RooAbsData::Tree);
  RooDataSet * tree_data = new RooDataSet("tree_with_weights", "Tree with sWeights", data, *data->get());
  auto tree = tree_data->tree();
  TFile treefile("fitter_exampletree_with_sweights.root", "recreate");
  tree->Write();

}
