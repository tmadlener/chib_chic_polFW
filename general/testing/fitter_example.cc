#include "Fitter.h"
#include "FitAnalyser.h"
#include "RooWorkspace.h"

#include <iostream>


int main() {
  Fitter f;

  f.SetInputData("/afs/hephy.at/work/j/jnecker/data/bug_study/jpsi_2016.root", "data");
  f.SetOutfile("fitter_testoutput.root");
  f.SetWorkspaceName("test_workspace_pt_18_20");
  f.AddBinVariable("jpsiPt", 18, 20);
  f.SetModel({
    "CBShape::jpsi_l(jpsiMass,mu[2.9,3.3],sigma[0.01,0.2],a_l[0.5,10],n_l[1,10])",
    "CBShape::jpsi_r(jpsiMass,mu,sigma,a_r[-10,-0.5],n_r[1,10])",
    "SUM::jpsi(frac_jpsi_left[0,1]*jpsi_l,jpsi_r)",
    "Exponential::background(jpsiMass,l[-2,-0.001])",
    "SUM::model(frac_jpsi[0,1]*jpsi,background)"
  },
    "model", "jpsiMass", 3, 3.2);
  f.Fit(8,false);
  

  // Another bin
  f.SetWorkspaceName("test_workspace_pt_20_21");
  f.AddBinVariable("jpsiPt", 20,21);
  f.Fit(8, false);

  FitAnalyser a(f);
  a.SetWorkspaceName("test_workspace_pt_18_20");

  bool ok = false;
  for (const auto &v : { "mu", "sigma", "frac_jpsi_left" }) {
    double val = a.GetVariableValue(v, ok);
    if (ok) std::cout << v << ": " << a.GetVariableValue(v, ok) << std::endl; 
    else std::cout << "Could not read value for variable '" << v <<"'." << std::endl;
  }
  a.GetWorkspace()->Print();

  a.PlotFitResult("fitresult_example.pdf");

}