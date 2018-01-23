#include "Fitter.h"
#include "FitAnalyser.h"
#include "TROOT.h"

#include <map>
#include <iostream>

// TODO:
//  - a class or function that creates workspacename and so on from fitvar, bins, ...

int main() {
  

  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // START OF PARAMETERS

  //
  // Parameters that should come from the cmd line in future
  //

  bool force_refit = false;
  bool force_file_recreation = false;

  std::string p_dimuon_fitvar = "dimuon_mass";
  std::pair<double, double> p_dimuon_fitrange = { 8.7, 11.3 };
  std::string p_chi_fitvar = "chi_mass";
  std::pair<double, double> p_chi_fitrange = { 9.72, 10.2 };

  // Binning - arbitrary number of bin variables
  std::map<std::string, std::pair<double, double> > p_binning;
  p_binning["dimuon_pt"] = { 10, 70 };

  // Preselected flat tuple
  std::string p_inputdata = "first_chib_preselection.root";
  std::string p_inputtree = "data";

  // TODO: parameter for 1P->1S, or mu+/-3sigma


  // Workspace name and output files are handled through an external class
  // that produces filenames just from 

  // Get workspace and output file names 
  // TODO: automatise the folder, file and workspace name generation with a file handler
  std::string h_outputfile = "first_chib_fitting.root"; //input: binning and fitvar
  std::string h_dimuonwsname = "ws_dimuon_pt_10_70"; // input: binning and fitvar
  std::string h_chiwsname = h_dimuonwsname+"_chib1P1S"; // input: dimuonwsname 1P->1S or mass regions, chib or chic


  // Dimuon models

  //Upsilon nS
  std::string upsilon_model_name = "model";
  std::string ups1s_1 = "CBShape::ups1s_1(" + p_dimuon_fitvar + ", mu1s[9.5,9.3,9.7],sigma1s[0.075,0.01,0.25],a_l[0.5,5],5)";
  std::string ups1s_2 = "CBShape::ups1s_2(" + p_dimuon_fitvar + ", mu1s,sigma1s,a_r[-5,-0.5],5)";
  std::string ups1s = "SUM::ups1s(tail_left[0.5,0.,1.]*ups1s_1,ups1s_2)";
  std::string ups2s = "CBShape::ups2s(" + p_dimuon_fitvar + ", expr('mu1s*1.0595',mu1s),sigma2s[0.07,0.01,0.2],a2s[0.5,5],5)"; //factor in expr is mass_pdg_ups1s/mass_pdg_ups2s
  std::string ups3s = "CBShape::ups3s(" + p_dimuon_fitvar + ",  expr('mu1s*1.0946',mu1s),sigma3s[0.07,0.01,0.2],a3s[0.5,5],5)"; //factor in expr is mass_pdg_ups1s/mass_pdg_ups3s
  std::string ups_bg = "Exponential::background(" + p_dimuon_fitvar + ",l[-0.14,-0.9,-0.001])";
  std::string ups_model = "SUM::model("
    "frac1s[0.3,0.,1]*ups1s,"
    "frac2s[0.1,0.,1]*ups2s,"
    "frac3s[0.06,0.,1]*ups3s,"
    "background)";
  std::vector<std::string > upsilon_model = { ups1s_1,ups1s_2, ups1s, ups2s, ups3s, ups_bg, ups_model };


  // Chi models

  // Chib{1,2} 1P
  std::string chib1P_modelname = "model";
  std::string chib1_1P = "CBShape::chib1(" + p_chi_fitvar + ", mu1[9.86,9.902],sigma1[0.001,0.1],a1[0.5,10],n1[1,5])";
  std::string chib2_1P = "CBShape::chib2(" + p_chi_fitvar + ", mu2[9.902,9.94],sigma2[0.001,0.1],a2[0.5,10],n2[1,5])";
  std::string chi1_delta = "expr::delta('" + p_chi_fitvar + " - q0'," + p_chi_fitvar + ",q0[9.68,9.0,15.])";
  std::string chi1_bg = "CEXPR::background('(TMath::Sign( -1, delta) + 1 ) * 0.5 * "
    "pow( fabs(delta), alpha) * exp(beta * delta)',delta, alpha[1.,0.1,15], beta[-1.,-10.,-0.1])";
  std::string chi1_model = "SUM::model("
    "frac1[0.01,0.,1.]*chib1,"
    "frac2[0.01,0.,1.]*chib2,"
    "background)";
  std::vector<std::string > chib1P_model = { chib1_1P, chib2_1P, chi1_delta, chi1_bg, chi1_model };

  // END OF PARAMETERS
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  // Dimuon fit

  Fitter dimuon_fitter;
  dimuon_fitter.SetInputData(p_inputdata, p_inputtree);
  dimuon_fitter.SetOutfile(h_outputfile, force_file_recreation);
  dimuon_fitter.SetWorkspaceName(h_dimuonwsname, force_refit);
  dimuon_fitter.AddVariable(p_chi_fitvar);

  for (const auto &bin : p_binning) {
    const auto &bn = bin.first;
    const auto &border = bin.second;
    dimuon_fitter.AddBinVariable(bn, border.first, border.second);
  }

  dimuon_fitter.SetModel(upsilon_model, upsilon_model_name, p_dimuon_fitvar, p_dimuon_fitrange.first, p_dimuon_fitrange.second);
  dimuon_fitter.Fit(8, true);

  std::pair<double, double> dimuon_cut = { 0,0 };
  {
    FitAnalyser dimuon_analyser(dimuon_fitter);

    std::string dimuon_plot_name(h_outputfile);
    dimuon_plot_name.replace(dimuon_plot_name.find(".root"), 5, "_dimuon_fitresult.pdf");
    dimuon_analyser.PlotFitResult(dimuon_plot_name);


    // TODO: get mu and sigma depending if 1S, 2S or 3S is wanted
    bool ok = false;
    auto mu = dimuon_analyser.GetVariableValue("mu1s", ok);
    if (!ok) std::cout << "Could not get mu from dimuon fit." << std::endl;
    auto sigma = dimuon_analyser.GetVariableValue("sigma1s", ok);
    if (!ok) std::cout << "Could not get sigma from dimuon fit." << std::endl;

    dimuon_cut = { mu - 3 * sigma, mu + 3 * sigma };
  }
  
  // Chi fit

  Fitter chi_fitter;
  chi_fitter.SetInputData(h_outputfile, h_dimuonwsname, Fitter::dataset_name);
  chi_fitter.SetOutfile(h_outputfile);
  chi_fitter.SetWorkspaceName(h_chiwsname, force_refit);


  // Just add binning variables, the binning itself should already be done correctly
  // Mainly to see them on the FitAnalysers plot
  for (const auto &bin : p_binning) dimuon_fitter.AddVariable(bin.first);

  // Cut depending on nS of dimuon fit
  chi_fitter.AddBinVariable(p_dimuon_fitvar, dimuon_cut.first, dimuon_cut.second);


  chi_fitter.SetModel(chib1P_model, chib1P_modelname, p_chi_fitvar, p_chi_fitrange.first, p_chi_fitrange.second);
  chi_fitter.Fit(8, true);
  
  FitAnalyser chi_analyser(chi_fitter);

  std::string chi_plot_name(h_outputfile);
  chi_plot_name.replace(chi_plot_name.find(".root"), 5, "_chib1P_fitresult.pdf"); // TODO: change name to 1P,2P or 3P automatically
  chi_analyser.PlotFitResult(chi_plot_name);
  
  // TODO: Create flat root tuple for polarisation analysis code
}