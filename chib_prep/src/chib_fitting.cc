#include "Fitter.h"
#include "FitAnalyser.h"
#include "ChiOrganizer.h"

#include "TROOT.h"

#include <map>
#include <iostream>
#include <fstream>

using strvec = std::vector<std::string>;


int main() {

  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // PARAMETERS


  //
  // Parameters from command line
  //
  // Overall configurations are in the config.json
  // Command-line argument are only used for the parameters that may change for each bin 
  //
  // TODO: implement ArgParser for below parameter
  //

  // Mandatory
  std::map<std::string, std::pair<double, double> > p_binvars_min_max; // one entry per bin var: e.g. dimuon_rap(0,0.6), dimuon_pt(10,20)
  p_binvars_min_max["dimuon_pt"] = { 15, 50 };

  // Optional
  std::string p_model_folder = ""; // If empty the current workingdir is used
  std::string p_configfile = ""; // if empty it is looked for config.json in the model folder

  auto c_force_refit = false;
  bool c_force_file_recreation = false;
  

  //
  // Parameters from config.json
  //

  // MIGHTDO: Implement getter functions for below parameters in ChiOrganizer, e.g. GetModel(); instead of corg.GetConfigParam<strvec>("dimuon_model", ok);
  ChiOrganizer corg(p_model_folder, p_configfile);
  bool ok = false;

  // Input data file and workspace
  std::string c_inputdata = corg.GetConfigParam<std::string>("input_data_file", ok);
  if (!ok) return 1;
  std::string c_inputtree = corg.GetConfigParam<std::string>("input_data_tree", ok);
  if (!ok) return 1;
    
  // Dimuon model specifications
  auto c_dimuon_modelname = corg.GetConfigParam<std::string>("dimuon_modelname", ok);
  if (!ok) return 1;
  auto c_dimuon_fitvar = corg.GetConfigParam<std::string>("dimuon_fitvar", ok);
  if (!ok) return 1;
  double c_dimuon_fitrange_min = corg.GetConfigParam<double>(strvec{ "dimuon_fitrange", "min" }, ok);
  if (!ok) return 1;
  double c_dimuon_fitrange_max = corg.GetConfigParam<double>(strvec{ "dimuon_fitrange", "max" }, ok);
  if (!ok) return 1;
  auto c_dimuon_model = corg.GetConfigParam<strvec>("dimuon_model", ok);
  if (!ok) return 1;

  // Dimuon cuts before chi fit
  std::string c_dimuon_cut_formula_min = corg.GetConfigParam<std::string>(strvec{ "dimuon_cut_formula", "min" }, ok);
  if (!ok) return 1;
  std::string c_dimuon_cut_formula_max = corg.GetConfigParam<std::string>(strvec{ "dimuon_cut_formula", "max" }, ok);
  if (!ok) return 1;

  // Chi model specifications
  auto c_chi_modelname = corg.GetConfigParam<std::string>("chi_modelname", ok);
  if (!ok) return 1;
  std::string c_chi_fitvar = corg.GetConfigParam<std::string>("chi_fitvar", ok);
  if (!ok) return 1;
  double c_chi_fitrange_min = corg.GetConfigParam<double>(strvec{ "chi_fitrange", "min" }, ok);
  if (!ok) return 1;
  double c_chi_fitrange_max = corg.GetConfigParam<double>(strvec{ "chi_fitrange", "max" }, ok);
  if (!ok) return 1;
  auto c_chi_model = corg.GetConfigParam<strvec>("chi_model", ok);
  if (!ok) return 1;
  auto c_extended_chi_fit = corg.GetConfigParam<bool>("extendedmll_chifit", false);


  // END OF PARAMETERS
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////


  // TODO: parameter for 1P->1S, or mu+/-3sigma

  // Get workspace and output file names
  std::string h_outputfile = corg.FileName(c_dimuon_fitvar, c_dimuon_fitrange_min, c_dimuon_fitrange_max, p_binvars_min_max); 
  std::string h_dimuonwsname = corg.WorkspaceName(c_dimuon_fitvar, c_dimuon_fitrange_min, c_dimuon_fitrange_max, p_binvars_min_max);
  std::string h_chiwsname = h_dimuonwsname + "_chib1P1S"; // TODO: dimuonwsname 1P->1S or mass regions, chib or chic
  std::string h_dimuonplotname = corg.FileName(c_dimuon_fitvar, c_dimuon_fitrange_min, c_dimuon_fitrange_max, p_binvars_min_max, "_dimuonfit.pdf");
  std::string h_chiplotname = corg.FileName(c_dimuon_fitvar, c_dimuon_fitrange_min, c_dimuon_fitrange_max, p_binvars_min_max, "_chifit1P1S.pdf");


  // Dimuon fit

  Fitter dimuon_fitter;
  dimuon_fitter.SetInputData(c_inputdata, c_inputtree);
  dimuon_fitter.SetOutfile(h_outputfile, c_force_file_recreation);
  dimuon_fitter.SetWorkspaceName(h_dimuonwsname, c_force_refit);
  dimuon_fitter.AddVariable(c_chi_fitvar);

  for (const auto &bin : p_binvars_min_max) {
    const auto &bn = bin.first;
    const auto &border = bin.second;
    dimuon_fitter.AddBinVariable(bn, border.first, border.second);
  }

  dimuon_fitter.SetModel(c_dimuon_model, c_dimuon_modelname, c_dimuon_fitvar, c_dimuon_fitrange_min, c_dimuon_fitrange_max);
  dimuon_fitter.Fit(8, true);

  std::pair<double, double> dimuon_cut{ 0,0 };
  {
    FitAnalyser dimuon_analyser(dimuon_fitter);

    std::string dimuon_plot_name(h_dimuonplotname);
    dimuon_analyser.PlotFitResult(dimuon_plot_name);
    bool okmin = false;
    bool okmax = false;

    dimuon_cut = { dimuon_analyser.EvaluateFormula(c_dimuon_cut_formula_min, okmin), dimuon_analyser.EvaluateFormula(c_dimuon_cut_formula_max, okmax) };
    if (!okmin || !okmax) {
      std::cout << "Problems evaluating dimuon mass cut for chi fit. \n\t"
        << c_dimuon_cut_formula_min << "\n\t"
        << c_dimuon_cut_formula_max << "\n"
        << "NOT continuing with the chi fit." << std::endl;
      return 1;
    }
  }

  // Chi fit
  Fitter chi_fitter;
  chi_fitter.SetInputData(h_outputfile, h_dimuonwsname, Fitter::dataset_name);
  chi_fitter.SetOutfile(h_outputfile);
  chi_fitter.SetWorkspaceName(h_chiwsname, c_force_refit);


  // Just add binning variables, the binning itself is already applied when the dataset for the dimuon is created
  // Mainly to see the variable ranges on the FitAnalysers plot
  for (const auto &bin : p_binvars_min_max) chi_fitter.AddVariable(bin.first);

  // Cut depending on nS of dimuon fit
  chi_fitter.AddBinVariable(c_dimuon_fitvar, dimuon_cut.first, dimuon_cut.second);


  chi_fitter.SetModel(c_chi_model, c_chi_modelname, c_chi_fitvar, c_chi_fitrange_min, c_chi_fitrange_max);
  chi_fitter.Fit(8, true, c_extended_chi_fit);

  FitAnalyser chi_analyser(chi_fitter);

  std::string chi_plot_name(h_outputfile);
  chi_analyser.PlotFitResult(h_chiplotname);

}