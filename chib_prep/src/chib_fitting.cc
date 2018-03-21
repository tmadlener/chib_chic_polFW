#include "Fitter.h"
#include "FitAnalyser.h"
#include "ChiOrganizer.h"

#include "ArgParser.h"

#include "TROOT.h"
#include "TUUID.h"
#include "TTree.h"
#include "RooDataSet.h"
#include "TList.h"

#include <map>
#include <iostream>
#include <fstream>

using strvec = std::vector<std::string>;

int main(int argc, char **argv) {


  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // PARAMETERS


  //
  // Parameters from command line
  //
  // Overall configurations are in the config.json
  // Command-line argument are only used for the parameters that may change for each bin 
  //

  ArgParser parser(argc, argv);

  // TODO: eventually out file argument for testing purposes

  // Mandatory
  auto binvarname = parser.getOptionVal < std::string>("--binvar");
  auto bin_min = parser.getOptionVal < double>("--binmin");
  auto bin_max = parser.getOptionVal < double>("--binmax");

  std::map<std::string, std::pair<double, double> > p_binvars_min_max; // one entry per bin var: e.g. dimuon_rap(0,0.6), dimuon_pt(10,20)
  if (!binvarname.empty()) p_binvars_min_max[binvarname] = { bin_min, bin_max };

  // Optional
  auto p_model_folder = parser.getOptionVal < std::string>("--folder", ""); // If empty the current workingdir is used
  auto p_configfile = parser.getOptionVal < std::string>("--config", "");// if empty it is looked for config.json in the model folder
  auto p_outfile = parser.getOptionVal < std::string>("--outfile", "");// mainly for testing purposes

  auto c_force_file_recreation = parser.getOptionVal <bool>("--forceFileRecreation", false);
  auto c_force_refit_dimuon = parser.getOptionVal <bool>("--forceDimuonRefit", false);
  auto c_force_refit_chi = parser.getOptionVal <bool>("--forceChiRefit", false);


  //
  // Parameters from config.json
  //

  // MIGHTDO: Implement getter functions for below parameters in ChiOrganizer, e.g. GetModel(); instead of corg.GetConfigParam<strvec>("dimuon_model", ok);
  ChiOrganizer corg(p_model_folder, p_configfile);
  bool ok = false;

  // Input data file and workspace
  auto c_inputdata = corg.GetConfigParam<std::string>("input_data_file", ok);
  if (!ok) return 1;
  auto c_inputtree = corg.GetConfigParam<std::string>("input_data_tree", "data");

  // Additional cuts
  auto c_additional_cutvars = corg.GetVariableList("cut_variables");
  for (auto &c : c_additional_cutvars) std::cout << c.first << ": " << c.second.first << ", " << c.second.second << std::endl;

  // Dimuon model specifications
  auto c_dimuon_modelname = corg.GetConfigParam<std::string>("dimuon_modelname", ok);
  if (!ok) return 1;
  auto c_dimuon_fitvar = corg.GetConfigParam<std::string>(strvec{ "dimuon_fitvar","name" }, ok);
  if (!ok) return 1;
  double c_dimuon_fitrange_min = corg.GetConfigParam<double>(strvec{ "dimuon_fitvar", "min" }, ok);
  if (!ok) return 1;
  double c_dimuon_fitrange_max = corg.GetConfigParam<double>(strvec{ "dimuon_fitvar", "max" }, ok);
  if (!ok) return 1;
  auto c_dimuon_model = corg.GetConfigParam<strvec>("dimuon_model", ok);
  if (!ok) return 1;

  // Dimuon cuts before chi fit
  auto c_dimuon_cut_formula_min = corg.GetConfigParam<std::string>(strvec{ "dimuon_cut_formula", "min" }, ok);
  if (!ok) return 1;
  auto c_dimuon_cut_formula_max = corg.GetConfigParam<std::string>(strvec{ "dimuon_cut_formula", "max" }, ok);
  if (!ok) return 1;

  // Chi model specifications
  auto c_chi_modelname = corg.GetConfigParam<std::string>("chi_modelname", ok);
  if (!ok) return 1;
  auto c_chi_fitvar = corg.GetConfigParam<std::string>(strvec{ "chi_fitvar","name" }, ok);
  if (!ok) return 1;
  double c_chi_fitrange_min = corg.GetConfigParam<double>(strvec{ "chi_fitvar", "min" }, ok);
  if (!ok) return 1;
  double c_chi_fitrange_max = corg.GetConfigParam<double>(strvec{ "chi_fitvar", "max" }, ok);
  if (!ok) return 1;
  auto c_chi_model = corg.GetConfigParam<strvec>("chi_model", ok);
  if (!ok) return 1;
  auto c_extended_chi_fit = corg.GetConfigParam<bool>("extendedmll_chifit", false);
  auto c_sweight_yields = corg.GetConfigParam<strvec>("sweight_yields", strvec());

  // END OF PARAMETERS
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////


  // TODO: parameter for 1P->1S, or mu+/-3sigma

  // Get workspace and output file names
  std::string out_file_base = corg.FileName(c_dimuon_fitvar, c_dimuon_fitrange_min, c_dimuon_fitrange_max, p_binvars_min_max, "");
  if (!p_outfile.empty()) {
    std::regex e(R"(\.root\b)");
    out_file_base = std::regex_replace(p_outfile, e, "");
  }

  std::string h_outputfile = out_file_base + ".root";
  std::string h_dimuonwsname = corg.WorkspaceName(c_dimuon_fitvar, c_dimuon_fitrange_min, c_dimuon_fitrange_max, p_binvars_min_max);
  // TODO: dimuonwsname 1P->1S or mass regions, chib or chic
  std::string h_chiwsname = h_dimuonwsname + "_chib1P1S";
  std::string h_dimuonplotname = out_file_base + "_dimuonfit.pdf";
  std::string h_chiplotname = out_file_base + "_chifit1P1S.pdf";
  std::string h_resultfilename = out_file_base + "_chifit1P1S_outputtree.root";

  std::string DataID;
  {
    TFile in(c_inputdata.c_str(), "read");
    if (!in.IsZombie()) {
      TNamed *id = nullptr;
      in.GetObject("DataID", id);
      if (id) DataID = id->GetTitle();
      delete id;
    }
  }

  // Check if ID of data file changed, if yes force a refit
  bool data_FileID_changed = false;

  if (!c_force_file_recreation) {
    TFile out(h_outputfile.c_str(), "read");
    if (!out.IsZombie()) {
      TNamed *id = nullptr;
      std::string id_out;
      out.GetObject("DataID", id);
      if (id) id_out = id->GetTitle();
      delete id;

      // Check ids, only force recreation of file if both ids are available to not accidently delete fitresults
      if ((!id_out.empty() && !DataID.empty()) && DataID != id_out) {
        std::cout << "chib_fitting: DataID changed, forcing complete refit!" << std::endl;
        data_FileID_changed = true;
      }
    }
  }

  // Dimuon fit

  Fitter dimuon_fitter;
  dimuon_fitter.SetInputData(c_inputdata, c_inputtree);
  dimuon_fitter.SetOutfile(h_outputfile, c_force_file_recreation || data_FileID_changed);
  dimuon_fitter.SetWorkspaceName(h_dimuonwsname, c_force_refit_dimuon);
  dimuon_fitter.AddVariable(c_chi_fitvar);

  // Additional cuts
  for (auto &cut : c_additional_cutvars) dimuon_fitter.AddBinVariable(cut.first, cut.second.first, cut.second.second);

  // Needed for later matching of weights to events:
  dimuon_fitter.AddVariable("EntryID_low");
  dimuon_fitter.AddVariable("EntryID_high");

  for (const auto &bin : p_binvars_min_max) {
    const auto &bn = bin.first;
    const auto &border = bin.second;
    dimuon_fitter.AddBinVariable(bn, border.first, border.second);
  }

  dimuon_fitter.SetModel(c_dimuon_model, c_dimuon_modelname, c_dimuon_fitvar, c_dimuon_fitrange_min, c_dimuon_fitrange_max);
  dimuon_fitter.Fit(8, true);

  // Add input data DataID to workspace file
  {
    TFile f(h_outputfile.c_str(), "update");
    TFile in(c_inputdata.c_str(), "read");
    if (!f.IsZombie() && !in.IsZombie()) {
      // Add file and tree names
      TList file_list;
      file_list.SetOwner(true);
      file_list.Add(new TNamed("InputDataFile", c_inputdata));
      file_list.Add(new TNamed("InputDataTree", c_inputtree));
      file_list.Add(new TNamed("OutputDataFile", h_resultfilename));
      file_list.Add(new TNamed("OutputDataTree", "data"));
      TNamed *id = nullptr; in.GetObject("DataID", id);
      file_list.Add(id);
      f.cd();
      file_list.Write(0, TObject::kWriteDelete);

      // Add yield variable names
      if (!c_sweight_yields.empty()) {
        TList yield_list;
        yield_list.SetOwner(true); // now the list is responsible for deleting its objects
        yield_list.SetName("sWeight_yield_names");
        for (auto &y : c_sweight_yields) yield_list.Add(new TObjString(y.c_str()));
        f.cd();
        yield_list.Write(0, TObject::kWriteDelete | TObject::kSingleKey);
      }
    }
  }

  std::pair<double, double> dimuon_cut{ 0,0 };
  {
    FitAnalyser dimuon_analyser(dimuon_fitter);

    std::string dimuon_plot_name(h_dimuonplotname);
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

    dimuon_analyser.PlotFitResult(dimuon_plot_name, { dimuon_cut.first, dimuon_cut.second });
  }

  // Chi fit
  Fitter chi_fitter;
  chi_fitter.SetInputData(h_outputfile, h_dimuonwsname, Fitter::dataset_name);
  chi_fitter.SetOutfile(h_outputfile);
  chi_fitter.SetWorkspaceName(h_chiwsname, c_force_refit_chi);

  // Needed for later matching of weights to events:
  chi_fitter.AddVariable("EntryID_low");
  chi_fitter.AddVariable("EntryID_high");

  // Just add binning variables, the binning itself is already applied when the dataset for the dimuon is created
  // Mainly to see the variable ranges on the FitAnalysers plot
  for (const auto &bin : p_binvars_min_max) chi_fitter.AddVariable(bin.first);

  // Cut depending on nS of dimuon fit
  chi_fitter.AddBinVariable(c_dimuon_fitvar, dimuon_cut.first, dimuon_cut.second);


  chi_fitter.SetModel(c_chi_model, c_chi_modelname, c_chi_fitvar, c_chi_fitrange_min, c_chi_fitrange_max);
  chi_fitter.Fit(8, true, c_extended_chi_fit);

  // Create output file with tree containing sWeights

  FitAnalyser chi_analyser(chi_fitter);

  std::string chi_plot_name(h_outputfile);
  chi_analyser.PlotFitResult(h_chiplotname);

  chi_analyser.AddSWeights(c_sweight_yields);

  //Create file containing tree with the sweights
  auto sdata = chi_analyser.GetDataset();

  RooAbsData::setDefaultStorageType(RooAbsData::Tree);
  RooDataSet * sweight_data = new RooDataSet("data", "data with sWeights", sdata, *sdata->get());
  {
    TFile f(h_resultfilename.c_str(), "recreate");
    auto stree = sweight_data->tree();
    stree->Write(0, TObject::kWriteDelete);
  }
  // Now add original EntryID to TTree (maybe move this code to chib_output.cc)
  {
    TFile f(h_resultfilename.c_str(), "update");
    TTree* t(nullptr);
    if (!f.IsZombie() && (t = (TTree*)f.Get("data"))) {

      Long64_t EntryID = 0,
        long_e_high = 0,
        lowmask = 0xFFFFFFFF;
      Double_t e_low = 0, e_high = 0;

      t->SetBranchAddress("EntryID_low", &e_low);
      t->SetBranchAddress("EntryID_high", &e_high);
      auto idbranch = t->Branch("EntryID", &EntryID);

      for (Long64_t e = 0, s = t->GetEntries(); e < s; ++e) {
        t->GetEntry(e);
        long_e_high = e_high;
        long_e_high <<= 32;
        EntryID = e_low;
        EntryID &= lowmask;
        EntryID |= long_e_high;
        idbranch->Fill();
      }

      // Write to file
      f.cd();
      t->Write(0, TObject::kWriteDelete);
      TNamed fid("DataID", DataID.c_str());
      fid.Write(0, TObject::kWriteDelete);
    }
  }


}