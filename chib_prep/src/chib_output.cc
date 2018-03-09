#include "TreeMerger.h"
#include "ChiOrganizer.h"

#include "ArgParser.h"

#include "TTree.h"
#include "TFile.h"

#include <iostream>


class ChibSWeightMerger : public TreeMerger
{
public:
  using csr = const std::string &;
  ChibSWeightMerger(csr input_file, csr sweight_file, csr input_tree, csr sweight_tree, csr outfile) :
    TreeMerger(sweight_file, input_file, sweight_tree, input_tree, outfile, "data", "EntryID") {}

  void SetYields(const std::vector<std::string> & yield_names) { m_yield_names = yield_names; }
  void SetDoubleBranches(const std::vector<std::string> & double_branches_names) { m_double_branches_names = double_branches_names; } // Branches of other types has to be added 'manually'

private:
  std::vector<std::string> m_yield_names;
  std::vector<std::string> m_double_branches_names;
  virtual ChibSWeightMerger* clone() const override { return new ChibSWeightMerger(*this); }
  virtual void setup_branches() override
  {

    // Set branches from input data

    // Setup double branches
    double_vars.reserve(m_double_branches_names.size());
    for (const auto &d : m_double_branches_names)
    {
      double_vars.emplace_back(0);
      m_secondary_in_tree->SetBranchAddress(d.c_str(), &double_vars.back());
      m_out_tree->Branch(d.c_str(), &double_vars.back());
    }

    // Setup sWeight branches
    yield_vars.reserve(m_yield_names.size() * 2);
    for (const auto & y : m_yield_names) {

      auto ysw = y + "_sw";
      yield_vars.emplace_back(0);
      m_in_tree->SetBranchAddress(ysw.c_str(), &yield_vars.back());
      m_out_tree->Branch(ysw.c_str(), &yield_vars.back());

      auto Ly = "L_" + y;
      yield_vars.emplace_back(0);
      m_in_tree->SetBranchAddress(Ly.c_str(), &yield_vars.back());
      m_out_tree->Branch(Ly.c_str(), &yield_vars.back());

    }

  }

  // Branches to merge
  Double_t dimuon_mass;
  Double_t dimuon_pt;

  Double_t chi_mass_rf1S;
  Double_t chi_pt_rf1S;

  Double_t cosTh_HX;
  Double_t phi_HX;
  Double_t cosTh_PX;
  Double_t phi_PX;
  Double_t cosTh_CS;
  Double_t phi_CS;
  Double_t cosAlpha;

  std::vector<Double_t> yield_vars;
  std::vector<Double_t> double_vars;
};

int main(int argc, char **argv) {

  ArgParser parser(argc, argv);

  auto binvarname = parser.getOptionVal < std::string>("--binvar");
  auto bin_min = parser.getOptionVal < double>("--binmin");
  auto bin_max = parser.getOptionVal < double>("--binmax");
  auto model_folder = parser.getOptionVal < std::string>("--folder", "");
  auto configfile = parser.getOptionVal < std::string>("--config", "");

  auto fitresult_filename = parser.getOptionVal<std::string>("--infile", "");
  auto outfile = parser.getOptionVal<std::string>("--outfile", "");
  auto intree = parser.getOptionVal < std::string>("--intree", "data");
  auto outtree = parser.getOptionVal<std::string>("--outtree", "data");
  auto nEvents = parser.getOptionVal<Long64_t>("--events", -1);
  auto nThreads = parser.getOptionVal<Long64_t>("--threads", 8);
  
  ChiOrganizer org(model_folder, configfile);

  std::map<std::string,std::pair<double,double> > binvars;
  binvars[binvarname] = { bin_min,bin_max };
  if (fitresult_filename.empty()) fitresult_filename = org.FileName(binvars);
  if (outfile.empty()) outfile = org.FileName(binvars, "_sWeightsData.root");

  std::string inputfile, sweightfile, inputtree, sweighttree; // To be read from workspace file
  std::vector<std::string> yield_names;

  // Read file and yield names from workspace file
  {
    TFile file(fitresult_filename.c_str(), "read");
    if (file.IsZombie()) return -1;
    TNamed *tmp = nullptr;
    std::vector<std::pair<std::string, std::string *> > to_read{
      {"InputDataFile", &inputfile },
      { "InputDataTree", &inputtree },
      { "OutputDataFile", &sweightfile },
      { "OutputDataTree", &sweighttree }
    };
    for (auto &p : to_read) {
      file.GetObject(p.first.c_str(), tmp);
      if (!tmp) {
        std::cout << "chib_output: TNamed::" << p.first << " not found in fitresult file '" << fitresult_filename << "'." << std::endl;
        return -1;
      }
      *(p.second) = tmp->GetTitle();
    }
    delete tmp;

    // Read yield names
    TList *list;
    file.GetObject("sWeight_yield_names", list);
    for (auto y : *list) {
      yield_names.push_back(((TObjString*)y)->String().Data());
    }

  }


  // Check DataID to be sure that EntryIDs correspond
  {
    std::string DataID_fitresult;
    std::string DataID_input;
    TNamed *id = nullptr;
    TFile fr(fitresult_filename.c_str());
    TFile fi(inputfile.c_str());
    if (!fr.IsZombie() && !fi.IsZombie()) {
      fr.GetObject("DataID", id);
      if (id) DataID_fitresult = id->GetTitle();
      fi.GetObject("DataID", id);
      if (id) DataID_input = id->GetTitle();
    }
    delete id;

    if (DataID_fitresult.empty() || DataID_input.empty() || DataID_fitresult != DataID_input) {
      std::cout << "chib_output: DataIDs from Input and from FitResult file do not match or are empty, stopping." << std::endl;
      return 1;
    }
  }


  // Merge input data and sWeights data by EntryID

  ChibSWeightMerger merger(inputfile, sweightfile, inputtree, sweighttree, outfile);
  merger.SetYields(yield_names);
  merger.SetDoubleBranches({ "dimuon_mass", "dimuon_pt", "chi_mass_rf1S", "chi_pt_rf1S",
    "costh_HX", "phi_HX", "costh_CS", "phi_CS", "costh_PX", "phi_PX", "cosalpha_HX" });
  merger.loop(nEvents, nThreads);

}
