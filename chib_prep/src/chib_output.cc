#include "TreeMerger.h"

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

int main()
{
  std::string fitresult_filename = "fitresults-dimuon_mass_8p7_11p1-dimuon_pt_20_50.root";
  std::string outfile = "fitresults-dimuon_mass_8p7_11p1-dimuon_pt_20_50_data_with_sWeights.root";

  std::string inputfile, sweightfile, inputtree, sweighttree; // To be read from workspace file
  std::vector<std::string> yield_names;

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

  // Merge input data and sWeights data by EntryID
  ChibSWeightMerger merger(inputfile, sweightfile, inputtree, sweighttree, outfile);
  merger.SetYields(yield_names);
  merger.SetDoubleBranches({ "dimuon_mass", "dimuon_pt", "chi_mass_rf1S", "chi_pt_rf1S",
    "cosTh_HX", "phi_HX", "cosTh_CS", "phi_CS", "cosTh_PX", "phi_PX", "cosAlpha_HX" });
  merger.loop(-1, 1);

}
