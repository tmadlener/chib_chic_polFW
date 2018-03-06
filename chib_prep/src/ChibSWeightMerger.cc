#include "TreeMerger.h"
#include "FitAnalyser.h"

#include "TTree.h"
#include "TFile.h"

#include <iostream>

class ChibSWeightMerger : public TreeMerger
{
public:
  using TreeMerger::TreeMerger;

private:
  virtual ChibSWeightMerger* clone() const override { return new ChibSWeightMerger(*this); }
  virtual void setup_branches() override
  {
    // Choose here which branch from which tree
    m_in_tree->SetBranchAddress("jpsiMass", &jpsiMass);
    m_secondary_in_tree->SetBranchAddress("jpsiPt", &jpsiPt);

    // Set branches for out tree
    m_out_tree->Branch("jpsiMass", &jpsiMass);
    m_out_tree->Branch("jpsiPt", &jpsiPt);
  }

  // Branches to merge
  Double_t jpsiMass;
  Double_t jpsiPt;

};


int main()
{
  std::string workspace_file = "";
  std::string model_name = "";
  std::string fitvar_name = "";
  std::string workspace_name = "";
  std::string file1;

  {
    TFile in(workspace_file.c_str(), "read");
    TNamed *datafilename = nullptr; in.GetObject("InputDataFile", datafilename);
    file1 = datafilename->GetTitle();
    delete datafilename;

    FitAnalyser a(workspace_file, model_name, fitvar_name, workspace_name);
    auto ds = a.GetDataset();

    TFile out =


  }

  // Workspaces and data TTree has to contain same id branches run && event.

  std::string file2 = ""; // First create file containing tree with weights and IDs
  std::string outfile = "";

  ChibSWeightMerger m(file1, file2, "data", "data", outfile, "data", "EntryID");
  m.loop(-1, 8);
}