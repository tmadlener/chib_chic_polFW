#include "TreeMerger.h"
#include "TTree.h"

class SimpleTreeMerger : public TreeMerger
{
public:
  using TreeMerger::TreeMerger;

private:
  virtual SimpleTreeMerger* clone() const override { return new SimpleTreeMerger(*this); } // In this example we can use the default copy constructor
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
  SimpleTreeMerger m("merger_testfile1.root", "merger_testfile2.root", "data", "data", "merger_example_output.root", "data", "entry_id");
  m.loop(-1,5);
}