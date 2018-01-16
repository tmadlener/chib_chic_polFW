#include "ParallelTreeLooper.h"
#include "TObject.h"
#include "TTree.h"

////////////////////////////////////////////////////////
// Example for implementation of a ParallelTreeLooper //
////////////////////////////////////////////////////////

class ExampleLooper : public ParallelTreeLooper {
public:
  using ParallelTreeLooper::ParallelTreeLooper;
  ExampleLooper *clone() const override 
  { 
    return new ExampleLooper(*this);
  }
  virtual ~ExampleLooper() { }

  bool fill_and_cut_variables() override
  {
    doubledJpsiMass = 2 * jpsiMass;
  }

private:
  Double_t jpsiMass;
  Double_t doubledJpsiMass;

  void init() override {
    if (hasError()) return;
    m_in_tree->SetBranchAddress("jpsiMass", &jpsiMass);
    m_out_tree->Branch("jpsiMass", &jpsiMass);
    m_out_tree->Branch("doubledJpsiMass", &doubledJpsiMass);
  }
};

int main(int argc, char** argv) 
{
  ExampleLooper l({ "/afs/hephy.at/work/j/jnecker/data/bug_study/jpsi_2016.root" },"data","paralooper_test_output.root","data");
  l.loop(1000000, 4);
}