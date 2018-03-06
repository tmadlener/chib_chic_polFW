#include "ParallelTreeLooper.h"

#include "TObject.h"
#include "TTree.h"
#include "TLorentzVector.h"

#include <thread>
///////////////////////////////////////////////////
// Example for implementing a ParallelTreeLooper //
///////////////////////////////////////////////////


class ExampleLooper : public ParallelTreeLooper {
public:
  using ParallelTreeLooper::ParallelTreeLooper;
  ExampleLooper *clone() const override 
  { 
    return new ExampleLooper(*this);
  }
  virtual ~ExampleLooper()
  {
    delete jpsi;
  }

  bool fill_and_cut_variables() override
  {
    doubledJpsiMass = 2 * jpsiMass;
    jpsiPt = jpsi->Pt();
    return true;
  }

private:

  Double_t jpsiMass = 0;
  Double_t doubledJpsiMass = 0;
  Double_t jpsiPt = 0;
  TLorentzVector* jpsi = nullptr;

  void init() override {
    if (hasError()) return;

    // set in branches
    m_in_tree->SetBranchAddress("jpsiMass", &jpsiMass);
    m_in_tree->SetBranchAddress("JpsiP", &jpsi);

    // set out branches
    m_out_tree->Branch("jpsiMass", &jpsiMass);
    //m_out_tree->Branch("doubledJpsiMass", &doubledJpsiMass);
    //m_out_tree->Branch("jpsiPt", &jpsiPt);


  }
};


int main(int argc, char** argv) 
{
  ExampleLooper l({ "/afs/hephy.at/work/j/jnecker/data/bug_study/jpsi_2016.root" },"data","paralooper_testfile.root","data");
  l.loop(100000, 4);
}