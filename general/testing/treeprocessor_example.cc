#include "TreeProcessor.h"
#include <cmath>

class SimpleTreeCopier : public TreeProcessor
{
public:
  using TreeProcessor::TreeProcessor;

private:
  // The following three functions always have to be implemented

  virtual bool fill_and_cut_variables() override { return true; } // No selection in this example
  virtual ParallelTreeLooper* clone() const override { return new SimpleTreeCopier(*this); } // In this example we can use the default copy constructor
  virtual void setup_new_branches() override { } // No new branches in this example

};

class SimpleBranchAdder : public TreeProcessor
{
public:
  using TreeProcessor::TreeProcessor;

private:
  virtual bool fill_and_cut_variables() override {

    // Always get the branches used for new branches or cuts in this function in that way
    static thread_local auto & cosTh = get_branch<Double_t>("cosTh_HX");
    static thread_local auto & jpsi = get_branch<TLorentzVector*>("JpsiP");

    // Fill the new branches
    Theta = acos(cosTh);
    jpsiPhi = jpsi->Phi();

    // Also cuts can be defined here, simple return false if the event should be rejected
    return true; 
  }

  virtual ParallelTreeLooper* clone() const override {
    // The class has to be cloned for parallel processing,
    // so if you want to use parallel processing and
    // the default copy constructor won't work for your implementation 
    // you can handle this stuff here before returning the new instance.

    return new SimpleBranchAdder(*this);
  }

  virtual void setup_new_branches() override {
    // New Branches has to be added here

    m_out_tree->Branch("Theta_HX", &Theta);
    m_out_tree->Branch("jpsiPhi", &jpsiPhi);
  }

  // Variables for new branches
  Double_t Theta = 0;
  Double_t jpsiPhi = 0;

};

int main()
{
  SimpleTreeCopier cp({ "/afs/hephy.at/work/j/jnecker/data/bug_study/jpsi_2016.root" }, "treecopier_testoutput.root", "data");
  cp.AddBranchesToCopy({ "jpsiMass", "cosTh_HX" });
  cp.process(100000,4);

  SimpleBranchAdder ba({ "/afs/hephy.at/work/j/jnecker/data/bug_study/jpsi_2016.root" }, "branchadder_testoutput.root", "data");
  ba.AddBranchesToCopy({ "jpsiMass", "cosTh_HX", "JpsiP" });
  ba.process(100000,4);
}