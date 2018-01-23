#include "TreeProcessor.h"

#include <algorithm>
#include <iostream>

#include "TTree.h"
#include "TLeaf.h"
#include "TLorentzVector.h"

void TreeProcessor::AddBranchesToCopy(const std::vector<std::string>&branchnames)
{
  // Already existing branches are set to true, that they are copied
  for (const auto &bn : branchnames) branch_names[bn] = true;
}

void TreeProcessor::AddBranchesNeededOnlyAsInput(const std::vector<std::string>& branchnames)
{
  // Check if branch already exists to not overwrite 'true'
  for (const auto &bn : branchnames)
    if (branch_names.find(bn) == branch_names.end()) branch_names[bn] = false;
}

void TreeProcessor::init()
{
  if (hasError()) return;

  // Activating only needed branches
  m_in_tree->SetBranchStatus("*", 0);
  for (const auto &b : branch_names) m_in_tree->SetBranchStatus(b.first.c_str(), 1);

  // Create variables of correct type for each branchname
  fill_branch_maps();

  // Set in and out branch addresses
  set_branch_adresses();

  // Setup new out branches
  setup_new_branches();

}

void TreeProcessor::fill_branch_maps()
{

  // Create BranchHolders
  double_branches = new BranchHolder<Double_t>(m_in_tree, m_out_tree);
  int_branches = new BranchHolder<Int_t>(m_in_tree, m_out_tree);
  long64_branches = new BranchHolder<Long64_t>(m_in_tree, m_out_tree);
  bool_branches = new BranchHolder<Bool_t>(m_in_tree, m_out_tree);
  uint_branches = new BranchHolder<UInt_t>(m_in_tree, m_out_tree);
  float_branches = new BranchHolder<Float_t>(m_in_tree, m_out_tree);
  lorentz_branches = new BranchHolder<TLorentzVector*>(m_in_tree, m_out_tree);
  vector_branches = new BranchHolder<TVector3*>(m_in_tree, m_out_tree);

  for (const auto &bi : branch_names) {

    // Get type
    auto b = dynamic_cast<TBranch*>(m_in_tree->GetBranch(bi.first.c_str()));
    if (!b) {
      std::cout << "Branch '" << bi.first << "' not found in '" << m_in_tree->GetName() << "'." << std::endl;
      continue;
    }
    std::string typestr = b->GetClassName();
    if (typestr.empty()) {
      auto l = b->GetLeaf(bi.first.c_str());
      if (l) typestr = l->GetTypeName();
    }

    // Add to correct BranchHolder
    if (typestr == "Double_t") double_branches->addBranch(bi.first, bi.second);
    else if (typestr == "Int_t") int_branches->addBranch(bi.first, bi.second);
    else if (typestr == "Long64_t") long64_branches->addBranch(bi.first, bi.second);
    else if (typestr == "Bool_t")  bool_branches->addBranch(bi.first, bi.second);
    else if (typestr == "UInt_t") uint_branches->addBranch(bi.first, bi.second);
    else if (typestr == "Float_t") float_branches->addBranch(bi.first, bi.second);
    else if (typestr == "TLorentzVector") lorentz_branches->addBranch(bi.first, bi.second);
    else if (typestr == "TVector3") vector_branches->addBranch(bi.first, bi.second);
    else std::cout << "Type '" << typestr << "' of branch '" << bi.first << "' is not implemented." << std::endl;
  }
}

void TreeProcessor::set_branch_adresses()
{
  double_branches->setBranches();
  int_branches->setBranches();
  long64_branches->setBranches();
  bool_branches->setBranches();
  uint_branches->setBranches();
  float_branches->setBranches();
  lorentz_branches->setBranches();
  vector_branches->setBranches();
}

TreeProcessor::~TreeProcessor()
{
  // Because TTree does not clean the objects we have to do this

  // Citation from : https://root.cern.ch/doc/master/classTTree.html
  // Whether the pointer is set to zero or not, the ownership of the object is not taken over by the TTree. 
  // I.e. even though an object will be allocated by TTree::Branch if the pointer p_object is zero, 
  // the object will not be deleted when the TTree is deleted.

  for (auto & b : lorentz_branches->branches) delete b.second.first;
  for (auto &b : vector_branches->branches) delete b.second.first;

  delete double_branches;
  delete int_branches;
  delete long64_branches;
  delete bool_branches;
  delete uint_branches;
  delete float_branches;
  delete lorentz_branches;

}