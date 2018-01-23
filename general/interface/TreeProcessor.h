#ifndef TREEPROCESSOR_JN_H
#define TREEPROCESSOR_JN_H

#include "ParallelTreeLooper.h"

#include <memory>
#include <map>
#include <iostream>

#include "TTree.h"
#include "TLorentzVector.h"

class TreeProcessor;

/////////////////////////////////////////////////////////////////////////////////
// 
//  BRANCHHOLDER j.n. 2018
//  -----------------------------------------------------------------------------
//  A Container to handle the TTree branches.
//  Only for internal use with TreeProcessor
//  Usage:
//   1. Add ALL branch names
//   2. Call setBranches() to connect them to the tree
//   3. To get a reference to a variable call getBranchReference(<branchname>);
//
/////////////////////////////////////////////////////////////////////////////////

template<class T>
class BranchHolder {

  friend class TreeProcessor;


private:
  void addBranch(const std::string & branchname, bool copy = true) 
  {
    // No checks here, assuming that all branches are unique and fine
    branches[branchname] = { default_val,copy };
  }

  T & getBranchReference(const std::string & branchname) 
  {
    if (branches.find(branchname) != branches.end()) return branches[branchname].first;
    std::cout <<"Branch '" << branchname << "' not found. Check spelling and type." << std::endl;
    return default_val;
  }

  void setBranches() 
  {
    for (auto & db : branches) {
      auto & name = db.first;
      auto & var = db.second.first;
      auto copy = db.second.second;
      set_branch(name, var, copy);
    }
  }

  BranchHolder(TTree* in, TTree* out) :
    m_in(in),
    m_out(out) 
  {  }

  BranchHolder(const BranchHolder &) = delete;
  BranchHolder & operator= (const BranchHolder &) = delete;

  TTree* m_in = nullptr;
  TTree* m_out = nullptr;
  std::map < std::string, std::pair< T, bool> > branches;

  void set_branch(const std::string & branchname, T &var, bool copy) {

    m_in->SetBranchAddress(branchname.c_str(), &var);
    if (copy) m_out->Branch(branchname.c_str(), &var);
  }

  T default_val = 0; // should be ok as default value for all types at the moment
};


/////////////////////////////////////////////////////////////////////////////////
// 
//  TREEPROCESSOR j.n. 2018
//  -----------------------------------------------------------------------------                                                                         
//  A class for creating a new TTree from an existing one 
//  with the following options:
//      - copy branches 1:1 with just the branch names as input,
//      - make some selection based on branch values and
//      - add and fill new branches.      
//
/////////////////////////////////////////////////////////////////////////////////

class TreeProcessor : public ParallelTreeLooper {
public:
  TreeProcessor(const std::vector<std::string> & infilenames, const std::string & outfilename, const std::string &treename) :
    TreeProcessor(infilenames, outfilename, treename, treename) { }
  TreeProcessor(const std::vector<std::string> & infilenames, const std::string & outfilename, const std::string &intreename, const std::string &outtreename) :
    ParallelTreeLooper(infilenames, intreename, outfilename, outtreename) { }

  virtual ~TreeProcessor();

  void AddBranchesToCopy(const std::vector<std::string> &);
  void AddBranchesNeededOnlyAsInput(const std::vector<std::string> &);
  void process(long long nEvents = -1, int nThreads = 4) { loop(nEvents, nThreads); } // just for naming reasons

protected:
  // Variables storing the returned reference has to be defined with:
  // static thread_local auto &  <variablename> = get_branch<typename>(<branchname>);
  // This function should be called in fill_and_cut_variables() to get the variables needed for filling new entries
  template <typename T>
  const T & get_branch(const std::string & branchname);

private:
  void init() override;

  //////////////////////////////////////////////////////////
  // HAVE TO BE IMPLEMENTED:

  //from ParallelTreeLooper
  //virtual bool fill_and_cut_variables();
  //virtual ParallelTreeLooper* clone() const;
  virtual void setup_new_branches() = 0; // Variables used in fill_and_cut_variables for adding new variables

  // END TO IMPLEMENT
  //////////////////////////////////////////////////////////

  BranchHolder<Double_t> * double_branches = nullptr;
  BranchHolder<Int_t > * int_branches = nullptr;
  BranchHolder<Long64_t > * long64_branches = nullptr;
  BranchHolder<Bool_t > * bool_branches = nullptr;
  BranchHolder<UInt_t > * uint_branches = nullptr;
  BranchHolder<Float_t > * float_branches = nullptr;
  BranchHolder<TLorentzVector * > * lorentz_branches = nullptr;
  BranchHolder<TVector3 * > * vector_branches = nullptr;

  std::map < std::string, bool> branch_names; //true means copy that variable, false means only read that variable
  
  void fill_branch_maps();
  void set_branch_adresses();

};

// Getter functions for each implemented type
template<>
inline const Double_t & TreeProcessor::get_branch<Double_t>(const std::string & branchname)
{
  return double_branches->getBranchReference(branchname);
}

template<>
inline  TLorentzVector * const & TreeProcessor::get_branch<TLorentzVector*>(const std::string & branchname)
{
  return lorentz_branches->getBranchReference(branchname);
}

template<>
inline  TVector3 * const & TreeProcessor::get_branch<TVector3*>(const std::string & branchname)
{
  return vector_branches->getBranchReference(branchname);
}

template<>
inline  Int_t const & TreeProcessor::get_branch<Int_t>(const std::string & branchname)
{
  return int_branches->getBranchReference(branchname);
}
template<>
inline  Long64_t const & TreeProcessor::get_branch<Long64_t>(const std::string & branchname)
{
  return long64_branches->getBranchReference(branchname);
}
template<>
inline  Bool_t const & TreeProcessor::get_branch<Bool_t>(const std::string & branchname)
{
  return bool_branches->getBranchReference(branchname);
}
template<>
inline  UInt_t const & TreeProcessor::get_branch<UInt_t>(const std::string & branchname)
{
  return uint_branches->getBranchReference(branchname);
}
template<>
inline  Float_t const & TreeProcessor::get_branch<Float_t>(const std::string & branchname)
{
  return float_branches->getBranchReference(branchname);
}

#endif //TREEPROCESSOR_JN_H