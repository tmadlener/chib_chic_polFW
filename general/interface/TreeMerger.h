#ifndef TREEMERGER_JN_H
#define TREEMERGER_JN_H

/////////////////////////////////////////////////////////////////////////////////
// 
//  TREEMERGER j.n. 2018
//  -----------------------------------------------------------------------------
//  A class for merging two TTrees with different branches and different
//  number of entries through an ID.
//
//  Only entries that exist in both trees are written to the 
//  output tree.
//
//  A TreeMerger impementation has to set the in and out branches 
//  in the function setup_branches().
//
/////////////////////////////////////////////////////////////////////////////////

#include "ParallelTreeLooper.h"

#include "RtypesCore.h"

class TreeMerger : public ParallelTreeLooper {

public:
  TreeMerger(const std::string & infilename1, const std::string infilename2, const std::string &intreename1, const std::string &intreename2,
    const std::string & outfilename, const std::string &outtreename, const std::string & major_id, const std::string & minor_id = "");
  
  // TODO: implement constructor with trees
  //TreeMerger(TTree* tree1, TTree* tree2, TTree* outtree, const std::string & major_id, const std::string & minor_id = "");

  virtual ~TreeMerger();

  /////////////////////////////
  // HAVE TO BE IMPLEMENTED: //
  /////////////////////////////
  // virtual ParallelTreeLooper* clone() const = 0;
  virtual void setup_branches() = 0; // Setup all branches (in and out) here
  /////////////////////////////

protected:
  TTree * m_secondary_in_tree = nullptr;

private:
  const std::string secondary_filename;
  const std::string secondary_treename;
  TFile * secondary_file = nullptr;
  // It seems that the TTreeIndex documentation is OUTDATED and that the "Int_t << 31" storage has been removed
  // and now two Long64_t ids are stored instead.
  const std::string major_id_name;
  const std::string minor_id_name;
  Long64_t major_id = 0;
  Long64_t minor_id = 0;
  Long64_t sec_major_id = 0; // stored for check
  Long64_t sec_minor_id = 0;

  void init() override;
  bool fill_and_cut_variables() override;

};
#endif //TREEMERGER_JN_H
