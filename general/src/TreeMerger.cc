#include "TreeMerger.h"

#include "TTree.h"
#include "TFile.h"

#include <iostream>

TreeMerger::TreeMerger(const std::string & infilename1, const std::string infilename2, const std::string & intreename1, const std::string & intreename2, 
  const std::string & outfilename, const std::string & outtreename,const std::string & major_id, const std::string & minor_id) :
  ParallelTreeLooper({ infilename1 }, intreename1, outfilename, outtreename),
  secondary_filename(infilename2),
  secondary_treename(intreename2),
  major_id_name(major_id),
  minor_id_name(minor_id)
{ }

TreeMerger::~TreeMerger()
{
  delete secondary_file;
}

void TreeMerger::init()
{
  if (hasError()) return;

  // Open secondary file
  secondary_file = TFile::Open(secondary_filename.c_str(), "read");
  if (!secondary_file) {
    std::cout << "TreeMerger: Could not open file '" << secondary_filename << "'." << std::endl;
    return;
  }
  m_secondary_in_tree = dynamic_cast<TTree*>(secondary_file->Get(secondary_treename.c_str()));
  if (!m_secondary_in_tree) {
    std::cout << "TreeMerger: No tree '" << secondary_treename << "' in file '" << secondary_filename << "'." << std::endl;
    return;
  }
  
  // Major id
  m_in_tree->SetBranchAddress(major_id_name.c_str(), &major_id);
  m_secondary_in_tree->SetBranchAddress(major_id_name.c_str(), &sec_major_id);
  m_out_tree->Branch(major_id_name.c_str(), &major_id);
  m_out_tree->Branch((major_id_name+"_check").c_str(), &sec_major_id);

  // Minor id
  if (!minor_id_name.empty()) {
    m_in_tree->SetBranchAddress(minor_id_name.c_str(), &minor_id);
    m_secondary_in_tree->SetBranchAddress(minor_id_name.c_str(), &sec_minor_id);
    m_out_tree->Branch(minor_id_name.c_str(), &minor_id);
    m_out_tree->Branch((minor_id_name + "_check").c_str(), &sec_minor_id);
  }

  // Build index
  m_in_tree->BuildIndex(major_id_name.c_str(), minor_id_name.empty() ? "0" : minor_id_name.c_str());
  m_secondary_in_tree->BuildIndex(major_id_name.c_str(), minor_id_name.empty() ? "0" : minor_id_name.c_str());

  setup_branches();
}

bool TreeMerger::fill_and_cut_variables()
{
  if (!m_secondary_in_tree) return false;

  // check if secondary tree has an entry with the index 
  auto entry = m_secondary_in_tree->GetEntryNumberWithIndex(major_id, minor_id);
  auto bytes = m_secondary_in_tree->GetEntry(entry);
  
  return  bytes > 0;
};