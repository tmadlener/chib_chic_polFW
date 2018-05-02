#include "interface/ChicTuplingMC.h"

#include "general/interface/ArgParser.h"
#include "general/interface/TTreeLooper.h"
#include "general/interface/root_utils.h"

#include "TFile.h"
#include "TTree.h"

#include <string>
#include <iostream>


#if !(defined(__CINT__) or defined (__CLING__))
int main(int argc, char *argv[])
{
  ArgParser parser(argc, argv);
  const auto mcfile = parser.getOptionVal<std::string>("--datafile");
  const auto outfile = parser.getOptionVal<std::string>("--outfile");
  const auto treename = parser.getOptionVal<std::string>("--treename", "rootuple/chicTree");
  const auto isGenOnly = parser.getOptionVal<bool>("--isGenOnly", false);


  auto* infile = checkOpenFile(mcfile);
  auto* intree = checkGetFromFile<TTree>(infile, treename);

  auto* ofile = new TFile(outfile.c_str(), "recreate");
  auto* tupleTree = new TTree("chic_mc_tuple", "tupled chic mc data");
  tupleTree->SetDirectory(ofile);

  if (!isGenOnly) {
    TTreeLooper<ChicMCInputEvent, ChicMCTupleEvent> treeLooper(intree, tupleTree);
    treeLooper.loop(chicMCTupling);
  } else {
    TTreeLooper<GenMCChicInputEvent, GenMCChicOutputEvent> treeLooper(intree, tupleTree);
    treeLooper.loop(chicGenMCTupling);
  }


  ofile->Write();
  ofile->Close();

  infile->Close();

  return 0;
}

#endif
