#include "ChicBasicTupling.h"
#include "ChicInputEvent.h"
#include "ChicTupleEvent.h"

#include "general/interface/TTreeLooper.h"
#include "general/interface/ArgParser.h"
#include "general/interface/root_utils.h"

#include <string>
#include <iostream>

#if !(defined(__CINT__) or defined(__CLING__))
int main(int argc, char *argv[])
{
  ArgParser parser(argc, argv);
  const auto inFileName = parser.getOptionVal<std::string>("--infile");
  const auto outFileName = parser.getOptionVal<std::string>("--outfile");
  const auto inTreeName = parser.getOptionVal<std::string>("--intree", "rootuple/chicTree");
  const auto outTreeName = parser.getOptionVal<std::string>("--outtree", "chic_tuple");

  auto* inFile = checkOpenFile(inFileName);
  auto* inTree = checkGetFromFile<TTree>(inFile, inTreeName);

  auto* outFile = new TFile(outFileName.c_str(), "recreate");
  auto* outTree = new TTree(outTreeName.c_str(), "chic raw data as flat ntuple");
  outTree->SetDirectory(outFile);

  TTreeLooper<ChicBasicTuplingInEvent, ChicBasicTuplingOutEvent> treeLooper(inTree, outTree);
  treeLooper.loop(chicBasicTupling, -1);

  outFile->Write();
  outFile->Close();

  inFile->Close();


  return 0;
}
#endif
