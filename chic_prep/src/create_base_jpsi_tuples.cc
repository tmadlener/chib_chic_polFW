#include "JpsiBasicTupling.h"

#include "general/interface/ArgParser.h"
#include "general/interface/TTreeLooper.h"
#include "general/interface/root_utils.h"

#include "TObject.h"

#include <string>
#include <vector>

#if !(defined(__CLING__) or defined(__CINT__))
int main(int argc, char *argv[])
{
  ArgParser parser(argc, argv);
  const auto inputFiles = parser.getOptionVal<std::vector<std::string>>("--files");
  const auto outfile = parser.getOptionVal<std::string>("--outfile");
  const auto intree = parser.getOptionVal<std::string>("--intree", "data");
  const auto outtree = parser.getOptionVal<std::string>("--outtree", "jpsi_tuple");
  const auto isMC = parser.getOptionVal<bool>("--mc", false);
  const auto nEvents = parser.getOptionVal<int>("--nevents", -1);

  auto *inTree = createTChain(inputFiles, intree);
  auto *outFile = new TFile(outfile.c_str(), "recreate");
  auto *outTree = new TTree(outtree.c_str(), "jpsi raw data as flat ntuple");
  outTree->SetDirectory(outFile);


  TTreeLooper<JpsiBasicTuplingInEvent, JpsiBasicTuplingOutEvent> treeLooper(inTree, outTree);
  if (isMC) {
    treeLooper.loop(BasicJpsiMCTupling, nEvents);
  } else {
    treeLooper.loop(BasicJpsiTupling, nEvents);
  }

  outFile->Write("", TObject::kWriteDelete);
  outFile->Close();

  return 0;
}

#endif
