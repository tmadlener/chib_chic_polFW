#include "ChicTupleEvent_new.h"

#include "general/interface/ArgParser.h"
#include "general/interface/TTreeLooper_new.h"
#include "general/interface/root_utils.h"

#include "TFile.h"
#include "TTree.h"

#include <string>

#if !(defined(__CINT__) or defined(__CLING__))
int main(int argc, char *argv[])
{
  const auto parser = ArgParser(argc, argv);
  const auto inFileName = parser.getOptionVal<std::string>("--infile");
  const auto outFileName = parser.getOptionVal<std::string>("--outfile");
  const auto inTreeName = parser.getOptionVal<std::string>("--intree", "data");
  const auto outTreeName = parser.getOptionVal<std::string>("--outtree", "jpsi_tuple");
  const auto gen = parser.getOptionVal<bool>("--gen", false);
  const auto maxEvents = parser.getOptionVal<int>("--maxEvents", -1);

  auto *inFile = checkOpenFile(inFileName);
  auto *inTree = checkGetFromFile<TTree>(inFile, inTreeName);

  auto *outFile = new TFile(outFileName.c_str(), "recreate");
  auto *tupleTree = new TTree(outTreeName.c_str(), "tupled jpsi data");
  tupleTree->SetDirectory(outFile);

  ChicTupleEvent event; // can also be used for jpsi tupling

  event.Add( // integer variables
            {{"HLT_Dimuon8_Jpsi_v3", "HLT_Dimuon8_Jpsi_v3"},
                {"HLT_Dimuon8_Jpsi_v4", "HLT_Dimuon8_Jpsi_v4"},
                  {"HLT_Dimuon8_Jpsi_v5", "HLT_Dimuon8_Jpsi_v5"},
                    {"HLT_Dimuon8_Jpsi_v6", "HLT_Dimuon8_Jpsi_v6"},
                      {"HLT_Dimuon8_Jpsi_v7", "HLT_Dimuon8_Jpsi_v7"},
                        {"HLT_Dimuon10_Jpsi_v3", "HLT_Dimuon10_Jpsi_v3"},
                          {"HLT_Dimuon10_Jpsi_v4", "HLT_Dimuon10_Jpsi_v4"},
                              {"HLT_Dimuon10_Jpsi_v6", "HLT_Dimuon10_Jpsi_v6"}},

            // double vars
            {{"Jpsict", "Jpsict"}, {"JpsictErr", "JpsictErr"}},
            // TLorentzVector vars
            {{"JpsiP", "Jpsi"}},
            // single muon vars with costh and phi calc
            {std::make_tuple("muPosP", "muNegP", "")}
             );

  if (gen) {
    event.Add({{"HLT_Dimuon10_Jpsi_v5", "HLT_Dimuon10_Jpsi_v5"}},
              {},
              {{"JpsiP_Gen", "gen_Jpsi"}},
              {std::make_tuple("muPosP_Gen", "muNegP_Gen", "gen_")});
  }

  TTreeLooper<ChicTupleEvent> treeLooper(inTree, tupleTree, event);
  treeLooper.loop(event, maxEvents);

  outFile->Write("", TObject::kWriteDelete);
  outFile->Close();

  inFile->Close();

  return 0;

}


#endif
