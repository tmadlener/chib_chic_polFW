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
  const auto inTreeName = parser.getOptionVal<std::string>("--intree", "rootuple/chicTree");
  const auto outTreeName = parser.getOptionVal<std::string>("--outtree", "chic_tuple");
  const auto gen = parser.getOptionVal<bool>("--gen", false);
  const auto genOnly = parser.getOptionVal<bool>("--genonly", false);
  const auto maxEvents = parser.getOptionVal<int>("--nEvents", -1);

  auto *inFile = checkOpenFile(inFileName);
  auto *inTree = checkGetFromFile<TTree>(inFile, inTreeName);

  auto *outFile = new TFile(outFileName.c_str(), "recreate");
  auto *tupleTree = new TTree(outTreeName.c_str(), "tupled chic data");
  tupleTree->SetDirectory(outFile);

  ChicTupleEvent event;

  // Declare the basic event content (valid for data and reco MC)
  if (!genOnly) {
    event.Add({{"trigger", "trigger"}}, // integer variables
              // double vars
              {{"ctpv", "Jpsict"}, {"conv_vertex", "convRadius"}, {"dz", "gammaDz"}, {"probFit1S", "vtxProb"},
               {"ctpv_error", "JpsictErr"}}, // double vars
              {{"dimuon_p4", "Jpsi"}, {"rf1S_chi_p4", "chic"}, {"photon_p4", "photon"},
               {"chi_p4", "mumugamma"}}, // TLorentzVector vars
              {std::make_tuple("muonP_p4", "muonN_p4", "")}); // single muon vars with costh and phi calc
              // {std::make_tuple("muonP_p4", "muonM_p4", "")}); // 2017
  }

  if (gen || genOnly) {
    event.Add({{"chic_pdgId", "pdgId"}},
              {},
              {{"gen_chic_p4", "gen_chic"}, {"gen_jpsi_p4", "gen_Jpsi"}, {"gen_photon_p4", "gen_photon"}},
              {std::make_tuple("gen_muonP_p4", "gen_muonM_p4", "gen_")});
  }

  TTreeLooper<ChicTupleEvent> treeLooper(inTree, tupleTree, event);
  treeLooper.loop(event, maxEvents);

  outFile->Write("", TObject::kWriteDelete);
  outFile->Close();

  inFile->Close();

  return 0;
}

#endif
