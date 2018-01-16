#include "Region.h"
#include "BkgWeightCalcChic.h" // keep this before ChicTuplingWithWeights (forward declarations)
#include "ChicTuplingWithWeights.h"
#include "ChicTupleEvent.h"
#include "ChicInputEvent.h"
#include "JpsiMassSelector.h"

#include "config/commonVar.h"

#include "general/interface/TTreeLooper.h"

#include "general/interface/ArgParser.h"
#include "general/interface/root_utils.h"
#include "general/interface/misc_utils.h"

#include "RooWorkspace.h"

#include <string>
#include <iostream>
#include <memory>

std::unique_ptr<DiMuonSelector> getDimuonSelector(const bool rejSeagulls, const bool rejCowboys)
{
  if (rejCowboys) {
    return std::unique_ptr<CowboySelector>(new CowboySelector());
  }
  if (rejSeagulls) {
    return std::unique_ptr<SeagullSelector>(new SeagullSelector());
  }

  return std::unique_ptr<AllSelector>(new AllSelector());
}

#if !(defined(__CINT__) or defined(__CLING__))
int main(int argc, char *argv[])
{
  ArgParser parser(argc, argv);
  const auto fitFileName = parser.getOptionVal<std::string>("--fitfile");
  const auto dataFileName = parser.getOptionVal<std::string>("--datafile");
  const auto outputFileName = parser.getOptionVal<std::string>("--outfile");
  const auto dataTreeName = parser.getOptionVal<std::string>("--tree", "selectedData");
  const auto outTreeName = parser.getOptionVal<std::string>("--otree", "chic_tuple");
  const auto rejCowboys = parser.getOptionVal<bool>("--rejCowboys", false);
  const auto rejSeagulls = parser.getOptionVal<bool>("--rejSeagulls", false);

  // do this check first to not do any unnecessary work
  const int ptBin = getPtBinFromFile(fitFileName);
  if (ptBin < 0) {
    std::cerr << "Can't deduce pt bin from\'" << fitFileName << "\'.\n";
    return 1;
  }

  auto* fitFile = checkOpenFile(fitFileName);
  auto* ws = checkGetFromFile<RooWorkspace>(fitFile, "ws_masslifetime");

  const auto binname = getBinFromFile(fitFileName);
  const auto snapname = "snapshot_" + binname;

  auto ltRegions = calcLifeTimeRegionsJpsi(ws, "data_" + binname + "_SR");
  std::cout << "Lifetime regions:\n" << ltRegions.PR << "\n" << ltRegions.NP << "\n";

  auto mRegions = calcMassRegions(ws, snapname);
  std::cout << "Mass regions:\n"
            << mRegions.LSB << "\n" << mRegions.RSB << "\n"
            << mRegions.SR1 << "\n" << mRegions.SR2 << "\n";

  auto bkgWeights = calculate2DWeights(ws, mRegions, ltRegions, snapname);

  std::cout << "Weights:" << bkgWeights.wChic1 << ", " << bkgWeights.wChic2 << "\n";

  JpsiMassSelector jpsiSelector(ws, 3.0);

  fitFile->Close();

  auto* dataFile = checkOpenFile(dataFileName);
  auto* dataTree = checkGetFromFile<TTree>(dataFile, dataTreeName);

  auto* outFile = new TFile(outputFileName.c_str(), "recreate");
  auto* tupleTree = new TTree(outTreeName.c_str(), "tupled chic data with 2 dim mass weights");
  tupleTree->SetDirectory(outFile);

  const double ptMin = config::ptBinBorders[ptBin - 1];
  const double ptMax = config::ptBinBorders[ptBin];

  const auto dimuonSelector = getDimuonSelector(rejSeagulls, rejCowboys);

  TTreeLooper<ChicInputEvent<>, ChicTupleEvent<>> treeLooper(dataTree, tupleTree);
  auto tuplingFunc = [&] (const ChicInputEvent<>& ie, ChicTupleEvent<>& e) {
    return chicTuplingWith2DWeights(ie, e, mRegions, ltRegions, bkgWeights, ptMin, ptMax,
                                    config::maxAbsRap, jpsiSelector, *dimuonSelector);
  };

  treeLooper.loop(tuplingFunc);

  outFile->Write();
  outFile->Close();

  dataFile->Close();

  return 0;
}

#endif
