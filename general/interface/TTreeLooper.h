#ifndef CHIBCHICPOLFW_GENERAL_TTREELOOPER_H__
#define CHIBCHICPOLFW_GENERAL_TTREELOOPER_H__

#include "root_utils.h"

#include "progress.h"

#include "TFile.h"
#include "TTree.h"

#include <string>
#include <iostream>
#include <chrono>

template<typename InEventT, typename OutEventT>
class TTreeLooper {
public:
  TTreeLooper() = delete; /** do not want a default constructor */

  /** ctor. */
  TTreeLooper(TTree* inTree, TTree* outTree);

  template<typename CondF, PrintStyle PS = PrintStyle::ProgressBar>
  void loop(CondF cond, const long int maxEvents = -1);

private:

  InEventT m_inEvent;

  OutEventT m_outEvent;

  TTree* m_inTree{nullptr};

  // TFile* m_outFile{nullptr};

  TTree* m_outTree{nullptr};
};

template<typename InEventT, typename OutEventT>
TTreeLooper<InEventT, OutEventT>::TTreeLooper(TTree* inTree, TTree* outTree) :
  m_inTree(inTree), m_outTree(outTree)
{
  m_inEvent.Init(m_inTree);

  // m_outFile = new TFile(outFile.c_str(), "recreate");
  // m_outTree = new TTree(outTree.c_str(), outTree.c_str());

  // m_outTree->SetDirectory(outFile); // just to be sure that it does not get a memory resident TTree

  m_outEvent.Create(m_outTree);
}

template<typename InEventT, typename OutEventT>
template<typename CondF, PrintStyle PS>
void TTreeLooper<InEventT, OutEventT>::loop(CondF cond, const long int maxEvents)
{
  // check first how many events we want to process and correct for a possible input error, where more events
  // then present are requested
  const long int nInputEvents = m_inTree->GetEntries();
  const size_t nEvents = (maxEvents < 0 || maxEvents > nInputEvents) ? nInputEvents : maxEvents;

  size_t count{};

  std::cout << "Looping over " << nEvents << "\n";
  auto startTime = std::chrono::high_resolution_clock::now();
  for (size_t i = 0; i < nEvents; ++i) {
    if (!checkGetEntry(m_inTree, i)) continue;
    // if (!(i % 2)) continue; // skip every second event

    if (cond(m_inEvent, m_outEvent)) {
      m_outTree->Fill();
      count++;
    }
    printProgress<PS>(i, nEvents - 1, startTime); // -1 to reach 100 %
  }

  std::cout << "number of reconstructed events: " << count << " of a total of " << nEvents << " events\n";

  // m_outFile->Write();
  // m_outFile->Close();
}


#endif
