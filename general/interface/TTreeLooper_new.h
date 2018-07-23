#ifndef CHIBCHICPOLFW_GENERAL_TTREELOOPER_H__new // TODO: remove new after done
#define CHIBCHICPOLFW_GENERAL_TTREELOOPER_H__new // TODO: remove new after done

#include "root_utils.h"
#include "progress.h"

#include "TTree.h"

#include <iostream>

template<typename EventT>
class TTreeLooper {
public:
  TTreeLooper() = delete;

  TTreeLooper(TTree* inTree, TTree *outTree, EventT& event);

  template<PrintStyle PS=PrintStyle::ProgressBar>
  void loop(EventT& event, const long int maxEvents=-1);

private:
  TTree *m_inTree{nullptr};
  TTree *m_outTree{nullptr};
};


template<typename EventT>
TTreeLooper<EventT>::TTreeLooper(TTree* inTree, TTree *outTree, EventT& event) :
  m_inTree(inTree), m_outTree(outTree)
{
  event.Init(m_inTree);
  event.Create(m_outTree);
}

template<typename EventT>
template<PrintStyle PS>
void TTreeLooper<EventT>::loop(EventT& event, const long int maxEvents)
{
  // check first how many events we want to process and correct for a possible input error, where more events
  // then present are requested
  const long int nInputEvents = m_inTree->GetEntries();
  const size_t nEvents = (maxEvents < 0 || maxEvents > nInputEvents) ? nInputEvents : maxEvents;

  size_t count{};


  auto startTime = ProgressClock::now();
  for (size_t i = 0; i < nEvents; ++i) {
    if (!checkGetEntry(m_inTree, i)) continue;
    event.Convert();
    if (event.Accept()) {
      m_outTree->Fill();
      count++;
    }

    printProgress<PS>(i, nEvents - 1, startTime);// -1 to reach 100 %
  }

  std::cout << "Converted " << count << " events from a total of " << nEvents << "\n";
}

#endif
