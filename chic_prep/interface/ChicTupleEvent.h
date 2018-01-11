#ifndef CHIBCHICPOLFW_CHICPREP_CHICTUPLEEVENT_H__
#define CHIBCHICPOLFW_CHICPREP_CHICTUPLEEVENT_H__

#include "ChicInputEvent.h"

#include "TTree.h"


class ChicTupleEvent {
public:
  ChicTupleEvent() = default;

  ChicTupleEvent(const ChicInputEvent &inEvent);

  void Create(TTree *t);

  double costh_HX{};
  double phi_HX{};
  double costh_PX{};
  double phi_PX{};
  double costh_CS{};
  double phi_CS{};

  double chicPt{};
  double chicMass{};
  double chicRap{};

  double Jpsict{};

  double wChic1{};
  double wChic2{};
};

void ChicTupleEvent::Create(TTree *t)
{
  t->Branch("costh_HX", &costh_HX);
  t->Branch("phi_HX", &phi_HX);
  t->Branch("costh_PX", &costh_PX);
  t->Branch("phi_PX", &phi_PX);
  t->Branch("costh_CS", &costh_CS);
  t->Branch("phi_CS", &phi_CS);

  t->Branch("chicPt", &chicPt);
  t->Branch("chicMass", &chicMass);
  t->Branch("chicRap", &chicRap);

  t->Branch("Jpsict", &Jpsict);

  t->Branch("wChic1", &wChic1);
  t->Branch("wChic2", &wChic2);
}

#endif
