#ifndef CHIBCHICPOLFW_CHICPREP_CHICTUPLINGMC_H__
#define CHIBCHICPOLFW_CHICPREP_CHICTUPLINGMC_H__

#include "general/interface/calcAngles.h"

#include "ChicInputEvent.h"
#include "ChicTupleEvent.h"
#include "ChicTuplingInfoMC.h"

// convenience typedefs
using ChicMCInputEvent = ChicInputEvent<MCAddInfoIn, MCDefaultNames>;
using ChicMCTupleEvent = ChicTupleEvent<MCAddInfoOut>;

bool chicMCTupling(const ChicMCInputEvent& inEvent, ChicMCTupleEvent& event)
{
  // weed out obviously mis-reconstructed events
  if (inEvent.chic().Pt() == 0) return false;

  event.Jpsict = inEvent.Jpsict;
  event.chicPt = inEvent.chic().Pt();
  event.chicMass = inEvent.chic().M();
  event.chicRap = inEvent.chic().Rapidity();

  event.wChic2 = inEvent.info().gen_chic->M() > 3.53;
  event.wChic1 = !event.wChic2;

  const auto anglesHX = calcAnglesInFrame(inEvent.muN(), inEvent.muP(), RefFrame::HX);
  event.costh_HX = anglesHX.costh;
  event.phi_HX = anglesHX.phi;
  event.cosalpha_HX = anglesHX.cosalpha;

  const auto anglesCS = calcAnglesInFrame(inEvent.muN(), inEvent.muP(), RefFrame::CS);
  event.costh_CS = anglesCS.costh;
  event.phi_CS = anglesCS.phi;
  event.cosalpha_CS = anglesCS.cosalpha;

  const auto anglesPX = calcAnglesInFrame(inEvent.muN(), inEvent.muP(), RefFrame::PX);
  event.costh_PX = anglesPX.costh;
  event.phi_PX = anglesPX.phi;
  event.cosalpha_PX = anglesPX.cosalpha;

  const auto genAnglesHX = calcAnglesInFrame(inEvent.muN(), inEvent.muP(), RefFrame::HX);
  event.info().costh_HX = genAnglesHX.costh;
  event.info().phi_HX = genAnglesHX.phi;
  event.info().cosalpha_HX = genAnglesHX.cosalpha;

  const auto genAnglesCS = calcAnglesInFrame(inEvent.muN(), inEvent.muP(), RefFrame::CS);
  event.info().costh_CS = genAnglesCS.costh;
  event.info().phi_CS = genAnglesCS.phi;
  event.info().cosalpha_CS = genAnglesCS.cosalpha;

  const auto genAnglesPX = calcAnglesInFrame(inEvent.muN(), inEvent.muP(), RefFrame::PX);
  event.info().costh_PX = genAnglesPX.costh;
  event.info().phi_PX = genAnglesPX.phi;
  event.info().cosalpha_PX = genAnglesPX.cosalpha;

  event.info().trigger = inEvent.info().trigger;
  event.info().vtxProb = inEvent.info().vtxProb;

  event.info().muP_pt = inEvent.muP().Pt();
  event.info().muN_pt = inEvent.muN().Pt();
  event.info().muP_eta = inEvent.muP().Eta();
  event.info().muN_eta = inEvent.muN().Eta();

  event.info().jpsiPt = inEvent.jpsi().Pt();
  event.info().jpsiRap = inEvent.jpsi().Rapidity();

  return true; // simply accept all events
}

#endif
