#ifndef CHIBCHICPOLFW_CHICPREP_CHICTUPLINGWITHWEIGHTS_H__
#define CHIBCHICPOLFW_CHICPREP_CHICTUPLINGWITHWEIGHTS_H__

#include "general/interface/misc_utils.h"

#include "ChicInputEvent.h"
#include "ChicTupleEvent.h"
#include "JpsiMassSelector.h"

#include "general/interface/calcAngles.h"
#include "general/interface/misc_utils.h"

struct MassRegions;
struct LifeTimeRegions;
struct BkgWeights;
struct BkgWeights2D;

inline bool containedInMass(const double mass, const MassRegions &mr)
{
  // std::cout << mr.SR1 << " " << mr.SR2 << " " << mr.LSB << " " << mr.RSB << " " << mass << "\n";
  return mr.SR1.contains(mass) || mr.SR2.contains(mass) || mr.LSB.contains(mass) || mr.RSB.contains(mass);
}

inline bool containedInLT(const double ctau, const LifeTimeRegions& ltr)
{
  // std::cout << ltr.NP << " " << ltr.PR << " " << ctau << "\n";
  return ltr.NP.contains(ctau) || ltr.PR.contains(ctau);
}


struct DiMuonSelector {
  virtual bool operator()(const TLorentzVector& muP, const TLorentzVector& muN) const = 0;
};

struct AllSelector : public DiMuonSelector {
  virtual bool operator()(const TLorentzVector&, const TLorentzVector&) const override
  {
    return true;
  }
};

struct NegativeSign {
  bool operator()(const double v) const { return v < 0; }
};

struct PositiveSign {
  bool operator()(const double v) const { return v > 0; }
};

template<typename SignSelector>
struct RelativePhiSelector {
  bool operator()(const TLorentzVector& muP, const TLorentzVector& muN) const
  {
    const double deltaPhi = reduceRange(muN.Phi() - muP.Phi());
    return selector(deltaPhi);
  }
  const SignSelector selector;
};

struct CowboySelector : public DiMuonSelector {
  virtual bool operator()(const TLorentzVector& muP, const TLorentzVector& muN) const override
  {
    return RelativePhiSelector<PositiveSign>()(muP, muN);
  }
};

struct SeagullSelector : public DiMuonSelector {
  virtual bool operator()(const TLorentzVector& muP, const TLorentzVector& muN) const override
  {
    return RelativePhiSelector<NegativeSign>()(muP, muN);
  }
};

bool chicTuplingWithWeights(const ChicInputEvent &inEvent, ChicTupleEvent &event,
                            const MassRegions &mr, const LifeTimeRegions &ltr,
                            const BkgWeights &weights)
{
  event.chicMass = inEvent.chic().M();
  event.Jpsict = inEvent.Jpsict;
  if (!(containedInMass(event.chicMass, mr) && containedInLT(event.Jpsict, ltr))) return false;

  event.chicPt = inEvent.chic().Pt();
  event.chicRap = inEvent.chic().Rapidity();

  const auto anglesHX = calcAnglesInFrame(inEvent.muN(), inEvent.muP(), RefFrame::HX);
  event.costh_HX = anglesHX.costh;
  event.phi_HX = anglesHX.phi;
  event.cosalpha_HX = anglesHX.cosalpha;

  const auto anglesPX = calcAnglesInFrame(inEvent.muN(), inEvent.muP(), RefFrame::PX);
  event.costh_PX = anglesPX.costh;
  event.phi_PX = anglesPX.phi;
  event.cosalpha_PX = anglesPX.cosalpha;

  const auto anglesCS = calcAnglesInFrame(inEvent.muN(), inEvent.muP(), RefFrame::CS);
  event.costh_CS = anglesCS.costh;
  event.phi_CS = anglesCS.phi;
  event.cosalpha_CS = anglesCS.cosalpha;

  // std::cout << "====================\n";

  double wChic1 = !mr.SR2.contains(event.chicMass);
  double wChic2 = !mr.SR1.contains(event.chicMass);

  // std::cout << wChic1 << " " << wChic2 << "\n";

  if (mr.LSB.contains(event.chicMass) || mr.RSB.contains(event.chicMass)) {
    wChic1 *= weights.wMBkg1;
    wChic2 *= weights.wMBkg2;
  }
  // std::cout << wChic1 << " " << wChic2 << "\n";

  if (ltr.NP.contains(event.Jpsict)) {
    wChic1 *= weights.wNP;
    wChic2 *= weights.wNP;
  }
  // std::cout << wChic1 << " " << wChic2 << "\n";

  if (!(mr.SR1.contains(event.chicMass) && ltr.PR.contains(event.Jpsict))) wChic1 *= -1;
  if (!(mr.SR2.contains(event.chicMass) && ltr.PR.contains(event.Jpsict))) wChic2 *= -1;
  // std::cout << wChic1 << " " << wChic2 << "\n";


  event.wChic1 = wChic1;
  event.wChic2 = wChic2;

  // std::cout << "====================\n";


  return true;
}


bool chicTuplingWithMassWeights(const ChicInputEvent &inEvent, ChicTupleEvent &event,
                                const MassRegions& mr, const Region<Boundary::TwoSided> &ltr,
                                const std::pair<double, double> weights,
                                const double ptMin, const double ptMax,
                                const JpsiMassSelector& jpsiMass)
{

  event.chicPt = inEvent.chic().Pt();
  event.chicRap = inEvent.chic().Rapidity();

  // NOTE: tmadlener, 04.08.17: hardcoded 1 pt bin into here to have some results
  if (event.chicPt < ptMin || event.chicPt > ptMax) return false;
  // if (event.chicRap > 1.2) return false; // NOTE: hardcoded

  if (!jpsiMass.contains(inEvent.jpsi().M(), inEvent.jpsi().Rapidity())) return false;

  event.chicMass = inEvent.chic().M();
  event.Jpsict = inEvent.Jpsict;

  // check if event falls into analysis range
  if (!containedInMass(event.chicMass, mr) || !ltr.contains(event.Jpsict)) return false;

  const auto anglesHX = calcAnglesInFrame(inEvent.muN(), inEvent.muP(), RefFrame::HX);
  event.costh_HX = anglesHX.costh;
  event.phi_HX = anglesHX.phi;
  event.cosalpha_HX = anglesHX.cosalpha;

  const auto anglesPX = calcAnglesInFrame(inEvent.muN(), inEvent.muP(), RefFrame::PX);
  event.costh_PX = anglesPX.costh;
  event.phi_PX = anglesPX.phi;
  event.cosalpha_PX = anglesPX.cosalpha;

  const auto anglesCS = calcAnglesInFrame(inEvent.muN(), inEvent.muP(), RefFrame::CS);
  event.costh_CS = anglesCS.costh;
  event.phi_CS = anglesCS.phi;
  event.cosalpha_CS = anglesCS.cosalpha;

  double wChic1 = mr.SR1.contains(event.chicMass) ? 1 : weights.first;
  if (mr.SR2.contains(event.chicMass)) wChic1 = 0;

  double wChic2 = mr.SR2.contains(event.chicMass) ? 1 : weights.second;
  if (mr.SR1.contains(event.chicMass)) wChic2 = 0;

  event.wChic1 = wChic1;
  event.wChic2 = wChic2;

  return true;
}

template<typename DiMuSelector>
bool chicTuplingWith2DWeights(const ChicInputEvent& inEvent, ChicTupleEvent& event,
                              const MassRegions& mr, const LifeTimeRegions& ltr,
                              const BkgWeights2D& weights,
                              const double ptMin, const double ptMax, const double absRapMax,
                              const JpsiMassSelector& jpsiMass, const DiMuSelector& dimuSelector)
{
  // check kinematic region
  event.chicPt = inEvent.chic().Pt();
  event.chicRap = inEvent.chic().Rapidity();
  if (event.chicPt < ptMin || event.chicPt > ptMax) return false;
  if (std::abs(event.chicRap) > absRapMax) return false;

  if (!jpsiMass.contains(inEvent.jpsi().M(), inEvent.jpsi().Rapidity())) return false;

  if (!dimuSelector(inEvent.muP(), inEvent.muN())) return false;

  // check if event is in analysis range
  event.chicMass = inEvent.chic().M();
  event.Jpsict = inEvent.Jpsict;
  if (!containedInMass(event.chicMass, mr) || !containedInLT(event.Jpsict, ltr)) return false;

  const auto anglesHX = calcAnglesInFrame(inEvent.muN(), inEvent.muP(), RefFrame::HX);
  event.costh_HX = anglesHX.costh;
  event.phi_HX = anglesHX.phi;
  event.cosalpha_HX = anglesHX.cosalpha;

  const auto anglesPX = calcAnglesInFrame(inEvent.muN(), inEvent.muP(), RefFrame::PX);
  event.costh_PX = anglesPX.costh;
  event.phi_PX = anglesPX.phi;
  event.cosalpha_PX = anglesPX.cosalpha;

  const auto anglesCS = calcAnglesInFrame(inEvent.muN(), inEvent.muP(), RefFrame::CS);
  event.costh_CS = anglesCS.costh;
  event.phi_CS = anglesCS.phi;
  event.cosalpha_CS = anglesCS.cosalpha;

  const bool massSR1 = mr.SR1.contains(event.chicMass);
  const bool massSR2 = mr.SR2.contains(event.chicMass);
  const bool massSB = mr.LSB.contains(event.chicMass) || mr.RSB.contains(event.chicMass);
  const bool PR = ltr.PR.contains(event.Jpsict);
  const bool NP = ltr.NP.contains(event.Jpsict);

  const bool SR1 = massSR1 && PR;
  const bool SR2 = massSR2 && PR;
  const bool bkg1 = massSB || (massSR1 && NP);
  const bool bkg2 = massSB || (massSR2 && NP);

  event.wChic1 = SR1 + bkg1 * weights.wChic1;
  event.wChic2 = SR2 + bkg2 * weights.wChic2;

  return true;
}

#endif
