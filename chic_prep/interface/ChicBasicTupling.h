#ifndef CHIBCHICPOLFW_CHICPREP_CHICBASICTUPLING_H__
#define CHIBCHICPOLFW_CHICPREP_CHICBASICTUPLING_H__

#include "ChicInputEvent.h"
#include "ChicTupleEvent.h"

#include "general/interface/calcAngles.h"

struct TTree;
struct TLorentzVector;

struct RawDefaultNames {
  static constexpr auto jpsi = "dimuon_p4";
  static constexpr auto chic = "rf1S_chi_p4";
  static constexpr auto lepP = "muonP_p4";
  static constexpr auto lepN = "muonN_p4";
  static constexpr auto Jpsict = "ctpv";
};

struct BasicAddInfo {
  double JpsiPt;
  double JpsiRap;
  double JpsiMass;
  double JpsictErr;
  double vtxProb;

  double muP_pt;
  double muP_eta;
  double muN_pt;
  double muN_eta;

  double convRadius;
  double gammaEta;
  double gammaDz;

  double photonPt;
  double photonEta;

  int trigger;

  TLorentzVector* gamma{nullptr};

  void Create(TTree* t)
  {
    t->Branch("JpsiPt", &JpsiPt);
    t->Branch("JpsiRap", &JpsiRap);
    t->Branch("JpsiMass", &JpsiMass);
    t->Branch("JpsictErr", &JpsictErr);
    t->Branch("vtxProb", &vtxProb);

    t->Branch("muP_pt", &muP_pt);
    t->Branch("muP_eta", &muP_eta);
    t->Branch("muN_pt", &muN_pt);
    t->Branch("muN_eta", &muN_eta);

    t->Branch("convRadius", &convRadius);
    t->Branch("gammaEta", &gammaEta);
    t->Branch("gammaDz", &gammaDz);

    t->Branch("photonPt", &photonPt);
    t->Branch("photonEta", &photonEta);

    t->Branch("trigger", &trigger);
  }

  void Init(TTree* t)
  {
    t->SetBranchAddress("probFit1S", &vtxProb);
    t->SetBranchAddress("ctpv_error", &JpsictErr);
    t->SetBranchAddress("dz", &gammaDz);
    t->SetBranchAddress("conv_vertex", &convRadius);
    t->SetBranchAddress("rf1S_photon_p4", &gamma);
    t->SetBranchAddress("trigger", &trigger);
  }
};

using ChicBasicTuplingOutEvent = ChicTupleEvent<BasicAddInfo>;
using ChicBasicTuplingInEvent = ChicInputEvent<BasicAddInfo, RawDefaultNames>;

bool chicBasicTupling(const ChicBasicTuplingInEvent& inEvent, ChicBasicTuplingOutEvent& event)
{
  event.Jpsict = inEvent.Jpsict;
  event.chicPt = inEvent.chic().Pt();
  event.chicMass = inEvent.chic().M();
  event.chicRap = inEvent.chic().Rapidity();

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

  event.info().muP_pt = inEvent.muP().Pt();
  event.info().muP_eta = inEvent.muP().Eta();
  event.info().muN_pt = inEvent.muN().Pt();
  event.info().muN_eta = inEvent.muN().Eta();

  event.info().JpsiPt = inEvent.jpsi().Pt();
  event.info().JpsiRap = inEvent.jpsi().Rapidity();
  event.info().JpsiMass = inEvent.jpsi().M();
  event.info().JpsictErr = inEvent.info().JpsictErr;

  event.info().vtxProb = inEvent.info().vtxProb;

  event.info().gammaEta = inEvent.info().gamma->Eta();
  event.info().convRadius = inEvent.info().convRadius;
  event.info().gammaDz = inEvent.info().gammaDz;

  event.info().photonPt = inEvent.info().gamma->Pt();
  event.info().photonEta = inEvent.info().gamma->Eta();

  event.info().trigger = inEvent.info().trigger;

  return true;
}


#endif
