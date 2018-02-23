#ifndef CHIBCHICPOLFW_CHICPREP_JPSIBASICTUPLING_H__
#define CHIBCHICPOLFW_CHICPREP_JPSIBASICTUPLING_H__


#include "ChicInputEvent.h"
#include "ChicTupleEvent.h"

#include "general/interface/calcAngles.h"

#include "TTree.h"

#include <array>
#include <vector>



// default branch names for 2012 Jpsi
struct DefaultJpsiBranches {
  static constexpr auto jpsi = "JpsiP";
  static constexpr auto chic = "There is no chic branch in J/psi data";
  static constexpr auto lepP = "muPosP";
  static constexpr auto lepN = "muNegP";
  static constexpr auto Jpsict = "Jpsict";
};

std::array<const char*, 9> triggerBranches = {
  "HLT_Dimuon8_Jpsi_v3",
  "HLT_Dimuon8_Jpsi_v4",
  "HLT_Dimuon8_Jpsi_v5",
  "HLT_Dimuon8_Jpsi_v6",
  "HLT_Dimuon8_Jpsi_v7",
  "HLT_Dimuon10_Jpsi_v3",
  "HLT_Dimuon10_Jpsi_v4",
  "HLT_Dimuon10_Jpsi_v5",
  "HLT_Dimuon10_Jpsi_v6"
};

struct BasicJpsiAddInfo {
  double JpsiPt;
  double JpsiRap;
  double JpsiMass;
  double JpsictErr;
  double vtxProb;

  double muN_pt;
  double muN_eta;
  double muP_pt;
  double muP_eta;

  int trigger;

  std::vector<int> trigger_branches;

  void Init(TTree *t)
  {
    t->SetBranchAddress("JpsictErr", &JpsictErr);
    t->SetBranchAddress("JpsiVprob", &vtxProb);

    trigger_branches.reserve(triggerBranches.size());
    for (size_t i = 0; i < triggerBranches.size(); ++i) {
      trigger_branches.push_back(-1);
      t->SetBranchAddress(triggerBranches[i], &trigger_branches[i]);
    }
  }

  void Create(TTree *t)
  {
    t->Branch("JpsiPt", &JpsiPt);
    t->Branch("JpsiRap", &JpsiRap);
    t->Branch("JpsiMass", &JpsiMass);
    t->Branch("JpsictErr", &JpsictErr);
    t->Branch("vtxProb", &vtxProb);

    t->Branch("muN_pt", &muN_pt);
    t->Branch("muN_eta", &muN_eta);
    t->Branch("muP_pt", &muP_pt);
    t->Branch("muP_eta", &muP_eta);
    t->Branch("trigger", &trigger);
  }
};

using JpsiBasicTuplingInEvent = ChicInputEvent<BasicJpsiAddInfo, DefaultJpsiBranches>;
using JpsiBasicTuplingOutEvent = ChicTupleEvent<BasicJpsiAddInfo>;

int triggerDecision(const std::vector<int>& triggers)
{
  for (const auto t : triggers) {
    if (t == 1) return 1;
  }
  return 0;
}

bool BasicJpsiTupling(const JpsiBasicTuplingInEvent& inEvent, JpsiBasicTuplingOutEvent& event)
{
  event.info().trigger = triggerDecision(inEvent.info().trigger_branches);
  event.info().vtxProb = inEvent.info().vtxProb;
  event.info().JpsictErr = inEvent.info().JpsictErr;
  event.Jpsict = inEvent.Jpsict;

  event.info().JpsiPt = inEvent.jpsi().Pt();
  event.info().JpsiRap = inEvent.jpsi().Rapidity();
  event.info().JpsiMass = inEvent.jpsi().M();

  event.info().muP_pt = inEvent.muP().Pt();
  event.info().muP_eta = inEvent.muP().Eta();
  event.info().muN_pt = inEvent.muN().Pt();
  event.info().muN_eta = inEvent.muN().Eta();

  const auto angles_HX  = calcAnglesInFrame(inEvent.muN(), inEvent.muP(), RefFrame::HX);
  event.costh_HX = angles_HX.costh;
  event.phi_HX = angles_HX.phi;
  event.cosalpha_HX = angles_HX.cosalpha;

  const auto angles_PX  = calcAnglesInFrame(inEvent.muN(), inEvent.muP(), RefFrame::PX);
  event.costh_PX = angles_PX.costh;
  event.phi_PX = angles_PX.phi;
  event.cosalpha_PX = angles_PX.cosalpha;

  const auto angles_CS  = calcAnglesInFrame(inEvent.muN(), inEvent.muP(), RefFrame::CS);
  event.costh_CS = angles_CS.costh;
  event.phi_CS = angles_CS.phi;
  event.cosalpha_CS = angles_CS.cosalpha;

  return true;
}


#endif
