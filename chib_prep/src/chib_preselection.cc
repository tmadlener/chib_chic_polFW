#include "chib_preselection.h"

#include "calcAngles.h"
#include "misc_utils.h"

#include <string>
#include <cmath>

ChibPreselection::ChibPreselection(const std::vector<std::string>& infilenames, const std::string & outfilename, const std::string & intreename, const std::string & outtreename) :
  TreeProcessor(infilenames, outfilename, intreename, outtreename)
{
  std::string trigger_base = "HLT_Dimuon8_Upsilon_Barrel_v";
  for (int i = 1; i < 8; ++i)  triggers.push_back(trigger_base + std::to_string(i));

  AddBranchesToCopy({ "chi_p4","dimuon_p4", "photon_p4", "dz", "vProb" });
  AddBranchesNeededOnlyAsInput({ "muonP_p4","muonN_p4" });
  AddBranchesNeededOnlyAsInput(triggers);
}

bool ChibPreselection::fill_and_cut_variables()
{
  // Get variables needed for new branches or cuts
  static thread_local auto & chi = get_branch<TLorentzVector*>("chi_p4");
  static thread_local auto & dimuon = get_branch<TLorentzVector*>("dimuon_p4");
  static thread_local auto & muPos = get_branch<TLorentzVector*>("muonP_p4");
  static thread_local auto & muNeg = get_branch<TLorentzVector*>("muonN_p4");

  static thread_local bool get_trigger = true;
  static thread_local std::vector<std::reference_wrapper<const UInt_t> > trigs;
  if (get_trigger) {
    get_trigger = false;
    trigs.reserve(triggers.size());
    for (const auto &t : triggers) {
      //std::cout << t << std::endl;
      trigs.push_back(get_branch<UInt_t>(t));
    }
  }


  // CUTS - TODO: config file for cuts
  // TODO: what about rf1S, rf2S,...

  // Trigger - TODO: only matched?
  bool trig_passed = false;
  for (const auto &t : trigs) if (t > 0) trig_passed = true;
  if (!trig_passed) return false;

  // Dimuon Mass
  dimuon_mass = dimuon->M();
  if (dimuon_mass < 8.5 || dimuon_mass > 11.5) return false;

  // Pt
  dimuon_pt = dimuon->Pt();
  if (dimuon_pt < 8 || dimuon_pt > 999) return false;

  // Rapidity
  dimuon_rap = dimuon->Rapidity();
  if (dimuon_rap < -1.2 || dimuon_rap > 1.2) return false;

  // Cowboys/Seagulls
  delta_phi = reduceRange(muPos->Phi() - muNeg->Phi());

  // Dimuon Vertex Probability




  // NEW BRANCHES (some are already filled in the CUTS section)

  chi_mass = chi->M();
  chi_pt = chi->Pt();
  chi_rap = chi->Rapidity();
    

  eta_mupos = muPos->Eta();
  pt_mupos = muPos->Pt();
  phi_mupos = muPos->Phi();
  eta_muneg = muNeg->Eta();
  pt_muneg = muNeg->Pt();
  phi_muneg = muNeg->Phi();

  //TODO: calcAnglesInFrame has a beam energy of 4000 GeV hardcoded, is this a problem
  const auto angles_hx = calcAnglesInFrame(*muNeg, *muPos, RefFrame::HX);
  const auto angles_px = calcAnglesInFrame(*muNeg, *muPos, RefFrame::PX);
  const auto angles_cs = calcAnglesInFrame(*muNeg, *muPos, RefFrame::CS);

  costh_hx = angles_hx.costh;
  phi_hx = angles_hx.phi;
  cosalpha_hx = angles_hx.cosalpha;

  costh_px = angles_px.costh;
  phi_px = angles_px.phi;
  cosalpha_px = angles_px.cosalpha;

  costh_cs = angles_cs.costh;
  phi_cs = angles_cs.phi;
  cosalpha_cs = angles_cs.cosalpha;

  return true;
}

void ChibPreselection::setup_new_branches()
{

  m_out_tree->Branch("dimuon_mass", &dimuon_mass);
  m_out_tree->Branch("dimuon_pt", &dimuon_pt);
  m_out_tree->Branch("dimuon_rap", &dimuon_rap);

  m_out_tree->Branch("chi_mass", &chi_mass);
  m_out_tree->Branch("chi_pt", &chi_pt);
  m_out_tree->Branch("chi_rap", &chi_rap);


  m_out_tree->Branch("cosTh_HX", &costh_hx);
  m_out_tree->Branch("phi_HX", &phi_hx);
  m_out_tree->Branch("cosAlpha_HX", &cosalpha_hx);

  m_out_tree->Branch("cosTh_PX", &costh_px);
  m_out_tree->Branch("phi_PX", &phi_px);
  m_out_tree->Branch("cosAlpha_PX", &cosalpha_px);

  m_out_tree->Branch("cosTh_CS", &costh_cs);
  m_out_tree->Branch("phi_CS", &phi_cs);
  m_out_tree->Branch("cosAlpha_CS", &cosalpha_cs);

  m_out_tree->Branch("delta_phi", &delta_phi);

  m_out_tree->Branch("eta_mupos", &eta_mupos);
  m_out_tree->Branch("pt_mupos", &pt_mupos);
  m_out_tree->Branch("phi_mupos", &phi_mupos);
  m_out_tree->Branch("eta_muneg", &eta_muneg);
  m_out_tree->Branch("pt_muneg", &pt_muneg);
  m_out_tree->Branch("phi_muneg", &phi_muneg);


}


int main() {
  ChibPreselection preselection({ "/afs/hephy.at/work/j/jnecker/data/full/chib_2016.root" }, "first_chib_preselection.root", "rootuple/chiTree", "data");
  preselection.process(100000, 8);
}