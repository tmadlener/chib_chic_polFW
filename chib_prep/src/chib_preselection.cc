#include "chib_preselection.h"

#include "calcAngles.h"
#include "misc_utils.h"

#include <string>
#include <cmath>

ChibPreselection::ChibPreselection(const std::vector<std::string>& infilenames, const std::string & outfilename, const std::string & intreename, const std::string & outtreename) :
  TreeProcessor(infilenames, outfilename, intreename, outtreename)
{
  for (const auto &tb : {
    "HLT_Dimuon8_Upsilon_Barrel_v"
  }) {
    for (int i = 1; i < 8; ++i)  triggers.push_back(tb + std::to_string(i));
  }

  var_collections = {
    {"chi_p4","dimuon_p4", "photon_p4", "muonP_p4","muonN_p4"},
    { "rf1S_chi_p4", /*"rf1S_dimuon_p4", "rf1S_muonP_p4", "rf1S_muonN_p4" */},
    { "rf2S_chi_p4", /*"rf2S_dimuon_p4", "rf2S_muonP_p4", "rf2S_muonN_p4"*/ },
    { "rf3S_chi_p4", /*"rf3S_dimuon_p4", "rf3S_muonP_p4", "rf3S_muonN_p4"*/ }
  };

  // Add branches to read
  for (const auto &c : var_collections) AddBranchesNeededOnlyAsInput(c);
  AddBranchesNeededOnlyAsInput(triggers);

  // Add Branches to copy
  AddBranchesToCopy({ "dz", "vProb" , "q_value", "conversionflag", "probFit1S", "rf1S_rank", "probFit2S", "rf2S_rank", "probFit3S", "rf3S_rank", });


}

bool ChibPreselection::fill_and_cut_variables()
{

  // Get variables needed for new branches or cuts
  static thread_local auto & chi = get_branch<TLorentzVector*>("chi_p4");
  static thread_local auto & dimuon = get_branch<TLorentzVector*>("dimuon_p4");
  static thread_local auto & muPos = get_branch<TLorentzVector*>("muonP_p4");
  static thread_local auto & muNeg = get_branch<TLorentzVector*>("muonN_p4");

  static thread_local auto & chi_rf1s = get_branch<TLorentzVector*>("rf1S_chi_p4");
  //static thread_local auto & dimuon_rf1s = get_branch<TLorentzVector*>("rf1S_dimuon_p4");
  //static thread_local auto & muPos_rf1s = get_branch<TLorentzVector*>("rf1S_muonP_p4");
  //static thread_local auto & muNeg_rf1s = get_branch<TLorentzVector*>("rf1S_muonN_p4");

  static thread_local auto & chi_rf2s = get_branch<TLorentzVector*>("rf2S_chi_p4");
  //static thread_local auto & dimuon_rf2s = get_branch<TLorentzVector*>("rf2S_dimuon_p4");
  //static thread_local auto & muPos_rf2s = get_branch<TLorentzVector*>("rf2S_muonP_p4");
  //static thread_local auto & muNeg_rf2s = get_branch<TLorentzVector*>("rf2S_muonN_p4");

  static thread_local auto & chi_rf3s = get_branch<TLorentzVector*>("rf3S_chi_p4");
  //static thread_local auto & dimuon_rf3s = get_branch<TLorentzVector*>("rf3S_dimuon_p4");
  //static thread_local auto & muPos_rf3s = get_branch<TLorentzVector*>("rf3S_muonP_p4");
  //static thread_local auto & muNeg_rf3s = get_branch<TLorentzVector*>("rf3S_muonN_p4");

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


  // CUTS - TODO: LOGGING
  
  //
  // Dimuon Cuts
  //

  // Trigger Cut
  bool trig_passed = false;
  for (const auto &t : trigs) if (t > 1 /*>1: only matched*/) { trig_passed = true;  break; }
  if (!trig_passed) return false;

  // Single Muon Cut
  if (!accept_muon(muPos) || !accept_muon(muNeg)) return false;

  // Dimuon Mass
  if (dimuon->M() < 8.5 || dimuon->M() > 11.5) return false;

  // Pt
  if (dimuon->Pt() < 8 || dimuon->Pt() > 999) return false;

  // Rapidity
  if (dimuon->Rapidity() < -1.2 || dimuon->Rapidity() > 1.2) return false;

  // Cowboys/Seagulls

  // Dimuon Vertex Probability  (vProb)

  // Dz

  //
  // Chi cuts?
  //


  // NEW BRANCHES

  auto dimuon_id = make_lorentz_flat(out_vars, chi, dimuon);
  make_mu_things(out_vars, muPos, muNeg, dimuon_id);

  // Only chi
  make_lorentz_flat(out_vars_rf1S, chi_rf1s);
  make_lorentz_flat(out_vars_rf2S, chi_rf2s);
  make_lorentz_flat(out_vars_rf3S, chi_rf3s);

  return true;
}

void ChibPreselection::setup_new_branches()
{
  setup_collections("", out_vars);
  setup_collections("_rf1S", out_vars_rf1S, true);
  setup_collections("_rf2S", out_vars_rf2S, true);
  setup_collections("_rf3S", out_vars_rf3S, true);
}

int ChibPreselection::make_lorentz_flat(std::vector<Double_t>& vars, TLorentzVector* chi, TLorentzVector * dimuon)
{
  int id = -1;

  vars.at(++id) = chi->M();
  vars.at(++id) = chi->Pt();
  vars.at(++id) = chi->Rapidity();

  if (dimuon) {
    vars.at(++id) = dimuon->M();
    vars.at(++id) = dimuon->Pt();
    vars.at(++id) = dimuon->Rapidity();
  }

  return id;
}

void ChibPreselection::make_mu_things(std::vector<Double_t>& vars, TLorentzVector * muPos, TLorentzVector * muNeg,int start)
{
  int id = start;

  //TODO: calcAnglesInFrame has a beam energy of 4000 GeV hardcoded, is this a problem?
  const auto angles_hx = calcAnglesInFrame(*muNeg, *muPos, RefFrame::HX);
  const auto angles_px = calcAnglesInFrame(*muNeg, *muPos, RefFrame::PX);
  const auto angles_cs = calcAnglesInFrame(*muNeg, *muPos, RefFrame::CS);

  vars.at(++id) = angles_hx.costh;
  vars.at(++id) = angles_hx.phi;
  vars.at(++id) = angles_hx.cosalpha;

  vars.at(++id) = angles_px.costh;
  vars.at(++id) = angles_px.phi;
  vars.at(++id) = angles_px.cosalpha;

  vars.at(++id) = angles_cs.costh;
  vars.at(++id) = angles_cs.phi;
  vars.at(++id) = angles_cs.cosalpha;

  vars.at(++id) = reduceRange(muPos->Phi() - muNeg->Phi());

  vars.at(++id) = muPos->Eta();
  vars.at(++id) = muPos->Pt();
  vars.at(++id) = muPos->Phi();

  vars.at(++id) = muNeg->Eta();
  vars.at(++id) = muNeg->Pt();
  vars.at(++id) = muNeg->Phi();

}

void ChibPreselection::setup_collections(const std::string & varsuffix, std::vector<Double_t>& vars, bool only_chi)
{
  // First make sure that there are enough elements in the vector
  for (size_t i = 0; i < 32; ++i) vars.push_back(0);

  int id = -1;

  //The ids defined here has to be used in make_mu_things and make_lorentz_flat functions
  m_out_tree->Branch(("chi_mass" + varsuffix).c_str(), &vars.at(++id));
  m_out_tree->Branch(("chi_pt" + varsuffix).c_str(), &vars.at(++id));
  m_out_tree->Branch(("chi_rap" + varsuffix).c_str(), &vars.at(++id));

  if (!only_chi) {

    m_out_tree->Branch(("dimuon_mass" + varsuffix).c_str(), &vars.at(++id));
    m_out_tree->Branch(("dimuon_pt" + varsuffix).c_str(), &vars.at(++id));
    m_out_tree->Branch(("dimuon_rap" + varsuffix).c_str(), &vars.at(++id));

    // id == 5

    m_out_tree->Branch(("cosTh_HX" + varsuffix).c_str(), &vars.at(++id));
    m_out_tree->Branch(("phi_HX" + varsuffix).c_str(), &vars.at(++id));
    m_out_tree->Branch(("cosAlpha_HX" + varsuffix).c_str(), &vars.at(++id));

    m_out_tree->Branch(("cosTh_PX" + varsuffix).c_str(), &vars.at(++id));
    m_out_tree->Branch(("phi_PX" + varsuffix).c_str(), &vars.at(++id));
    m_out_tree->Branch(("cosAlpha_PX" + varsuffix).c_str(), &vars.at(++id));

    m_out_tree->Branch(("cosTh_CS" + varsuffix).c_str(), &vars.at(++id));
    m_out_tree->Branch(("phi_CS" + varsuffix).c_str(), &vars.at(++id));
    m_out_tree->Branch(("cosAlpha_CS" + varsuffix).c_str(), &vars.at(++id));
    // id == 14
    m_out_tree->Branch(("delta_phi" + varsuffix).c_str(), &vars.at(++id));

    m_out_tree->Branch(("eta_mupos" + varsuffix).c_str(), &vars.at(++id));
    m_out_tree->Branch(("pt_mupos" + varsuffix).c_str(), &vars.at(++id));
    m_out_tree->Branch(("phi_mupos" + varsuffix).c_str(), &vars.at(++id));

    m_out_tree->Branch(("eta_muneg" + varsuffix).c_str(), &vars.at(++id));
    m_out_tree->Branch(("pt_muneg" + varsuffix).c_str(), &vars.at(++id));
    m_out_tree->Branch(("phi_muneg" + varsuffix).c_str(), &vars.at(++id));
    // id == 21
  }
}

bool ChibPreselection::accept_muon(const TLorentzVector * mu)
{
  // TODO: rapidity dependent muon cut

  if (mu->Pt() < 3.) return false;
  return true;
}


int main() {
  ChibPreselection preselection({ "/afs/hephy.at/work/j/jnecker/data/full/chib_2016.root" }, "preselected_data.root", "rootuple/chiTree", "data");
  preselection.process(-1, 8);
}