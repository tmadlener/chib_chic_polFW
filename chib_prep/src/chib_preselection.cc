#include "chib_preselection.h"

#include "utils.h"

#include "calcAngles.h"
#include "misc_utils.h"
#include "ArgParser.h"

#include "TFile.h"
#include "TUUID.h"

#include <string>
#include <cmath>

Long64_t ChibPreselection::entry_idx = 0;
std::mutex ChibPreselection::entry_mtx;

ChibPreselection::ChibPreselection(const std::vector<std::string>& infilenames, const std::string & outfilename, const std::string & intreename, const std::string & outtreename, int year) :
  TreeProcessor(infilenames, outfilename, intreename, outtreename),
  m_year(year)
{
  if (m_year == 2016) {
    for (const auto &tb : {
      "HLT_Dimuon8_Upsilon_Barrel_v"
    }) {
      for (int i = 1; i < 8; ++i)  triggers.push_back(tb + std::to_string(i));
    }

    var_collections = {
      {"chi_p4","dimuon_p4", "photon_p4", "muonP_p4","muonN_p4"},
      { "rf1S_chi_p4", "rf1S_photon_p4"/*"rf1S_dimuon_p4", "rf1S_muonP_p4", "rf1S_muonN_p4" */},
      { "rf2S_chi_p4", "rf2S_photon_p4"/*"rf2S_dimuon_p4", "rf2S_muonP_p4", "rf2S_muonN_p4"*/ },
      { "rf3S_chi_p4", "rf3S_photon_p4"/*"rf3S_dimuon_p4", "rf3S_muonP_p4", "rf3S_muonN_p4"*/ }
    };
    
    // Add trigger to branches to read
    for (const auto &c : var_collections) AddBranchesNeededOnlyAsInput(c);
    AddBranchesNeededOnlyAsInput(triggers);

    // Add Branches to copy
    AddBranchesToCopy({ "dz", "vProb" , "q_value", "conversionflag", "probFit1S", "rf1S_rank", "probFit2S", "rf2S_rank", "probFit3S", "rf3S_rank", "pi0rejected", "pi0_abs_mass", "numPrimaryVertices" });

    min_pt = 8;
    m_trigminval = 0; // bei rerco hat das mit matching nicht funktioniert, 1 only matched else 0
    muonNname = "muonN_p4";

  }

  if (m_year == 2017) {
    triggers.emplace_back("dimuon10ups_trigger");
    triggers.emplace_back("dimuon12ups_trigger");

    AddBranchesToCopy(triggers);

    AddBranchesNeededOnlyAsInput({ "chi_p4","dimuon_p4", "photon_p4", "muonP_p4", "muonM_p4", "rf1S_chi_p4", "rf2S_chi_p4", "rf3S_chi_p4" }); 
    
    AddBranchesToCopy({ "dz", "photon_flags", "probFit1S", "rf1S_rank", "probFit2S", "probFit3S", "numPrimaryVertices" });
    
    min_pt = 10;
    m_trigminval = 0;
    muonNname = "muonM_p4";
  }


}

bool ChibPreselection::fill_and_cut_variables()
{

  // Get variables needed for new branches or cuts
  static thread_local auto & chi = get_branch<TLorentzVector*>("chi_p4");
  static thread_local auto & dimuon = get_branch<TLorentzVector*>("dimuon_p4");
  static thread_local auto & muPos = get_branch<TLorentzVector*>("muonP_p4");
  static thread_local auto & muNeg = get_branch<TLorentzVector*>(muonNname.c_str());

  static thread_local auto & photon = get_branch<TLorentzVector*>("photon_p4");

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
  for (const auto &t : trigs) if (t > m_trigminval /*>1: 2016 this means only matched*/) { trig_passed = true;  break; }
  if (!trig_passed) return false;

   // Single Muon Cut
  if (!accept_muon(muPos) || !accept_muon(muNeg)) return false;

  // Dimuon Mass
  if (dimuon->M() < 8.5 || dimuon->M() > 11.5) return false;

  // Pt
  if (dimuon->Pt() < min_pt || dimuon->Pt() > 999) return false;

  // Rapidity
  if (dimuon->Rapidity() < -1.2 || dimuon->Rapidity() > 1.2) return false;

  // Cowboys/Seagulls

  // Dimuon Vertex Probability  (vProb)
  if (m_year == 2016) {
    static thread_local auto & vProb = get_branch<double>("vProb");
    if (vProb < 0.01) return false;
  }

  // Photon Eta, Rap

  // Dz

  //
  // Chi cuts?
  //


  // NEW BRANCHES

  auto dimuon_id = make_lorentz_flat(out_vars, chi, dimuon);
  make_mu_things(out_vars, muPos, muNeg, dimuon_id);

  if (m_year == 2016) {
    static thread_local auto & rf1SPhoton = get_branch<TLorentzVector*>("rf1S_photon_p4");
    static thread_local auto & rf2SPhoton = get_branch<TLorentzVector*>("rf2S_photon_p4");
    static thread_local auto & rf3SPhoton = get_branch<TLorentzVector*>("rf2S_photon_p4");
    fill_photon_vars({ photon, rf1SPhoton, rf2SPhoton, rf3SPhoton });
  }

  if (m_year == 2017) fill_photon_vars({ photon});

  // Only chi
  make_lorentz_flat(out_vars_rf1S, chi_rf1s);
  make_lorentz_flat(out_vars_rf2S, chi_rf2s);
  make_lorentz_flat(out_vars_rf3S, chi_rf3s);
  {
    std::lock_guard<std::mutex> lock(entry_mtx);
    // This copying to a local variable is important, 
    // because when TTree::Fill is called the global ID probably already has another value,
    // because it is no longer locked at that point.
    local_entry_idx = ++entry_idx;
  }
  // RooDataset cannot store a Long_64t, so split it up in two Int_t:
  entry_idx_high = (local_entry_idx >> 32);
  entry_idx_low = local_entry_idx;
  //-----------------------------------------------------
  // REBUILD EntryID:
  //-----------------------------------------------------
  // Long64_t lowmask = 0xFFFFFFFF;
  // Long64_t high_long = EntryID_high;
  // high_long <<= 32;
  // Long64_t EntryID_rebuild = EntryID_low & lowmask;
  // EntryID_rebuild |= high_long;
  //-----------------------------------------------------

  return true;
}

void ChibPreselection::setup_new_branches()
{
  // Add id for each entry
  m_out_tree->Branch("EntryID", &local_entry_idx);
  m_out_tree->Branch("EntryID_high", &entry_idx_high);
  m_out_tree->Branch("EntryID_low", &entry_idx_low);

  // Add file uid
  get_outfile()->cd();
  // NB: running the TreeLooper parallel, every worker will have another uuid, but at the end the one of the first worker is used.
  TNamed id_s("DataID", TUUID().AsString());
  id_s.Write(0, TObject::kWriteDelete);

  setup_collections("", out_vars);
  setup_collections("_rf1S", out_vars_rf1S, true);
  setup_collections("_rf2S", out_vars_rf2S, true);
  setup_collections("_rf3S", out_vars_rf3S, true);

  setup_photon_vars();

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

void ChibPreselection::make_mu_things(std::vector<Double_t>& vars, TLorentzVector * muPos, TLorentzVector * muNeg, int start)
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

    m_out_tree->Branch(("costh_HX" + varsuffix).c_str(), &vars.at(++id));
    m_out_tree->Branch(("phi_HX" + varsuffix).c_str(), &vars.at(++id));
    m_out_tree->Branch(("cosalpha_HX" + varsuffix).c_str(), &vars.at(++id));

    m_out_tree->Branch(("costh_PX" + varsuffix).c_str(), &vars.at(++id));
    m_out_tree->Branch(("phi_PX" + varsuffix).c_str(), &vars.at(++id));
    m_out_tree->Branch(("cosalpha_PX" + varsuffix).c_str(), &vars.at(++id));

    m_out_tree->Branch(("costh_CS" + varsuffix).c_str(), &vars.at(++id));
    m_out_tree->Branch(("phi_CS" + varsuffix).c_str(), &vars.at(++id));
    m_out_tree->Branch(("cosalpha_CS" + varsuffix).c_str(), &vars.at(++id));
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

void ChibPreselection::fill_photon_vars(std::vector<TLorentzVector*> photon_vecs) // has to correspond to setup_photon_vars
{
  int id = -1;
  for (auto vec : photon_vecs) {
    photon_vars.at(++id) = vec->Rapidity();
    photon_vars.at(++id) = vec->Eta();
    photon_vars.at(++id) = vec->E();
  }
}

void ChibPreselection::setup_photon_vars() // has to correspond to fill_photon_vars
{
  for (size_t i = 0; i < 32; ++i) photon_vars.push_back(0);
  int id = -1;
  std::vector<std::string> suffixes;
  if (m_year == 2016) suffixes = { "", "_rf1S", "_rf2S", "_rf3S" };
  if (m_year == 2017) suffixes = { "" };

  for (auto &suf : suffixes) {
    m_out_tree->Branch(("photon_rap" + suf).c_str(), &photon_vars.at(++id));
    m_out_tree->Branch(("photon_eta" + suf).c_str(), &photon_vars.at(++id));
    m_out_tree->Branch(("photon_energy" + suf).c_str(), &photon_vars.at(++id));
  }
}

bool ChibPreselection::accept_muon(const TLorentzVector * mu)
{
  // TODO: rapidity dependent muon cut

  if (mu->Pt() < 3.) return false;
  return true;
}


int main(int argc, char **argv) {

  ArgParser parser(argc, argv);
  auto infiles = parser.getOptionVal<std::vector<std::string> >("--infiles");
  auto outfile = parser.getOptionVal<std::string>("--outfile");
  auto intree = parser.getOptionVal < std::string>("--intree");
  auto outtree = parser.getOptionVal<std::string>("--outtree", "data");
  auto nEvents = parser.getOptionVal<Long64_t>("--events", -1);
  auto nThreads = parser.getOptionVal<Long64_t>("--threads", 8);
  auto recreate = parser.getOptionVal<bool>("--recreate", false);
  auto year = parser.getOptionVal<int>("--year", 2016);

  if (!recreate && file_exists(outfile)) {
    std::cout << "File '" << outfile << "' exists already, to force recreation use the option '--recreate true'." << std::endl;
    return -1;
  }

  ChibPreselection preselection(infiles, outfile, intree, outtree, year);
  preselection.process(nEvents, nThreads);

}