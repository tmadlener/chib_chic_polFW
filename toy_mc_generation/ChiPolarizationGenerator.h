#ifndef CHI_POL_GEN_JN_H
#define CHI_POL_GEN_JN_H

// SIMULATION OF POLARIZED CHI events in the HELICITY HX frame 

#include <utility>
#include <string>
#include <chrono>
#include <memory>

#include "TF1.h"
#include "TRotation.h"
#include "TVector3.h"
#include "TLorentzVector.h"

#include "calcAngles_defs.h"

class TTree;

#define GENRAPIDITY 1
#define DONT_WRITE_LORENTZ // to reduce size of files

struct Mass {
  double central = 0;
  double width = 0;
};

class ChiPolarizationGenerator final
{
private:
  ChiPolarizationGenerator(const ChiPolarizationGenerator &) = delete;
  void operator=(const ChiPolarizationGenerator &) = delete;

  TLorentzVector* generate_chi(); // in center-of-mass (CM) frame
  std::pair< Angles, Angles > generate_dimuon_angles(); // first:{cosTh, phi} in chi restframe, second:{costh_lep_hx,phi_lep_hx} in chi restframe

  // Generate Dimuon, Gamma in CM frame:
  // rot2xyz in the chi restframe, boost2cm from chi restframe and cosTh and phi in chi restframe
  TLorentzVector* generate_dimuon(TRotation *rot2xyz, TVector3* boost2cm, double cosTh, double phi);
  TLorentzVector* generate_gamma(TRotation *rot2xyz, TVector3* boost2cm, double cosTh, double phi);
  std::pair<TLorentzVector*, TLorentzVector*> generate_leptons(TRotation *rot2xyz, TVector3* boost2cm, double cosTh, double phi); //first:lepPos,second:lepNeg

  TRotation get_rotation(const TLorentzVector &chi, bool for_lepton = false);

  bool accept_event(TLorentzVector* chi, TLorentzVector* dimuon, TLorentzVector* muPos, TLorentzVector* muNeg, TLorentzVector* gamma);

  double tmp_dimuon_M = 0;
  double tmp_chi_M = 0;
  TLorentzVector tmp_dimuon_in_chiframe{};

  // Constants:  [GeV]
  const Mass PdgMassPsi{ 3.097, 9.29e-5 };
  const std::array<Mass, 3> PdgMassChic{
      Mass{3.415, 0.0105}, // chi0
      Mass{3.511, 0.00088}, // chi1
      Mass{3.556, 0.002}, // chi2
  };

  const Mass PdgMassUps1S{ 9.4603, 5.402e-5 };
  const std::array<Mass, 3> PdgMassChib_1P{
      Mass{9.85944, 0}, // chib0
      Mass{9.89278, 0}, // chib1
      Mass{9.91221, 0}, // chib2
  };

  const double mass_proton{ 0.9382720 };
  const double mass_muon{ 0.10566 };

  const double p_beam = 6500;

  const double Ebeam = std::sqrt(p_beam * p_beam + mass_proton * mass_proton);
  const TLorentzVector beam1{};
  const TLorentzVector beam2{};

  ////////////
  // Settings:
  ////////////

  // some example parameter sets:
  // * chic1 unpolarized: R = 2/3, R2 = 0
  // * chic1 with lambdatheta = +1 (maximum positive): R = 0, R2 = 0
  // * chic1 with lambdatheta = -1/3 (maximum negative): R = 1, R2 = 0
  // * chic2 unpolarized: R = 2/5, R2 = 2/5
  // * chic2 with lambdatheta = +1 (maximum positive): R = 0, R2 = 1
  // * chic2 with lambdatheta = -3/5 (maximum negative): R = 0, R2 = 0
  // Simulation follows https://arxiv.org/abs/1103.4882
  double R = 2. / 3.; // chi(j=1) helicity R = |b(m=-1)|^2 + |b(m=+1)|^2
  double R2 = 0; // chi(j=1) helicity R2 = |b(m=-2)|^2 + |b(m=+2)|^2
  // where R+R2+|b(m=0)|^2 = 1

  RefFrame chi_polarization_frame = RefFrame::HX;

  ULong_t rnd_seed = 0; // if 0: "... a TUUID is generated and used to fill the first 8 integers of the seed array ..."
                        //       "... the seed is guaranteed to be unique in space and time ..." (from ROOT Doc TRandom3)

  std::string outfilename;
  Mass chi_mass = PdgMassChic[1];
  Mass dimuon_mass = PdgMassPsi;
  double dimuon_mass_pdg = PdgMassPsi.central;

  double min_pT = 5; // chi
  double max_pT = 60; // chi
  double min_rap = 0; // Absolute Rapidity
  double max_rap = 1.3;
  int chi_state = 1;
  std::unique_ptr<TF1> pTM_distr;
  std::unique_ptr<TF1> rap_distr;
  std::unique_ptr<TF1> photonCrystalBall; // smearing function for photons (fit to MC + some tuning)

  // Selections
  bool apply_selections = true;
  double min_photon_pt = 0.4;
  double max_photon_abseta = 1.5;
  bool apply_loose_muon_selection = true;
  double min_dimuon_pt = 8;
  double max_dimuon_pt = 50;
  double min_abs_dimuon_rap = 0;
  double max_abs_dimuon_rap = 1.2;

  void setup_distributions();
  
  ////////////
  void print_settings();

  bool chi_width_is_zero = false;
  bool dimuon_width_is_zero = false;

  void fill_branches();
  void setup_branches(TTree *t);

  Long64_t accepted_events{};

  // LorentzVector

  TLorentzVector m_chi{};
  TLorentzVector m_dimuon{};
  TLorentzVector m_gamma{};
  TLorentzVector m_muPos{};
  TLorentzVector m_muNeg{};
  TLorentzVector m_smearedLepP{};
  TLorentzVector m_smearedLepN{};
  TLorentzVector m_smearedGamma{};
  TLorentzVector m_smearedJpsi{};
  TLorentzVector m_smearedChi{};
  TLorentzVector m_halfSmearedChi{};


  // BRANCHES
  //

  Long64_t event_id{};

  TLorentzVector *chi{};
  TLorentzVector *dimuon{};
  TLorentzVector *gamma{};
  TLorentzVector* muPos{};
  TLorentzVector* muNeg{};

  double costh_dimuon_in_chi_restframe{};
  double phi_dimuon_in_chi_restframe{};
  double costh_lepton_in_dimuon_restframe{};
  double phi_lepton_in_dimuon_restframe{};
  double costh_HX{};
  double phi_HX{};
  double costh_CS{};
  double phi_CS{};
  double costh_PX{};
  double phi_PX{};

  TLorentzVector* smearedLepP{};
  TLorentzVector* smearedLepN{};
  TLorentzVector* smearedGamma{};
  TLorentzVector* smearedJpsi{};
  TLorentzVector* smearedChi{};
  TLorentzVector* halfSmearedChi{};

  double costh_HX_sm{};
  double phi_HX_sm{};
  double costh_CS_sm{};
  double phi_CS_sm{};
  double costh_PX_sm{};
  double phi_PX_sm{};


  // compatibility branches:

  double pT_chi{};
  double pT{};
  double pL_chi{};

  double y_chi{};
  double y{};
  double Mchi{};
  double Mpsi{};


  double pT_gamma{};
  double pL_gamma{};
  double y_gamma{};
  double pT_lepP{};
  double eta_lepP{};
  double pT_lepN{};
  double eta_lepN{};

  // angle of psi direction in chic rest frame, wrt to chosen chic polarization axis
  double cosTH_psi{};
  // angles of dilepton direction in the psi rest frame, wrt the psi direction in the chic rest frame
  // (axis definitions as in Fig 1b of PRD 83, 096001 (2011))
  double costh_chihe{};
  double phi_chihe{};
  // psi decay angles in the helicity frame
  double costh_he{};
  double phi_he{};
  // psi decay angles in the CS frame
  double costh_cs{};
  double phi_cs{};

  // smeared variables with "_sm" postfix
  double pT_chi_sm{};
  double y_chi_sm{};
  double M_chi_sm{};
  double qM_chi_sm{};
  double pT_gamma_sm{};
  double y_gamma_sm{};
  double eta_gamma_sm{};
  double pT_jpsi_sm{};
  double y_jpsi_sm{};
  double M_jpsi_sm{};
  double pT_lepP_sm{};
  double eta_lepP_sm{};
  double pT_lepN_sm{};
  double eta_lepN_sm{};
  double Mchic{};

  //TODO: Efficiency branches {lepP_eff_sm, lepN_eff_sm, gamma_eff_sm}

  //
  // END BRANCHES

  // Performance
  std::chrono::nanoseconds total_time_chigen{};
  ULong64_t total_number_chigen{};
  std::chrono::nanoseconds total_time_dimuongen{};
  ULong64_t total_number_dimuongen{};
  std::chrono::nanoseconds total_time_anglesgen{};
  ULong64_t total_number_anglesgen{};
  std::chrono::nanoseconds total_time_smearing{};
  ULong64_t total_number_smearing{};
  std::chrono::nanoseconds total_time_fillbranches{};
  ULong64_t total_number_fillbranches{};
  std::chrono::nanoseconds total_time_selection{};
  ULong64_t total_number_selection{};

  void print_performance();

public:
  ChiPolarizationGenerator(const std::string &filename) :
    outfilename(filename),
    beam1(0, 0, -p_beam, Ebeam),
    beam2(0., 0., p_beam, Ebeam)
  {
    //TODO: either ups,psi chooser or mass parameter in constructor
    // and either 1,2,3 chooser or chi mass parameter in constructor:
    chi_mass = PdgMassChic[chi_state];
    dimuon_mass = PdgMassPsi;

  }

  ~ChiPolarizationGenerator() {
    print_performance();
  }

  void generate(ULong64_t n = 123456);
  void setChiHelicityFractions(double R_1, double R_2 = 0)
  {
    if (R_1 > 1 || R_1 < 0 || (chi_state > 1 && (R_1 + R_2 > 1 || R_2 > 1 || R_2 < 0))) {
      std::cout << "% % % ERROR: R1[" << R_1 << "] and R2[" << R_2 << "] not changed,\n\t they have to fulfill R1+R2<1 and 0<R1<1, 0<R2<1\n"
        "\t (current values: " << R << ", " << R2 << ")" << std::endl;
      return;
    }
    R = R_1;
    R2 = R_2;
  }

  void setSeed(ULong_t seed = 0) { rnd_seed = seed; }

  void setChib(int chib_state = 1)
  {
    chi_width_is_zero = false;
    dimuon_width_is_zero = false;
    if (chib_state < 0 || chib_state > 2) {
      std::cout << "% % % ERROR: NO VALID CHIB STATE (has to be 0, 1 or 2): " << chib_state << '\n'
        << "\tSetting it to 1." << std::endl;
      chib_state = 1;
    }
    chi_state = chib_state;
    chi_mass = PdgMassChib_1P[chi_state];
    dimuon_mass = PdgMassUps1S;
    dimuon_mass_pdg = PdgMassUps1S.central;
    if (chi_mass.width == 0) {
      std::cout << "WARNING: chi_mass width is zero" << std::endl;
      chi_width_is_zero = true;
    }
    if (dimuon_mass.width == 0) {
      std::cout << "WARNING: dimuon_mass width is zero" << std::endl;
      dimuon_width_is_zero = true;
    }
}
  void applySelections(bool apply = true) { apply_selections = apply; }

  void setKinematics(double pt_min, double pt_max, double absrap_min, double absrap_max) {
    if (pt_min < 0 || pt_max < pt_min || absrap_min < 0 || absrap_max < absrap_min) {
      std::cout << "% % % ERROR: Impossible pT or rapidity values, nothing gets changed:\n"
        "\t0 < " << pt_min << " < pT < " << pt_max << "\n"
        "\t0 < " << absrap_min << " < rapidity < " << absrap_max << "\n"
        "current values: pT[" << min_pT << ", " << max_pT << "], rap[" << min_rap << ", " << max_rap << "]" << std::endl;
      return;
    }

    min_pT = pt_min;
    max_pT = pt_max;
    min_rap = absrap_min;
    max_rap = absrap_max;
  }

};

#endif