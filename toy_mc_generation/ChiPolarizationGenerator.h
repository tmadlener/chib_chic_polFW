#ifndef CHI_POL_GEN_JN_H
#define CHI_POL_GEN_JN_H

//TODO: optimize generation

#include <utility>
#include <string>
#include <chrono>

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

class ChiPolarizationGenerator
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

  // pseudo-globals
  double tmp_dimuon_M = 0;
  double tmp_chi_M = 0;
  TLorentzVector tmp_dimuon_in_chiframe;

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
  const TLorentzVector beam1;
  const TLorentzVector beam2;

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
  double R{ 2. / 3. }; // chi(j=1) helicity R = |b(m=-1)|^2 + |b(m=+1)|^2
  double R2{ 0 }; // chi(j=1) helicity R2 = |b(m=-2)|^2 + |b(m=+2)|^2
  // where R+R2+|b(m=0)|^2 = 0 has to fulfilled

  RefFrame chi_polarization_frame = RefFrame::HX;

  ULong_t rnd_seed = 0; // if 0: "... a TUUID is generated and used to fill the first 8 integers of the seed array ..."
                        //       "... the seed is guaranteed to be unique in space and time ..." (from ROOT Doc TRandom3)

  std::string outfilename;
  Mass chi_mass = PdgMassChic[1];
  Mass dimuon_mass = PdgMassPsi;
  double dimuon_mass_pdg = PdgMassPsi.central;

  double min_pT = 7; //chi
  double max_pT = 30; // chi
  double min_rap = 0; // Absolute Rapidity
  double max_rap = 1.3;
  int chi_state = 1;
  TF1* pT_distr = nullptr;
  TF1* rap_distr = nullptr;
  TF1* photonCrystalBall = nullptr; // smearing function for photons (fit to MC + some tuning)

  const int pT_distr_npx = 30; //The integral of the function is computed at fNpx points.
  //  If the function has sharp peaks, you should increase the number of points(SetNpx) 
  //  such that the peak is correctly tabulated at several points.) from ROOT TF1::GetRandom() Documentation

  ////////////
  void print_settings();

  bool chi_width_is_zero = false;
  bool dimuon_width_is_zero = false;

  void fill_branches();
  void setup_branches(TTree *t);

  // BRANCHES
  //

  Long64_t event_id = 0;

  TLorentzVector *chi = nullptr;
  TLorentzVector *dimuon = nullptr;
  TLorentzVector *gamma = nullptr;
  TLorentzVector* muPos = nullptr;
  TLorentzVector* muNeg = nullptr;

  double costh_dimuon_in_chi_restframe = 0;
  double phi_dimuon_in_chi_restframe = 0;
  double costh_lepton_in_dimuon_restframe = 0;
  double phi_lepton_in_dimuon_restframe = 0;
  double costh_HX = 0;
  double phi_HX = 0;
  double costh_CS = 0;
  double phi_CS = 0;

  TLorentzVector* smearedLepP = nullptr;
  TLorentzVector* smearedLepN = nullptr;
  TLorentzVector* smearedGamma = nullptr;
  TLorentzVector* smearedJpsi = nullptr;
  TLorentzVector* smearedChi = nullptr;
  TLorentzVector* halfSmearedChi = nullptr;

  double costh_HX_sm = 0;
  double phi_HX_sm = 0;
  double costh_CS_sm = 0;
  double phi_CS_sm = 0;


  // compatibility branches:

  double pT_chi = 0;
  double pT = 0;
  double pL_chi = 0;

  double y_chi = 0;
  double y = 0;
  double Mchi = 0;
  double Mpsi = 0;


  double pT_gamma = 0;
  double pL_gamma = 0;
  double y_gamma = 0;
  double pT_lepP = 0;
  double eta_lepP = 0;
  double pT_lepN = 0;
  double eta_lepN = 0;

  // angle of psi direction in chic rest frame, wrt to chosen chic polarization axis
  double cosTH_psi = 0;
  // angles of dilepton direction in the psi rest frame, wrt the psi direction in the chic rest frame
  // (axis definitions as in Fig 1b of PRD 83, 096001 (2011))
  double costh_chihe = 0;
  double phi_chihe = 0;
  // psi decay angles in the helicity frame
  double costh_he = 0;
  double phi_he = 0;
  // psi decay angles in the CS frame
  double costh_cs = 0;
  double phi_cs = 0;

  // smeared variables with "_sm" postfix
  double pT_chi_sm = 0;
  double y_chi_sm = 0;
  double M_chi_sm = 0;
  double qM_chi_sm = 0;
  double pT_gamma_sm = 0;
  double y_gamma_sm = 0;
  double eta_gamma_sm = 0;
  double pT_jpsi_sm = 0;
  double y_jpsi_sm = 0;
  double M_jpsi_sm = 0;
  double pT_lepP_sm = 0;
  double eta_lepP_sm = 0;
  double pT_lepN_sm = 0;
  double eta_lepN_sm = 0;
  double Mchic = 0;

  //TODO: Efficiency branches {lepP_eff_sm, lepN_eff_sm, gamma_eff_sm}

  //
  // END BRANCHES

  // Performance
  std::chrono::nanoseconds total_time_chigen{};
  ULong64_t total_number_chigen = 0;
  std::chrono::nanoseconds total_time_dimuongen{};
  ULong64_t total_number_dimuongen = 0;
  std::chrono::nanoseconds total_time_anglesgen{};
  ULong64_t total_number_anglesgen = 0;
  std::chrono::nanoseconds total_time_smearing{};
  ULong64_t total_number_smearing = 0;
  std::chrono::nanoseconds total_time_fillbranches{};
  ULong64_t total_number_fillbranches = 0;

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

    rap_distr = new TF1("rap_distr", [](double* x, double* p) { return 1; }, min_rap, max_rap, 0);
    //rap_distr->Npx(5000);

    pT_distr = new TF1("pT_distr", [/*&mass=chi_mass.central*/](double* x, double* p) {
      // The time consuming part is the SetParameter, that forces the integral to be updated for EACH event:
      // https://root-forum.cern.ch/t/faster-update-for-tf1/20991

      // beta: CHECK HERE FUNCTION AND PARAMETER VALUES: USE THOSE OF GLOBAL FIT (considering that this is a pT distribution, not a pT/M distribution)
      constexpr double beta = 3.39924; // same as in MC generation from Alberto, (was 3.45)
      constexpr double gamma = 0.635858; // same as in MC generation from Alberto, (was 0.73)
      constexpr double A = 1. / (beta - 2.) / gamma;
      const double pT_over_chimass = x[0]/p[0]; // pT = x[0], chimass = p[0]

      return x[0] * pow(1. + A * pT_over_chimass*pT_over_chimass, -beta);
    }, min_pT, max_pT, 1);
    pT_distr->SetNpx(pT_distr_npx); // for performance reasons

    // smearing function for photons (fit to MC + some tuning)
    photonCrystalBall = new TF1("photonCrystalBall", "ROOT::Math::crystalball_pdf(x[0], [2], [3], [1], [0])", -1.5, 1.5);
    photonCrystalBall->FixParameter(0, 0);
    photonCrystalBall->FixParameter(1, 1.7e-2);
    photonCrystalBall->FixParameter(2, 0.82); // alpha
    photonCrystalBall->FixParameter(3, 1.9); // N
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