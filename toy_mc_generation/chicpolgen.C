#include "smearing.h"
#include "efficiencies.h"
#include "select.h"

#include "../general/interface/calcAngles.h"

#include "Riostream.h"
#include "TSystem.h"
#include "TString.h"
#include "TROOT.h"
#include "TTree.h"
#include "TMath.h"
#include "TRandom3.h"
#include "TVector3.h"
#include "TLorentzVector.h"
#include "TRotation.h"
#include "TF1.h"
#include "TFile.h"

#include <iostream>
#include <string>
#include <memory>

#define GENRAPIDITY 1 // set to zero to generate the chic eta instead of the chic rapidity

/** * Configuration and settings for the generation
 *
 * Including some default values
 */
struct gen_config {
  int n_events{3000000};

  // some example parameter sets:
  // * chic1 unpolarized: R = 2/3, R2 = 0
  // * chic1 with lambdatheta = +1 (maximum positive): R = 0, R2 = 0
  // * chic1 with lambdatheta = -1/3 (maximum negative): R = 1, R2 = 0
  // * chic2 unpolarized: R = 2/5, R2 = 2/5
  // * chic2 with lambdatheta = +1 (maximum positive): R = 0, R2 = 1
  // * chic2 with lambdatheta = -3/5 (maximum negative): R = 0, R2 = 0
  int chic_state{1};               //  0 = chi_c0,  1 = chi_c1,  2 = chi_c2
  double R{2./3.};                 // fraction of chic helicity 1: 0 < R  < 1
  double R2{0};                   // fraction of chic helicity 2: 0 < R2 < 1

  double pbeam{4000.}; // (this is actually irrelevant, as long as pbeam >> proton mass)
  double y_min{0.0}; // min abs rapidity of the chi
  double y_max{1.3}; // max abs rapidity of the chi

  double pT_min{7.0}; // min pt of the chi
  double pT_max{23.0}; // max pt of the chi

  bool CSframeIsNatural{false};   // generate chic polarization in the CS frame (true)
                                  // or in the HX frame (false)

  size_t n_accepted{0}; // number of events that need to be accepted before the generation stops.
                      // Set to 0, to never stop (i.e. simply run until n_events) is reached
                      // This only takes effect if the number of accepted events is reached
                      // before the number of generated events (n_events) is exhausted, so that this
                      // can only be used to stop early, not to guarantee a certain number of
                      // accepted events

  std::string genfile{"chicpolgen.root"}; // name of the output file
  // To not produce efficiency branches leave the efficiency file names empty
  std::string muonEffs{""}; // file name from where the muon efficiencies should be loaded
  std::string photonEffs{""}; // file name from where the photon efficiencies should be loaded

  void print(bool verbose) const {
    std::cout << "--------------------------------------------------\n";
    std::cout << "gen_config settings used for generation:\n"
              << "n_events = " << n_events << '\n'
              << "chic_state = " << chic_state << '\n'
              << "R (fraction of helicity 1) = " << R << '\n'
              << "R2 (fraction of helicity 2) = " << R2 << '\n'
              << "min pT = " << pT_min << " GeV\n"
              << "max pT = " << pT_max << " GeV\n";
    if (!muonEffs.empty()) {
      std::cout << "muon efficiency file: " << muonEffs << '\n';
    }
    if (!photonEffs.empty()) {
      std::cout << "photon efficiency file: " << photonEffs << '\n';
    }
    if (verbose) {
      std::cout << "generate in the CS frame: " << std::boolalpha << CSframeIsNatural << '\n'
                << "min abs rapidity = " << y_min << '\n'
                << "max abs rapidity = " << y_max << '\n'
                << "beam energy = " << pbeam << " GeV\n";
    }
    std::cout << "output file: " << genfile << '\n';
    std::cout << "--------------------------------------------------" << std::endl;
  }
};

/**
 * Struct holding the width and the central value for a given state
 */
struct Mass {
  double central;
  double width;
};

const struct MassSettings {
  static constexpr Mass MpsiPDG{3.097, 9.29e-5}; // J/psi mass and width
  static constexpr std::array<Mass, 3> MchiPDG{
    Mass{3.415, 0.0105}, // chic0
      Mass{3.511, 0.00088}, // chic1
        Mass{3.556, 0.002}, // chic2
          };

  static constexpr double Mprot{0.9382720};
  static constexpr double Mlepton{0.10566}; // muon GeV
} GenMassSettings;

// need to define the array, so the linker can see it
// declaration and initialization are handled by the definition of the class above
constexpr std::array<Mass, 3> MassSettings::MchiPDG;
constexpr Mass MassSettings::MpsiPDG;
constexpr double MassSettings::Mprot;

// global state variables necessary to hold the chi and psi mass in a given event
// NOTE: when refactoring make these parameters that get passed around
double Mchi;
double Mpsi;

const double PIG = TMath::Pi();

// tmadlener, 08.06.2018: no longer necessary setup for smearing according to MC distributions
// constexpr auto residualMapFileName = "res_maps.root"; // The file containing the smearing maps for photons and muons
// constexpr auto photonMapName = "photon_rel_res_map"; // The name of the photon smearing map in the file
// constexpr auto muonXYMapName = "muon_xy_rel_res_map"; // The name of the muonXY smearing map in the file
// constexpr auto muonZMapName = "muon_z_rel_res_map"; // The name of the muonZ smearing map in the file

template<typename Func>
double getMass(Func distribution)
{
  return distribution();
}

double func_rap_gen(double* /*x*/, double* /*par*/)
{
  return   1.;
}

double func_pT_gen(double* x, double* /*par*/)
{
  // const double beta = 3.45;  //  CHECK HERE FUNCTION AND PARAMETER VALUES: USE THOSE OF GLOBAL FIT (considering that this is a pT distribution, not a pT/M distribution)
  const double beta = 3.39924;  // same as in MC generation from Alberto
  // const double gamma = 0.73;
  const double gamma = 0.635858; // same as in MC generation from Alberto
  double pT = x[0];
  return pT * pow( 1. + 1./(beta - 2.) * pT/Mchi*pT/Mchi / gamma, -beta  );
}


void chicpolgen(const gen_config& config = gen_config{}){
  gROOT->SetBatch();
  config.print(true);
  // translate configuration into const variables
  const double Ebeam = std::sqrt(config.pbeam * config.pbeam + GenMassSettings.Mprot * GenMassSettings.Mprot);
  const TLorentzVector targ(0., 0. , -config.pbeam, Ebeam); // "targ" = second beam
  const TLorentzVector beam(0., 0. , config.pbeam, Ebeam);

  const double chic_state = config.chic_state;
  const double R = config.R;
  const double R2 = config.R2;
  const double n_events = config.n_events;

  const double y_min = config.y_min; const double y_max = config.y_max;
  const double pT_min = config.pT_min; const double pT_max = config.pT_max;

  const double Mlepton = GenMassSettings.Mlepton;
  const double MpsiPDG = GenMassSettings.MpsiPDG.central;

  const bool check_accept = config.n_accepted > 0;

  EffProv *muonEffs = nullptr;
  EffProv *photonEffs = nullptr;

  if (!config.muonEffs.empty()) {
    std::cout << "Reading muon efficiencies\n";
    // non-parametrized single muon efficiencies
    // muonEffs = new EfficiencyProvider<TGraphAsymmErrors>(config.muonEffs, "muon_eff_pt", RangeFromGraph{});
    // parametrized single muon efficiencies
    muonEffs = new EfficiencyProvider<TF1>(config.muonEffs, "muon_eff_pt", RangeFromFit{});
  }
  if (!config.photonEffs.empty()) {
    std::cout << "Reading photon efficiencies\n";
    // using a fixed range for the photon efficiencies since they are parametrized in the region from 400 MeV to 7 GeV
    // photonEffs = new EfficiencyProvider<TF1>(config.photonEffs, "photon_eff_pt", FixedRange<TF1>{0.4, 7});
    photonEffs = new EfficiencyProvider<TF1>(config.photonEffs, "photon_eff_pt", RangeFromFit{});
  }

  // Selectors to act on the smeared variables
  // TODO: Make this configurable in some way
  const auto jpsiSelector = std::make_unique<PtRangeAbsRapiditySelector>(Range{8, 20}, 1.25);
  const auto muonSelector = std::make_unique<LooseMuonSelector>();
  const auto photonSelector = std::make_unique<MinPtMaxEtaSelector>(0.41, 1.5);

  // To accept all events without filters
  // const auto jpsiSelector = std::make_unique<AllSelector>();
  // const auto muonSelector = std::make_unique<AllSelector>();
  // const auto photonSelector = std::make_unique<AllSelector>();

  delete gRandom;
  gRandom = new TRandom3(0);

  // create the functions from which the chic and the jpsi masses are drawn
  const auto chicMassDist = std::bind(&TRandom3::Gaus, std::ref(gRandom),
                                      GenMassSettings.MchiPDG[config.chic_state].central,
                                      GenMassSettings.MchiPDG[config.chic_state].width);

  const auto jpsiMassDist = std::bind(&TRandom3::Gaus, std::ref(gRandom),
                                      MpsiPDG,
                                      GenMassSettings.MpsiPDG.width);

  // tmadlener 08.06.2018: For generating according to smeared distributions
  // auto *modelFile = TFile::Open("./mass_distributions_data.root");
  // if (chic_state == 1) {
  //   chicMass = static_cast<TF1*>(modelFile->Get("chic1_mass"));
  // }
  // if (chic_state == 2) {
  //   chicMass = static_cast<TF1*>(modelFile->Get("chic2_mass"));
  // }

  // jpsiMass = static_cast<TF1*>(modelFile->Get("jpsi_mass"));


  TF1* pT_distr = new TF1("pT_distr",func_pT_gen,pT_min,pT_max,0);
  TF1* rap_distr = new TF1("rap_distr",func_rap_gen,y_min,y_max,0);


  TFile* hfile = new TFile( config.genfile.c_str(), "RECREATE", "chicpolgen");

  TTree* tr = new TTree("tr", "tr");

  double pT_chi;        tr->Branch( "gen_chicPt",         &pT_chi);
  double pT;            tr->Branch( "gen_JpsiPt",             &pT);
#if GENRAPIDITY == 1
  double pL_chi;     //   tr->Branch( "pL_chi",         &pL_chi,         "pL_chi/D" );
#endif
  // double pL;         //   tr->Branch( "pL",             &pL,             "pL/D" );
  double y_chi;         tr->Branch( "gen_chicRap",          &y_chi);
  double y;             tr->Branch( "gen_JpsiRap",              &y);

  tr->Branch("gen_chicMass", &Mchi);
  tr->Branch("gen_JpsiMass", &Mpsi);

  double pT_gamma;      tr->Branch( "gen_photonPt",       &pT_gamma);
  double pL_gamma;      tr->Branch( "gen_photonPl",       &pL_gamma);
  double y_gamma;       tr->Branch( "gen_photonEta",        &y_gamma);

  double pT_lepP;       tr->Branch( "gen_muPPt",        &pT_lepP);
  double eta_lepP;      tr->Branch( "gen_muPEta",       &eta_lepP);

  double pT_lepN;       tr->Branch( "gen_muNPt",        &pT_lepN);
  double eta_lepN;      tr->Branch( "gen_muNEta",       &eta_lepN);

  // int inAcc0;            tr->Branch( "inAcc0",          &inAcc0,          "inAcc0/I"    );
  // int inAcc1;            tr->Branch( "inAcc1",          &inAcc1,          "inAcc1/I"    );

// angle of psi direction in chic rest frame, wrt to chosen chic polarization axis
  double cosTH_psi;     tr->Branch( "cosTH_psi",      &cosTH_psi,      "cosTH_psi/D" );

// angles of dilepton direction in the psi rest frame, wrt the psi direction in the chic rest frame
// (axis definitions as in Fig 1b of PRD 83, 096001 (2011))
  double costh_chihe;   tr->Branch( "costh_chihe",    &costh_chihe,    "costh_chihe/D" );
  double phi_chihe;     tr->Branch( "phi_chihe",      &phi_chihe,      "phi_chihe/D" );

// psi decay angles in the helicity frame
  double costh_he;      tr->Branch( "gen_costh_HX",       &costh_he);
  double phi_he;        tr->Branch( "gen_phi_HX",         &phi_he);

// psi decay angles in the CS frame
  double costh_cs;      tr->Branch( "gen_costh_CS",       &costh_cs);
  double phi_cs;        tr->Branch( "gen_phi_CS",         &phi_cs);


  // double M_gamma; tr->Branch("M_gamma", &M_gamma);
  // double qM_chi; tr->Branch("qM_chi", &qM_chi);

  // smeared variables with "_sm" postfix
  double pT_chi_sm;     tr->Branch("chicPt", &pT_chi_sm);
  double y_chi_sm;     tr->Branch("chicRap", &y_chi_sm);
  double M_chi_sm;      tr->Branch("mumugammaMass", &M_chi_sm);
  double qM_chi_sm;      tr->Branch("chicMass", &qM_chi_sm);

  double pT_gamma_sm;     tr->Branch("photonPt", &pT_gamma_sm);
  double y_gamma_sm;     tr->Branch("photonEta", &y_gamma_sm);
  double eta_gamma_sm;     tr->Branch("eta_gamma_sm", &eta_gamma_sm);

  double pT_jpsi_sm;     tr->Branch("JpsiPt", &pT_jpsi_sm);
  double y_jpsi_sm;     tr->Branch("JpsiRap", &y_jpsi_sm);
  double M_jpsi_sm;      tr->Branch("JpsiMass", &M_jpsi_sm);

  double pT_lepP_sm;     tr->Branch("muPPt", &pT_lepP_sm);
  double eta_lepP_sm;     tr->Branch("muPEta", &eta_lepP_sm);

  double pT_lepN_sm;     tr->Branch("muNPt", &pT_lepN_sm);
  double eta_lepN_sm;     tr->Branch("muNEta", &eta_lepN_sm);

  double Mchic;    tr->Branch("Q_value_gen", &Mchic);

  double costh_HX_sm; tr->Branch("costh_HX", &costh_HX_sm);
  double phi_HX_sm; tr->Branch("phi_HX", &phi_HX_sm);

  double costh_CS_sm; tr->Branch("costh_CS", &costh_CS_sm);
  double phi_CS_sm; tr->Branch("phi_CS", &phi_CS_sm);

  // double ca_gamma_jpsi;     tr->Branch("ca_gamma_jpsi", &ca_gamma_jpsi);
  // double ca_mu_mu;     tr->Branch("ca_mu_mu", &ca_mu_mu);
  // double ca_sm_gamma_jpsi;     tr->Branch("ca_sm_gamma_jpsi", &ca_sm_gamma_jpsi);
  // double ca_sm_mu_mu;     tr->Branch("ca_sm_mu_mu", &ca_sm_mu_mu);

  // double lepP_eff, lepN_eff;
  double lepP_eff_sm, lepN_eff_sm;
  // double gamma_eff;
  double gamma_eff_sm;

  if (!config.muonEffs.empty()) {
    // tr->Branch("lepP_eff", &lepP_eff);
    // tr->Branch("lepN_eff", &lepN_eff);
    tr->Branch("lepP_eff_sm", &lepP_eff_sm);
    tr->Branch("lepN_eff_sm", &lepN_eff_sm);
  }

  if (!config.photonEffs.empty()) {
    // tr->Branch("gamma_eff", &gamma_eff);
    tr->Branch("gamma_eff_sm", &gamma_eff_sm);
  }

  // smearing initialization (for smearing according to MC)
  // auto *smearingFile = TFile::Open(residualMapFileName);
  // auto *photonResMap = static_cast<TH2D*>(smearingFile->Get(photonMapName));
  // const auto photonSmearing = SmearingProvider(photonResMap);

  // auto *muonXYResMap = static_cast<TH2D*>(smearingFile->Get(muonXYMapName));
  // const auto muonXYSmearing = SmearingProvider(muonXYResMap);

  // auto *muonZResMap = static_cast<TH2D*>(smearingFile->Get(muonZMapName));
  // const auto muonZSmearing = SmearingProvider(muonZResMap);


  // smearing function for photons (fit to MC + some tuning)
  TF1* photonCrystalBall = new TF1("photonCrystalBall", "ROOT::Math::crystalball_pdf(x[0], [2], [3], [1], [0])", -1.5, 1.5);
  photonCrystalBall->FixParameter(0, 0);
  photonCrystalBall->FixParameter(1, 1.7e-2);
  photonCrystalBall->FixParameter(2, 0.82); // alpha
  photonCrystalBall->FixParameter(3, 1.9); // N

  photonCrystalBall->Draw();

  const int n_step = n_events/50;
  std::cout << '\n';
  std::cout << "------------------------------------------------------------" << '\n';
  std::cout << "Progress: ";

  size_t accepted = 0;
  int i_event = 0;
/////////////////// CYCLE OF EVENTS ////////////////////////
  for(; i_event < n_events; i_event++){

    if ((i_event + 1) % n_step == 0) {
      std::cout << "X";  std::cout.flush();
    }

  // generation of chic in the CMS of the proton-proton event

    // M:
    // generate the chic mass and the jpsi mass such that they are correlated similar to what they are in data:
    // 1) get the chic mass and the J/psi mass according to the fitted data distributions
    // 2) in data the chic mass is actually the mass from the KVF but we assume that we can treat it as a Q-value
    //    (M_mumugamma - M_mumu + M_jpsi), so we reverse it to get to the chic mass we actually want to generate
    Mchic = getMass(chicMassDist);
    Mpsi = getMass(jpsiMassDist);
    Mchi = Mchic - MpsiPDG + Mpsi;

    // Have to differentiate between generating the rapidity according to the rapidity distribution or
    // generating eta according to the rapidity distribution
    TLorentzVector chi;
    // Phi:
    const double Phi_chi   = 2. * PIG * gRandom->Rndm(); // needed for any of the two cases

    // pT:

    pT_chi = pT_distr->GetRandom();

#if GENRAPIDITY == 1
    // pL:
    double rap_sign = gRandom->Uniform(-1., 1.); rap_sign /= fabs(rap_sign);
    y_chi = rap_distr->GetRandom() * rap_sign;

    double mT = sqrt( Mchi*Mchi + pT_chi*pT_chi );
    double pL1 = 0.5 *mT * exp(y_chi);
    double pL2 = - 0.5 *mT * exp(-y_chi);
    pL_chi = pL1 + pL2;

    chi.SetXYZM( pT_chi * cos(Phi_chi) , pT_chi * sin(Phi_chi), pL_chi, Mchi );
#else
    const bool eta_sign = gRandom->Uniform(-1., 1) > 0;
    const double genEta = rap_distr->GetRandom() * (1 - 2 * eta_sign);

    chi.SetPtEtaPhiM(pT_chi, genEta, Phi_chi, Mchi);
    y_chi = chi.Rapidity();
#endif

  // generation of full angular distribution


    const double angdistr_max = 0.02;


    double angdistr_rnd;
    // initialize this to the max double value, so that a printout appears if an invalid setting for
    // chic_state has been chosen
    // NOTE: this also makes the macro go into an infinite loop
    double angdistr = std::numeric_limits<double>::max();


    double sinTH_psi  = 100.;
    double PHI_psi = 2. * PIG * gRandom->Rndm();
    double sinth_chihe = 100.;
    double cosphi_chihe = 100.;


    do {
         cosTH_psi = -1. + 2. * gRandom->Rndm();
              // direction of the PSI in the CHI rest frame (wrt to a reference frame, HE or CS, chosen afterwards)
              // PHI_psi is the second coordinate, generated outside the loop
              // because the global angular decay distribution does not depend on it.

         costh_chihe = -1. + 2. * gRandom->Rndm();  // direction of the lepton in the PSI rest frame
         phi_chihe   = 360. * gRandom->Rndm();      // (wrt the PSI direction seen from the CHI rest frame)


         double cosTH2_psi = cosTH_psi*cosTH_psi;
         double cosTH4_psi = cosTH2_psi*cosTH2_psi;
         double costh2_chihe = costh_chihe*costh_chihe;
         double sinth2_chihe = 1 - costh2_chihe;

         sinTH_psi   = sqrt( 1. -   cosTH2_psi );
         sinth_chihe = sqrt( sinth2_chihe );

         double sin2TH_psi   = 2.*sinTH_psi*cosTH_psi;
         double sin2th_chihe = 2.*sinth_chihe*costh_chihe;

         cosphi_chihe = cos( phi_chihe * PIG/180. );
         double cos2phi_chihe = 2.*cosphi_chihe*cosphi_chihe -1.;



      // chic_0 angular distribution
         if ( chic_state == 0 ) {


           angdistr = 1. + costh2_chihe;

           angdistr *= 3. / ( 64.* PIG*PIG );


         }



      // chic_1 angular distribution
         if ( chic_state == 1 ) {

         //  double a2 = gRandom->Gaus( -0.006, 0.013 );  // from Crystal Ball and E835 measurements
           double a2 = 0.;

           double a1 = sqrt( 1. - a2*a2 );   // (a1 taken to be positive)

           double A0 = sqrt(1./2.) * ( a1 + a2 );
           double A1 = sqrt(1./2.) * ( a1 - a2 );

           double k1 = A1*A1 + 1./2.* R * ( A0*A0 - A1*A1 );
           double k2 = ( 1. - 3./2.* R ) * ( A0*A0 - A1*A1 );
           double k3 = -A1*A1 + 1./2.* R;
           double k4 = 1. - 3./2.* R;
           double k5 = 1./4.* A1*A0 * ( 3.* R - 2. );

           angdistr = k1 + k2 * cosTH2_psi + ( k3 + k4 * cosTH2_psi ) * costh2_chihe
                         + k5 * sin2TH_psi * sin2th_chihe * cosphi_chihe ;

           angdistr *= 9. / ( 64.* PIG*PIG );



         }


      // chic_2 angular distribution
         else if ( chic_state == 2 ) {

        //   double a3 = gRandom->Gaus( 0.01, 0.04 );  // from E760 and E835 measurements
           double a3 = 0.;

        //   double a2 = gRandom->Gaus( -0.13, 0.04 );  // from Crystal Ball, E760 and E835 measurements
           double a2 = 0.;

           double a1 = sqrt( 1. - a2*a2 - a3*a3 );  // (a1 taken to be positive)

           double A0 = sqrt(1./10.)*a1 + sqrt(1./2.)*a2 + sqrt(2./5.)*a3;
           double A1 = sqrt(3./10.)*a1 + sqrt(1./6.)*a2 - sqrt(8./15.)*a3;
           double A2 = sqrt(3./5.)*a1  - sqrt(1./3.)*a2 + sqrt(1./15.)*a3;


/*
           double k1 = 1./8.* ( 2.* A0*A0 + 3.* A2*A2 - R * ( 2.* A0*A0 - 4.* A1*A1 + A2*A2 ) );
           double k2 = 3./4.* ( -2.* A0*A0 + 4.* A1*A1 - A2*A2 + R * ( 4.* A0*A0 -6.* A1*A1 + A2*A2 ) );
           double k3 = 1./8.* ( 6.* A0*A0 -8.* A1*A1 + A2*A2 ) * ( 3. -5.* R );
           double k4 = 1./8.* ( 2.* A0*A0 +3.* A2*A2 -R * ( 2.* A0*A0 + 4.* A1*A1 + A2*A2 ) );
           double k5 = 3./4.* ( -2.* A0*A0 -4.* A1*A1 - A2*A2 +R * ( 4.* A0*A0 + 6.* A1*A1 + A2*A2 ));
           double k6 = 1./8.* ( 6.* A0*A0 + 8.* A1*A1 + A2*A2 ) * ( 3. - 5.* R );
           double k7 = sqrt(6.)/4.* ( R - 1. ) * A0*A2;
           double k8 = sqrt(6.)/4.* ( 4. - 6.* R) * A0*A2;
           double k9 = sqrt(6.)/4.* ( 5.* R - 3. ) * A0*A2;
           double k10 = sqrt(3.)/4.* ( A0*A1 + sqrt(3./2.)* A1*A2 - R * ( 2.* A0*A1 + sqrt(3./2.)* A1*A2 ) );
           double k11 = 1./(4.*sqrt(3.))* ( 5.* R - 3. ) * ( 3.* A0*A1 + sqrt(3./2.)* A1*A2 );
*/

           double k1_0 = 1./4.* A0*A0 + 3./8.* A2*A2;
           double k1_1 = 1./2.* A1*A1 + 1./4.* A2*A2;
           double k1_2 = 3./8.* A0*A0 + 1./2.* A1*A1 + 1./16.* A2*A2;
           double k2_0 = -3./2.* A0*A0 + 3.* A1*A1 - 3./4.* A2*A2;
           double k2_1 = 3./2.* A0*A0 - 3./2.* A1*A1;
           double k2_2 = -3./4.* A0*A0 + 3./8.* A2*A2;
           double k3_0 = 9./4.* A0*A0 - 3.* A1*A1 + 3./8.* A2*A2;
           double k3_1 = -3./2.* A0*A0 + 2.* A1*A1 - 1./4.* A2*A2;
           double k3_2 = 3./8.* A0*A0 - 1./2.* A1*A1 + 1./16.* A2*A2;
           double k4_0 = 1./4.* A0*A0 + 3./8.* A2*A2;
           double k4_1 = -1./2.* A1*A1 + 1./4.* A2*A2;
           double k4_2 = 3./8.* A0*A0 - 1./2.* A1*A1 + 1./16.* A2*A2;
           double k5_0 = -3./2.* A0*A0 - 3.*A1*A1 - 3./4.* A2*A2;
           double k5_1 = 3./2.* A0*A0 + 3./2.* A1*A1;
           double k5_2 = -3./4.* A0*A0 + 3./8.* A2*A2;
           double k6_0 = 9./4.* A0*A0 + 3.* A1*A1 + 3./8.* A2*A2;
           double k6_1 = -3./2.* A0*A0 - 2.*A1*A1 - 1./4.* A2*A2;
           double k6_2 = 3./8.* A0*A0 + 1./2.* A1*A1 + 1./16.* A2*A2;
           double k7_0 = -sqrt(6.)/4.* A0*A2;
           double k7_1 = 0.;
           double k7_2 = sqrt(6.)/8.* A0*A2;
           double k8_0 = sqrt(6.)* A0*A2;
           double k8_1 = -sqrt(6.)/2. *A0*A2;
           double k8_2 = 0.;
           double k9_0 = -3.* sqrt(6.)/4.* A0*A2;
           double k9_1 = sqrt(6.)/2.* A0*A2;
           double k9_2 = -sqrt(6.)/8.* A0*A2;
           double k10_0 = sqrt(3.)/4.* A0*A1 + 3.*sqrt(2.)/8.* A1*A2;
           double k10_1 = -sqrt(3.)/4.* A0*A1;
           double k10_2 = sqrt(3.)/8.* A0*A1 - 3.*sqrt(2.)/16.* A1*A2;
           double k11_0 = -3.*sqrt(3.)/4.* A0*A1 - 3.*sqrt(2.)/8.*A1*A2;
           double k11_1 = sqrt(3.)/2.* A0*A1 + sqrt(2.)/4.* A1*A2;
           double k11_2 = -sqrt(3.)/8.* A0*A1 - sqrt(2.)/16.* A1*A2;

           double R1 = R; double R0=1.-R1-R2;

           double k1 = R0*k1_0 +R1*k1_1 +R2*k1_2;
           double k2 = R0*k2_0 +R1*k2_1 +R2*k2_2;
           double k3 = R0*k3_0 +R1*k3_1 +R2*k3_2;
           double k4 = R0*k4_0 +R1*k4_1 +R2*k4_2;
           double k5 = R0*k5_0 +R1*k5_1 +R2*k5_2;
           double k6 = R0*k6_0 +R1*k6_1 +R2*k6_2;
           double k7 = R0*k7_0 +R1*k7_1 +R2*k7_2;
           double k8 = R0*k8_0 +R1*k8_1 +R2*k8_2;
           double k9 = R0*k9_0 +R1*k9_1 +R2*k9_2;
           double k10= R0*k10_0+R1*k10_1+R2*k10_2;
           double k11= R0*k11_0+R1*k11_1+R2*k11_2;

           angdistr = k1 + k2 * cosTH2_psi + k3 * cosTH4_psi + ( k4 + k5 * cosTH2_psi + k6 * cosTH4_psi ) * costh2_chihe
                         + ( k7 + k8 * cosTH2_psi + k9 * cosTH4_psi ) * sinth2_chihe * cos2phi_chihe
                         + ( k10 + k11 * cosTH2_psi )* sin2TH_psi * sin2th_chihe * cosphi_chihe ;

           angdistr *= 15. / ( 64.* PIG*PIG );
         }

         if (angdistr > angdistr_max) { std::cout << "PASSED LIMIT" << std::endl; }

         angdistr_rnd = angdistr_max * gRandom->Rndm();

    } while ( angdistr_rnd > angdistr );



 // psi 4-momentum in the chi rest frame, wrt the chosen chi_c polarization axes:

    double p_psi_chi = 0.5 * ( Mchi*Mchi - Mpsi*Mpsi ) / Mchi;

    TLorentzVector psi_chi;
    psi_chi.SetXYZM( p_psi_chi * sinTH_psi * cos(PHI_psi),
                     p_psi_chi * sinTH_psi * sin(PHI_psi),
                     p_psi_chi * cosTH_psi,
                     Mpsi );


 // gamma 4-momentum in the chi rest frame, wrt the chosen chi_c polarization axes:

    double p_gamma_chi = p_psi_chi;

    TLorentzVector gamma_chi;
    gamma_chi.SetXYZM( -p_gamma_chi * sinTH_psi * cos(PHI_psi),
                       -p_gamma_chi * sinTH_psi * sin(PHI_psi),
                       -p_gamma_chi * cosTH_psi,
                       0. );

 // calculate psi 4-momentum in the chi rest frame, wrt the xyz axes:
 // need to rotate from the "CS" or "HX" system of axes to the "xyz" system of axes

   // calculate reference directions in the chic rest frame:

    TVector3 cm_to_chi = -chi.BoostVector();
    TVector3 chi_to_cm = chi.BoostVector();

    TLorentzVector targ_chi = targ;
    targ_chi.Boost(cm_to_chi);         // target in the chi rest frame
    TLorentzVector beam_chi = beam;
    beam_chi.Boost(cm_to_chi);         // beam in the chi rest frame

    TVector3 beam_direction_chi     = beam_chi.Vect().Unit();
    TVector3 targ_direction_chi     = targ_chi.Vect().Unit();
    TVector3 chi_direction          = chi.Vect().Unit();
    TVector3 beam_targ_bisec_chi    = (beam_direction_chi - targ_direction_chi).Unit();

   // all polarization frames have the same Y axis = the normal to the plane formed by
   // the directions of the colliding hadrons

    TVector3 Yaxis = ( beam_direction_chi.Cross( targ_direction_chi ) ).Unit();

   // transform (rotation) psi momentum components from polarization axis system
   // to the system with x,y,z axes as in the laboratory

    TVector3 ChiPolAxis = beam_targ_bisec_chi;            // definition of the polarization axis: CS frame
    if ( !config.CSframeIsNatural ) ChiPolAxis = chi_direction;  // or helicity frame

    TVector3 oldZaxis = ChiPolAxis;
    TVector3 oldYaxis = Yaxis;
    TVector3 oldXaxis = oldYaxis.Cross(oldZaxis);

    TRotation rotation;
    rotation.RotateAxes(oldXaxis, oldYaxis, oldZaxis);
                     // transforms coordinates from the "old" frame to the "xyz" frame


    psi_chi.Transform(rotation);
                     // psi in the chi rest frame
                     // wrt to the xyz axes

    gamma_chi.Transform(rotation);


    // const auto test_chi = psi_chi + gamma_chi;
    // std::cout << Mchi << " " << test_chi.M() << "\n";


 // boost psi from the chic rest frame into the proton-proton CM frame:

    TLorentzVector psi = psi_chi;
    psi.Boost(chi_to_cm);


 // kinematics of the psi (measured in the p-p CM):

    pT = psi.Perp();
    // pL = psi.Pz();
    y  = psi.Rapidity();


 // boost gamma from the chic rest frame into the proton-proton CM frame:

    TLorentzVector gamma = gamma_chi;
    gamma.Boost(chi_to_cm);


 // kinematics of the gamma (measured in the p-p CM):

    pT_gamma = gamma.Perp();
    pL_gamma = gamma.Pz();
    y_gamma  = gamma.Rapidity();


 // lepton 4-momentum in the psi rest frame, wrt the "natural" polarization axes (z = psi direction in the chic rest frame):

    double p_lepton_psi = sqrt( 0.25*Mpsi*Mpsi - Mlepton*Mlepton );

    TLorentzVector lepton_psi;

    double sinphi_chihe = sin( phi_chihe * PIG/180. );

    lepton_psi.SetXYZM( p_lepton_psi * sinth_chihe * cosphi_chihe,
                        p_lepton_psi * sinth_chihe * sinphi_chihe,
                        p_lepton_psi * costh_chihe,
                        Mlepton );


    // now "psi" is the psi in the proton-proton CM frame;
    // find reference directions in the psi rest frame
    // for lepton decay distribution:


    TVector3 cm_to_psi = -psi.BoostVector();
    TVector3 psi_to_cm = psi.BoostVector();

    TLorentzVector targ_psi = targ;
    targ_psi.Boost(cm_to_psi);         // target in the psi rest frame
    TLorentzVector beam_psi = beam;
    beam_psi.Boost(cm_to_psi);         // beam in the psi rest frame


    TVector3 beam_direction_psi     = beam_psi.Vect().Unit();
    TVector3 targ_direction_psi     = targ_psi.Vect().Unit();
    TVector3 psi_direction          = psi.Vect().Unit();      // psi as seen in the CM system
    TVector3 psi_direction_chi      = psi_chi.Vect().Unit();     // psi as seen in the chi rest frame!
    TVector3 beam_targ_bisect_psi   = ( beam_direction_psi - targ_direction_psi ).Unit();


    // also need chi direction in the psi rest frame: the y axis of the
    // "natural polarization frame" is defined as the axis perpendicular to the directions
    // of psi and of the chic polarization axis


    Yaxis = ( ChiPolAxis.Cross( psi_direction_chi ) ).Unit();


    // transform (rotation) lepton momentum components from generation frame
    // to the xyz coordinate system

    oldZaxis = psi_direction_chi;
    oldYaxis = Yaxis;
    oldXaxis = oldYaxis.Cross(oldZaxis);

    rotation.SetToIdentity();
    rotation.RotateAxes(oldXaxis, oldYaxis, oldZaxis);
                     // transforms coordinates from the "old" frame to the "xyz" frame

    lepton_psi.Transform(rotation);



  // now calculate angles in the observable frames of the inclusive psi detection

  // the y axis is common to CS and helicity:

    TVector3 newYaxis = ( beam_direction_psi.Cross( targ_direction_psi ) ).Unit();


 /////////////////////////////////////////////////////////////////////
 // CS frame: polar axis = bisector of beam and arget directions

    TVector3 newZaxis = beam_targ_bisect_psi;
    TVector3 newXaxis = newYaxis.Cross(newZaxis);

    rotation.SetToIdentity();
    rotation.RotateAxes(newXaxis,newYaxis,newZaxis);
    rotation.Invert(); // transforms coordinates from the "xyz" frame to the new frame

    TVector3 lepton_psi_rotated = lepton_psi.Vect();

    lepton_psi_rotated.Transform(rotation);

    costh_cs = lepton_psi_rotated.CosTheta();

    phi_cs = lepton_psi_rotated.Phi() * 180. / PIG;
    // if ( phi_cs < 0. ) phi_cs = 360. + phi_cs;


 /////////////////////////////////////////////////////////////////////
 // HELICITY frame: polar axis = psi direction in the CM system

    newZaxis = psi_direction;
    newXaxis = newYaxis.Cross(newZaxis);

    rotation.SetToIdentity();
    rotation.RotateAxes(newXaxis,newYaxis,newZaxis);
    rotation.Invert(); // transforms coordinates from the "xyz" frame to the new frame

    lepton_psi_rotated = lepton_psi.Vect();

    lepton_psi_rotated.Transform(rotation);

    costh_he = lepton_psi_rotated.CosTheta();

    phi_he = lepton_psi_rotated.Phi() * 180. / PIG;
    // if ( phi_he < 0. ) phi_he = 360. + phi_he;



  // leptons in the laboratory: using the above-defined
  // TVector3 psi_to_cm = psi.BoostVector();

    TLorentzVector lepP = lepton_psi;
    lepP.Boost(psi_to_cm);
    pT_lepP = lepP.Perp();
    eta_lepP = lepP.PseudoRapidity();

    TLorentzVector lepN;
    lepN.SetPxPyPzE(-lepton_psi.Px(),-lepton_psi.Py(),-lepton_psi.Pz(),lepton_psi.E());
    lepN.Boost(psi_to_cm);
    pT_lepN = lepN.Perp();
    eta_lepN = lepN.PseudoRapidity();


  // accepted events:

    // inAcc0 = pT > pT_psi_min && pT < pT_psi_max && fabs(psi.Rapidity()) < 1.0 &&   // some basic acceptance cuts
    //          pT_lepP > 4.0  &&  fabs(eta_lepP) < 1.4  &&
    //          pT_lepN > 4.0  &&  fabs(eta_lepN) < 1.4 ;

    // inAcc1 = pT_gamma > 1.0;  // some further cuts

  //   // obtaining smeared four-momenta for the muons and photons and calculating the smeared variables for the chi and j/psi from there
    const auto smearedLepP = smearParticleGaus(lepP, 0, 0.03);
    const auto smearedLepN = smearParticleGaus(lepN, 0, 0.03);
    const auto smearedGamma = smearParticleTF1(gamma, photonCrystalBall);

    const auto smearedJpsi = smearedLepP + smearedLepN;
    const auto smearedChi = smearedJpsi + smearedGamma;

    // apply the selectors (as soon as possible in this case)
    if (!(jpsiSelector->accept(smearedJpsi) &&
         muonSelector->accept(smearedLepP) && muonSelector->accept(smearedLepN) &&
         photonSelector->accept(smearedGamma))) {
      continue;
    }

    // const TLorentzVector halfSmearedChi = smearedJpsi + gamma;

    // const TLorentzVector fullSmearedChi = smearedLepP + smearedLepN + smearedGamma;

    pT_chi_sm = smearedChi.Pt();
    y_chi_sm = smearedChi.Rapidity();
    M_chi_sm = smearedChi.M();

    pT_jpsi_sm = smearedJpsi.Pt();
    y_jpsi_sm = smearedJpsi.Rapidity();
    M_jpsi_sm = smearedJpsi.M();

    qM_chi_sm = M_chi_sm - M_jpsi_sm + MpsiPDG;

    pT_gamma_sm = smearedGamma.Pt();
    y_gamma_sm = smearedGamma.Rapidity();

    pT_lepP_sm = smearedLepP.Pt();
    eta_lepP_sm = smearedLepP.Eta();
    pT_lepN_sm = smearedLepN.Pt();
    eta_lepN_sm = smearedLepN.Eta();

    //  filling of the ntuple:

    // qM_chi = chi.M() - psi.M() + MpsiPDG;
    // M_gamma = gamma.M();

    // ca_gamma_jpsi = TMath::Cos(gamma.Angle(psi.Vect()));
    // ca_sm_gamma_jpsi = TMath::Cos(smearedGamma.Angle(smearedJpsi.Vect()));

    // ca_mu_mu = TMath::Cos(lepP.Angle(lepN.Vect()));
    // ca_sm_mu_mu = TMath::Cos(smearedLepP.Angle(smearedLepN.Vect()));

    const auto angles_HX = calcAnglesInFrame(smearedLepN, smearedLepP, RefFrame::HX);
    costh_HX_sm = angles_HX.costh;
    phi_HX_sm = angles_HX.phi;

    const auto angles_CS = calcAnglesInFrame(smearedLepN, smearedLepP, RefFrame::CS);
    costh_CS_sm = angles_CS.costh;
    phi_CS_sm = angles_CS.phi;

    // add the desired efficiencies
    if (muonEffs) {
      // lepP_eff = muonEffs->Eval(pT_lepP, eta_lepP);
      // lepN_eff = muonEffs->Eval(pT_lepN, eta_lepN);
      lepP_eff_sm = muonEffs->Eval(pT_lepP_sm, eta_lepP_sm);
      lepN_eff_sm = muonEffs->Eval(pT_lepN_sm, eta_lepN_sm);
    }
    if (photonEffs) {
      // gamma_eff = photonEffs->Eval(pT_gamma, y_gamma);
      gamma_eff_sm = photonEffs->Eval(pT_gamma_sm, y_gamma_sm);
    }

    tr->Fill();
    accepted++;

    if (check_accept && accepted >= config.n_accepted) {
      break;
    }
  } // end of external loop (generated events)

  std::cout << '\n' << '\n';

  std::cout << "accepted " << accepted <<" of " << i_event << " generated events" << std::endl;

  hfile->Write("", TObject::kWriteDelete);

} // end of main
