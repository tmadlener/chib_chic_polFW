#include "chipolgen.h"
#include "smearing.h"
#include "efficiencies.h"

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

void ChiPolarizationGenerator::generate() 
{
    if (!muonEffsFile.empty()) {
        std::cout << "Reading muon efficiencies\n";
        // non-parametrized single muon efficiencies
        // muonEffs = new EfficiencyProvider<TGraphAsymmErrors>(config.muonEffs, "muon_eff_pt", RangeFromGraph{});
        // parametrized single muon efficiencies
        muonEffs = new EfficiencyProvider<TF1>(config.muonEffs, "muon_eff_pt", RangeFromFit{});
    }
    if (!photonEffsFile.empty()) {
        std::cout << "Reading photon efficiencies\n";
        // using a fixed range for the photon efficiencies since they are parametrized in the region from 400 MeV to 7 GeV
        // photonEffs = new EfficiencyProvider<TF1>(config.photonEffs, "photon_eff_pt", FixedRange<TF1>{0.4, 7});
        photonEffs = new EfficiencyProvider<TF1>(config.photonEffs, "photon_eff_pt", RangeFromFit{});
    }

    pT_distr = new TF1("pT_distr",func_pT_gen,pT_min,pT_max,0);
    rap_distr = new TF1("rap_distr",func_rap_gen,y_min,y_max,0);

    setup_file();
    // reset gRandom
    // loop over events and fill tree variables
}

void ChiPolarizationGenerator::fill_random()
{
    // generate chi and dimuon mass
    chi_mass = gRandom->Gaus( PdgMassChi[chi_state].central, PdgMassChi[chi_state].width);
    dimuon_mass = gRandom->Gaus( PdgMassDimuon.central, PdgMassDimuon.width );

    // generate polarization
    auto pol = random_polarization();
    costh_hx = pol.first;
    phi_hx = pol.second;
    // translate to cs frame

    // generate chi pT and rapidity distributions




}

void ChiPolarizationGenerator::setup_file()
{
  TFile* hfile = new TFile( config.genfile.c_str(), "RECREATE", "chicpolgen");

  TTree* tr = new TTree("tr", "tr");

  tr->Branch( "gen_chicPt",         &brn.pT_chi);
  tr->Branch( "gen_JpsiPt",             &brn.pT);


  tr->Branch( "gen_chicRap",          &brn.y_chi);
  tr->Branch( "gen_JpsiRap",              &brn.y);

  tr->Branch("gen_chicMass", &brn.Mchi);
  tr->Branch("gen_JpsiMass", &brn.Mpsi);

  tr->Branch( "gen_photonPt",       &brn.pT_gamma);
  tr->Branch( "gen_photonPl",       &brn.pL_gamma);
  tr->Branch( "gen_photonEta",        &brn.y_gamma);

  tr->Branch( "gen_muPPt",        &brn.pT_lepP);
  tr->Branch( "gen_muPEta",       &brn.eta_lepP);

  tr->Branch( "gen_muNPt",        &brn.pT_lepN);
  tr->Branch( "gen_muNEta",       &brn.eta_lepN);

 // angle of psi direction in chic rest frame, wrt to chosen chic polarization axis
  tr->Branch( "cosTH_psi",      &brn.cosTH_psi,      "cosTH_psi/D" );

 // angles of dilepton direction in the psi rest frame, wrt the psi direction in the chic rest frame
 // (axis definitions as in Fig 1b of PRD 83, 096001 (2011))
  tr->Branch( "costh_chihe",    &brn.costh_chihe,    "costh_chihe/D" );
  tr->Branch( "phi_chihe",      &brn.phi_chihe,      "phi_chihe/D" );

 // psi decay angles in the helicity frame
  tr->Branch( "gen_costh_HX",       &brn.costh_he);
  tr->Branch( "gen_phi_HX",         &brn.phi_he);

 // psi decay angles in the CS frame
  tr->Branch( "gen_costh_CS",       &brn.costh_cs);
  tr->Branch( "gen_phi_CS",         &brn.phi_cs);

  // smeared variables with "_sm" postfix
  tr->Branch("chicPt", &brn.pT_chi_sm);
  tr->Branch("chicRap", &brn.y_chi_sm);
  tr->Branch("mumugammaMass", &brn.M_chi_sm);
  tr->Branch("chicMass", &brn.qM_chi_sm);

  tr->Branch("photonPt", &brn.pT_gamma_sm);
  tr->Branch("photonEta", &brn.y_gamma_sm);
  tr->Branch("eta_gamma_sm", &brn.eta_gamma_sm);

  tr->Branch("JpsiPt", &brn.pT_jpsi_sm);
  tr->Branch("JpsiRap", &brn.y_jpsi_sm);
  tr->Branch("JpsiMass", &brn.M_jpsi_sm);

  tr->Branch("muPPt", &brn.pT_lepP_sm);
  tr->Branch("muPEta", &brn.eta_lepP_sm);

  tr->Branch("muNPt", &brn.pT_lepN_sm);
  tr->Branch("muNEta", &brn.eta_lepN_sm);

  tr->Branch("Q_value_gen", &brn.Mchic);

  tr->Branch("costh_HX", &brn.costh_HX_sm);
  tr->Branch("phi_HX", &brn.phi_HX_sm);

  tr->Branch("costh_CS", &brn.costh_CS_sm);
  tr->Branch("phi_CS", &brn.phi_CS_sm);

  if (muonEffs) {
    tr->Branch("lepP_eff", &brn.lepP_eff);
    tr->Branch("lepN_eff", &brn.lepN_eff);
    tr->Branch("lepP_eff_sm", &brn.lepP_eff_sm);
    tr->Branch("lepN_eff_sm", &brn.lepN_eff_sm);
  }

  if (photonEffs) {
    tr->Branch("gamma_eff", &brn.gamma_eff);
    tr->Branch("gamma_eff_sm", &brn.gamma_eff_sm);
  }

}


int main(int argc, char const *argv[])
{
    ChiPolarizationGenerator chib_polgen;
    chib_polgen.generate();
    
    return 0;
}
std::pair<double, double> ChiPolarizationGenerateor::random_polarization()
{
    
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
              // direction of the dimuon in the CHI rest frame (wrt to a reference frame, HE or CS, chosen afterwards)
              // PHI_dimuon is the second coordinate, generated outside the loop
              // because the global angular decay distribution does not depend on it.

         costh_chihe = -1. + 2. * gRandom->Rndm();  // direction of the lepton in the dimuon rest frame
         phi_chihe   = 360. * gRandom->Rndm();      // (wrt the dimuon direction seen from the CHI rest frame)


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
                         + ( k10 + k11 * cosTk1 + k2 * cosTH2_psi + k3 * cosTH4_psi + ( k4 + k5 * cosTH2_psi + k6 * cosTH4_psi ) * costh2_chihe
                         + ( k7 + k8 * cosTH2_psi + k9 * cosTH4_psi ) * sinth2_chihe * cos2phi_chihe
                         + ( k10 + k11 * cosTH2_psi )* sin2TH_psi * sin2th_chihe * cosphi_chihe ;

           angdistr *= 15. / ( 64.* PIG*PIG );
         }

         if (angdistr > angdistr_max) { std::cout << "PASSED LIMIT" << std::endl; }

         angdistr_rnd = angdistr_max * gRandom->Rndm();

    } while ( angdistr_rnd > angdistr );

    return {sinth_ch}
}