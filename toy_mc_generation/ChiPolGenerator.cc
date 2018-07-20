#include "ChiPolGenerator.h"

#include "smearing.h"
#include "efficiencies.h"

#include "Riostream.h"
#include "TSystem.h"
#include "TString.h"
#include "TROOT.h"
#include "TTree.h"
#include "TMath.h"
#include "TRandom3.h"
#include "TFile.h"

#include <iostream>
#include <string>

#include "../general/interface/calcAngles.h"

const double PI = TMath::Pi();
const double PIover180 = PI/180.;
template <typename T> int sgn(T val) {
    return (T(0) < val) - (val < T(0));
}

TLorentzVector* ChiPolarizationGenerator::generate_chi()
{
    // Generation of chi in the CMS of the proton-proton event

    // M:
    // generate the chi mass and the dimuon mass such that they are correlated similar to what they are in data:
    // 1) get the chi mass and the dimuon mass according to the fitted data distributions
    // 2) in data the chi mass is actually the mass from the KVF but we assume that we can treat it as a Q-value
    //    (M_mumugamma - M_mumu + M_pdg_dimuon), so we reverse it to get to the chic mass we actually want to generate
    
    double tmp_chi_mass = gRandom->Gaus( chi_mass.central, chi_mass.width);    
    tmp_dimuon_M = gRandom->Gaus( dimuon_mass.central, dimuon_mass.width);

    tmp_chi_M = tmp_chi_mass - dimuon_mass_pdg + tmp_dimuon_M;

    // pT
    pT_distr->SetParameter(0, tmp_chi_M);
    double chi_pT = pT_distr->GetRandom();

    // pL
    int rap_sign = sgn<>(gRandom->Uniform(-1., 1.));
    double chi_y = rap_distr->GetRandom() * rap_sign;
    double chi_pL = 0.5 * sqrt( tmp_chi_M*tmp_chi_M + chi_pT*chi_pT ) * (exp(chi_y) - exp(-chi_y));

    double chi_phi   = 2. * PI * gRandom->Rndm();

    TLorentzVector* chi = new TLorentzVector();
    chi->SetXYZM( chi_pT * cos(chi_phi) , chi_pT * sin(chi_phi), chi_pL, tmp_chi_M );

#ifdef TESTING
    std::cout << chi_pT << ", " << tmp_chi_mass << ", " << tmp_dimuon_M << ", " << tmp_chi_M << std::endl;
#endif
    return chi;
}

std::pair< Angles, Angles > ChiPolarizationGenerator::generate_dimuon_angles()
{
  // Generation of full angular distribution
    const double angdistr_max = 0.02;
    double angdistr_rnd;

    // initialize this to the max double value, so that a printout appears if an invalid setting for
    // chi_state has been chosen
    // NOTE: this also makes the macro go into an infinite loop
    double angdistr = std::numeric_limits<double>::max();

    double sinTH_psi  = 100.;
    double cosTH_psi  = 100.;
    double PHI_psi = 2. * PI * gRandom->Rndm();
    double sinth_chihe = 100.;
    double cosphi_chihe = 100.;

    double costh_chihe = 100.;
    double phi_chihe = 100.;

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

        cosphi_chihe = cos( phi_chihe * PI/180. );
        double cos2phi_chihe = 2.*cosphi_chihe*cosphi_chihe -1.;

        
        if ( chi_state == 0 ) { // chic_0 angular distribution

            angdistr = 1. + costh2_chihe;
            angdistr *= 3. / ( 64.* PI*PI );

        } else if ( chi_state == 1 ) { // chic_1 angular distribution

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
            angdistr *= 9. / ( 64.* PI*PI );

        } else if ( chi_state == 2 ) { // chic_2 angular distribution

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
            angdistr *= 15. / ( 64.* PI*PI );
        }

         if (angdistr > angdistr_max) { std::cout << "PASSED LIMIT" << std::endl; }

         angdistr_rnd = angdistr_max * gRandom->Rndm();

    } while ( angdistr_rnd > angdistr );

    return {{cosTH_psi, PHI_psi,0},{costh_chihe, phi_chihe,0}};
}

TLorentzVector* ChiPolarizationGenerator::generate_dimuon(TRotation *rot2xyz, TVector3* boost2cm, double cosTh, double phi)
{
    // Dimuon 4-momentum in the chi rest frame, wrt the chosen chi polarization axes
    double p_dimuon_chi = 0.5 * ( tmp_chi_M*tmp_chi_M - tmp_dimuon_M*tmp_dimuon_M ) / tmp_chi_M;
    double sinTh = sqrt(1 - cosTh*cosTh);
    
    tmp_dimuon_in_chiframe.SetXYZM( p_dimuon_chi * sinTh * cos(phi),
                     p_dimuon_chi * sinTh * sin(phi),
                     p_dimuon_chi * cosTh,
                     tmp_dimuon_M );
                     
    // Dimuon 4-momentum in the chi rest frame, wrt the xyz axes:  
    tmp_dimuon_in_chiframe.Transform(*rot2xyz);
    
    // Dimuon in the proton-proton CM frame:
    TLorentzVector *dimuon = new TLorentzVector(tmp_dimuon_in_chiframe);
    dimuon->Boost(*boost2cm);

    return dimuon;
}

std::pair< TLorentzVector*, TLorentzVector*> ChiPolarizationGenerator::generate_leptons(TRotation *rot2xyz, TVector3* boost2cm, double cosTh, double phi)
{
    double p_lepton_dimuon = sqrt( 0.25 * tmp_dimuon_M * tmp_dimuon_M - mass_muon * mass_muon );
    TLorentzVector lepton_dimuon;
    double sinth = sqrt(1-cosTh*cosTh);
    double phi_r = phi * PIover180;

    lepton_dimuon.SetXYZM( p_lepton_dimuon * sinth * cos(phi_r),
                        p_lepton_dimuon * sinth * sin(phi_r),
                        p_lepton_dimuon * cosTh,
                        mass_muon );

    // transforms coordinates from the "old" frame to the "xyz" frame
    lepton_dimuon.Transform(*rot2xyz);

    // leptons in the laboratory frame
    TLorentzVector* lepP = new TLorentzVector(lepton_dimuon);
    lepP->Boost(*boost2cm);   
    TLorentzVector* lepN = new TLorentzVector(lepton_dimuon);
    lepN->SetPxPyPzE(-lepton_dimuon.Px(), -lepton_dimuon.Py(), -lepton_dimuon.Pz(), lepton_dimuon.E());
    lepN->Boost(*boost2cm);

    return {lepP, lepN};
}

TLorentzVector* ChiPolarizationGenerator::generate_gamma(TRotation *rot2xyz, TVector3* boost2cm, double cosTh, double phi)
{ 
    // gamma 4-momentum in the chi rest frame, wrt the chosen chi polarization axes:
    double p_gamma_chi = -0.5 * ( tmp_chi_M * tmp_chi_M - tmp_dimuon_M * tmp_dimuon_M ) / tmp_chi_M;
    double sinTh = sqrt(1 - cosTh*cosTh);

    TLorentzVector gamma_chi;
    gamma_chi.SetXYZM( p_gamma_chi * sinTh * cos(phi),
                       p_gamma_chi * sinTh * sin(phi),
                       p_gamma_chi * cosTh,
                       0. );


    gamma_chi.Transform(*rot2xyz);

 // boost gamma from the chic rest frame into the proton-proton CM frame:
    TLorentzVector *gamma = new TLorentzVector(gamma_chi);
    gamma->Boost(*boost2cm);
    
    return gamma;
}


TLorentzVector* ChiPolarizationGenerator::smeared_chi(TLorentzVector* chi)
{
}

TLorentzVector* ChiPolarizationGenerator::smeared_dimuon(TLorentzVector* psi)
{
}

TLorentzVector* ChiPolarizationGenerator::smeared_gamma(TLorentzVector* gamma)
{
}


TRotation ChiPolarizationGenerator::get_rotation(const TLorentzVector &vec, bool for_lepton)
{    
    //vec.Print();
    // need to rotate from the "CS" or "HX" system of axes to the "xyz" system of axes

   // calculate reference directions in the chic rest frame:

    TVector3 cm_to_chi = -vec.BoostVector();

    TLorentzVector targ_chi = beam1;
    targ_chi.Boost(cm_to_chi);         // target in the chi rest frame
    TLorentzVector beam_chi = beam2;
    beam_chi.Boost(cm_to_chi);         // beam in the chi rest frame

    TVector3 beam_direction_chi     = beam_chi.Vect().Unit();
    TVector3 targ_direction_chi     = targ_chi.Vect().Unit();
    TVector3 chi_direction          = vec.Vect().Unit();
    TVector3 psi_direction_chi      = tmp_dimuon_in_chiframe.Vect().Unit(); //LEPTON specific  // psi as seen in the chi rest frame!
    TVector3 beam_targ_bisec_chi    = (beam_direction_chi - targ_direction_chi).Unit();

   // all polarization frames have the same Y axis = the normal to the plane formed by
   // the directions of the colliding hadrons  


   // rotation from polarization axis system to the system with x,y,z axes as in the laboratory

    TVector3 ChiPolAxis;    
    switch (chi_polarization_frame)
    {
        case RefFrame::CS:
            ChiPolAxis = beam_targ_bisec_chi;
            break;
        case RefFrame::HX:
        default:
            ChiPolAxis = chi_direction;
    }


    TVector3 Yaxis = ( beam_direction_chi.Cross( targ_direction_chi ) ).Unit();
    if (for_lepton) Yaxis = ( ChiPolAxis.Cross( psi_direction_chi ) ).Unit(); //LEPTON specific 

    TVector3 oldZaxis = ChiPolAxis;
    if (for_lepton) oldZaxis = psi_direction_chi; //LEPTON specific 
    TVector3 oldYaxis = Yaxis;
    TVector3 oldXaxis = oldYaxis.Cross(oldZaxis);

    TRotation rotation;
    rotation.RotateAxes(oldXaxis, oldYaxis, oldZaxis);
    // transforms coordinates from the "old" frame to the "xyz" frame

    return rotation;    
}


void ChiPolarizationGenerator::generate(const ULong64_t n)
{
    TFile f(outfilename.c_str(), "RECREATE");
    auto t = new TTree("toy_mc","toy_mc");
    
    Long64_t event_id = 0; 
    t->Branch("event_id", &event_id);

    TLorentzVector *chi = nullptr;
    t->Branch( "gen_chi", &chi);
    TLorentzVector *dimuon = nullptr;
    t->Branch( "gen_dimuon", &dimuon);
    TLorentzVector *gamma = nullptr;
    t->Branch( "gen_gamma", &gamma);
    double costh_dimuon_in_chi_restframe = 0;
    t->Branch( "costh_dimuon_in_chi_restframe", &costh_dimuon_in_chi_restframe);
    double phi_dimuon_in_chi_restframe = 0;
    t->Branch( "phi_dimuon_in_chi_restframe", &phi_dimuon_in_chi_restframe);
    double costh_lepton_in_dimuon_restframe = 0;
    t->Branch( "costh_lepton_in_dimuon_restframe", &costh_lepton_in_dimuon_restframe);
    double phi_lepton_in_dimuon_restframe = 0;
    t->Branch( "phi_lepton_in_dimuon_restframe", &phi_lepton_in_dimuon_restframe);
    TLorentzVector* muPos = nullptr;
    t->Branch("gen_muPos", &muPos);
    TLorentzVector* muNeg = nullptr;
    t->Branch("gen_muNeg", &muNeg);
    double costh_HX = 0;
    t->Branch("costh_HX", &costh_HX);
    double phi_HX = 0;
    t->Branch("phi_HX", &phi_HX);
    double costh_CS = 0;
    t->Branch("costh_CS", &costh_CS);
    double phi_CS = 0;
    t->Branch("phi_CS", &phi_CS);
    TLorentzVector* smearedLepP = nullptr;
    t->Branch("smearedLepP", &smearedLepP);
    TLorentzVector* smearedLepN = nullptr;
    t->Branch("smearedLepN", &smearedLepN);
    TLorentzVector* smearedGamma = nullptr;
    t->Branch("smearedGamma", &smearedGamma);
    TLorentzVector* smearedJpsi = nullptr;
    t->Branch("smearedJpsi", &smearedJpsi);
    TLorentzVector* smearedChi = nullptr;
    t->Branch("smearedChi", &smearedChi);
    TLorentzVector* halfSmearedChi = nullptr;
    t->Branch("halfSmearedChi", &halfSmearedChi);
    double costh_HX_sm = 0;
    t->Branch("costh_HX_sm", &costh_HX_sm);
    double phi_HX_sm = 0;
    t->Branch("phi_HX_sm", &phi_HX_sm);
    double costh_CS_sm = 0;
    t->Branch("costh_CS_sm", &costh_CS_sm);
    double phi_CS_sm = 0;
    t->Branch("phi_CS_sm", &phi_CS_sm);

    delete gRandom;
    gRandom = new TRandom3(rnd_seed);
    
    std::cout << "seed(actual seed): " << rnd_seed <<"(" <<gRandom->GetSeed() <<"), first Rndm() val: " << gRandom->Rndm() << "\n\n";
    
    event_id = 0;
    
    for (ULong64_t i = 0; i<n; ++i) {
        ++event_id;

        chi = generate_chi();
        //chi->Print();

        auto dimuon_angles_in_chi_restframe = generate_dimuon_angles();
        costh_dimuon_in_chi_restframe = dimuon_angles_in_chi_restframe.first.costh;
        phi_dimuon_in_chi_restframe = dimuon_angles_in_chi_restframe.first.phi;
        costh_lepton_in_dimuon_restframe = dimuon_angles_in_chi_restframe.second.costh;  // direction of the lepton in the PSI rest frame
        phi_lepton_in_dimuon_restframe = dimuon_angles_in_chi_restframe.second.phi;       // (wrt the PSI direction seen from the CHI rest frame)
        
        auto chi_rotation = get_rotation(*chi);
        auto chi_boost = chi->BoostVector();
        dimuon = generate_dimuon(&chi_rotation, &chi_boost, costh_dimuon_in_chi_restframe, phi_dimuon_in_chi_restframe);
        gamma = generate_gamma(&chi_rotation, &chi_boost, costh_dimuon_in_chi_restframe, phi_dimuon_in_chi_restframe);
        
        auto dimuon_boost = dimuon->BoostVector();
        auto dimuon_rotation = get_rotation(*dimuon, true);
        auto leptons = generate_leptons(&dimuon_rotation, &dimuon_boost, costh_lepton_in_dimuon_restframe, phi_lepton_in_dimuon_restframe);
        muPos = leptons.first;
        muNeg = leptons.second;
        auto angles_HX = calcAnglesInFrame(*muNeg, *muPos, RefFrame::HX);
        costh_HX = angles_HX.costh;
        phi_HX = angles_HX.phi;
        auto angles_CS = calcAnglesInFrame(*muNeg, *muPos, RefFrame::CS);
        costh_CS = angles_CS.costh;
        phi_CS = angles_CS.phi;

#ifdef TESTING
        muPos->Print();
        muNeg->Print();
        std::cout << "HX:" << costh_HX <<"/" << phi_HX
                 << " CS:" << costh_CS <<"/" << phi_CS << std::endl;
#endif

        // accepted events:

        // inAcc0 = pT > pT_psi_min && pT < pT_psi_max && fabs(psi.Rapidity()) < 1.0 &&   // some basic acceptance cuts
        //          pT_lepP > 4.0  &&  fabs(eta_lepP) < 1.4  &&
        //          pT_lepN > 4.0  &&  fabs(eta_lepN) < 1.4 ;

        // inAcc1 = pT_gamma > 1.0;  // some further cuts

        // obtaining smeared four-momenta for the muons and photons and calculating the smeared variables for the chi and j/psi from there
        smearedLepP = new TLorentzVector(smearParticleGaus(*muPos, 0, 0.03));
        smearedLepN = new TLorentzVector(smearParticleGaus(*muNeg, 0, 0.03));
        smearedGamma = new TLorentzVector(smearParticleTF1(*gamma, photonCrystalBall));

        smearedJpsi = new TLorentzVector((*smearedLepP) + (*smearedLepN));
        smearedChi = new TLorentzVector((*smearedJpsi) + (*smearedGamma));
        halfSmearedChi = new TLorentzVector((*smearedJpsi) + (*gamma));

        auto angles_sm_HX = calcAnglesInFrame(*smearedLepN, *smearedLepP, RefFrame::HX);
        costh_HX_sm = angles_sm_HX.costh;
        phi_HX_sm = angles_sm_HX.phi;

        auto angles_sm_CS = calcAnglesInFrame(*smearedLepN, *smearedLepP, RefFrame::CS);
        costh_CS_sm = angles_sm_CS.costh;
        phi_CS_sm = angles_sm_CS.phi;

#ifdef THIS_IS_STILL_TODO
        // add the desired efficiencies
        if (muonEffs) {
            lepP_eff = muonEffs->Eval(pT_lepP, eta_lepP);
            lepN_eff = muonEffs->Eval(pT_lepN, eta_lepN);
            lepP_eff_sm = muonEffs->Eval(pT_lepP_sm, eta_lepP_sm);
            lepN_eff_sm = muonEffs->Eval(pT_lepN_sm, eta_lepN_sm);
        }
        if (photonEffs) {
            gamma_eff = photonEffs->Eval(pT_gamma, y_gamma);
            gamma_eff_sm = photonEffs->Eval(pT_gamma_sm, y_gamma_sm);
        }
#endif

        t->Fill();
    }
    t->Write();
    f.Close();

    std::cout << std::flush;
    delete chi;
    delete dimuon;
    delete gamma;
    delete muPos;
    delete muNeg;
    delete smearedLepP;
    delete smearedLepN;
    delete smearedGamma;
    delete smearedJpsi;
    delete smearedChi;
    delete halfSmearedChi;
}
