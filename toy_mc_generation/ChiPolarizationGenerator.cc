#include "ChiPolarizationGenerator.h"

#include <iostream>
#include <iomanip>
#include <string>
#include <sstream>

#include "Riostream.h"
#include "TSystem.h"
#include "TString.h"
#include "TROOT.h"
#include "TTree.h"
#include "TMath.h"
#include "TRandom3.h"
#include "TFile.h"

#include "smearing.h"
#include "efficiencies.h"
#include "calcAngles.h"
#include "utils.h"



const double PI = TMath::Pi();
const double PIover180 = PI / 180.;
template <typename T> int sgn(T val) {
  return (T(0) < val) - (val < T(0));
}

TLorentzVector* ChiPolarizationGenerator::generate_chi()
{
  addTime time(total_time_chigen);
  ++total_number_chigen;
  // Generation of chi in the CMS of the proton-proton event

  // M:
  // generate the chi mass and the dimuon mass such that they are correlated similar to what they are in data:
  // 1) get the chi mass and the dimuon mass according to the fitted data distributions
  // 2) in data the chi mass is actually the mass from the KVF but we assume that we can treat it as a Q-value
  //    (M_mumugamma - M_mumu + M_pdg_dimuon), so we reverse it to get to the chic mass we actually want to generate

  double tmp_chi_mass;
  {
    //stopwatch s("random mass");
    tmp_chi_mass = (chi_width_is_zero) ? chi_mass.central : gRandom->Gaus(chi_mass.central, chi_mass.width);
    tmp_dimuon_M = (dimuon_width_is_zero) ? dimuon_mass.central : gRandom->Gaus(dimuon_mass.central, dimuon_mass.width);
  }
  tmp_chi_M = tmp_chi_mass - dimuon_mass_pdg + tmp_dimuon_M;

  //The following line is for compatibility reasons:
  Mchic = tmp_chi_mass;
  //

  // pT
  double chi_pT;
  {
    //stopwatch s("random pT");
    pT_distr->SetParameter(0, tmp_chi_M);
    chi_pT = pT_distr->GetRandom();

  }

  // Phi
  double chi_phi;
  {
    //stopwatch s("random chi");
    chi_phi = 2. * PI * gRandom->Rndm();
  }
  TLorentzVector* chi = new TLorentzVector();

  double rnd_rap;
  {
    //stopwatch s("random rap");
    const int rnd_sign = sgn<>(gRandom->Uniform(-1., 1.));
    rnd_rap = rap_distr->GetRandom() * rnd_sign; // this value is also used as random eta
  }
#if GENRAPIDITY == 1
  // pL
  double chi_pL = 0.5 * sqrt(tmp_chi_M*tmp_chi_M + chi_pT*chi_pT) * (exp(rnd_rap) - exp(-rnd_rap));
  chi->SetXYZM(chi_pT * cos(chi_phi), chi_pT * sin(chi_phi), chi_pL, tmp_chi_M);
#else
  chi->SetPtEtaPhiM(chi_pT, rnd_rap, chi_phi, tmp_chi_M);
#endif

#ifdef TESTING
  std::cout << chi_pT << ", " << tmp_chi_mass << ", " << tmp_dimuon_M << ", " << tmp_chi_M << std::endl;
#endif

  return chi;
}

std::pair< Angles, Angles > ChiPolarizationGenerator::generate_dimuon_angles()
{
  addTime time(total_time_anglesgen);
  ++total_number_anglesgen;

  // Generation of full angular distribution
  const double angdistr_max = 0.02;
  double angdistr_rnd;

  // initialize this to the max double value, so that a printout appears if an invalid setting for
  // chi_state has been chosen
  // NOTE: this also makes the macro go into an infinite loop
  double angdistr = std::numeric_limits<double>::max();

  double sinTH_psi = 100.;
  double cosTH_dimuon = 100.;
  double PHI_psi = 2. * PI * gRandom->Rndm();
  double sinth_chihe = 100.;
  double cosphi_chihe = 100.;

  double costh_chihx = 100.;
  double phi_chihe = 100.;

  do {
    cosTH_dimuon = gRandom->Uniform(-1, 1);// -1. + 2. * gRandom->Rndm();
    // direction of the PSI in the CHI rest frame (wrt to a reference frame, HE or CS, chosen afterwards)
    // PHI_psi is the second coordinate, generated outside the loop
    // because the global angular decay distribution does not depend on it.

    costh_chihx = gRandom->Uniform(-1, 1); //-1. + 2. * gRandom->Rndm();  // direction of the lepton in the PSI rest frame
    phi_chihe = 360. * gRandom->Rndm();      // (wrt the PSI direction seen from the CHI rest frame)

    double cosTH2_psi = cosTH_dimuon*cosTH_dimuon;
    double cosTH4_psi = cosTH2_psi*cosTH2_psi;
    double costh2_chihe = costh_chihx*costh_chihx;
    double sinth2_chihe = 1 - costh2_chihe;

    sinTH_psi = sqrt(1. - cosTH2_psi);
    sinth_chihe = sqrt(sinth2_chihe);

    double sin2TH_psi = 2.*sinTH_psi*cosTH_dimuon;
    double sin2th_chihe = 2.*sinth_chihe*costh_chihx;

    cosphi_chihe = cos(phi_chihe * PI / 180.);
    double cos2phi_chihe = 2.*cosphi_chihe*cosphi_chihe - 1.;


    if (chi_state == 0) { // chic_0 angular distribution

      angdistr = 1. + costh2_chihe;
      angdistr *= 3. / (64.* PI*PI);

    }
    else if (chi_state == 1) { // chic_1 angular distribution

     //  double a2 = gRandom->Gaus( -0.006, 0.013 );  // from Crystal Ball and E835 measurements
      double a2 = 0.;
      double a1 = sqrt(1. - a2*a2);   // (a1 taken to be positive)
      double A0 = sqrt(1. / 2.) * (a1 + a2);
      double A1 = sqrt(1. / 2.) * (a1 - a2);
      double k1 = A1*A1 + 1. / 2.* R * (A0*A0 - A1*A1);
      double k2 = (1. - 3. / 2.* R) * (A0*A0 - A1*A1);
      double k3 = -A1*A1 + 1. / 2.* R;
      double k4 = 1. - 3. / 2.* R;
      double k5 = 1. / 4.* A1*A0 * (3.* R - 2.);
      angdistr = k1 + k2 * cosTH2_psi + (k3 + k4 * cosTH2_psi) * costh2_chihe
        + k5 * sin2TH_psi * sin2th_chihe * cosphi_chihe;
      angdistr *= 9. / (64.* PI*PI);

    }
    else if (chi_state == 2) { // chic_2 angular distribution

 //   double a3 = gRandom->Gaus( 0.01, 0.04 );  // from E760 and E835 measurements
      double a3 = 0.;
      //   double a2 = gRandom->Gaus( -0.13, 0.04 );  // from Crystal Ball, E760 and E835 measurements
      double a2 = 0.;
      double a1 = sqrt(1. - a2*a2 - a3*a3);  // (a1 taken to be positive)
      double A0 = sqrt(1. / 10.)*a1 + sqrt(1. / 2.)*a2 + sqrt(2. / 5.)*a3;
      double A1 = sqrt(3. / 10.)*a1 + sqrt(1. / 6.)*a2 - sqrt(8. / 15.)*a3;
      double A2 = sqrt(3. / 5.)*a1 - sqrt(1. / 3.)*a2 + sqrt(1. / 15.)*a3;
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
      double k1_0 = 1. / 4.* A0*A0 + 3. / 8.* A2*A2;
      double k1_1 = 1. / 2.* A1*A1 + 1. / 4.* A2*A2;
      double k1_2 = 3. / 8.* A0*A0 + 1. / 2.* A1*A1 + 1. / 16.* A2*A2;
      double k2_0 = -3. / 2.* A0*A0 + 3.* A1*A1 - 3. / 4.* A2*A2;
      double k2_1 = 3. / 2.* A0*A0 - 3. / 2.* A1*A1;
      double k2_2 = -3. / 4.* A0*A0 + 3. / 8.* A2*A2;
      double k3_0 = 9. / 4.* A0*A0 - 3.* A1*A1 + 3. / 8.* A2*A2;
      double k3_1 = -3. / 2.* A0*A0 + 2.* A1*A1 - 1. / 4.* A2*A2;
      double k3_2 = 3. / 8.* A0*A0 - 1. / 2.* A1*A1 + 1. / 16.* A2*A2;
      double k4_0 = 1. / 4.* A0*A0 + 3. / 8.* A2*A2;
      double k4_1 = -1. / 2.* A1*A1 + 1. / 4.* A2*A2;
      double k4_2 = 3. / 8.* A0*A0 - 1. / 2.* A1*A1 + 1. / 16.* A2*A2;
      double k5_0 = -3. / 2.* A0*A0 - 3.*A1*A1 - 3. / 4.* A2*A2;
      double k5_1 = 3. / 2.* A0*A0 + 3. / 2.* A1*A1;
      double k5_2 = -3. / 4.* A0*A0 + 3. / 8.* A2*A2;
      double k6_0 = 9. / 4.* A0*A0 + 3.* A1*A1 + 3. / 8.* A2*A2;
      double k6_1 = -3. / 2.* A0*A0 - 2.*A1*A1 - 1. / 4.* A2*A2;
      double k6_2 = 3. / 8.* A0*A0 + 1. / 2.* A1*A1 + 1. / 16.* A2*A2;
      double k7_0 = -sqrt(6.) / 4.* A0*A2;
      double k7_1 = 0.;
      double k7_2 = sqrt(6.) / 8.* A0*A2;
      double k8_0 = sqrt(6.)* A0*A2;
      double k8_1 = -sqrt(6.) / 2. *A0*A2;
      double k8_2 = 0.;
      double k9_0 = -3.* sqrt(6.) / 4.* A0*A2;
      double k9_1 = sqrt(6.) / 2.* A0*A2;
      double k9_2 = -sqrt(6.) / 8.* A0*A2;
      double k10_0 = sqrt(3.) / 4.* A0*A1 + 3.*sqrt(2.) / 8.* A1*A2;
      double k10_1 = -sqrt(3.) / 4.* A0*A1;
      double k10_2 = sqrt(3.) / 8.* A0*A1 - 3.*sqrt(2.) / 16.* A1*A2;
      double k11_0 = -3.*sqrt(3.) / 4.* A0*A1 - 3.*sqrt(2.) / 8.*A1*A2;
      double k11_1 = sqrt(3.) / 2.* A0*A1 + sqrt(2.) / 4.* A1*A2;
      double k11_2 = -sqrt(3.) / 8.* A0*A1 - sqrt(2.) / 16.* A1*A2;
      double R1 = R; double R0 = 1. - R1 - R2;
      double k1 = R0*k1_0 + R1*k1_1 + R2*k1_2;
      double k2 = R0*k2_0 + R1*k2_1 + R2*k2_2;
      double k3 = R0*k3_0 + R1*k3_1 + R2*k3_2;
      double k4 = R0*k4_0 + R1*k4_1 + R2*k4_2;
      double k5 = R0*k5_0 + R1*k5_1 + R2*k5_2;
      double k6 = R0*k6_0 + R1*k6_1 + R2*k6_2;
      double k7 = R0*k7_0 + R1*k7_1 + R2*k7_2;
      double k8 = R0*k8_0 + R1*k8_1 + R2*k8_2;
      double k9 = R0*k9_0 + R1*k9_1 + R2*k9_2;
      double k10 = R0*k10_0 + R1*k10_1 + R2*k10_2;
      double k11 = R0*k11_0 + R1*k11_1 + R2*k11_2;
      angdistr = k1 + k2 * cosTH2_psi + k3 * cosTH4_psi + (k4 + k5 * cosTH2_psi + k6 * cosTH4_psi) * costh2_chihe
        + (k7 + k8 * cosTH2_psi + k9 * cosTH4_psi) * sinth2_chihe * cos2phi_chihe
        + (k10 + k11 * cosTH2_psi)* sin2TH_psi * sin2th_chihe * cosphi_chihe;
      angdistr *= 15. / (64.* PI*PI);
    }

    if (angdistr > angdistr_max) { std::cout << "PASSED LIMIT" << std::endl; }

    angdistr_rnd = angdistr_max * gRandom->Rndm();

  } while (angdistr_rnd > angdistr);

  return{ {cosTH_dimuon, PHI_psi,0},{costh_chihx, phi_chihe,0} };
}

TLorentzVector* ChiPolarizationGenerator::generate_dimuon(TRotation *rot2xyz, TVector3* boost2cm, double cosTh, double phi)
{
  addTime time(total_time_dimuongen);
  ++total_number_dimuongen;
  // Dimuon 4-momentum in the chi rest frame, wrt the chosen chi polarization axes
  double p_dimuon_chi = 0.5 * (tmp_chi_M*tmp_chi_M - tmp_dimuon_M*tmp_dimuon_M) / tmp_chi_M;
  double sinTh = sqrt(1 - cosTh*cosTh);

  tmp_dimuon_in_chiframe.SetXYZM(p_dimuon_chi * sinTh * cos(phi),
    p_dimuon_chi * sinTh * sin(phi),
    p_dimuon_chi * cosTh,
    tmp_dimuon_M);

  // Dimuon 4-momentum in the chi rest frame, wrt the xyz axes:  
  tmp_dimuon_in_chiframe.Transform(*rot2xyz);

  // Dimuon in the proton-proton CM frame:
  TLorentzVector *dimuon = new TLorentzVector(tmp_dimuon_in_chiframe);
  dimuon->Boost(*boost2cm);

  return dimuon;
}

std::pair< TLorentzVector*, TLorentzVector*> ChiPolarizationGenerator::generate_leptons(TRotation *rot2xyz, TVector3* boost2cm, double cosTh, double phi)
{
  double p_lepton_dimuon = sqrt(0.25 * tmp_dimuon_M * tmp_dimuon_M - mass_muon * mass_muon);
  TLorentzVector lepton_dimuon;
  double sinth = sqrt(1 - cosTh*cosTh);
  double phi_r = phi * PIover180;

  lepton_dimuon.SetXYZM(p_lepton_dimuon * sinth * cos(phi_r),
    p_lepton_dimuon * sinth * sin(phi_r),
    p_lepton_dimuon * cosTh,
    mass_muon);

  // transforms coordinates from the "old" frame to the "xyz" frame
  lepton_dimuon.Transform(*rot2xyz);

  // leptons in the laboratory frame
  TLorentzVector* lepP = new TLorentzVector(lepton_dimuon);
  lepP->Boost(*boost2cm);
  TLorentzVector* lepN = new TLorentzVector(lepton_dimuon);
  lepN->SetPxPyPzE(-lepton_dimuon.Px(), -lepton_dimuon.Py(), -lepton_dimuon.Pz(), lepton_dimuon.E());
  lepN->Boost(*boost2cm);

  return{ lepP, lepN };
}

TLorentzVector* ChiPolarizationGenerator::generate_gamma(TRotation *rot2xyz, TVector3* boost2cm, double cosTh, double phi)
{
  // gamma 4-momentum in the chi rest frame, wrt the chosen chi polarization axes:
  double p_gamma_chi = -0.5 * (tmp_chi_M * tmp_chi_M - tmp_dimuon_M * tmp_dimuon_M) / tmp_chi_M;
  double sinTh = sqrt(1 - cosTh*cosTh);

  TLorentzVector gamma_chi;
  gamma_chi.SetXYZM(p_gamma_chi * sinTh * cos(phi),
    p_gamma_chi * sinTh * sin(phi),
    p_gamma_chi * cosTh,
    0.);


  gamma_chi.Transform(*rot2xyz);

  // boost gamma from the chic rest frame into the proton-proton CM frame:
  TLorentzVector *gamma = new TLorentzVector(gamma_chi);
  gamma->Boost(*boost2cm);

  return gamma;
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

  TVector3 beam_direction_chi = beam_chi.Vect().Unit();
  TVector3 targ_direction_chi = targ_chi.Vect().Unit();
  TVector3 chi_direction = vec.Vect().Unit();
  TVector3 psi_direction_chi = tmp_dimuon_in_chiframe.Vect().Unit(); //LEPTON specific  // psi as seen in the chi rest frame!
  TVector3 beam_targ_bisec_chi = (beam_direction_chi - targ_direction_chi).Unit();

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


  TVector3 Yaxis = (beam_direction_chi.Cross(targ_direction_chi)).Unit();
  if (for_lepton) Yaxis = (ChiPolAxis.Cross(psi_direction_chi)).Unit(); //LEPTON specific 

  TVector3 oldZaxis = ChiPolAxis;
  if (for_lepton) oldZaxis = psi_direction_chi; //LEPTON specific 
  TVector3 oldYaxis = Yaxis;
  TVector3 oldXaxis = oldYaxis.Cross(oldZaxis);

  TRotation rotation;
  rotation.RotateAxes(oldXaxis, oldYaxis, oldZaxis);
  // transforms coordinates from the "old" frame to the "xyz" frame

  return rotation;
}

void ChiPolarizationGenerator::generate(const ULong64_t nevents)
{
  TFile f(outfilename.c_str(), "RECREATE");
  auto t = new TTree("toy_mc", "toy_mc");

  setup_branches(t);

  delete gRandom;
  gRandom = new TRandom3(rnd_seed); // important to use TRandom3


  event_id = 0;

  std::cout << "STARTING the GENERATION of " << nevents << " events...\n";
  std::cout << "SEED[actual seed]: " << rnd_seed << "[" << gRandom->GetSeed() << "], first Rndm() val: " << gRandom->Rndm() << "\n";
  print_settings();

  for (ULong64_t i = 0; i < nevents; ++i) {
    ++event_id;

    chi = generate_chi();
#ifdef TESTING
    chi->Print();
#endif

    auto dimuon_angles_in_chi_restframe = generate_dimuon_angles();
    costh_dimuon_in_chi_restframe = dimuon_angles_in_chi_restframe.first.costh;
    phi_dimuon_in_chi_restframe = dimuon_angles_in_chi_restframe.first.phi;

    // direction of the lepton in the DIMUON rest frame
    // (wrt the PSI direction seen from the CHI rest frame)
    costh_lepton_in_dimuon_restframe = dimuon_angles_in_chi_restframe.second.costh;
    phi_lepton_in_dimuon_restframe = dimuon_angles_in_chi_restframe.second.phi;

    auto chi_rotation = get_rotation(*chi);
    auto chi_boost = chi->BoostVector();
    dimuon = generate_dimuon(&chi_rotation, &chi_boost, costh_dimuon_in_chi_restframe, phi_dimuon_in_chi_restframe);
    gamma = generate_gamma(&chi_rotation, &chi_boost, costh_dimuon_in_chi_restframe, phi_dimuon_in_chi_restframe);

    auto dimuon_boost = dimuon->BoostVector();
    auto dimuon_rotation = get_rotation(*dimuon, true);
    auto leptons = generate_leptons(&dimuon_rotation, &dimuon_boost, costh_lepton_in_dimuon_restframe, phi_lepton_in_dimuon_restframe);
    muPos = leptons.first;
    muNeg = leptons.second;


#ifdef TESTING
    muPos->Print();
    muNeg->Print();
    std::cout << "HX:" << costh_HX << "/" << phi_HX
      << " CS:" << costh_CS << "/" << phi_CS << std::endl;
#endif

    {
      addTime time(total_time_smearing);
      ++total_number_smearing;
      // obtaining smeared four-momenta for the muons and photons and 
      // calculating the smeared variables for the chi and dimuon from there
      smearedLepP = new TLorentzVector(smearParticleGaus(*muPos, 0, 0.03));
      smearedLepN = new TLorentzVector(smearParticleGaus(*muNeg, 0, 0.03));
      smearedGamma = new TLorentzVector(smearParticleTF1(*gamma, photonCrystalBall));

      smearedJpsi = new TLorentzVector((*smearedLepP) + (*smearedLepN));
      smearedChi = new TLorentzVector((*smearedJpsi) + (*smearedGamma));
      halfSmearedChi = new TLorentzVector((*smearedJpsi) + (*gamma));
    }

    // TODO: efficiencies
        //// add the desired efficiencies
        //if (muonEffs) {
        //  lepP_eff = muonEffs->Eval(pT_lepP, eta_lepP);
        //  lepN_eff = muonEffs->Eval(pT_lepN, eta_lepN);
        //  lepP_eff_sm = muonEffs->Eval(pT_lepP_sm, eta_lepP_sm);
        //  lepN_eff_sm = muonEffs->Eval(pT_lepN_sm, eta_lepN_sm);
        //}
        //if (photonEffs) {
        //  gamma_eff = photonEffs->Eval(pT_gamma, y_gamma);
        //  gamma_eff_sm = photonEffs->Eval(pT_gamma_sm, y_gamma_sm);
        //}

    fill_branches();
    t->Fill();
  }
  t->Write();
  f.Close();

  std::cout << "FINISHED GENERATION\n"
    "Written output to file " << outfilename << std::endl;

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

void ChiPolarizationGenerator::fill_branches()
{

  addTime time(total_time_fillbranches);
  ++total_number_fillbranches;

  auto angles_HX = calcAnglesInFrame(*muNeg, *muPos, RefFrame::HX);
  costh_HX = angles_HX.costh;
  phi_HX = angles_HX.phi;
  auto angles_CS = calcAnglesInFrame(*muNeg, *muPos, RefFrame::CS);
  costh_CS = angles_CS.costh;
  phi_CS = angles_CS.phi;

  auto angles_sm_HX = calcAnglesInFrame(*smearedLepN, *smearedLepP, RefFrame::HX);
  costh_HX_sm = angles_sm_HX.costh;
  phi_HX_sm = angles_sm_HX.phi;
  auto angles_sm_CS = calcAnglesInFrame(*smearedLepN, *smearedLepP, RefFrame::CS);
  costh_CS_sm = angles_sm_CS.costh;
  phi_CS_sm = angles_sm_CS.phi;

  // compatibility branches

  pT_chi = chi->Pt();
  pT = dimuon->Pt();
  pL_chi = chi->Pz();

  y_chi = chi->Rapidity();
  y = dimuon->Rapidity();
  Mchi = chi->M();
  Mpsi = dimuon->M();

  pT_gamma = gamma->Pt();
  pL_gamma = gamma->Pz();
  y_gamma = gamma->Rapidity();
  pT_lepP = muPos->Pt();
  eta_lepP = muPos->Eta();
  pT_lepN = muNeg->Pt();
  eta_lepN = muNeg->Eta();

  // angle of psi direction in chic rest frame, wrt to chosen chic polarization axis
  cosTH_psi = costh_dimuon_in_chi_restframe;

  // angles of dilepton direction in the psi rest frame, wrt the psi direction in the chic rest frame
  // (axis definitions as in Fig 1b of PRD 83, 096001 (2011))
  costh_chihe = costh_lepton_in_dimuon_restframe;
  phi_chihe = phi_lepton_in_dimuon_restframe;

  // psi decay angles in the helicity frame
  costh_he = costh_HX;
  phi_he = phi_HX;

  // psi decay angles in the CS frame
  costh_cs = costh_CS;
  phi_cs = phi_CS;

  // smeared variables with "_sm" postfix

  pT_chi_sm = smearedChi->Pt();
  y_chi_sm = smearedChi->Rapidity();
  M_chi_sm = smearedChi->M();

  pT_jpsi_sm = smearedJpsi->Pt();
  y_jpsi_sm = smearedJpsi->Rapidity();
  M_jpsi_sm = smearedJpsi->M();

  qM_chi_sm = M_chi_sm - M_jpsi_sm + dimuon_mass_pdg;

  pT_gamma_sm = smearedGamma->Pt();
  y_gamma_sm = smearedGamma->Rapidity();
  eta_gamma_sm = smearedGamma->Eta();

  pT_lepP_sm = smearedLepP->Pt();
  eta_lepP_sm = smearedLepP->Eta();
  pT_lepN_sm = smearedLepN->Pt();
  eta_lepN_sm = smearedLepN->Eta();

  //Mchic is filled in generate_chi()


  //TODO: Efficiency branches {lepP_eff_sm, lepN_eff_sm, gamma_eff_sm}



}

void ChiPolarizationGenerator::setup_branches(TTree *t)
{
  t->Branch("event_id", &event_id);

#ifndef DONT_WRITE_LORENTZ
  t->Branch("gen_chi", &chi);
  t->Branch("gen_dimuon", &dimuon);
  t->Branch("gen_gamma", &gamma);
  t->Branch("gen_muPos", &muPos);
  t->Branch("gen_muNeg", &muNeg);
#endif

  t->Branch("costh_dimuon_in_chi_restframe", &costh_dimuon_in_chi_restframe);
  t->Branch("phi_dimuon_in_chi_restframe", &phi_dimuon_in_chi_restframe);
  t->Branch("costh_lepton_in_dimuon_restframe", &costh_lepton_in_dimuon_restframe);
  t->Branch("phi_lepton_in_dimuon_restframe", &phi_lepton_in_dimuon_restframe);
  t->Branch("costh_HX", &costh_HX);
  t->Branch("phi_HX", &phi_HX);
  t->Branch("costh_CS", &costh_CS);
  t->Branch("phi_CS", &phi_CS);

#ifndef DONT_WRITE_LORENTZ
  t->Branch("smearedLepP", &smearedLepP);
  t->Branch("smearedLepN", &smearedLepN);
  t->Branch("smearedGamma", &smearedGamma);
  t->Branch("smearedDimuon", &smearedJpsi);
  t->Branch("smearedChi", &smearedChi);
  t->Branch("halfSmearedChi", &halfSmearedChi);
#endif

  t->Branch("costh_HX_sm", &costh_HX_sm);
  t->Branch("phi_HX_sm", &phi_HX_sm);
  t->Branch("costh_CS_sm", &costh_CS_sm);
  t->Branch("phi_CS_sm", &phi_CS_sm);

  // compatibility branches:

  t->Branch("gen_chicPt", &pT_chi);
  t->Branch("gen_JpsiPt", &pT);
  t->Branch("pL_chi", &pL_chi);

  t->Branch("gen_chicRap", &y_chi);
  t->Branch("gen_JpsiRap", &y);

  t->Branch("gen_chicMass", &Mchi);
  t->Branch("gen_JpsiMass", &Mpsi);

  t->Branch("gen_photonPt", &pT_gamma);
  t->Branch("gen_photonPl", &pL_gamma);
  t->Branch("gen_photonEta", &y_gamma);

  t->Branch("gen_muPPt", &pT_lepP);
  t->Branch("gen_muPEta", &eta_lepP);

  t->Branch("gen_muNPt", &pT_lepN);
  t->Branch("gen_muNEta", &eta_lepN);


  // angle of psi direction in chic rest frame, wrt to chosen chic polarization axis
  t->Branch("cosTH_psi", &cosTH_psi, "cosTH_psi/D");

  // angles of dilepton direction in the psi rest frame, wrt the psi direction in the chic rest frame
  // (axis definitions as in Fig 1b of PRD 83, 096001 (2011))
  t->Branch("costh_chihe", &costh_chihe, "costh_chihe/D");
  t->Branch("phi_chihe", &phi_chihe, "phi_chihe/D");

  // psi decay angles in the helicity frame
  t->Branch("gen_costh_HX", &costh_he);
  t->Branch("gen_phi_HX", &phi_he);

  // psi decay angles in the CS frame
  t->Branch("gen_costh_CS", &costh_cs);
  t->Branch("gen_phi_CS", &phi_cs);


  // smeared variables with "_sm" postfix
  t->Branch("chicPt", &pT_chi_sm);
  t->Branch("chicRap", &y_chi_sm);
  t->Branch("mumugammaMass", &M_chi_sm);
  t->Branch("chicMass", &qM_chi_sm);

  t->Branch("photonPt", &pT_gamma_sm);
  t->Branch("photonEta", &y_gamma_sm);
  t->Branch("eta_gamma_sm", &eta_gamma_sm);

  t->Branch("JpsiPt", &pT_jpsi_sm);
  t->Branch("JpsiRap", &y_jpsi_sm);
  t->Branch("JpsiMass", &M_jpsi_sm);

  t->Branch("muPPt", &pT_lepP_sm);
  t->Branch("muPEta", &eta_lepP_sm);

  t->Branch("muNPt", &pT_lepN_sm);
  t->Branch("muNEta", &eta_lepN_sm);

  t->Branch("Q_value_gen", &Mchic);

  //TODO: efficiency branches
}

void ChiPolarizationGenerator::print_settings()
{
  std::stringstream ss;
  ss << std::string(50, '-') << "\n"
    << "SETTINGS:\n"
    << "j(chi) = " << chi_state << "\n";

  if (chi_state > 0) ss << "R1 = " << R << "\n";
  if (chi_state > 1) ss << "R2 = " << R2 << "\n";

  ss << min_pT << " < pT(chi) < " << max_pT << "\n" <<
    min_rap << " < |rap(chi)| < " << max_rap << "\n"
    "mass(chi) = " << chi_mass.central << "\n"
    "width(chi) = " << chi_mass.width << "\n"
    "mass(dimuon) = " << dimuon_mass.central << "\n"
    "width(dimuon) = " << dimuon_mass.width << "\n"
    "pdg_mass(dimuon) = " << dimuon_mass_pdg << "\n"
    "Npx of pT_distr = " << pT_distr->GetNpx() << "\n"
    "output file: " << outfilename << " (will be overwritten)\n"
    << std::string(50, '-');

  std::cout << ss.rdbuf() << std::endl;

}


void ChiPolarizationGenerator::print_performance()
{
  const int b = 25;

  auto total_time = total_time_chigen +
    total_time_dimuongen +
    total_time_anglesgen +
    total_time_smearing +
    total_time_fillbranches;

  auto format = [b, total_time](std::string what, std::chrono::nanoseconds time, ULong64_t count) {
    std::stringstream ss;
    std::stringstream stime;
    auto duration_ms = time.count() / 1.e6;
    auto percentoftotaltime = double(total_time.count()) / time.count() * 100.;
    unsigned long long mins = duration_ms / 1000. / 60.;
    unsigned long long secs = (duration_ms - mins * 1000 * 60) / 1000.;
    unsigned long long ms = (duration_ms - secs * 1000 - mins * 1000 * 60);
    stime << mins << ":" << std::setw(2) << std::setfill('0') << secs << ',' << std::setw(3) << std::setfill('0') << ms;
    ss << std::setw(b) << what << std::setw(b) << stime.str()
      << std::setw(b) << count << std::setw(b) << duration_ms / double(count)
      << std::setw(b) << std::setprecision(2) << percentoftotaltime;
    return ss.str();
  };


  std::cout << "\nPerformance:\n"
    << std::setw(b) << "what" << std::setw(b) << "total time [mm:ss,ms]" << std::setw(b) << "how often" << std::setw(b) << "time per event [ms]" << std::setw(b) << "% of total time\n"
    << std::string(b * 5, '-') << "\n"
    << format("chi", total_time_chigen, total_number_chigen) << "\n"
    << format("dimuon", total_time_dimuongen, total_number_dimuongen) << "\n"
    << format("angles", total_time_anglesgen, total_number_anglesgen) << "\n"
    << format("smearing", total_time_smearing, total_number_smearing) << "\n"
    << format("fillbranches", total_time_fillbranches, total_number_fillbranches) << "\n"
    << std::string(b * 5, '-') << "\n"
    << format("total time", total_time, event_id) << "\n\n\n"
    << std::flush;
}
