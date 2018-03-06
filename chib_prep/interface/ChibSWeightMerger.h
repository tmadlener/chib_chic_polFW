#include "TreeMerger.h"
#include "FitAnalyser.h"

#include "TTree.h"
#include "TFile.h"
#include "TLorentzVector.h"

#include <iostream>

class ChibSWeightMerger : public TreeMerger
{
public:
  using TreeMerger::TreeMerger;

private:
  virtual ChibSWeightMerger* clone() const override { return new ChibSWeightMerger(*this); }
  virtual void setup_branches() override
  {
    // Choose here which branch from which tree
    m_in_tree->SetBranchAddress("dimuon_mass", &dimuon_mass);
    m_in_tree->SetBranchAddress("dimuon_pt", &dimuon_pt);
    m_in_tree->SetBranchAddress("chi_mass_rf1S", &chi_mass_rf1S);
    m_in_tree->SetBranchAddress("chi_pt_rf1S", &chi_pt_rf1S);
    m_in_tree->SetBranchAddress("cosTh_HX", &cosTh_HX);
    m_in_tree->SetBranchAddress("phi_HX", &phi_HX);
    m_in_tree->SetBranchAddress("cosTh_PX", &cosTh_PX);
    m_in_tree->SetBranchAddress("phi_PX", &phi_PX);
    m_in_tree->SetBranchAddress("cosTh_CS", &cosTh_CS);
    m_in_tree->SetBranchAddress("phi_CS", &phi_CS);
    m_in_tree->SetBranchAddress("cosAlpha", &cosAlpha);

    m_secondary_in_tree->SetBranchAddress("N_chib1_sw", &N_chib1_sw);
    m_secondary_in_tree->SetBranchAddress("N_chib2_sw", &N_chib2_sw);
    m_secondary_in_tree->SetBranchAddress("N_bkg_sw", &N_bkg_sw);
    m_secondary_in_tree->SetBranchAddress("L_N_chib1", &L_N_chib1);
    m_secondary_in_tree->SetBranchAddress("L_N_chib2", &L_N_chib2);
    m_secondary_in_tree->SetBranchAddress("L_N_bkg", &L_N_bkg);

    // Set branches for out tree
    m_out_tree->Branch("dimuon_mass", &dimuon_mass);
    m_out_tree->Branch("dimuon_pt", &dimuon_pt);
    m_out_tree->Branch("chi_mass_rf1S", &chi_mass_rf1S);
    m_out_tree->Branch("chi_pt_rf1S", &chi_pt_rf1S);
    m_out_tree->Branch("cosTh_HX", &cosTh_HX);
    m_out_tree->Branch("phi_HX", &phi_HX);
    m_out_tree->Branch("cosTh_PX", &cosTh_PX);
    m_out_tree->Branch("phi_PX", &phi_PX);
    m_out_tree->Branch("cosTh_CS", &cosTh_CS);
    m_out_tree->Branch("phi_CS", &phi_CS);
    m_out_tree->Branch("cosAlpha", &cosAlpha);
    m_out_tree->Branch("N_chib1_sw", &N_chib1_sw);
    m_out_tree->Branch("N_chib2_sw", &N_chib2_sw);
    m_out_tree->Branch("N_bkg_sw", &N_bkg_sw);
    m_out_tree->Branch("L_N_chib1", &L_N_chib1);
    m_out_tree->Branch("L_N_chib2", &L_N_chib2);
    m_out_tree->Branch("L_N_bkg", &L_N_bkg);
  }

  // Branches to merge
  Double_t dimuon_mass;
  Double_t dimuon_pt;

  Double_t chi_mass_rf1S;
  Double_t chi_pt_rf1S;

  Double_t cosTh_HX;
  Double_t phi_HX;
  Double_t cosTh_PX;
  Double_t phi_PX;
  Double_t cosTh_CS;
  Double_t phi_CS;
  Double_t cosAlpha;

  Double_t N_chib1_sw;
  Double_t N_chib2_sw;
  Double_t N_bkg_sw;
  Double_t L_N_chib1;
  Double_t L_N_chib2;
  Double_t L_N_bkg;


};