#ifndef CHIBCHICPOLFW_CHICPREP_CHICTUPLINGINFOMC_H__
#define CHIBCHICPOLFW_CHICPREP_CHICTUPLINGINFOMC_H__

class TTree;
class TLorentzVector;

struct MCDefaultNames {
  static constexpr auto jpsi = "rf1S_dimuon_p4";
  static constexpr auto chic = "rf1S_chi_p4";
  static constexpr auto lepP = "rf1S_muonP_p4";
  static constexpr auto lepN = "rf1S_muonN_p4";
  static constexpr auto Jpsict = "ctpv";
};


struct MCAddInfoIn {
  void Init(TTree* t)
  {
    t->SetBranchAddress("gen_jpsi_p4", &gen_jpsi);
    t->SetBranchAddress("gen_chic_p4", &gen_chic);
    t->SetBranchAddress("gen_muonP_p4", &gen_muP);
    t->SetBranchAddress("gen_muonM_p4", &gen_muN);
    t->SetBranchAddress("trigger", &trigger);
    t->SetBranchAddress("probFit1S", &vtxProb);
    t->SetBranchAddress("gen_photon_p4", &gen_photon);
    t->SetBranchAddress("photon_p4", &photon); // using non-refitted for the moment
  }

  TLorentzVector* gen_jpsi{nullptr};
  TLorentzVector* gen_chic{nullptr};
  TLorentzVector* gen_muP{nullptr};
  TLorentzVector* gen_muN{nullptr};

  TLorentzVector* gen_photon{nullptr};
  TLorentzVector* photon{nullptr};

  int trigger;
  double vtxProb;
};


struct MCAddInfoOut {
  void Create(TTree* t)
  {
    t->Branch("gen_costh_HX", &costh_HX);
    t->Branch("gen_phi_HX", &phi_HX);
    t->Branch("gen_cosalpha_HX", &cosalpha_HX);

    t->Branch("gen_costh_PX", &costh_PX);
    t->Branch("gen_phi_PX", &phi_PX);
    t->Branch("gen_cosalpha_PX", &cosalpha_PX);

    t->Branch("gen_costh_CS", &costh_CS);
    t->Branch("gen_phi_CS", &phi_CS);
    t->Branch("gen_cosalpha_CS", &cosalpha_CS);

    t->Branch("muP_pt", &muP_pt);
    t->Branch("muN_pt", &muN_pt);
    t->Branch("muP_eta", &muP_eta);
    t->Branch("muN_eta", &muN_eta);

    t->Branch("trigger", &trigger);
    t->Branch("vtxProb", &vtxProb);

    t->Branch("JpsiPt", &jpsiPt);
    t->Branch("JpsiRap", &jpsiRap);

    t->Branch("photonPt", &photonPt);
    t->Branch("photonRap", &photonRap);
    t->Branch("photonEta", &photonEta);

    t->Branch("gen_photonPt", &gen_photonPt);
    t->Branch("gen_photonRap", &gen_photonRap);
    t->Branch("gen_photonEta", &gen_photonEta);

  }


  double costh_HX;
  double phi_HX;
  double cosalpha_HX;

  double costh_PX;
  double phi_PX;
  double cosalpha_PX;

  double costh_CS;
  double phi_CS;
  double cosalpha_CS;

  double muP_pt;
  double muN_pt;
  double muP_eta;
  double muN_eta;

  int trigger;
  double vtxProb;

  double jpsiPt;
  double jpsiRap;

  double photonPt;
  double photonRap;
  double photonEta;

  double gen_photonPt;
  double gen_photonRap;
  double gen_photonEta;
};


struct GenMCDefaultNames {
  static constexpr auto jpsi = "gen_jpsi_p4";
  static constexpr auto chic = "gen_chic_p4";
  static constexpr auto lepN = "gen_muonM_p4"; // No this is not a type here
  static constexpr auto lepP = "gen_muonP_p4";
  static constexpr auto Jpsict = "No lifetime information for gen MC";
};

struct GenMCAddInfo {
  void Init(TTree *t) {
    t->SetBranchAddress("gen_photon_p4", &gen_photon);
  }

  void Create(TTree *t) {
    t->Branch("muP_pt", &muP_pt);
    t->Branch("muN_pt", &muN_pt);
    t->Branch("muP_eta", &muP_eta);
    t->Branch("muN_eta", &muN_eta);

    t->Branch("JpsiPt", &jpsiPt);
    t->Branch("JpsiRap", &jpsiRap);

    t->Branch("photonPt", &photonPt);
    t->Branch("photonRap", &photonRap);
    t->Branch("photonEta", &photonEta);
  }

  TLorentzVector *gen_photon{nullptr};

  double muP_pt;
  double muN_pt;
  double muP_eta;
  double muN_eta;

  double jpsiPt;
  double jpsiRap;

  double photonPt;
  double photonRap;
  double photonEta;
};

#endif
