#ifndef CHIB_POL_GEN_JN_H
#define CHIB_POL_GEN_JN_H

struct Mass {
  double central;
  double width;
};

struct gen_branches
{
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


    // double M_gamma = 0; tr->Branch("M_gamma", &M_gamma) = 0;
    // double qM_chi = 0; tr->Branch("qM_chi", &qM_chi) = 0;

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

    double costh_HX_sm = 0;
    double phi_HX_sm = 0;

    double costh_CS_sm = 0;
    double phi_CS_sm = 0;

    double lepP_eff = 0;
    double lepN_eff = 0;
    double lepP_eff_sm = 0;
    double lepN_eff_sm = 0;
    double gamma_eff = 0;
    double gamma_eff_sm = 0;
}

class TF1;
class EffProv;

class ChiPolarizationGenerator
{
private:
    void setup_file();
    void fill_random();
    std::pair<double, double> random_polarization(); // first:costh_hx, second:phi_hx
    TF1* pT_distr = nullptr;
    TF1* rap_distr = nullptr;
    gen_branches brn;

// masses
    const Mass PdgMassDimuon{3.097, 9.29e-5};
    const std::array<Mass, 3> PdgMassChi{
        Mass{3.415, 0.0105}, // chi0
        Mass{3.511, 0.00088}, // chi1
        Mass{3.556, 0.002}, // chi2
        };

    const double Mprot{0.9382720};
    const double Mlepton{0.10566}; // muon GeV
    
    // Settings:
    const double PdgMassChi = 0;
    const double PdgWidthChi = 0;
    const double PdgMassDimuon = 0;
    const double PdgWidthDimuon = 0;


    EffProv *muonEffs = nullptr;
    EffProv *photonEffs = nullptr;

    Long64_t n_events{3000000};

    // some example parameter sets:
    // * chic1 unpolarized: R = 2/3, R2 = 0
    // * chic1 with lambdatheta = +1 (maximum positive): R = 0, R2 = 0
    // * chic1 with lambdatheta = -1/3 (maximum negative): R = 1, R2 = 0
    // * chic2 unpolarized: R = 2/5, R2 = 2/5
    // * chic2 with lambdatheta = +1 (maximum positive): R = 0, R2 = 1
    // * chic2 with lambdatheta = -3/5 (maximum negative): R = 0, R2 = 0
    int chi_state{1};               //  0 = chi_c0,  1 = chi_c1,  2 = chi_c2
    double R{2./3.};                 // fraction of chic helicity 1: 0 < R  < 1
    double R2{0};                   // fraction of chic helicity 2: 0 < R2 < 1

    double pbeam{4000.}; // (this is actually irrelevant, as long as pbeam >> proton mass)
    double y_min{0.0}; // min abs rapidity of the chi
    double y_max{1.3}; // max abs rapidity of the chi

    double pT_min{7.0}; // min pt of the chi
    double pT_max{23.0}; // max pt of the chi

    bool CSframeIsNatural{false};   // generate chic polarization in the CS frame (true)
                                    // or in the HX frame (false)

    std::string genfile{"chicpolgen.root"}; // name of the output file
    // To not produce efficiency branches leave the efficiency file names empty
    std::string muonEffsFile{""}; // file name from where the muon efficiencies should be loaded
    std::string photonEffsFile{""}; // file name from where the photon efficiencies should be loaded




public:
    void generate();

};

#endif