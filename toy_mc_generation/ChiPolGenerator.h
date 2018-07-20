#ifndef CHI_POL_GEN_JN_H
#define CHI_POL_GEN_JN_H

#include <utility>
#include <string>

#include "TF1.h"
#include "TRotation.h"
#include "TVector3.h"
#include "TLorentzVector.h"

#include "../general/interface/calcAngles_defs.h"

struct Mass{
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

    TLorentzVector* smeared_chi(TLorentzVector* chi);
    TLorentzVector* smeared_dimuon(TLorentzVector* psi);
    TLorentzVector* smeared_gamma(TLorentzVector* gamma);

    TRotation get_rotation(const TLorentzVector &chi, bool for_lepton = false);

    // pseudo-globals
    double tmp_dimuon_M = 0;
    double tmp_chi_M = 0;
    TLorentzVector tmp_dimuon_in_chiframe;

    // Constants:  [GeV]
    const Mass PdgMassPsi{3.097, 9.29e-5};
    const std::array<Mass, 3> PdgMassChic{
        Mass{3.415, 0.0105}, // chi0
        Mass{3.511, 0.00088}, // chi1
        Mass{3.556, 0.002}, // chi2
        };

    const Mass PdgMassUps1S{9.4603, 5.402e-5};
    const std::array<Mass, 3> PdgMassChib_1P{
        Mass{9.85944, 3e-4}, // chib0
        Mass{9.89278, 3e-6}, // chib1
        Mass{9.91221, 6e-6}, // chib2
        };

    const double mass_proton{0.9382720};
    const double mass_muon{0.10566};

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
    double R{2./3.};                 // fraction of chic helicity 1: 0 < R  < 1
    double R2{0};                   // fraction of chic helicity 2: 0 < R2 < 1

    RefFrame chi_polarization_frame = RefFrame::HX;

    ULong_t rnd_seed = 0; // if 0: "... a TUUID is generated and used to fill the first 8 integers of the seed array ..."
                          //       "... the seed is guaranteed to be unique in space and time ..." (from ROOT Doc TRandom3)

    std::string outfilename;
    Mass chi_mass = PdgMassChic[1];
    Mass dimuon_mass = PdgMassPsi;
    double dimuon_mass_pdg = PdgMassPsi.central;

    double min_pT = 7;
    double max_pT = 30;
    double min_rap = 0; // Absolute Rapidity
    double max_rap = 1.3;
    int chi_state = 1;
    TF1* pT_distr = nullptr;
    TF1* rap_distr = nullptr;
    TF1* photonCrystalBall = nullptr; // smearing function for photons (fit to MC + some tuning)

    ////////////

 public:
    ChiPolarizationGenerator(const std::string &filename):
        outfilename(filename),
        beam1(0, 0, -p_beam, Ebeam),
        beam2(0., 0. , p_beam, Ebeam)
       {
            //TODO: either ups,psi chooser or mass parameter in constructor
            // and either 1,2,3 chooser or chi mass parameter in constructor:
            chi_mass = PdgMassChic[chi_state];
            dimuon_mass = PdgMassPsi;

            rap_distr = new TF1("pT_distr", [](double* x, double* p){ return 1; }, min_rap, max_rap, 0);
            //rap_distr->Npx(5000);
            pT_distr = new TF1("rap_distr", [](double* x, double* p){                
                double chimass = p[0];
                double pT = x[0];

                // const double beta = 3.45;  //  CHECK HERE FUNCTION AND PARAMETER VALUES: USE THOSE OF GLOBAL FIT (considering that this is a pT distribution, not a pT/M distribution)
                const double beta = 3.39924;  // same as in MC generation from Alberto
                // const double gamma = 0.73;
                const double gamma = 0.635858; // same as in MC generation from Alberto        
                
                return pT * pow( 1. + 1./(beta - 2.) * pT/chimass*pT/chimass / gamma, -beta  );
                }, min_pT, max_pT, 1);
            //pT_distr->Npx(5000);
            
            // smearing function for photons (fit to MC + some tuning)
            photonCrystalBall = new TF1("photonCrystalBall", "ROOT::Math::crystalball_pdf(x[0], [2], [3], [1], [0])", -1.5, 1.5);
            photonCrystalBall->FixParameter(0, 0);
            photonCrystalBall->FixParameter(1, 1.7e-2);
            photonCrystalBall->FixParameter(2, 0.82); // alpha
            photonCrystalBall->FixParameter(3, 1.9); // N
        }

    ~ChiPolarizationGenerator(){}

    void generate(ULong64_t n = 1000000);
    void setChiHelicityFraction1(double hel_frac) { R = hel_frac;}
    void setChiHelicityFraction2(double hel_frac) { R2 = hel_frac;}
    void setSeed(ULong_t seed = 0) {rnd_seed=seed;}
    void setChib(int chib_state) { 
        std::cout << "TODO: CHIB WIDTHs ARE NOT VALID!" << std::endl;
        if (chib_state < 0 || chib_state > 2) {
            std::cout << "WARNING: NON VALID CHIB STATE: " << chib_state << '\n'
            << "\tSetting it to 1." << std::endl;
            chib_state = 1;
        }
        chi_state = chib_state;
        chi_mass = PdgMassChib_1P[chi_state];
        dimuon_mass = PdgMassUps1S;
        dimuon_mass_pdg = PdgMassUps1S.central;
    }

};

#endif