#include "../general/interface/ArgParser.h"

#include "chicpolgen.C"

#include <string>


#if !(defined(__CINT__) or defined(__CLING__))
int main(int argc, char *argv[])
{
  const auto parser = ArgParser(argc, argv);
  const auto genFileName = parser.getOptionVal<std::string>("--genfile", "chicpolgen.root");
  const auto helicity1 = parser.getOptionVal<double>("--helicity1", 2./3.);
  const auto helicity2 = parser.getOptionVal<double>("--helicity2", 0);
  const auto state = parser.getOptionVal<int>("--state", 1);
  const auto nevents = parser.getOptionVal<int>("--nevents", 3000000);
  // chic kinematics at generation
  const auto ptmin = parser.getOptionVal<double>("--ptmin", 7.0);
  const auto ptmax = parser.getOptionVal<double>("--ptmax", 30.0);
  const auto ymin = parser.getOptionVal<double>("--ymin", 0);
  const auto ymax = parser.getOptionVal<double>("--ymax", 1.3);
  // natural pol frame
  const auto CSframe = parser.getOptionVal<bool>("--CSframe", false);

  const auto naccept = parser.getOptionVal<size_t>("--naccept", 0);

  const auto muonEffFile = parser.getOptionVal<std::string>("--muonEffs", "");
  const auto photonEffFile = parser.getOptionVal<std::string>("--photonEffs", "");

  const auto storeBranches = parser.getOptionVal<std::vector<std::string>>("--storeBranches", {"all"});
  const auto storeHists = parser.getOptionVal<bool>("--storeHists", false);

  // j/psi kinematics at reconstruction
  const auto jpsiSel = parser.getOptionVal<bool>("--jpsiSel", false);
  const auto psiPtMin = parser.getOptionVal<double>("--psiPtMin", 8.0);
  const auto psiPtMax = parser.getOptionVal<double>("--psiPtMax", 20.0);
  const auto psiRapMax = parser.getOptionVal<double>("--psiRapMax", 1.2);
  const auto psiRapMin = parser.getOptionVal<double>("--psiRapMin", 0);

  // muon and photon selection
  const auto muonSel = parser.getOptionVal<bool>("--muonSel", false);
  const auto photonSel = parser.getOptionVal<bool>("--photonSel", false);


  // TODO: smearing configurable in the end
  // TODO: efficiencies configurable from json (including range finding)

  gen_config config;
  config.genfile = genFileName;
  config.R = helicity1;
  config.R2 = helicity2;
  config.chic_state = state;
  config.n_events = nevents;
  config.pT_min = ptmin;
  config.pT_max = ptmax;
  config.y_min = ymin;
  config.y_max = ymax;
  config.CSframeIsNatural = CSframe;

  config.n_accepted = naccept;

  config.muonEffs = muonEffFile;
  config.photonEffs = photonEffFile;

  sel_config sel_conf;
  sel_conf.jpsi_sel = jpsiSel;
  sel_conf.muon_sel = muonSel;
  sel_conf.photon_sel = photonSel;

  sel_conf.psiPtMin = psiPtMin;
  sel_conf.psiPtMax = psiPtMax;
  sel_conf.psiRapMax = psiRapMax;
  sel_conf.psiRapMin = psiRapMin;

  store_config store_conf;
  store_conf.storeBranches = storeBranches;
  store_conf.storeHists = storeHists;

  chicpolgen(config, sel_conf, store_conf);

  return 0;
}
#endif
