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
  const auto state = parser.getOptionVal<int>("--state", 0);
  const auto nevents = parser.getOptionVal<int>("--nevents", 3000000);
  const auto ptmin = parser.getOptionVal<double>("--ptmin", 7.0);
  const auto ptmax = parser.getOptionVal<double>("--ptmax", 21.0);
  const auto ymin = parser.getOptionVal<double>("--ymin", 0);
  const auto ymax = parser.getOptionVal<double>("--ymax", 1.3);
  const auto CSframe = parser.getOptionVal<bool>("--CSframe", false);

  // TODO: smearing configurable in the end

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

  chicpolgen(config);

  return 0;
}
#endif
