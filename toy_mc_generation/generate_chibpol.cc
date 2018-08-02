#include "ChiPolarizationGenerator.h"
#include "ArgParser.h"
#include "utils.h"
#include <string>

int main(int argc, char *argv[])
{
  const auto parser = ArgParser(argc, argv);

  const auto outfile = parser.getOptionVal< std::string >("--outfile", "test_mc.root");
  const auto nevents = parser.getOptionVal<int>("--nevents", 100000);
  const auto rnd_seed = parser.getOptionVal<ULong_t>("--seed", 0);
  const auto chib_state = parser.getOptionVal<int>("--chibstate", 1);
  const auto helicity1 = parser.getOptionVal<double>("--helicity1", 2. / 3.);
  const auto helicity2 = parser.getOptionVal<double>("--helicity2", 0);
  const auto ptmin = parser.getOptionVal<double>("--ptmin", 0);
  const auto ptmax = parser.getOptionVal<double>("--ptmax", 100);
  const auto rapmin = parser.getOptionVal<double>("--absrapmin", 0);
  const auto rapmax = parser.getOptionVal<double>("--absrapmax", 1.5);

  ChiPolarizationGenerator chib(outfile);
  chib.setChib(chib_state);
  chib.setSeed(rnd_seed);
  chib.setChiHelicityFractions(helicity1, helicity2);
  chib.setKinematics(ptmin, ptmax, rapmin, rapmax);

  {
    std::stringstream ss;
    ss << "Generation of " << nevents << " chib events";

    stopwatch stop(ss.str());
    chib.generate(nevents);
  }

  std::cout << std::endl;

  return 0;
}