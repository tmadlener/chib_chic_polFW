#include "ChiPolarizationGenerator.h"
#include "ArgParser.h"
#include "utils.h"
#include <string>
#include <iostream>

int main(int argc, char *argv[])
{
  const auto parser = ArgParser(argc, argv);

  const auto outfile = parser.getOptionVal< std::string >("--outfile", "test_mc.root");
  const auto nevents = parser.getOptionVal<int>("--nevents", 100000);
  const auto rnd_seed = parser.getOptionVal<ULong_t>("--seed", 0);
  const auto chib_state = parser.getOptionVal<int>("--chibstate", 1);
  const auto helicity1 = parser.getOptionVal<std::string>("--helicity1", "2 / 3");
  const auto helicity2 = parser.getOptionVal<std::string>("--helicity2", "0");
  const auto ptmin = parser.getOptionVal<double>("--ptmin", 0);
  const auto ptmax = parser.getOptionVal<double>("--ptmax", 100);
  const auto rapmin = parser.getOptionVal<double>("--absrapmin", 0);
  const auto rapmax = parser.getOptionVal<double>("--absrapmax", 1.5);
  const auto applysel = parser.getOptionVal<bool>("--applyselection", true);

  double hel1 = 0;
  double hel2 = 0;
  std::string delimiter = "/";
  auto del_pos1 = helicity1.find(delimiter);
  auto del_pos2 = helicity2.find(delimiter);
  if (del_pos1 == std::string::npos) hel1 = std::stod(helicity1);
  else hel1 = std::stod(helicity1.substr(0, del_pos1)) / std::stod(helicity1.substr(del_pos1+1, helicity1.size()));
  if (del_pos2 == std::string::npos) hel2 = std::stod(helicity2);
  else hel2 = std::stod(helicity2.substr(0, del_pos2)) / std::stod(helicity2.substr(del_pos2+1, helicity2.size()));

  std::cout << "R1: " << hel1 << ", R2: " << hel2 << std::endl;
  std::cout << (applysel ? "Applying selections" : "Selections NOT applied.") << std::endl;

  ChiPolarizationGenerator chib(outfile);
  chib.setChib(chib_state);
  chib.setSeed(rnd_seed);
  chib.setChiHelicityFractions(hel1, hel2);
  chib.setKinematics(ptmin, ptmax, rapmin, rapmax);
  chib.applySelections(applysel);


  {
    std::stringstream ss;
    ss << "Generation of " << nevents << " chib events";

    stopwatch stop(ss.str());
    chib.generate(nevents);
  }

  std::cout << std::endl;

  return 0;
}