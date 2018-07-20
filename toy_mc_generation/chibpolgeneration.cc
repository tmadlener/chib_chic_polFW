#include "ChiPolGenerator.h"
#include "../general/interface/ArgParser.h"
#include <string>

int main(int argc, char *argv[])
{
    const auto parser = ArgParser(argc, argv);
    
    const auto nevents = parser.getOptionVal<int>("--nevents", 100000);
    const auto rnd_seed = parser.getOptionVal<ULong_t>("--seed", 0);
    const auto helicity1 = parser.getOptionVal<double>("--helicity1", 2./3.);
    const auto helicity2 = parser.getOptionVal<double>("--helicity2", 0);

    ChiPolarizationGenerator chib("test_mc.root");
    chib.setChib(1);
    chib.setSeed(rnd_seed);
    
    chib.setChiHelicityFraction1(helicity1);
    chib.setChiHelicityFraction2(helicity2);

    chib.generate(nevents);

    return 0;
}
