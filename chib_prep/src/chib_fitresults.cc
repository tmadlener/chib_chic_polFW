#include "FitAnalyser.h"
#include "ChiOrganizer.h"
#include "ArgParser.h"
#include "Fitter.h"

int main(int argc, char **argv) {

  ArgParser parser(argc, argv);

  // Mandatory ... in the future
  auto binvarname = parser.getOptionVal < std::string>("--binvar","dimuon_pt");
  auto bin_min = parser.getOptionVal < double>("--binmin",8);
  auto bin_max = parser.getOptionVal < double>("--binmax",50);

  // Optional
  auto p_model_folder = parser.getOptionVal < std::string>("--folder", ""); // If empty the current workingdir is used
  auto p_configfile = parser.getOptionVal < std::string>("--config", "");// if empty it is looked for config.json in the model folder
  auto p_outfile = parser.getOptionVal < std::string>("--outfile", "");// mainly for testing purposes

  std::cout << "FitResults for " << bin_min << " < " << binvarname << " < " << bin_max << "." << std::endl;
  
  std::string infilename = "/afs/hephy.at/work/j/jnecker/data/chib_results/bin_8_50/fitresults-dimuon_mass_8p7_11p1-dimuon_pt_8_50.root";

  FitAnalyser a(infilename, "model", "chi_mass_rf1S", "ws-dimuon_mass_8p7_11p1-dimuon_pt_8_50_chib1P1S");
  a.CustomPlot("chi_fitresults.pdf");

  FitAnalyser b(infilename, "model", "dimuon_mass", "ws-dimuon_mass_8p7_11p1-dimuon_pt_8_50");
  b.CustomPlot("dimuon_fitresults.pdf");

  return 0;
}