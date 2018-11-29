#include "ChiOrganizer.h"
#include "ArgParser.h"

int main(int argc, char **argv) {

  std::streambuf* backup = std::cout.rdbuf(); // Store original cout buffer
  std::ostringstream tmpostream;  
  std::cout.rdbuf(tmpostream.rdbuf()); // Redirect cout buffer, that only filename is printed

  ArgParser parser(argc, argv);
  auto binvarname = parser.getOptionVal < std::string>("--binvar");
  auto bin_min = parser.getOptionVal < double>("--binmin");
  auto bin_max = parser.getOptionVal < double>("--binmax");
  auto configfile = parser.getOptionVal < std::string>("--config", "");
  auto basedir = parser.getOptionVal < std::string>("--basedir", "");
  auto extension = parser.getOptionVal < std::string>("--extension", ChiOrganizer::sweight_extension);

  std::map<std::string, std::pair<double, double> > bin;
  bin[binvarname] = { bin_min, bin_max };

  if (extension == "NONE") extension = "";

  ChiOrganizer o(basedir, configfile);
  auto result = o.FileName(bin, extension);

  // Restore original cout buffer and print result
  std::cout.rdbuf(backup);
  std::cout << result << std::endl;

}