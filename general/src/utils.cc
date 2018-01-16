#include "utils.h"

#include <fstream>

// Function implementations

bool file_exists(const std::string &fileName)
{
  std::ifstream infile(fileName);
  return infile.good();
}
