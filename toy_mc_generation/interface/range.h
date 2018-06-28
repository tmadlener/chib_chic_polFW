#ifndef CHIBCHICPOLFW_TOYMCGENERATION_RANGE_H__
#define CHIBCHICPOLFW_TOYMCGENERATION_RANGE_H__

#include <regex>

struct Range {
  double min;
  double max;
};

Range getEtaRange(const std::string& graphName)
{
  constexpr auto rgxstr = ".*_eta_([0-9]+)p?([0-9]*)_([0-9]+)p?([0-9]*).*";

  const std::regex rgx(rgxstr);
  std::smatch cm;
  if (std::regex_match(graphName, cm, rgx)) {
    const double min = atof((cm[1].str() + "." + cm[2].str()).c_str());
    const double max = atof((cm[3].str() + "." + cm[4].str()).c_str());
    return Range{min, max};
  }

  return Range{-1, -1};
}

#endif
