#ifndef CHIBCHICPOLFW_CHICPREP_CONFIG_COMMONVAR_H__
#define CHIBCHICPOLFW_CHICPREP_CONFIG_COMMONVAR_H__

#include <array>

namespace config {
  constexpr std::array<double,6> ptBinBorders = {10, 15, 20, 25, 30, 50};
  const double maxAbsRap = 2.4;
}
#endif
