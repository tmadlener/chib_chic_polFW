#ifndef CHIBCHICPOLFW_GENERAL_MISCUTILS_H__
#define CHIBCHICPOLFW_GENERAL_MISCUTILS_H__

#include <string>
#include <regex>
#include <iostream>

#include "RooRealVar.h"

/**
 * extract the 'rapX_ptY' part from the passed string (must not necessarily be a filename)
 */
std::string getBinFromFile(const std::string& filename)
{
  // declaring a separate (char *) string here to be able to print it later in case
  constexpr auto rgxstr = ".*_rap([0-9]+)_pt([0-9]+).*";
  std::regex rgx(rgxstr);
  std::smatch cm;
  if (std::regex_match(filename, cm, rgx)) {
    return "rap" + cm[1].str() + "_pt" + cm[2].str();
  }

  std::cerr << "Could't match " << rgxstr << " to " << filename << "\n";
  return "";
}

int getPtBinFromFile(const std::string& filename)
{
  constexpr auto rgxstr =".*(rap[0-9]+_)?pt([0-9]+).*";
  const std::regex rgx(rgxstr);
  std::smatch cm;

  if (std::regex_match(filename, cm, rgx)) {
    return std::stoi(cm[2].str());
  }

  return -1;
}

/** get the ending \'_rapX_ptY\' from the passed values. */
template<typename T>
std::string getBinString(const T iRap, const T iPt)
{
  std::stringstream str;
  str << "_rap" << iRap << "_pt" << iPt;
  return str.str();
}

/** check if the passed value is in ragne (min, max) */
template<typename T>
inline bool inRange(const T val, const T min, const T max)
{
  return (val >= min && val < max);
}

/**
 * check if the passed variable is in the (var->getMin(), var->getMax()).
 */
inline bool inRange(const double val, const RooRealVar* var)
{
  return inRange(val, var->getMin(), var->getMax());
}

/** get the bin of val in binning. return -1 if no bin can be found. */
template<typename T, typename C>
int getBin(const T val, const C& binning)
{
  for (size_t i = 0; i < binning.size() -1; ++i) {
    if (inRange(val, binning[i], binning[i+1])) return i;
  }
  return -1;
}

/**
 * Get the lower bound of the bin in the passed binning.
 * CAUTION: no boundary checks
 */
template<typename C>
inline typename C::value_type getBinMin(const size_t bin, const C& binning)
{
  return (!bin ? binning.front() : binning[bin - 1]);
}

/**
 * Get the upper bound of the bin in the passed binning.
 * CAUTION: no boundary checks
 */
template<typename C>
inline typename C::value_type getBinMax(const size_t bin, const C& binning)
{
  return (!bin ? binning.back() : binning[bin]);
}

/**
 * reduce the range of the passed variable into the range of [-pi, pi].
 */
template<typename T>
T reduceRange(T x)
{
  constexpr T pi = T(M_PI);
  constexpr T o2pi = 1.0 / (2*pi);

  if(std::abs(x) <= T(M_PI)) return x;

  auto n = std::round(x * o2pi);
  return x - n * 2 * pi;
}

/**
 * produce a cut-string that can be used in e.g. TTree::Draw() from the passed binning and bin.
 * If bin == 0, the full range is returned.
 * NOTE: there is no out-of-range check!
 */
template<typename T>
std::string getCutString(const size_t bin, const std::vector<T>& binning, const std::string& varname)
{
  // get the range from the binning, depending if the full range is wanted, or only one bin
  // const T min = !bin ? binning.front() : binning[bin - 1];
  // const T max = !bin ? binning.back() : binning[bin];
  const T min = getBinMin(bin, binning);
  const T max = getBinMax(bin, binning);

  std::stringstream cutstr;
  cutstr << "(" << varname << " >= " << min << " && " << varname << " <= " << max << ")";

  return cutstr.str();
}


/** create a vector with evenly spaced values between min and max (see MATLABs linspace). */
template<typename T>
std::vector<T> linspace(const T min, const T max, const size_t steps)
{
  std::vector<T> vals;
  vals.reserve(steps);
  const T step = (max - min) / (steps - 1);
  for (size_t i = 0; i < steps; ++i) {
    vals.push_back(min + i * step);
  }
  return vals;
}

/**
 * Find the (a) root of the given function using the secant method
 */
template<typename Func>
double findRoot(Func f, double x0, double x1, const double eps=1e-8, const size_t maxSteps=100)
{
  const double fx0 = f(x0);
  const double fx1 = f(x1);
  const double x2 = (x0 * fx1 - x1 * fx0) / (fx1 - fx0);

  if (std::abs(f(x2)) < eps) return x2;

  return findRoot(f, x1, x2, eps, maxSteps - 1);
}

#endif
