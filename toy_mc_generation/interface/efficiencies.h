#ifndef CHIBCHICPOLFW_TOY_MC_GENERATION_EFFICIENCIES_H__
#define CHIBCHICPOLFW_TOY_MC_GENERATION_EFFICIENCIES_H__

#include "range.h"
#include "../../general/interface/misc_utils.h"

#include "TFile.h"
#include "TIterator.h"
#include "TKey.h"
#include "TGraphAsymmErrors.h"
#include "TF1.h"

#include <vector>
#include <string>
#include <iostream>
#include <numeric>


template<typename EffType>
class Efficiency {
public:
  Efficiency(EffType* eff, double minEta, double maxEta, double ptMin, double ptMax) :
    m_eff(eff), m_eta(Range{minEta, maxEta}), m_pT(Range{ptMin, ptMax}) {}

  double Eval(double pT) const { return m_eff->Eval(pT); }

  bool etaContain(double absEta) const { return m_eta.min <= absEta && m_eta.max > absEta; }
  bool ptContain(double pt) const { return m_pT.min <= pt && m_pT.max > pt; }

  std::string const getName() const { return m_eff->GetName(); }

private:
  EffType *m_eff;

  Range m_eta;
  Range m_pT;
};


/**
 * Base class to derive from to be able to provide an Eval function
 *
 * Checked on godbolt.org and in principle this should be optimized out by the compiler
 */
struct EffProv {
  virtual double Eval(double, double) const = 0;
};


template<typename EffType>
class EfficiencyProvider : public EffProv {
public:
  template<typename RangeFinder>
  EfficiencyProvider(const std::string& effFile, const std::string& identifier, RangeFinder rangeFinder);

  double Eval(double pT, double eta) const override final;

private:
  TFile *m_effFile;
  std::vector<Efficiency<EffType>> m_efficiencies;
};

template<typename EffType>
template<typename RangeFinder>
EfficiencyProvider<EffType>::EfficiencyProvider(const std::string& effFile, const std::string& identifier, RangeFinder rangeFinder)
{
  m_effFile = TFile::Open(effFile.c_str());
  // NOTE: here we just assume that all objects matching the identifier can be correctly static_cast'ed
  m_efficiencies.reserve(m_effFile->GetListOfKeys()->GetSize());
  auto nextKey = TIter(m_effFile->GetListOfKeys());
  TKey* key = nullptr;
  while((key = static_cast<TKey*>(nextKey()))) {
    std::string const keyName = key->GetName();
    if (keyName.find(identifier) != std::string::npos) {
      // check if not already present. NOTE: this is necessary because the keys do not have to be unique in TFiles
      if (std::find_if(m_efficiencies.cbegin(), m_efficiencies.cend(),
                       [&keyName] (Efficiency<EffType> const& eff)
                       { return keyName == eff.getName(); }) == m_efficiencies.cend()) {
        auto *obj = static_cast<EffType*>(key->ReadObj());
        const auto etaRange = getEtaRange(obj->GetName());
        const auto ptRange = rangeFinder(obj);

        m_efficiencies.emplace_back(obj, etaRange.min, etaRange.max, ptRange.min, ptRange.max);
      }
    }
  }
}

template<typename EffType>
double EfficiencyProvider<EffType>::Eval(double pT, double eta) const
{
  // NOTE: this could be optimized by indexing directly to the efficiency but for now let's see if this is good enough
  for (const auto& eff : m_efficiencies) {
    if (eff.etaContain(std::abs(eta)) && eff.ptContain(pT)) {
      return eff.Eval(pT);
    }
  }

  return -1; // to clearly denote values that are not in the acceptance
}

template<typename EffType>
struct FixedRange {
  FixedRange(Range range) : m_range(range) {}
  FixedRange(double min, double max): m_range(Range{min, max}) {}

  Range operator()(EffType*) {
    return m_range;
  }
private:
  Range m_range;
};


struct RangeFromGraph {
  Range operator()(TGraphAsymmErrors *graph) {
    double x, y;
    graph->GetPoint(0, x, y);
    const double min = x - graph->GetErrorXlow(0);

    const int nPoints = graph->GetN();
    graph->GetPoint(nPoints - 1, x, y);
    const double max = x + graph->GetErrorXhigh(nPoints - 1);

    return Range{min, max};
  }
};

/**
 * For some of the fits the efficiency parametrization are TF1 that are built from more than 1 original function.
 * In that case TF1->GetMaximumX() and TF1->GetMinimumX() do not return the correct range, but it can be determined
 * by checking in which range the function is not 0
 * NOTE: assuming that this will be a contiguous range
 */
Range findRangeFromValues(const TF1 *func, const Range widestRange, const size_t nPoints=10000)
{
  double min = widestRange.min;
  double max = widestRange.max;

  double lastFuncVal = func->Eval(widestRange.min);
  for (const auto x : linspace<double>(widestRange.min, widestRange.max, nPoints)) {
    const double funcVal = func->Eval(x);
    if (funcVal != 0 && lastFuncVal == 0) {
      min = x;
    }
    if (lastFuncVal !=0 && funcVal == 0) {
      max = x;
    }

    lastFuncVal = funcVal;
  }

  return Range{min, max};
}


/**
 * Determine the range from a parametrization, such that the resulting efficiency will always be larger than or equal to 0
 * NOTE: this currently only checks for values smaller than 0 at the lower edge
 */
struct RangeFromFit {
  /** default constructor uses the range of the TF1 to check for a positive range */
  RangeFromFit() = default;
  /** constructor for using a proposal Range and check whether the function is non 0 over that proposal range */
  RangeFromFit(Range propRange) : m_findNonNullRange(true), m_proposalRange(propRange) {}

  Range operator()(TF1 *fit) {
    std::cout << "Finding range from fit: " << fit->GetName() << "\n";
    const auto funcRange = m_findNonNullRange ? findRangeFromValues(fit, m_proposalRange, 10000) : Range{fit->GetXmin(), fit->GetXmax()};

    const double minX = funcRange.min;
    const double maxX = funcRange.max;

    std::cout << "Searching for range in (" << minX << ", " << maxX << ")\n";

    double min = minX;

    // evaluate the function at 1000 evenly spaced points between min and max and see if any of them are below 0
    // if none are, assume the whole function is positive and simply return the range of the function
    // if there are points below zero do a simple Newton method to find the root and use that as lower limit
    constexpr size_t nPoints = 1000;
    double lastNegative = std::numeric_limits<double>::min();
    bool foundNegative = false;
    for (const auto x : linspace<double>(minX, maxX, nPoints)) {
      if (fit->Eval(x) < 0) {
        lastNegative = x;
        foundNegative = true;
      }
    }
    if (foundNegative) {
      const double defPositive = (maxX - minX) / nPoints + lastNegative; // has to be positive by definition
      std::cout << "Found negative values between " << lastNegative << " and " << defPositive << ". Now trying to find a more precise point\n";
      min = findRoot([fit](double x) { return fit->Eval(x); }, lastNegative, defPositive);
    }

    std::cout << "Final range: (" << min << "," << maxX << ")\n";
    return Range{min, maxX};
  }
private:
  bool m_findNonNullRange{false};
  Range m_proposalRange{0, 0};
};

#endif
