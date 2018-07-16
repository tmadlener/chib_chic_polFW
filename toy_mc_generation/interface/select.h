#ifndef CHIBCHICPOLFW_TOYMCGENERATION_SELECT_H__
#define CHIBCHICPOLFW_TOYMCGENERATION_SELECT_H__

#include "range.h"

#include "TLorentzVector.h"

#include <memory>
#include <string>

struct Selector {
  virtual bool accept(TLorentzVector const&) const = 0;
};

struct AllSelector : public Selector {
  virtual bool accept(TLorentzVector const&) const override final { return true; }
};

class MinPtSelector : public Selector {
  double minPt;
public:
  MinPtSelector(double MinPt) : minPt(MinPt) {}

  virtual bool accept(TLorentzVector const& p4) const override final { return p4.Pt() > minPt; }
};

class MinPtMaxEtaSelector : public Selector {
  double minPt;
  double maxEta;
public:
  MinPtMaxEtaSelector(double MinPt, double MaxEta) : minPt(MinPt), maxEta(MaxEta) {}

  virtual bool accept(TLorentzVector const& p4) const override final {
    return p4.Pt() > minPt && std::abs(p4.Eta()) < maxEta;
  }
};

class PtRangeAbsRapiditySelector : public Selector {
  Range ptRange;
  double maxRap;
public:
  PtRangeAbsRapiditySelector(Range PtRange, double MaxRap) : ptRange(PtRange), maxRap(MaxRap) {}

  virtual bool accept(TLorentzVector const& p4) const override final {
    return std::abs(p4.Rapidity()) < maxRap && p4.Pt() > ptRange.min && p4.Pt() < ptRange.max;
  }
};

class LooseMuonSelector : public Selector {
public:
  virtual bool accept(TLorentzVector const& p4) const override final {
    const double absEta = std::abs(p4.Eta());
    const double pT = p4.Pt();

    // pT > 3.5 GeV, if |eta| < 1.2
    // pT > 3.5 - (|eta| - 1.2) * 2.5 GeV, if 1.2 < |eta| < 1.6
    return ((absEta < 1.2 && pT > 3.5) ||
            (absEta > 1.2 && absEta < 1.6 && pT > 3.5 - (absEta - 1.2) * 2.5));
  }
};

#endif
