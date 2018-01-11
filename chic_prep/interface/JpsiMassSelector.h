#ifndef CHIBCHICPOLFW_CHICPREP_JPSIMASSSELECTOR_H__
#define CHIBCHICPOLFW_CHICPREP_JPSIMASSSELECTOR_H__

#include "general/interface/root_utils.h"

#include "RooWorkspace.h"

#include <string>

class JpsiMassSelector {
public:
  JpsiMassSelector() = delete;
  JpsiMassSelector(RooWorkspace* ws, const double nSigma) :
    m_nSigma(nSigma),
    m_p0(getVarVal(ws, m_p0Name)),
    m_pSig0(getVarVal(ws, m_pSig0Name)),
    m_pSig1(getVarVal(ws, m_pSig1Name)),
    m_pSig2(getVarVal(ws, m_pSig2Name))
  { }

  bool contains(const double mass, const double rap) const
  {
    return std::abs(mass - m_p0) < m_nSigma * rapSigma(std::abs(rap));
  }

  std::string getTFormulaExpr(const std::string& massName, const std::string& rapName) const
  {
    return "TMath::Abs(" + massName + " - " + std::to_string(m_p0) + ") < "
      + std::to_string(m_nSigma) + " * ("
      + std::to_string(m_pSig0) + " + "
      + std::to_string(m_pSig1) + " * TMath::Abs(" + rapName + ")" + " + "
      + std::to_string(m_pSig2) + " * " + rapName + " * " + rapName + ")";
  }

private:
  static constexpr auto m_p0Name = "CBmass_p0_jpsi";
  static constexpr auto m_pSig0Name = "CBsigma_p0_jpsi";
  static constexpr auto m_pSig1Name = "CBsigma_p1_jpsi";
  static constexpr auto m_pSig2Name = "CBsigma_p2_jpsi";

  const double m_nSigma;
  const double m_p0;
  const double m_pSig0;
  const double m_pSig1;
  const double m_pSig2;

  double rapSigma(const double absRap) const
  {
    return m_pSig0 + m_pSig1 * absRap + m_pSig2 * absRap * absRap;
  }
};

#endif
