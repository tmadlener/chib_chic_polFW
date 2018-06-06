#ifndef CHIB_CHIC_POLFW_TOYMCGENERATION_SMEARIN_H__
#define CHIB_CHIC_POLFW_TOYMCGENERATION_SMEARIN_H__

#include "../general/interface/root_utils.h"

#include "TH1D.h"
#include "TH2D.h"
#include "TFile.h"
#include "TLorentzVector.h"

#include <vector>
#include <iostream>
#include <string>

class SmearingProvider {
public:
  SmearingProvider(TH2D *residualsMap);

  double getSmearing(double varVal) const;

private:
  std::vector<TH1D*> m_varProjections;
  TH2D *m_residualsMap;
};

SmearingProvider::SmearingProvider(TH2D *residualMap) : m_residualsMap(residualMap)
{
  const int nVarBins = m_residualsMap->GetXaxis()->GetNbins();
  m_varProjections.reserve(nVarBins + 2);

  const std::string mapName = m_residualsMap->GetName();

  // copy the first bin into the underflow bin
  auto projection = m_residualsMap->ProjectionY((mapName + "bin_0").c_str(), 1, 1, "e");
  m_varProjections.push_back(projection);

  for (int iBin = 1; iBin <= nVarBins; ++iBin) {
    auto *binProjection = m_residualsMap->ProjectionY((mapName + "bin_" + std::to_string(iBin)).c_str(), iBin, iBin, "e");
    m_varProjections.push_back(binProjection);
  }

  // copy the last bin into the overflow bin
  projection = m_residualsMap->ProjectionY((mapName + "bin_" + std::to_string(nVarBins + 1)).c_str(), nVarBins, nVarBins, "e");
  m_varProjections.push_back(projection);

  std::cout << m_residualsMap->GetName() << "\n";
  for (auto * h: m_varProjections) std::cout << h << " " << h->GetName() << "\n";
}

double SmearingProvider::getSmearing(double varVal) const
{
  const int varBin = m_residualsMap->GetXaxis()->FindBin(varVal);
  if (varBin >= m_varProjections.size()) {
    std::cout << "bin for value: " << varVal << " " << varBin << " in " << m_residualsMap->GetName() << " | " << m_varProjections.size() << "\n";
  }
  return m_varProjections[varBin]->GetRandom();
}


TLorentzVector smearParticle(const TLorentzVector& particle, const SmearingProvider& xyzSmearer)
{
  const double spx = particle.Px() + particle.Px() * xyzSmearer.getSmearing(particle.Px());
  const double spy = particle.Py() + particle.Py() * xyzSmearer.getSmearing(particle.Py());
  const double spz = particle.Pz() + particle.Pz() * xyzSmearer.getSmearing(particle.Pz());

  TLorentzVector smeared;
  smeared.SetXYZM(spx, spy, spz, particle.M());
  return smeared;
}


TLorentzVector smearParticle(const TLorentzVector& particle, const SmearingProvider& xySmearer, const SmearingProvider& zSmearer)
{
  const double spx = particle.Px() + particle.Px() * xySmearer.getSmearing(particle.Px());
  const double spy = particle.Py() + particle.Py() * xySmearer.getSmearing(particle.Py());
  const double spz = particle.Pz() + particle.Pz() * zSmearer.getSmearing(particle.Pz());


  TLorentzVector smeared;
  smeared.SetXYZM(spx, spy, spz, particle.M());
  return smeared;
}

#endif
