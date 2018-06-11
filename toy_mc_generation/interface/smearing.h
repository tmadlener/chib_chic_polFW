#ifndef CHIB_CHIC_POLFW_TOYMCGENERATION_SMEARIN_H__
#define CHIB_CHIC_POLFW_TOYMCGENERATION_SMEARIN_H__

#include "../../general/interface/root_utils.h"

#include "TH1D.h"
#include "TH2D.h"
#include "TFile.h"
#include "TLorentzVector.h"
#include "TRandom.h"
#include "TF1.h"

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


TLorentzVector smearParticle(const TLorentzVector& particle, const double deltaP)
{
  const double smearFactor = (particle.P() + deltaP) / particle.P();

  const double spx = particle.Px() * smearFactor;
  const double spy = particle.Py() * smearFactor;
  const double spz = particle.Pz() * smearFactor;

  TLorentzVector smeared;
  smeared.SetXYZM(spx, spy, spz, particle.M());
  return smeared;
}

inline TLorentzVector smearParticleTF1(const TLorentzVector& particle, TF1* crystalBall)
{
  return smearParticle(particle, crystalBall->GetRandom() * particle.P());
}

inline TLorentzVector smearParticleGaus(const TLorentzVector& particle, const double mean, const double sigma)
{
  return smearParticle(particle, gRandom->Gaus(mean, sigma));
}

// NOTE: this is not really tested in this way, so make sure that MC smearings make sense in using them this way
inline TLorentzVector smearParticleMCSmearing(const TLorentzVector& particle, const SmearingProvider& smearingProvider)
{
  return smearParticle(particle, smearingProvider.getSmearing(particle.P()));
}

#endif
