#ifndef CHIBCHICPOLFW_TOY_MC_GENERATION_HISTOGRAM_H__
#define CHIBCHICPOLFW_TOY_MC_GENERATION_HISTOGRAM_H__

#include "THn.h"

#include <array>
#include <functional>
#include <string>


template<size_t N>
class StorageHistograms {

public:
  // do nothing on construction
  StorageHistograms() = default;

  ~StorageHistograms();

  // this one will initialize the histograms (as necessary)
  void Init(const std::string& basename, const bool sampling, const bool effs,
            const std::array<int, N>& bins, const std::array<double,N>& low, const std::array<double,N>& high);

  void fillGen();
  void fillAcc();
  void fillReco();
  void fillNoweight();

  template<typename F>
  void setSamplingWeightF(F getSamplingWeight) { m_samplingWeight = getSamplingWeight; }

  template<typename F>
  void setRecoEffWeightF(F getRecoEffWeight) { m_recoEfficiency = getRecoEffWeight; }

  template<typename F>
  void setHXFillF(F getHXValues) { m_fillFuncs[0] = getHXValues; }

  template<typename F>
  void setPXFillF(F getPXValues) { m_fillFuncs[1] = getPXValues; }

  template<typename F>
  void setCSFillF(F getCSValues) { m_fillFuncs[2] = getCSValues; }

private:
  std::array<THnD*, 3> m_genHists; // one per ref frame
  std::array<THnD*, 3> m_accHists; // one per ref frame
  std::array<THnD*, 3> m_recoHists; // one per ref frame

  std::array<THnD*, 9> m_noWeightHists; // one per histogram above

  std::array<std::function<std::array<double, N>()>, 3> m_fillFuncs; // one per frame
  std::function<double()> m_samplingWeight; // function to retrieve the samplingWeight
  std::function<double()> m_recoEfficiency; // function to retrieve the reco efficiency weight

  double m_sampleWCache{0};
  double m_recoWCache{0};
  bool m_sampling{false};
  bool m_effs{false}; // is it possible to fill the reco hists, if this is true, it is
  // use this to indicate if the whole class has been initialized. If not then all the calls
  // to fillXXX amount to checking this variable and doing nothing else
  bool m_filling{false};
  int n_called{0};
};


template<size_t N>
void StorageHistograms<N>::Init(const std::string& basename, const bool sampling, const bool effs,
                                const std::array<int, N>& bins, const std::array<double,N>& low, const std::array<double,N>& high)
{
  m_filling = true;
  m_sampling = sampling;
  m_effs = effs;

  std::cout << "----------------------------------------" << std::endl;

  // NOTE: Make sure this corresponds to the assignment of the fill functions
  const std::array<std::string, 3> frames = {"HX", "PX", "CS"};
  const std::string noweight = "noweight_";
  for (size_t i = 0; i < 3; ++i) {
    const std::string genName = basename + "_gen_" + frames[i];
    const std::string accName = basename + "_acc_" + frames[i];
    const std::string recoName = basename + "_reco_" + frames[i];

    m_genHists[i] = new THnD(genName.c_str(), ";cos#vartheta;#varphi", N, bins.data(), low.data(), high.data());
    m_accHists[i] = new THnD(accName.c_str(), ";cos#vartheta;#varphi", N, bins.data(), low.data(), high.data());
    m_recoHists[i] = new THnD(recoName.c_str(), ";cos#vartheta;#varphi", N, bins.data(), low.data(), high.data());

    if (m_sampling) {
      m_noWeightHists[i] = new THnD((noweight + genName).c_str(), ";cos#vartheta; #varphi", N, bins.data(), low.data(), high.data());
      m_noWeightHists[i+3] = new THnD((noweight + accName).c_str(), ";cos#vartheta; #varphi", N, bins.data(), low.data(), high.data());
      m_noWeightHists[i+6] = new THnD((noweight + recoName).c_str(), ";cos#vartheta; #varphi", N, bins.data(), low.data(), high.data());
    }
  }

  for (auto* h : m_genHists) h->Sumw2();
  for (auto* h : m_accHists) h->Sumw2();
  for (auto* h : m_recoHists) h->Sumw2();
  if (m_sampling) {
    for (auto* h : m_noWeightHists) h->Sumw2();
  }
}


template<size_t N>
StorageHistograms<N>::~StorageHistograms()
{
  if (m_filling) {
    for (auto* h : m_genHists) h->Write();
    for (auto* h : m_accHists) h->Write();
    for (auto* h : m_recoHists) h->Write();

    if (m_sampling) {
      for (auto* h : m_noWeightHists) h->Write();
    }
  }
}

template<size_t N>
void StorageHistograms<N>::fillGen()
{
  if (!m_filling) return;
  m_sampleWCache = m_samplingWeight();
  for (size_t i = 0; i < 3; ++i) {
    const auto vals = m_fillFuncs[i]();
    m_genHists[i]->Fill(vals.data(), m_sampleWCache);
    if (m_sampling) {
      m_noWeightHists[i]->Fill(vals.data());
    }
  }
}


template<size_t N>
void StorageHistograms<N>::fillAcc()
{
  if (!m_filling) return;
  m_sampleWCache = m_samplingWeight();
  for (size_t i = 0; i < 3; ++i) {
    const auto vals = m_fillFuncs[i]();
    m_accHists[i]->Fill(vals.data(), m_sampleWCache);
    if (m_sampling) {
      m_noWeightHists[3 + i]->Fill(vals.data());
    }
  }
}

template<size_t N>
void StorageHistograms<N>::fillReco()
{
  if (!m_filling || !m_effs) return;

  m_recoWCache = m_recoEfficiency();
  m_sampleWCache = m_samplingWeight();

  for (size_t i = 0; i < 3; ++i) {
    const auto vals = m_fillFuncs[i]();
    m_recoHists[i]->Fill(vals.data(), m_sampleWCache * m_recoWCache);
    if (m_sampling) {
      m_noWeightHists[6 + i]->Fill(vals.data(), m_recoWCache);
    }
  }
}


#endif
