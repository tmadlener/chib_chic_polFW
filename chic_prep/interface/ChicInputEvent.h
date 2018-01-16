#ifndef CHIBCHICPOLFW_CHICPREP_CHICINPUTEVENT_H__
#define CHIBCHICPOLFW_CHICPREP_CHICINPUTEVENT_H__

#include "TTree.h"
#include "TLorentzVector.h"

struct VoidInfo {
  void Init(TTree*) { /* No-OP */ }
  void Create(TTree*) { /* No-OP */ }
};

struct DefaultBranchNames  {
  static constexpr auto jpsi = "jpsi";
  static constexpr auto chic = "chic";
  static constexpr auto lepP = "lepP";
  static constexpr auto lepN = "lepN";
  static constexpr auto Jpsict = "Jpsict";
};

template<typename AdditionalInfo = VoidInfo,
         typename BranchNames = DefaultBranchNames>
class ChicInputEvent {
public:
  ChicInputEvent() = default;

  ~ChicInputEvent();

  void Init(TTree *t);

  const TLorentzVector &jpsi() const { return *m_jpsi; }
  const TLorentzVector &chic() const { return *m_chic; }
  const TLorentzVector &muP() const { return *m_muP; }
  const TLorentzVector &muN() const { return *m_muN; }

  const AdditionalInfo& info() const { return m_addInfo; }

  double Jpsict;
private:

  TLorentzVector *m_jpsi{nullptr};
  TLorentzVector *m_chic{nullptr};
  TLorentzVector *m_muP{nullptr};
  TLorentzVector *m_muN{nullptr};

  AdditionalInfo m_addInfo;
};


template<typename AI, typename BN>
ChicInputEvent<AI, BN>::~ChicInputEvent()
{
  delete m_jpsi;
  delete m_chic;
  delete m_muP;
  delete m_muN;
}

template<typename AI, typename BN>
void ChicInputEvent<AI, BN>::Init(TTree *t)
{
  t->SetBranchAddress(BN::jpsi, &m_jpsi);
  t->SetBranchAddress(BN::chic, &m_chic);
  t->SetBranchAddress(BN::lepP, &m_muP);
  t->SetBranchAddress(BN::lepN, &m_muN);

  t->SetBranchAddress(BN::Jpsict, &Jpsict);

  m_addInfo.Init(t);
}

#endif
