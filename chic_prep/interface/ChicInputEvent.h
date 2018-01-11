#ifndef CHIBCHICPOLFW_CHICPREP_CHICINPUTEVENT_H__
#define CHIBCHICPOLFW_CHICPREP_CHICINPUTEVENT_H__

#include "TTree.h"
#include "TLorentzVector.h"


class ChicInputEvent {
public:
  ChicInputEvent() = default;

  ~ChicInputEvent();

  void Init(TTree *t);

  const TLorentzVector &jpsi() const { return *m_jpsi; }
  const TLorentzVector &chic() const { return *m_chic; }
  const TLorentzVector &muP() const { return *m_muP; }
  const TLorentzVector &muN() const { return *m_muN; }

  double Jpsict;
private:

  TLorentzVector *m_jpsi{nullptr};
  TLorentzVector *m_chic{nullptr};
  TLorentzVector *m_muP{nullptr};
  TLorentzVector *m_muN{nullptr};
};


ChicInputEvent::~ChicInputEvent()
{
  delete m_jpsi;
  delete m_chic;
  delete m_muP;
  delete m_muN;
}

void ChicInputEvent::Init(TTree *t)
{
  t->SetBranchAddress("jpsi", &m_jpsi);
  t->SetBranchAddress("chic", &m_chic);
  t->SetBranchAddress("lepP", &m_muP);
  t->SetBranchAddress("lepN", &m_muN);

  t->SetBranchAddress("Jpsict", &Jpsict);
}

#endif
