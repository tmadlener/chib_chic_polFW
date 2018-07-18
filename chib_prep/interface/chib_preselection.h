#include "TreeProcessor.h"

class ChibPreselection : public TreeProcessor
{

public:
  ChibPreselection(const std::vector<std::string> & infilenames, const std::string & outfilename, const std::string &intreename, const std::string &outtreename, int year);
  virtual ~ChibPreselection() {};

private:
  int m_year = 0;
  double min_pt = 8;
  int m_trigminval = 1;
  std::string muonNname;

  virtual bool fill_and_cut_variables() override;
  virtual ParallelTreeLooper* clone() const override { return new ChibPreselection(*this); }
  virtual void setup_new_branches() override;

  std::vector<std::string> triggers;

  // Variables for new branches
 
  std::vector<std::vector<std::string>> var_collections;

  std::vector<Double_t> out_vars;
  std::vector<Double_t> out_vars_rf1S;
  std::vector<Double_t> out_vars_rf2S;
  std::vector<Double_t> out_vars_rf3S;
  std::vector<Double_t> photon_vars;

  static std::mutex entry_mtx;
  static Long64_t entry_idx;
  Long64_t local_entry_idx;

  // RooDataSet cannot handle a Long64_t, so split it into two Int_t:
  Int_t entry_idx_low = 0;
  Int_t entry_idx_high = 0;

  int make_lorentz_flat(std::vector<Double_t> & vars, TLorentzVector* chi, TLorentzVector* dimuon = nullptr);
  void make_mu_things(std::vector<Double_t> & vars, TLorentzVector* muPos, TLorentzVector *muNeg, int start);
  void setup_collections(const std::string & varsuffix, std::vector<Double_t> & vars, bool only_chi = false);
  void setup_photon_vars();
  void fill_photon_vars(std::vector<TLorentzVector*> photon_vecs);
  
  bool accept_muon(const TLorentzVector * mu);
  
};