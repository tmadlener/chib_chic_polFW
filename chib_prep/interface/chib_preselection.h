#include "TreeProcessor.h"

class ChibPreselection : public TreeProcessor
{

public:
  ChibPreselection(const std::vector<std::string> & infilenames, const std::string & outfilename, const std::string &intreename, const std::string &outtreename);

private:
  virtual bool fill_and_cut_variables() override;
  virtual ParallelTreeLooper* clone() const override { return new ChibPreselection(*this); }
  virtual void setup_new_branches() override;

  std::vector<std::string> triggers;

  // Variables for new branches

  Double_t dimuon_mass;
  Double_t chi_mass;
  Double_t dimuon_pt;
  Double_t chi_pt;
  Double_t dimuon_rap;
  Double_t chi_rap;

  Double_t costh_hx;
  Double_t phi_hx;
  Double_t cosalpha_hx;
  Double_t costh_px;
  Double_t phi_px;
  Double_t cosalpha_px;
  Double_t costh_cs;
  Double_t phi_cs;
  Double_t cosalpha_cs;

  Double_t delta_phi;

  Double_t eta_mupos;
  Double_t pt_mupos;
  Double_t phi_mupos;
  Double_t eta_muneg;
  Double_t pt_muneg;
  Double_t phi_muneg;
  
};