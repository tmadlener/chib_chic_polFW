#ifndef FITANALYSER_JN_H
#define FITANALYSER_JN_H


/////////////////////////////////////////////////////////////////////////////////
// 
//  FITANALYSER j.n. 2018
//  -----------------------------------------------------------------------------
//  A class to read out fitting results from a RooWorkspace
//  and to make some calculations and plotting.
//
/////////////////////////////////////////////////////////////////////////////////

#include "Fitter.h"

#include "TFile.h"
#include <string>


class RooRealVar;
class RooAbsData;
class RooFitResult;
class RooAbsPdf;

class FitAnalyser {
public:
  using csr = const std::string &;

  FitAnalyser(const Fitter &);

  void SetFileName(csr file_name) {
    m_filename = file_name;
  };
  void SetWorkspaceName(csr workspace_name) {
    m_wsname = workspace_name;
  }
  void SetFitVariableName(csr variable_name) {
    m_fitvarname = variable_name;
  }
  void SetModelName(csr model_name) {
    m_modelname = model_name;
  }
  
  RooWorkspace *GetWorkspace(bool load_snapshot = false);
  RooAbsData *GetDataset(csr dataset_name = "");
  RooRealVar *GetVariable(csr variable_name = "");
  RooFitResult *GetFitResult(csr result_name = "");
  RooAbsPdf *GetPdf(csr pdf_name = "");
  double GetVariableValue(csr variable_name, bool &ok);
  void PlotFitResult(csr output_file = "");
  
private:
  std::string m_filename;
  std::string m_wsname;
  std::string m_snapshotname = Fitter::snapshot_name;
  std::string m_datasetname = Fitter::dataset_name;
  std::string m_fitvarname;
  std::string m_modelname;

  TFile* open_file();
  std::unique_ptr<TFile> m_file;

  std::pair<double, double> get_chi2_ndof();

};

#endif // !FITANALYSER_JN_H
