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
#include <memory>


class RooRealVar;
class RooAbsData;
class RooFitResult;
class RooAbsPdf;
class RooWorkspace;
class TTree;

class FitAnalyser {
public:
  using csr = const std::string &;

  FitAnalyser(const Fitter &);
  FitAnalyser(csr filename, csr model_name, csr fit_variable_name, csr workspace_name = "", csr dataset_name = "", csr snapshot_name = "");

  void SetWorkspaceName(csr workspace_name);
  void SetFitVariableName(csr variable_name) { m_fitvarname = variable_name; }
  void SetModelName(csr model_name) { m_modelname = model_name; }
  void SetSnapshotName(csr snapshot_name) { m_snapshotname = snapshot_name; }
  void SetDatasetName(csr dataset_name) { m_datasetname = dataset_name; }

  RooWorkspace *GetWorkspace(bool load_snapshot = false);
  RooDataSet *GetDataset(csr dataset_name = "");
  RooRealVar *GetVariable(csr variable_name = "");
  RooFitResult *GetFitResult(csr result_name = "");
  RooAbsPdf *GetPdf(csr pdf_name = "");
  double GetVariableValue(csr variable_name, bool &ok);
  void PlotFitResult(csr output_file = "", const std::vector<double> & lines = std::vector<double>(), bool addResultBox = true);
  double EvaluateFormula(csr formula, bool &ok);
  double GetXForQuantile(csr pdf_name, double quantileq, bool &ok);
  //double GetIntegralValue(/*TODO: pdf_name, range*/);
  RooDataSet * AddSWeights(const std::vector<std::string> & yield_names);

private:
  const std::string m_filename;
  std::string m_wsname;
  std::string m_snapshotname = Fitter::snapshot_name;
  std::string m_datasetname = Fitter::dataset_name;
  std::string m_fitvarname;
  std::string m_modelname;
  RooWorkspace *m_workspace = nullptr;

  TFile* get_file();
  std::unique_ptr<TFile> m_file;

  std::pair<double, double> get_chi2_ndof();

};

#endif // !FITANALYSER_JN_H
