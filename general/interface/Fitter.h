#ifndef FITTER_JN_H
#define FITTER_JN_H

/////////////////////////////////////////////////////////////////////////////////
// 
//  FITTER j.n. 2018
//  -----------------------------------------------------------------------------
//  A helper class for setting up a fit for RooFit.
//
/////////////////////////////////////////////////////////////////////////////////

#include <string>
#include <vector>
#include <memory>
#include <algorithm>

class RooWorkspace;
class RooDataSet;
class TFile;
class FitAnalyser;

class Fitter
{

  friend FitAnalyser;

public:
  using csr = const std::string &;
  void SetWorkspaceName(csr workspace_name, bool force_refit = false) {
    m_wsname = workspace_name;
    m_force_refit = force_refit;
  }
  void SetOutfile(csr outfile_name, bool force_recreation = false) {
    m_ofname = outfile_name; 
    m_force_recreation = force_recreation;
  }

  void SetInputData(csr datafile_name, csr tree_name) { 
    m_in_filename = datafile_name;
    m_in_wsname = "";
    m_in_dsname = "";
    m_in_treename = tree_name;
  }
  
  void SetInputData(csr datafile_name, csr workspace_name, csr dataset_name) {
    m_in_filename = datafile_name;
    m_in_treename = "";
    m_in_wsname = workspace_name;
    m_in_dsname = dataset_name;

  }
  
  void SetModel(const std::vector<std::string>& factory_strings, csr model_name, csr fitvariable_name, double min, double max) {
    m_factory_strings = factory_strings;
    m_modelname = model_name;
    m_fitvarname = fitvariable_name;
    AddBinVariable(fitvariable_name, min, max);
  }

  void AddVariable(csr variable_name) {
    AddBinVariable(variable_name, 0, 0);
  }

  void AddBinVariable(csr variable_name, double min, double max) {

    // If min == max -> full range
    
    // If variable already exists it is updated with the given borders
    auto it = std::find(m_roo_vars.begin(), m_roo_vars.end(), variable_name);
    if (it != m_roo_vars.end()) {
      m_roovar_borders.at(std::distance(m_roo_vars.begin(), it)) = { min,max };
    }
    else {
      m_roo_vars.push_back(variable_name);
      m_roovar_borders.push_back({ min, max });
    }
  }

  void RemoveBinVariables(csr variable_name) {
    //TODO: implement
  }

  void Fit(int numCPUs = 8, bool enableMinos = true);

  virtual ~Fitter();

  static const std::string dataset_name;
  static const std::string snapshot_name;

  

private:
  std::string m_wsname;
  bool m_force_refit = false;
  std::string m_ofname;
  bool m_force_recreation = false;
  std::string m_in_filename;
  std::string m_in_treename;
  std::string m_in_wsname;
  std::string m_in_dsname;
  std::string m_fitvarname;
  std::string m_modelname;

  std::vector <std::string> m_roo_vars;
  std::vector < std::pair<double, double> > m_roovar_borders;
  std::vector <std::string> m_factory_strings;

  RooDataSet get_dataset(TFile* f, bool *ok = nullptr);
  bool setup_workspace(RooWorkspace*, RooDataSet */*, bool *ok = nullptr*/);
  bool params_ok();

};


#endif //FITTER_JN_H