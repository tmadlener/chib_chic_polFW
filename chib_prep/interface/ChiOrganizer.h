#ifndef CHIORGANIZER_JN_H
#define CHIORGANIZER_JN_H

/////////////////////////////////////////////////////////////////////////////////
// 
//  ChiOrganizer j.n. 2018
//  -----------------------------------------------------------------------------
//  A class to help to organize the file and workspace names as well as the 
//  configuration for the chi fitting procedure.
//
/////////////////////////////////////////////////////////////////////////////////


#include <string>
#include <map>
#include <fstream>
#include <iostream>

#include "json.hpp" //https://github.com/nlohmann/json
using json = nlohmann::json;

class ChiOrganizer final {

public:
  ChiOrganizer(const std::string &preselection_root_directory = ""); // Using current working directory as default
  ChiOrganizer(const std::string &preselection_root_directory, const std::string & config_file);

  std::string WorkspaceName(const std::string &fitvarname, double min, double max, const std::map<std::string, std::pair<double, double> > bin_varnames_borders);
  std::string FileName(const std::string &fitvarname, double min, double max, const std::map<std::string, std::pair<double, double> > bin_varnames_borders, const std::string extension = ".root");

  std::string WorkspaceName(const std::map<std::string, std::pair<double, double> > bin_varnames_borders);
  std::string FileName(const std::map<std::string, std::pair<double, double> > bin_varnames_borders, const std::string extension = ".root");
  std::map < std::string, std::pair<double, double> > GetVariableList(const std::string & listname);

  template<class T>
  T GetConfigParam(const std::string & valuename, bool &ok);
  template<class T>
  T GetConfigParam(const std::string & valuename, const T & defaultvalue);
  template<class T>
  T GetConfigParam(const std::vector<std::string> & valuename, bool &ok);
  template<class T>
  T GetConfigParam(const std::vector<std::string> &nested_names, const T & defaultvalue);


private:
  std::string make_id(const std::string &fitvarname, double min, double max, const std::map<std::string, std::pair<double, double> > bin_varnames_borders);
  std::string make_id(const std::map<std::string, std::pair<double, double> > bin_varnames_borders);
  const std::string var_separator = "-";
  const std::string range_separator = "_";
  const char decimal_point_character = 'p';
  const std::string workspace_base = "ws" + var_separator;
  const std::string filename_base = "fitresults" + var_separator;
  const std::string default_configfile = "config.json";

  const std::string m_rootdir;
  const std::string config_file;

  std::string prepare_directory(std::string dir);
  std::string check_file(std::string file, const std::string & default_file = "");

  bool is_ok = true;


  template<class T>
  T get_param(const std::vector<std::string> &nested_names, json &j, bool &ok);

};


template<class T>
T ChiOrganizer::GetConfigParam(const std::string & valuename, bool &ok)
{
  return GetConfigParam<T>(std::vector<std::string>{ valuename }, ok);
}

template<class T>
T ChiOrganizer::GetConfigParam(const std::string & valuename, const T & defaultvalue)
{
  return GetConfigParam<T>(std::vector<std::string>{ valuename }, defaultvalue);
}

template<class T>
inline T ChiOrganizer::GetConfigParam(const std::vector<std::string>& nested_names, bool & ok)
{
  ok = false;

  std::ifstream conf(config_file, std::ifstream::in);
  json j;
  conf >> j;

  return get_param<T>(nested_names, j, ok);

}

template<class T>
inline T ChiOrganizer::GetConfigParam(const std::vector<std::string>& nested_names, const T & defaultvalue)
{

  bool ok = false;
  T val = GetConfigParam<T>(nested_names, ok);
  if (ok) return val;

  return defaultvalue;

}

template<class T>
inline T ChiOrganizer::get_param(const std::vector<std::string>& nested_names, json &j, bool & ok)
{
  ok = false;
  T val;

  if (nested_names.empty()) return val;

  if (nested_names.size() == 1) {
    std::cout << "Reading parameter '" << nested_names.front() << "'" << std::endl;
    auto valuename = nested_names.front();
    try
    { 
      val = j[valuename].get<T>();
      ok = true;
    }
    catch (const std::exception&e)
    {
      std::cout << e.what() << std::endl;
      std::cout << "Could not get value for '" << valuename << "' from config file '" << config_file << "'." << std::endl;
    }
    return val;
  }

  return get_param<T>(std::vector<std::string>(nested_names.begin() + 1, nested_names.end()),j[nested_names.front()], ok);

}


#endif //CHIORGANIZER_JN_H