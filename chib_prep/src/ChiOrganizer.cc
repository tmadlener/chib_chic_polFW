#include "ChiOrganizer.h"
#include "utils.h"

#include <sstream>
#include <algorithm>
#include <map>

ChiOrganizer::ChiOrganizer(const std::string & rd) :
  m_rootdir(prepare_directory(rd)),
  config_file(check_file(m_rootdir + default_configfile))
{ }

ChiOrganizer::ChiOrganizer(const std::string & rd, const std::string & cf) :
  m_rootdir(prepare_directory(rd)),
  config_file(check_file(cf, m_rootdir + default_configfile))
{ }

std::string ChiOrganizer::WorkspaceName(const std::string & fitvarname, double min, double max, const std::map<std::string, std::pair<double, double>> bin_varnames_borders)
{
  return workspace_base + make_id(fitvarname, min, max, bin_varnames_borders);
}

std::string ChiOrganizer::FileName(const std::string & fitvarname, double min, double max, const std::map<std::string, std::pair<double, double>> bin_varnames_borders, const std::string extension)
{
  return m_rootdir + filename_base + make_id(fitvarname, min, max, bin_varnames_borders) + extension;
}

std::string ChiOrganizer::WorkspaceName(const std::map<std::string, std::pair<double, double>> bin_varnames_borders)
{
  return workspace_base + make_id(bin_varnames_borders);
}

std::string ChiOrganizer::FileName(const std::map<std::string, std::pair<double, double>> bin_varnames_borders, const std::string extension)
{
  return m_rootdir + filename_base + make_id(bin_varnames_borders) + extension;
}

std::map < std::string, std::pair<double, double> > ChiOrganizer::GetVariableList(const std::string & listname)
{
  std::map < std::string, std::pair<double, double> > varlist;
  bool ok = false;
  auto jsonlist = GetConfigParam< json >(listname, ok);
  if (ok) {
    for (json &j : jsonlist) {
      bool okname = false, okmin = false, okmax = false;
      auto name = get_param<std::string>(std::vector<std::string>{ "name" }, j, okname);
      auto min = get_param<double>(std::vector<std::string>{ "min" }, j, okmin);
      auto max = get_param<double>(std::vector<std::string>{ "max" }, j, okmax);
      if (okname && okmin && okmax) varlist[name] = { min,max };
    }
  }

  return varlist;
}

std::string ChiOrganizer::make_id(const std::string & fitvarname, double min, double max, const std::map<std::string, std::pair<double, double>> bin_varnames_borders)
{
  std::stringstream ss;
  ss << fitvarname << range_separator << min << range_separator << max;
  for (const auto &bin : bin_varnames_borders) ss << var_separator << bin.first << range_separator << bin.second.first << range_separator << bin.second.second;

  auto s = ss.str();
  std::replace(s.begin(), s.end(), '.', decimal_point_character);

  return s;
}

std::string ChiOrganizer::make_id(const std::map<std::string, std::pair<double, double>> bin_varnames_borders)
{
  bool ok = false;

  std::cout << "ChiOrganizer: Using default variables for id: dimuon_fitvar{name,min,max}" << std::endl;

  auto fitvar = GetConfigParam<std::string>(std::vector < std::string> { "dimuon_fitvar", "name" }, "NODIMUONFITVARNAME");
  double min = GetConfigParam<double>(std::vector < std::string> { "dimuon_fitvar", "min" }, 0);
  double max = GetConfigParam<double>(std::vector < std::string>{ "dimuon_fitvar", "max" }, 0);

  return make_id(fitvar, min, max, bin_varnames_borders);
}

std::string ChiOrganizer::prepare_directory(std::string dir)
{
  if (dir.empty()) dir = get_current_workingdir();
  if (dir.at(dir.size() - 1) != '/') dir += "/";
  if (is_directory) return dir;

  is_ok = false;
  return "";
}

std::string ChiOrganizer::check_file(std::string file, const std::string & default_val)
{
  if (file.empty()) file = default_val;

  if (!file_exists(file)) {
    std::cout << "File '" << file << "' not found." << std::endl;
    is_ok = false;
    return "";
  }

  return file;
}