#ifndef CHIBCHICPOLFW_UNBINNEDPOLFW_REGION_H__
#define CHIBCHICPOLFW_UNBINNEDPOLFW_REGION_H__

#include <limits>
#include <type_traits>
#include <ostream>
#include <string>
#include <sstream>

enum class Boundary {
  Low,
  High,
  TwoSided
};

template<Boundary B,
         class = typename std::enable_if<std::numeric_limits<double>::has_infinity>::type>
class Region;

template<Boundary B,
         class = typename std::enable_if<std::numeric_limits<double>::has_infinity>::type>
std::ostream& operator<<(std::ostream&, const Region<B>&);

template<Boundary B, class>
class Region {
public:
  Region(double min, double max, const std::string&& name = "");
  Region(double bound, const std::string&& name);

  double min() const { return m_min; }
  double max() const { return m_max; }
  const std::string& name() const { return m_name; }

  bool contains(const double val) const { return val > m_min && val < m_max; }

  std::string getCutStr() const { return getCutStr(m_name); }

  std::string getCutStr(const std::string&& name) const;

  friend std::ostream& operator<< <> (std::ostream& os, const Region& region);
private:
  double m_min{};
  double m_max{};
  std::string m_name{};
};

template<>
Region<Boundary::Low>::Region(double bound, const std::string&& name) :
  m_min(bound), m_max(std::numeric_limits<double>::infinity()), m_name(name) {}

template<>
Region<Boundary::High>::Region(double bound, const std::string&& name) :
  m_min(-std::numeric_limits<double>::infinity()), m_max(bound), m_name(name) {}

template<>
Region<Boundary::TwoSided>::Region(double min, double max, const std::string&& name) :
  m_min(min), m_max(max), m_name(name) {}


template<>
std::string Region<Boundary::TwoSided>::getCutStr(const std::string&& name) const
{
  std::stringstream sstr;
  sstr << "(" << name << " > " << m_min << " && " << name << " < " << m_max << ")";
  return sstr.str();
}

template<>
std::string Region<Boundary::Low>::getCutStr(const std::string&& name) const
{
  std::stringstream sstr;
  sstr << "(" << name << " > " << m_min << ")";
  return sstr.str();
}

template<>
std::string Region<Boundary::High>::getCutStr(const std::string&& name) const
{
  std::stringstream sstr;
  sstr << "(" << name << " < " << m_max << ")";
  return sstr.str();
}


template<>
std::ostream& operator<<(std::ostream& os, const Region<Boundary::TwoSided>& reg)
{
  os << reg.m_name << ": (" << reg.m_min <<  ", " << reg.m_max << ")";
  return os;
}

template<>
std::ostream& operator<<(std::ostream& os, const Region<Boundary::High>& reg)
{
  os << reg.m_name << ": (-infty, " << reg.m_min << ")";
  return os;
}

template<>
std::ostream& operator<<(std::ostream& os, const Region<Boundary::Low>& reg)
{
  os << reg.m_name << ": (" << reg.m_min << ", +infty)";
  return os;
}


#endif
