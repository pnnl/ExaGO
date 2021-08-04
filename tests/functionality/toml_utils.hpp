#pragma once
#include <toml.hpp>
#include <utils.hpp>
#include <version.hpp>
#include <exago_config.h>

/**
 * 
 * @brief Utilities for working with TOML files which describe test suites for
 * ExaGO's functionality tests.
 *
 * @author Asher Mancinelli <asher.mancinelli@pnnl.gov>
 *
 */

template<typename ValueType>
void set_if_found(ValueType& value, toml::value config, const std::string& key)
{
  if (config.contains(key) and config.count(key) > 0)
  {
    value = toml::find<ValueType>(config, key);
  }
}

void resolve_datafiles_path(std::string &path)
{
  std::vector<std::string> prefix;
  prefix.push_back("./");
  prefix.push_back( std::string(EXAGO_OPTIONS_DIR) + "/../");
  for (int i=0; i<prefix.size(); i++) {
    std::ifstream f{(prefix[i]+path).c_str()};
    if (f.is_open()) {
      path = prefix[i]+path;
      f.close();
      break;
    }
  }
}

/* For formatting reports as TOML */
void fmt_row(std::ostream& summary, int col_width, std::string key,
    std::string value)
{
  std::stringstream value_fmt;
  value_fmt << "'" << value << "'";

  summary
    << std::setw(col_width-1) << std::left << key << std::right << "="
    << std::setw(col_width-1) << std::left << value_fmt.str() << std::right
    << "\n";
}

template<typename T>
void fmt_row(std::ostream& summary, int col_width, std::string key, T value)
{
  summary
    << std::setw(col_width-1) << std::left << key << std::right << "="
    << std::setw(col_width-1) << std::left << value << std::right << "\n";
}

template<typename T>
void fmt_comment(std::ostream& summary, int col_width, std::string key, T value)
{
  summary
    << "#"
    << std::setw(col_width-2) << std::left << key << std::right << "="
    << std::setw(col_width-1) << std::left << value << std::right << "\n";
}

static std::string bool2str(bool b)
{
  return (b?"true":"false");
};
