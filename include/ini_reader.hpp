/*
 * Barcode
 * Copyright E.G.P. Bos
 *
 * Distributed under the terms of the MIT License.
 * The full license is in the file LICENSE, distributed with this software.
 */

#ifndef BARCODE_INI_READER_HPP
#define BARCODE_INI_READER_HPP

#include <map>
#include <string>
#include <sstream>

class parameter_inifile {
  std::map<std::string, std::string> parameters;
 public:
  parameter_inifile(const std::string& filename);

  template <typename T> T find(const std::string& key) {
    T value;
    std::stringstream ss(parameters[key]);
    ss >> std::boolalpha >> value;
    return value;
  }
};

#endif //BARCODE_INI_READER_HPP
