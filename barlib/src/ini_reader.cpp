/*
 * Barcode
 * Copyright E.G.P. Bos
 *
 * Distributed under the terms of the MIT License.
 * The full license is in the file LICENSE, distributed with this software.
 */

#include <fstream>
#include <algorithm>  // remove_if
#include <iostream>

#include <ini_reader.hpp>

parameter_inifile::parameter_inifile(const std::string& filename)  {
  std::ifstream cFile(filename);

  if (cFile.is_open()) {
    std::string line;
    while (getline(cFile, line)) {
      // remove spaces
      line.erase(std::remove_if(line.begin(), line.end(), isspace), line.end());
      // skip comment or empty lines
      if(line[0] == '#' || line.empty()) {
        continue;
      }
      // trim trailing comments
      auto trailing_comment_ix = line.find('#');
      if (trailing_comment_ix != std::string::npos) {
        line.erase(line.begin() + static_cast<std::ptrdiff_t>(trailing_comment_ix), line.end());
      }
      // split key and value and store in output map
      auto delimiterPos = line.find('=');
      auto key = line.substr(0, delimiterPos);
      auto value = line.substr(delimiterPos + 1);
      parameters[key] = value;
    }
  }
  else {
    std::cerr << "Couldn't open config file " << filename << " for reading.\n";
  }
}
