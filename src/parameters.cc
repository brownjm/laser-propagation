/*
Copyright (C) 2015 Jeffrey M Brown

This program is free software: you can redistribute it and/or modify it
under the terms of the GNU General Public License as published by the
Free Software Foundation, either version 3 of the License, or (at your
option) any later version.

This program is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
General Public License for more details.

You should have received a copy of the GNU General Public License
along with This program. If not, see <http://www.gnu.org/licenses/>.
*/

#include "parameters.h"

namespace Parameters {
// implementations

Parameters::Parameters() {}

Parameters::Parameters(const Parameters& p)
  :parameters(p.parameters), ss(p.ss.str()) {

}

  
Parameters::Parameters(const std::string& filename) {
  load(filename);
}

Parameters::~Parameters() {}

void Parameters::load(const std::string& filename) {
  // open file
  std::ifstream ifs(filename.c_str());
  if (!ifs) throw ParametersFileError("Cannot open input file: " + filename);
  
  // read lines from file
  std::string line;
  std::vector<std::string> text;
  while (!ifs.eof()) {
    getline(ifs, line);

    // remove comments
    std::size_t found = line.find('#');
    if (found != std::string::npos) {
      line = line.substr(0, found);
    }

    // remove leading and trailing whitespace
    line = trimWhitespace(line);
    if (line.length() > 0) {
      text.push_back(line); // add remaining line
    }
  }

  // iterate through text and build the dictionary
  std::string sectionName, key, value, fullPathKey;
  std::vector<std::string>::const_iterator it;
  for (it = text.begin(); it != text.end(); ++it) {
    line = *it;
    // check if line is a section header
    if (*line.begin() == '[' && *line.rbegin() == ']') {
      sectionName = trimWhitespace(line.substr(1, line.length()-2));
    }
    else { // line might contain: key = value
      std::size_t found = line.find("=");
      if (found == std::string::npos) {
	std::string msg("Under section '" + sectionName + "', malformed expression line: '" + line + "'\n");
	throw ParametersFileError(msg);
      }

      // split expression by '=' into key and value
      key = trimWhitespace(line.substr(0, found));
      value = trimWhitespace(line.substr(found+1));

      // if missing key or value
      if (key.length() == 0 || value.length() == 0) {
	std::string msg("Under section '" + sectionName + "', missing key or value: '" + key + "=" + value + "'\n");
	throw ParametersFileError(msg);
      }

      // store keys as 'section/key'
      fullPathKey = sectionName + "/" + key;
      parameters[fullPathKey] = value;
    }
  }
}

// save parameters map back to a configuration file format
void Parameters::save(const std::string& filename) {
  // open file
  std::ofstream ofs(filename.c_str());
  if (!ofs) throw ParametersFileError("Cannot open output file: " + filename);
  
  // iterate through the dictionary
  std::string fullPathKey, sectionName, key, value, currentSectionName;
  std::map<std::string, std::string>::const_iterator map_iter;
  for (map_iter=parameters.begin(); map_iter!=parameters.end(); map_iter++) {
    fullPathKey = map_iter->first;
    value = map_iter->second;
    
    // split by '/'
    int i = fullPathKey.find("/");
    sectionName = fullPathKey.substr(0, i);
    key = fullPathKey.substr(i+1);

    // check if we need to start a new section
    if (currentSectionName == sectionName) { // same section still
      ofs << key << " = " << value << "\n";
    }
    else { // create a new section
      ofs << "\n" << "[" << sectionName << "]\n";
      ofs << key << " = " << value << "\n";
      currentSectionName = sectionName; // update the section we are under
    }
  }
}







// print out all of the key/value pairs of parameters
void Parameters::print(std::ostream& os) {
  os << "*** Parameters ***\n";
  for (iterator it = begin(); it != end(); ++it) {
    os << it->first << ": " << it->second << "\n";
  }
  std::cout << "\n";
}


// get a map of the entries within a section
Parameters Parameters::getSection(const std::string& sectionName) {
  // iterate through the parameters map and find all entries with matching section name
  Parameters p;
  std::map<std::string, std::string>::const_iterator it;
  std::string fullPathKey, value, section, key;
  for (it = parameters.begin(); it != parameters.end(); ++it) {
    fullPathKey = it->first;
    value = it->second;
    int i = fullPathKey.find("/");
    section = fullPathKey.substr(0, i);
    key = fullPathKey.substr(i+1);
    if (section == sectionName) {
      p.set(key, value);
    }
  }
  return p;
}
  
// get a map of the entries within a section
std::map<std::string, std::string> Parameters::getSectionMap(const std::string& sectionName) {
  // iterate through the parameters map and find all entries with matching section name
  std::map<std::string, std::string> sectionMap;
  std::map<std::string, std::string>::const_iterator it;
  std::string fullPathKey, value, section, key;
  for (it = parameters.begin(); it != parameters.end(); ++it) {
    fullPathKey = it->first;
    value = it->second;
    int i = fullPathKey.find("/");
    section = fullPathKey.substr(0, i);
    key = fullPathKey.substr(i+1);
    if (section == sectionName) {
      sectionMap[key] = value;
    }
  }
  return sectionMap;
}


// trim leading and trailing whitespace
std::string Parameters::trimWhitespace(const std::string& str) {
  // find first and last characters that are not spaces
  int first = str.find_first_not_of(" "); // returns -1 if only spaces present
  int last = str.find_last_not_of(" ");
  
  // check if we found anything other than spaces
  if (first >= 0 && last >= 0) {
    return str.substr(first, last-first+1);
  }
  else { // string only contains spaces
    return std::string("");
  }
}

} // end namespace Parameters
