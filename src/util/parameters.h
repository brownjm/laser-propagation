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

#ifndef PARAMETERS_H_
#define PARAMETERS_H_

#include <string>
#include <map>
#include <vector>
#include <sstream>
#include <fstream>
#include <stdexcept>
#include <iostream>

namespace Parameters {

// general exception
struct ParametersFileError : std::runtime_error {
  ParametersFileError(const std::string& message)
    :std::runtime_error(message) {}
};

// exception to be thrown when a map is missing a key
struct KeyError : std::runtime_error {
  KeyError(const std::string& key)
    :std::runtime_error("Could not find key: '" + key + "'") {}
};


class Parameters {
public:
  Parameters();
  Parameters(const Parameters& p);
  Parameters(const std::string& filename); // loads file upon instantiation
  virtual ~Parameters();

  void load(const std::string& filename);  // add these parameters to current ones
  void save(const std::string& filename);  // save current parameters to a file
  
  // check if key exists
  bool key_exists(const std::string& key) {
    auto iter = parameters.find(key);
    return iter != parameters.end();
  }

  bool section_exists(const std::string& section_name) {
    auto section_map = getSectionMap(section_name);
    return !section_map.empty();
  }
  
  // get the value of a key
  template <class T>
  void get(const std::string& key, T& value);

  template <class T>
  T get(const std::string& key) {
    T value;
    get(key, value);
    return value;
  }
  
  // set a key to a new value
  template <class T>
  void set(const std::string& key, const T& value);

  // print out the current parameters
  void print(std::ostream& os);

  // get all the entries within a section
  Parameters getSection(const std::string& sectionName);
  std::map<std::string, std::string> getSectionMap(const std::string& sectionName);

  // iterator over all key-value pairs
  class iterator {
  public:
    iterator(std::map<std::string, std::string>::iterator it)
      :current(it) {}
    iterator& operator++() {++current; return *this;}
    iterator& operator--() {--current; return *this;}
    std::map<std::string, std::string>::iterator& operator->() {return current;}
    std::map<std::string, std::string>::value_type& operator*() {return *current;}
    bool operator==(const iterator& other) const {return current == other.current;}
    bool operator!=(const iterator& other) const {return current != other.current;}

  private:
    std::map<std::string, std::string>::iterator current;
  };
  iterator begin() {return iterator(parameters.begin());}
  iterator end()   {return iterator(parameters.end());}


private:
  std::map<std::string, std::string> parameters;
  std::stringstream ss; // used for type conversions
  
  // removes leading and trailing whitespace
  std::string trimWhitespace(const std::string& str);
};


template <class T>
void Parameters::get(const std::string& key, T& value) {
  std::map<std::string, std::string>::iterator iter = parameters.find(key);
  if (iter != parameters.end()) { // key exists
    ss << iter->second;    // send string value into stream
    ss >> value;           // read out value into proper type
    ss.str(std::string()); // set contents to an empty string
    ss.clear();            // clear any errors (most likely eof)
  }
  else throw KeyError(key);
}

  
template <class T>
void Parameters::set(const std::string& key, const T& value) {
  ss << value;                // read in value
  parameters[key] = ss.str(); // save string value for key
  ss.str(std::string());      // set contents to an empty string
  ss.clear();                 // clear any errors (most likely eof)
}
  
} // end namespace Parameters

#endif // PARAMETERS_H_
