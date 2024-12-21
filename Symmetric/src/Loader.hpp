#pragma once

#include <iostream>
#include <cmath>
#include <ios>
#include <sstream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <string>
#include <gsl/gsl_complex.h>

using std::vector;
using std::string;
using R = double;
using C = gsl_complex;

class CSVWriter {
public:
  CSVWriter(const char* absoluteOutputPrefix):
    AbsoluteOutputPrefix(absoluteOutputPrefix) {}
  CSVWriter(const std::string& absoluteOutputPrefix):
    AbsoluteOutputPrefix(absoluteOutputPrefix) {}

  string const AbsoluteOutputPrefix;

  void testFunction(vector<int> a) { return; }
  void testFunction(vector<R> a) { return; }

  bool WriteVector(vector<R>, string fileName);
  bool WriteVector(vector<C>, string fileName);

  bool WriteMatrix(vector<vector<R>>, string fileName);
  bool WriteMatrix(vector<vector<C>>, string fileName);

};


bool
WriteLn(vector<R> arr, std::ofstream& f) {
  if (arr.empty() || !f.is_open()) return false;
  f << std::scientific << std::setprecision(8) << arr[0];
  for (auto it = arr.begin()+1;
       it != arr.end(); it++) {
    f << "," << std::scientific << std::setprecision(8) << *(it);
  }
  f << std::endl;
  return true;
}

inline std::string 
CtoRI(C z){
  auto real = GSL_REAL(z);
  auto imag = GSL_IMAG(z);
  std::stringstream ss;
  ss << std::scientific << std::setprecision(8)
    << (real > 0 ? "+" : "-")
    << std::fabs(real)
    << (imag > 0 ? "+" : "-")
    << std::fabs(imag)
    << "i";
  return ss.str();
}

bool
WriteLn(vector<C> arr, std::ofstream& f) {
  if (arr.empty() || !f.is_open()) return false;
  f.precision(8);
  f << CtoRI(arr[0]);
  for (auto it = arr.begin()+1;
       it != arr.end(); it++) {
    f << "," << CtoRI(*(it));
  }
  f << std::endl;
  return true;
}

bool
CSVWriter::WriteVector(vector<R> arr, string fileName) {
  std::ofstream f (this->AbsoluteOutputPrefix + fileName);
  auto retv = WriteLn(arr, f);
  f.close();
  return retv;
}

bool
CSVWriter::WriteVector(vector<C> arr, string fileName) {
  std::ofstream f (this->AbsoluteOutputPrefix + fileName);
  auto retv = WriteLn(arr, f);
  f.close();
  return retv;
}

bool
CSVWriter::WriteMatrix(vector<vector<R>> arr, string fileName) {
  std::ofstream f (this->AbsoluteOutputPrefix + fileName);
  for (const auto& ln: arr) {
    if (!WriteLn(ln, f)) { f.close(); return false; }
  }
  return true;
}

bool
CSVWriter::WriteMatrix(vector<vector<C>> arr, string fileName) {
  std::ofstream f (this->AbsoluteOutputPrefix + fileName);
  for (const auto& ln: arr) {
    if (!WriteLn(ln, f)) { f.close(); return false; }
  }
  return true;
}
