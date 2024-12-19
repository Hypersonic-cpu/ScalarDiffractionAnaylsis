#pragma once

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

  bool WriteVector(vector<R>, string fileName);
  bool WriteVector(vector<C>, string fileName);

  bool WriteMatrix(vector<vector<R>>, string fileName);
  bool WriteMatrix(vector<vector<C>>, string fileName);

};

