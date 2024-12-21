#pragma once

#include <vector>
#include <gsl/gsl_complex.h>

using std::vector;

class ParamValue {
  static double BesselProdX(double x, void* degree);
public:
  /// Integrate x * J_n(x) in [xstrt, xend].
  static double 
  BesselIntegral(int n, double xstrt, double xend);

  /// (J_n(x)).*2 * x on [0, mu_m]
  /// mu_m is the m-th positive root of BesselJ_n.
  static double 
  BesselEigenSquare(int n, double mu_m);

  /// Vector length MUST BE 2N;
  static void 
  ThetaIntegral(vector< std::pair<double, double> > const & interval,
                vector<double>& v);

  /// Compute inegral (sin|cos nx).^2 on [0, 2pi]. n=0..N-1,
  /// Vector length MUST BE 2N;
  static void
  ThetaEigenSquare(vector<double>& v);

};

class FieldGenerate {
public:
  /// <param name="ret"> vec for theta = 0..pi </param>
  static void 
  FieldThetaDirArray(const double rho, const double L, const double z,
                                   const double WaveNum, 
                                   vector< double > const & theta, 
                                      vector< vector<double> > const & BesselRoots,
                                      vector< vector<double> > const & Params,
                                      vector<gsl_complex>& ret);
};
