#pragma once

#include <vector>
#include <algorithm>
#include <gsl/gsl_complex.h>

using std::vector;

class ParamValue {

    /// Integrate x * J_n(x) in [xstrt, xend].
    static inline double 
    BesselIntegral(int n, double xstrt, double xend); // TODO

    /// (Eigenfunction).*2 * x on [0, mu_m]
    /// mu_m is the m-th positive root of BesselJ_n.
    static inline double 
    BesselEigenSquare(int n, double mu_m);

    /// Integrate Theta_{+n} = cos nx or Theta_{-n} = sin dx.
    /// Vector length is 2n+1;
    static inline void 
    ThetaEigenIntegral(vector<double>& v);

    /// (Eigenfunction).^2 on [0, 2*pi];
    static inline void
    ThetaEigenSquare(vector<double>& v);

};

class FieldGenerate {
    /// <param name="ret"> vec for theta = 0..pi </param>
    static void 
    FieldThetaDirArray(const double rho, const double L, const double z,
                                      const double WaveNum, 
                                      vector< double > const & theta, 
                                      vector< vector<double> > const & BesselRoots,
                                      vector< vector<double> > const & Params,
                                      vector<gsl_complex>& ret);
};
