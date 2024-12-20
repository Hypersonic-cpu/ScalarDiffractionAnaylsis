#pragma once

#include <vector>

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

    /// <summary> Generate all terms of the seires. </summary>
    /// <param name="v"> vecOf</param>
    static void 
    TermsMatrixAtZ0 (vector<double> const& BesselRoots,
                     vector<double> const& ThetaParams,
                     vector< vector<double> > const& EigenParams,
                     vector< vector< vector<double> > >& v);
};