#include "ValueServe.h"

#include <vector>
#include <gsl/gsl_sf_bessel.h>

using std::vector;

inline double 
ParamValue::BesselIntegral(int n, double xstrt, double xend) {
    return 0;// TODO
}

inline double 
ParamValue::BesselEigenSquare(int n, double mu_m) {
    return 0; // TODO
}

inline void
ParamValue::ThetaEigenIntegral(vector<double>& v) {
        auto n = v.size() >> 1;
        for (auto i = 0; i <= n; ++i) {
            //TODO: Cosine.
        }
        for (auto i = n+1; i < 2*n+1; ++i) {
            //TODO: Sine.
        }
} 

inline void 
ParamValue::ThetaEigenSquare(vector<double>& v) {
        std::fill(v.begin(), v.end(), M_PI);
        v[0] = M_PI * 2;
}

void
FieldGenerate::TermsMatrixAtZ0 (vector<double> const& BesselRoots,
                                vector<double> const& ThetaParams,
                                vector< vector<double> > const& EigenParams,
                                vector< vector< vector<double> > >& v) {
    int a = 1;
    // TODO
}