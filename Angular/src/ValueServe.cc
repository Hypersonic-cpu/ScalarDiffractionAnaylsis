#include "ValueServe.h"

#include <algorithm>
#include <vector>
#include <print>
#include <cmath>
#include <numeric>
#include <limits>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_integration.h>

using std::vector;

double 
ParamValue::BesselProdX(double x, void* deg) {
  auto n = *(int *) deg;
  if (std::log10(std::fabs(x)) * n < -37) { return 0; }
  return x * gsl_sf_bessel_Jn(n, x);
}

double 
ParamValue::BesselIntegral(int n, double xstrt, double xend) {
  gsl_integration_workspace * w
    = gsl_integration_workspace_alloc (1000);
  double result, error;

  gsl_function F { &ParamValue::BesselProdX, &n };
  gsl_integration_qags (&F, xstrt, xend, 1e-4, 0, 1000,
                        w, &result, &error);
  return result;
}

/// Int x * Jn^2(x) over [0, m-th root (mu_m)].
double 
ParamValue::BesselEigenSquare(int n, double mu_m) {
  double besselDerivative = 0;
  if (n == 0) {
    besselDerivative = gsl_sf_bessel_J1(mu_m);
  } else {
    besselDerivative = 
      (gsl_sf_bessel_Jn(n+1, mu_m) - gsl_sf_bessel_Jn(n-1, mu_m))
      * 0.5;
  }
  return 0.5 * besselDerivative * besselDerivative;
}

/// length of v SHOULD be 2N.
void
ParamValue::ThetaIntegral(vector< std::pair<double, double> > const & interval,
                          vector<double>& v) {
  auto N = v.size() >> 1;
  auto cosInt = [](int n, double x) { if (n==0) return x; return std::sin(n*x) / n; };
  auto sinInt = [](int n, double x) { if (n==0) return 0.0; return -std::cos(n*x) / n; };
  for (auto i = 0; i < N; ++i) {
    v[i] = 0;
    for (auto [beg_, end_]: interval) {
      v[i] += cosInt(i, end_) - cosInt(i, beg_);
    }
  }
  for (auto i = 0; i < N; ++i) {
    v[i+N] = 0;
    for (auto [beg_, end_]: interval) {
      v[i+N] += sinInt(i, end_) - sinInt(i, beg_);
    }
  }
} 

void 
ParamValue::ThetaEigenSquare(vector<double>& v) {
  std::fill(v.begin(), v.end(), M_PI);
  v[0] = M_PI * 2;
}

void 
FieldGenerate::FieldThetaDirArray(const double rho, const double L, const double z,
           const double WaveNum, 
           vector< double > const & theta, 
           vector< vector<double> > const & BesselRoots,
           vector< vector<double> > const & Params,
           vector<gsl_complex>& ret) {
  /* Param
   * .first  - cos: 1,  cos 1x,   cos 2x,    cos(N-1)x,
   * .second - sin: 0,  sin 1x,   sin 2x,    sin(N-1)x.
   */
  auto N = Params.size() >> 1;
  auto M = BesselRoots[0].size();
  auto mat = vector< vector<gsl_complex> > (N);
  for (auto n = 0; n < N; ++n) {
    // std::println("n = {}", n);
    mat[n].resize(M);
    std::transform(BesselRoots[n].begin(), BesselRoots[n].end(),
                   mat[n].begin(),
                   [&n, &rho, &L, &z, &WaveNum] (double mu) 
                     { 
                       auto k = std::sqrt(WaveNum * WaveNum - std::pow(mu/L, 2));
                       auto xj = mu * rho / L;
                       if (std::log10(std::fabs(xj)) * n < -37) { return gsl_complex_rect(0,0); }
                       return gsl_complex_mul_real(
                                                   gsl_complex_polar(1, -k * z),
                                                   gsl_sf_bessel_Jn(n, mu * rho / L)
                                                  );
                     }
                  );
  }
  auto mul = vector<gsl_complex> (0);
  auto zero = gsl_complex_rect(0,0);
  mul.reserve(2*N);
  for (auto n = 0; n < N; ++n) {
    mul.emplace_back(
      std::inner_product(mat[n].begin(), mat[n].end(), 
                         Params[n].begin(), zero,
                         [](gsl_complex z_, gsl_complex y_) { return gsl_complex_add(z_, y_); },
                         [](gsl_complex z_, double r_) { return gsl_complex_mul_real(z_, r_); }
                        )
      );
  }
  for (auto n = 0; n < N; ++n) {
    mul.emplace_back(
      std::inner_product(mat[n].begin(), mat[n].end(), 
                         Params[n+N].begin(), zero,
                         [](gsl_complex z_, gsl_complex y_) { return gsl_complex_add(z_, y_); },
                         [](gsl_complex z_, double r_) { return gsl_complex_mul_real(z_, r_); }
                        )
      );
  }

  ret.clear();
  for (auto th: theta) {
    auto thetaVec = vector<double> (0);
    thetaVec.reserve(2*N);
    for (auto n = 0; n < N; ++n) { thetaVec.emplace_back(std::cos(n*th)); }
    for (auto n = 0; n < N; ++n) { thetaVec.emplace_back(std::sin(n*th)); }
    ret.emplace_back(
                     std::inner_product(mul.begin(), mul.end(),
                                        thetaVec.begin(), zero,
                                        [](gsl_complex z_, gsl_complex y_) { return gsl_complex_add(z_, y_); },
                                        [](gsl_complex z_, double r_) { return gsl_complex_mul_real(z_, r_); }
                                       )
                    );
  }
}

