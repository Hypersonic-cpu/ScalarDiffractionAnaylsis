#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cmath>
#include <chrono>
#include <cstdio>
#include <print>

#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>

#include "LoadData.h"
#include "ValueServe.h"

using std::vector;
using std::print;
using std::println;
using std::cout;
using std::cin;
using std::endl;

constexpr int Terms = 5000; // Bessel root num m in mu^(n)_m
constexpr int FuncN = 100;   // Bessel func num n in mu(n)_m
vector< vector<double> > BesselRoots (FuncN);
// vector< double > phi (Terms, 0);

constexpr int Grids = 5001; // Plotting R-dir Total Grids.
constexpr int Plots = 100;  // Plotting R-dir calculated Grids.
constexpr int Layers = 2;// 201; // Plotting Z-dir 
vector< double > dx;
vector< double > dz;
vector< vector<gsl_complex> > val (Layers);

double const RadiusB = 0.10e-3; // 100 um 
double const RadiusA = 0.04e-3;
double const RadiusL = 99e-3; // 50 mm
double const WaveLen = 1000.0e-9;
double const WaveNumK= 2.0 * M_PI / WaveLen;

void 
CalcRhoAtZ0(const vector<double>& dx, vector<gsl_complex>& ret, double z0) {
  const auto Len = dx.size();
  for (auto n = 0; n < Len; ++n) {
    auto x = dx[n];
    auto res = gsl_complex_rect(0,0);
    for (auto i = 0; i < Terms; ++i) {
      auto lambdaI = std::pow(BesselRoots[0][i] / RadiusL, 2);
      auto waveKI  = std::sqrt(WaveNumK * WaveNumK - lambdaI);
      auto mikz = gsl_complex_polar(1,-waveKI * z0);
      auto coef = phi[i] * gsl_sf_bessel_J0(BesselRoots[0][i] / RadiusL * x);
      res = gsl_complex_add(
                            res, 
                            gsl_complex_mul_real(gsl_complex_exp(mikz), coef)
                           );
    }
    ret[n] = res;
  }
}
/* VALUE ACCUMULATION END */ 

int
main (int argc, char* argv[])
{
  if (argc == 1) {
    print("{}", "Usage:\n"
          "  argv[1]\tZ-direction distance in SI unit\n"
          "  argv[2]\tPathPrefix,\n"
          "         \tthe file will be {OutputPrefix}-val|dz|dx.csv\n"
         );
    std::exit(0);
  }
  std::string const OutputCsvPrefix (argv[2]);
  double Z0 = std::stod(argv[1]);
  println("Z0\t{:.6f}\n", Z0);

  auto tStart = std::chrono::high_resolution_clock::now();

  /* Generate Roots of BesselJn */ 
  {
    for (auto nu = 0; nu < FuncN; ++nu) {
      BesselRoots[nu].clear();
      BesselRoots[nu].reserve(Terms+2);
      for (unsigned m = 1; m < Terms+1; ++m) {
        BesselRoots[nu].emplace_back(gsl_sf_bessel_zero_Jnu(nu, m));
      }
    }
  }

  auto tRootOver = std::chrono::high_resolution_clock::now();
  auto durationRoot = std::chrono::duration_cast<std::chrono::milliseconds>(tRootOver - tStart);
  println("Root calc over: {}us", durationRoot.count());

  /* Calc param Psi[n,m] */
  auto thetaInt = vector<double> (FuncN * 2 + 1);
  ParamValue::ThetaEigenIntegral(thetaInt);
  auto thetaEig = vector<double> (FuncN * 2 + 1);
  ParamValue::ThetaEigenSquare(thetaEig);

  auto besselInt = vector< vector<double> > (FuncN);
  for (auto n = 0; n < FuncN; ++n) {
    besselInt[n].clear();
    besselInt[n].reserve(Terms);
    for (auto m = 0; m < Terms; ++m) {
      besselInt[n].emplace_back(
                                std::pow(BesselRoots[n][m], X)
                                * ParamValue::BesselIntegral(n, RadiusB * xxx, RadiusB * xxx)
                               );
      // TODO  
    }
  }

  auto besselEig = vector< vector<double> > (FuncN);
  for (auto n = 0; n < FuncN; ++n) {
    besselInt[n].clear();
    besselInt[n].reserve(Terms);
    for (auto rt: BesselRoots[n]) {
      besselInt[n].emplace_back(
                                ParamValue::BesselEigenSquare(n, rt)
                                //TODO
                               );
    }
  }

  auto tPhiParam = std::chrono::high_resolution_clock::now();
  auto durationParam= std::chrono::duration_cast<std::chrono::milliseconds>(tPhiParam - tRootOver);
  println("Param calc done: {}us", durationParam.count());

  /* Grid & Layer Setup. */
  {
    dx.clear();
    dx.reserve(Plots);
    for (auto i = 0; i < Plots; ++i) {
      dx.push_back(RadiusL * i / (Grids-1));
    }
  }
  {
    dz.clear();
    dz.reserve(Layers);
    for (auto i = 0; i < Layers; ++i) {
      dz.push_back(Z0 * i / (Layers-1));
    }
  }

  println("Grids and Layers setup.");

  /* Calc Field */
  for (int i = 0; i < Layers; i++) {
    gsl_complex zero;
    GSL_SET_COMPLEX(&zero, 0, 0);
    val[i].assign(Terms, zero);
    auto z= Z0 * i / (Layers-1); 
    printf("\tCalc Layers %d\tat Z = %.8f\n", i, z);
    CalcRhoAtZ0(dx, val[i], z);
  }
  println("Layers Complete.");

  auto tEval = std::chrono::high_resolution_clock::now();

  auto writer = CSVWriter(OutputCsvPrefix);
  writer.WriteVector(dz, "dz-test.csv");
  writer.WriteVector(dx, "dx-test.csv");
  writer.WriteMatrix(val, "val-test.csv");

  auto durationEval = std::chrono::duration_cast<std::chrono::milliseconds>(tEval - tPhiParam);
  printf("Time used:\nRoot Calc\t%12lldms\nParam Calc\t%12lldms\nEval Vals\t%12lldms\n",
         durationRoot.count(),
         durationParam.count(),
         durationEval.count()
        );
  return 0;
}
