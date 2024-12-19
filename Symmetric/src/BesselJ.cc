#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cmath>
#include <chrono>
#include <cstdio>
#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>

#include "LoadData.h"

using std::vector;
using std::cout;
using std::cin;
using std::endl;

constexpr int N = 1;
constexpr int Terms = 100; // 100000;
std::string const OutputCsvPrefix = "/Users/kong/Documents/Proj/NumericalC/DataGen/CylindarTest/";
vector< vector<double> > BesselRoots (N);
vector< double > phi (Terms, 0);

constexpr int Grids = 5001;
constexpr int Plots = 100;
constexpr int Layers = 201;
vector< double > dx;
vector< double > dz;
vector< vector<gsl_complex> > val (Layers);

double const RadiusB = 0.10e-3; // 100 um 
double const RadiusA = 0.04e-3;
double const RadiusL = 99e-3; // 50 mm
double const WaveLen = 1000.0e-9;
double const WaveNumK= 2.0 * M_PI / WaveLen;

double PhiCalc (double root0M) {
  double J1a = gsl_sf_bessel_J1(root0M * RadiusA / RadiusL);
  double J1b = gsl_sf_bessel_J1(root0M * RadiusB / RadiusL);
  double J1  = gsl_sf_bessel_J1(root0M);
  double num = 
    (
      RadiusB / RadiusL * J1b
    - RadiusA / RadiusL * J1a
    ) / root0M;
  double den = std::pow(J1, 2) / 2; 
  return num / den;
}
    
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

int
main (int argc, char* argv[])
{
  double Z0 = std::stod(argv[1]);
  printf("Z0\t%.8f\n", Z0);

  auto tStart = std::chrono::high_resolution_clock::now();

  /* Read Roots of BesselJn */ 
  {
    // std::ifstream ifs (BesselRootTSVPath);
    // LoadBesselRoots(ifs, BesselRoots, N, Terms);
    // ifs.close();
    for (auto nu = 0; nu < N; ++nu) {
      BesselRoots[nu].clear();
      BesselRoots[nu].reserve(Terms+2);
      if (nu == 0) {
        for (unsigned m = 1; m < Terms+1; ++m) {
          BesselRoots[nu].push_back(gsl_sf_bessel_zero_J0(m));
        }
      }
    }
    // for (unsigned m = 0; m < 10; ++m) {
    //   printf("%.8f ", BesselRoots[0][m]);
    // }
    // printf("\n");
  }

  auto tRootOver = std::chrono::high_resolution_clock::now();

  for (auto i = 0; i < Terms; ++i) {
    phi[i] = PhiCalc(BesselRoots[0][i]);
  }

  auto tPhiParam = std::chrono::high_resolution_clock::now();

  printf("Phi[:]\t");
  for (auto i = 0; i < Terms; i += Terms / 25) {
    printf("%.3e ", phi[i]);
  }
  printf("\n");

  {
    dx.clear();
    dx.reserve(Plots);
    for (auto i = 0; i < Plots; i += 1) {
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

  printf("Grid Generate Complete.\n");

  for (int i = 0; i < Layers; i++) {
    gsl_complex zero;
    GSL_SET_COMPLEX(&zero, 0, 0);
    val[i].assign(Terms, zero);
    auto z= Z0 * i / (Layers-1); 
    printf("\tCalc Layers %d\tat Z = %.8f\n", i, z);
    CalcRhoAtZ0(dx, val[i], z);
  }
  printf("Layers Complete.\n");

  auto tEval = std::chrono::high_resolution_clock::now();

  auto writer = CSVWriter(OutputCsvPrefix);
  writer.WriteVector(dz, "dz-test.csv");
  writer.WriteVector(dx, "dx-test.csv");
  writer.WriteMatrix(val, "val-test.csv");
  
  auto durationRoot = std::chrono::duration_cast<std::chrono::milliseconds>(tRootOver - tStart);
  auto durationParam= std::chrono::duration_cast<std::chrono::milliseconds>(tPhiParam - tRootOver);
  auto durationEval = std::chrono::duration_cast<std::chrono::milliseconds>(tEval - tPhiParam);
  printf("Time used:\nRoot Calc\t%12lldms\nParam Calc\t%12lldms\nEval Vals\t%12lldms\n",
         durationRoot.count(),
         durationParam.count(),
         durationEval.count()
        );
  return 0;
}

