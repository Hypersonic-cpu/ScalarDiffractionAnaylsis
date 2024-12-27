#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cmath>
#include <chrono>
#include <cassert>
#include <print>

#include <omp.h>

#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>

#include "Loader.hpp"
#include "ValueServe.h"

using std::vector;
using std::print;
using std::println;
using std::cout;
using std::cin;
using std::endl;

constexpr int Terms = 5000; // Bessel root num m in mu^(n)_m
constexpr int FuncN = 121;   // Bessel func num n in mu(n)_m
vector< vector<double> > BesselRoots (FuncN);
vector< vector<double> > param (FuncN << 1);

constexpr int Grids = 4001; // Plotting R-dir Total Grids.
constexpr int Plots = 200;  // Plotting R-dir calculated Grids.
constexpr int Layers = 5;// 201; // Plotting Z-dir 
vector< double > dRho;
vector< double > dTheta;
vector< double > dz;

double const RadiusB = 0.10e-3; // 100 um 
double const RadiusA = 0.01e-3;
double const RadiusL = 99e-3; // 99 mm
double const Z0      = 10e-3; // 10 mm
double const WaveLen = 1000.0e-9;
double const WaveNumK= 2.0 * M_PI / WaveLen;

inline auto
PrintHelpExit() {
    print("{}", "Usage:\n"
          "  argv[1]\tGraph code,\n"
          "  argv[2]\tPathPrefix,\n"
          "         \tthe file will be ../DataGen/{preifx}-mtx|dz|dx.csv\n"
         );
    std::exit(0);
}

inline auto 
CalcBesselRoots() {
  for (auto nu = 0; nu < FuncN; ++nu) {
    BesselRoots[nu].clear();
    BesselRoots[nu].reserve(Terms+2);
    for (unsigned m = 1; m < Terms+1; ++m) {
      BesselRoots[nu].emplace_back(gsl_sf_bessel_zero_Jnu(nu, m));
    }
  }
}

auto const HalfWidth = M_PI * 2.5 / 180;
auto const IntervalSet = std::vector < std::vector<std::pair<R,R> > >
{
  {
    std::make_pair(0.0 * M_PI + HalfWidth, 1.0 * M_PI - HalfWidth),
    std::make_pair(1.0 * M_PI + HalfWidth, 2.0 * M_PI - HalfWidth),
  },
  {
    std::make_pair(0.0 * M_PI + HalfWidth, 0.5 * M_PI - HalfWidth),
    std::make_pair(0.5 * M_PI + HalfWidth, 1.0 * M_PI - HalfWidth),
    std::make_pair(1.0 * M_PI + HalfWidth, 1.5 * M_PI - HalfWidth),
    std::make_pair(1.5 * M_PI + HalfWidth, 2.0 * M_PI - HalfWidth),
  },
  {
    std::make_pair(0.0 * M_PI + HalfWidth, 0.75 * M_PI - HalfWidth),
    std::make_pair(0.75 * M_PI + HalfWidth, 1.25 * M_PI - HalfWidth),
    std::make_pair(1.25 * M_PI + HalfWidth, 2.0 * M_PI - HalfWidth),
  },
};

inline auto 
CalcParamAll(int graphCode) {
  auto thetaInt = vector<double> (FuncN * 2);
  ParamValue::ThetaIntegral(IntervalSet[graphCode], thetaInt);
  auto thetaEig = vector<double> (FuncN * 2);
  ParamValue::ThetaEigenSquare(thetaEig);
  auto besselInt = vector< vector<double> > (FuncN);
  for (auto n = 0; n < FuncN; ++n) {
    besselInt[n].clear();
    besselInt[n].reserve(Terms);
    for (auto rtnm: BesselRoots[n]) {
      besselInt[n].emplace_back(
                                ParamValue::BesselIntegral(n,
                                                                 RadiusA/RadiusL*rtnm, 
                                                                 RadiusB/RadiusL*rtnm)
                               / std::pow(rtnm, 2)
                               );
    }
  }

  auto besselEig = vector< vector<double> > (FuncN);
  for (auto n = 0; n < FuncN; ++n) {
    besselEig[n].clear();
    besselEig[n].reserve(Terms);
    for (auto rt: BesselRoots[n]) {
      besselEig[n].emplace_back(
                                ParamValue::BesselEigenSquare(n, rt) 
                               );
    }
  }

  /* Assign */
  for (auto n = 0; n < FuncN * 2; ++n) {
    auto th_ = thetaInt[n] / thetaEig[n];
    // println("{}, {}", n,  th_);
    auto usedN = n >= FuncN ? (n-FuncN) : n;
    param[n].resize(Terms);
    std::transform(
                   besselInt[usedN].begin(),
                   besselInt[usedN].end(),
                   besselEig[usedN].begin(), 
                   param[n].begin(), 
                   [&th_](auto const & int_, auto const & eig_)
                     { return int_ / eig_ * th_; }
                  );
  }
}

int
main (int argc, char* argv[])
{
  if (argc == 1) { PrintHelpExit(); }
  assert(argc == 3);
  std::string const OutputCsvPrefix (argv[2]);
  short const GraphCode (std::atoi(argv[1]));

  auto tStart = std::chrono::high_resolution_clock::now();

  /* Generate Roots of BesselJn */ 
  CalcBesselRoots();

  auto tRootOver = std::chrono::high_resolution_clock::now();
  auto durationRoot = std::chrono::duration_cast<std::chrono::milliseconds>(tRootOver - tStart);
  println("Root calc over: {}ms", durationRoot.count());

  /* Calc param Param[n,m] */
  CalcParamAll(GraphCode);

  auto tPhiParam = std::chrono::high_resolution_clock::now();
  auto durationParam= std::chrono::duration_cast<std::chrono::milliseconds>(tPhiParam - tRootOver);
  println("Param calc done: {}ms", durationParam.count());

  /* Grid & Layer Setup. */
  {
    dRho.resize(Plots);
    for (auto i = 0; i < Plots; ++i) {
      dRho[i] = RadiusL * i / (Grids-1);
    }
  }
  {
    dz.resize(Layers);
    for (auto i = 0; i < Layers; ++i) {
      dz[i] = Z0 * i / (Layers-1);
    }
  }
  {
    dTheta.resize(FuncN);
    for (auto i = 0; i < FuncN; ++i) {
      dTheta[i] = M_PI * i / (FuncN-1);
    }
  }
  println("Grids and Layers setup.");
  
  /* Calc Field */
#pragma omp parallel for 
  for (auto zIdx = 1; zIdx < Layers; ++zIdx) {
    print("[OMP  ] Z = {:.8g} by thread #{}", dz[zIdx], omp_get_thread_num());
    vector< vector<gsl_complex> > field;
    field.clear();
    field.resize(Plots);
    for (auto idxRho = 0; idxRho < Plots; ++idxRho) {
      if (idxRho % 10 == 0) { println("[RHO  ] processing #{}", idxRho); }
      FieldGenerate::FieldThetaDirArray(dRho[idxRho], RadiusL, dz[zIdx],
                                        WaveNumK,
                                        dTheta, 
                                        BesselRoots,
                                        param,
                                        field[idxRho]
                                       );
    }

    println("Field computation complete.");

    auto writeFile = OutputCsvPrefix + "L" + std::to_string(zIdx) + "-";
    auto writer = 
      CSVWriter("/Users/kong/Documents/Proj/NumericalC/DataGen/" + writeFile);
    assert(writer.WriteVector({dz[zIdx]}      , "dz.csv")         );
    assert(writer.WriteVector(dTheta  , "dTheta.csv")         );
    assert(writer.WriteVector(dRho    , "dRho.csv")       );
    assert(writer.WriteMatrix(field, "fld.csv")  );
    println("Write with prefix {} finished.", OutputCsvPrefix);

  }
  auto tEval = std::chrono::high_resolution_clock::now();
  auto durationEval = std::chrono::duration_cast<std::chrono::milliseconds>(tEval - tPhiParam);
  println("[Summary] time used:\nRoot Calc\t{:12}ms\nParam Calc\t{:12}ms\nEval Vals\t{:12}ms\n",
         durationRoot.count(),
         durationParam.count(),
         durationEval.count()
        );
  return 0;
}

