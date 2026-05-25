// multi_start.cc — Global optimization with random restarts
#include "TMinimization.h"

#include "Math/Functor.h"

#include <cmath>
#include <iostream>

// Rastrigin function (many local minima)
double rastrigin(const double* par) {
  double x = par[0], y = par[1];
  return 20.0 + x * x - 10.0 * std::cos(2.0 * M_PI * x) + y * y -
         10.0 * std::cos(2.0 * M_PI * y);
}

int main() {
  TMini<> mini(-1);
  mini.SetParameter(0, "x", 3.0, 0.1);
  mini.SetParLimits(0, -5.12, 5.12);
  mini.SetParameter(1, "y", 3.0, 0.1);
  mini.SetParLimits(1, -5.12, 5.12);

  ROOT::Math::Functor f(&rastrigin, 2);
  mini.SetFunction(f);

  MultiStartConfig config;
  config.n_starts = 20;
  config.seed = 42;
  config.par_ranges = {{-5.12, 5.12}, {-5.12, 5.12}};
  config.algo = Algorithm::Combined;

  FitResult result = mini.MultiStartFit(config);
  result.Print();

  std::cout << "\nGlobal minimum should be at (0, 0) = 0" << std::endl;
  return (std::abs(result.Parameter(0)) < 0.5 && std::abs(result.Parameter(1)) < 0.5) ? 0 : 1;
}
