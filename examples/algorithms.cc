// algorithms.cc — Migrad vs Simplex vs Combined comparison
#include "TMinimization.h"

#include "Math/Functor.h"

#include <cmath>
#include <iostream>

double rosenbrock(const double* par) {
  double x = par[0], y = par[1];
  return (1.0 - x) * (1.0 - x) + 100.0 * (y - x * x) * (y - x * x);
}

int main() {
  ROOT::Math::Functor f(&rosenbrock, 2);

  Algorithm algos[] = {Algorithm::Migrad, Algorithm::Simplex, Algorithm::Combined};
  const char* names[] = {"Migrad", "Simplex", "Combined"};

  for (int i = 0; i < 3; ++i) {
    TMini<> mini(-1);
    mini.SetParameter(0, "x", -1.5, 0.1);
    mini.SetParameter(1, "y", -1.5, 0.1);
    mini.SetFunction(f);

    mini.FitWith(algos[i], -1);
    auto result = mini.GetResult();

    std::cout << names[i] << ": x=" << result.Parameter(0) << " y=" << result.Parameter(1)
              << " min=" << result.min_value << " nfcn=" << result.nfcn << std::endl;
  }
  return 0;
}
