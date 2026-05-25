// iteration_trace.cc — TSV trace output (Minuit2 backend)
#include "TMinimization.h"

#include "Math/Functor.h"

#include <cmath>
#include <iostream>

double oscillation_chi2(const double* par) {
  double theta = par[0];
  double dm2 = par[1];
  double model = theta * std::sin(dm2 * 10.0);
  return (model - 0.3) * (model - 0.3) / (0.05 * 0.05);
}

int main() {
  TMini<> mini(2);
  mini.SetParameter(0, "sin2_2theta13", 0.08, 0.001);
  mini.SetParLimits(0, 0.0, 0.2);
  mini.SetParameter(1, "dm2", 2.5, 0.01);
  mini.SetParLimits(1, 2.0, 3.0);

  ROOT::Math::Functor f(&oscillation_chi2, 2);
  mini.SetFunction(f);

  mini.EnableIterationTrace("iteration_trace.tsv", 0, 1);
  mini.Fit(2);
  mini.DisableIterationTrace();

  std::cout << "Trace written to iteration_trace.tsv" << std::endl;
  std::cout << "Best: theta=" << mini.GetParameter(0) << " dm2=" << mini.GetParameter(1)
            << std::endl;
  return 0;
}
