// basic_fit.cc — Gaussian MLE, FitResult, ParameterHandle demo
#include "TMinimization.h"

#include "Math/Functor.h"

#include <cmath>
#include <iostream>
#include <vector>

static std::vector<double> data;

double chi2_func(const double* par) {
  double mu = par[0];
  double sigma = par[1];
  double chi2 = 0.0;
  for (double x : data) {
    double pull = (x - mu) / sigma;
    chi2 += pull * pull;
  }
  chi2 += 2.0 * static_cast<double>(data.size()) * std::log(sigma);
  return chi2;
}

int main() {
  // Generate pseudo-data: 100 samples around mu=5.0, sigma=1.5
  TRandom3 rng(42);
  data.resize(100);
  for (auto& d : data)
    d = rng.Gaus(5.0, 1.5);

  TMini<> mini(2, 1e-6);

  // Fluent parameter handles
  auto h_mu = mini.AddParam("mu", 4.0, 0.1).SetLimits(0.0, 10.0).SetGroup("physics");
  auto h_sigma = mini.AddParam("sigma", 1.0, 0.1).SetLimits(0.01, 10.0).SetGroup("physics");

  ROOT::Math::Functor f(&chi2_func, 2);
  mini.SetFunction(f);

  FitResult result = mini.FitAndGetResult();
  result.Print();

  std::cout << "\nmu handle: " << h_mu.Value() << " +/- " << h_mu.Error() << std::endl;
  std::cout << "sigma handle: " << h_sigma.Value() << " +/- " << h_sigma.Error() << std::endl;

  result.SaveJSON("basic_fit_result.json");
  return (result.status == FitStatus::Success) ? 0 : 1;
}
