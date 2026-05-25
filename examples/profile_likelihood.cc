// profile_likelihood.cc — Profile scan + confidence intervals
#include "TMinimization.h"

#include "Math/Functor.h"

#include <cmath>
#include <iostream>

double quadratic_chi2(const double* par) {
  double a = par[0] - 2.5;
  double b = par[1] - 7.0;
  return a * a / (0.3 * 0.3) + b * b / (0.8 * 0.8);
}

int main() {
  TMini<> mini(-1);
  mini.SetParameter(0, "alpha", 2.0, 0.01);
  mini.SetParameter(1, "beta", 6.5, 0.01);

  ROOT::Math::Functor f(&quadratic_chi2, 2);
  mini.SetFunction(f);
  mini.Fit();

  ProfileResult prof = mini.ProfileLikelihood(0, 50, 1.5, 3.5);
  prof.SaveTSV("profile_alpha.tsv");

  std::cout << "Best fit: " << prof.best_fit << std::endl;
  std::cout << "1-sigma: [" << prof.lower_1sigma << ", " << prof.upper_1sigma << "]" << std::endl;
  std::cout << "2-sigma: [" << prof.lower_2sigma << ", " << prof.upper_2sigma << "]" << std::endl;

  // MINOS comparison
  auto [err_lo, err_hi] = mini.ConfidenceInterval(0, 1.0);
  std::cout << "MINOS: " << err_lo << " +" << err_hi << std::endl;

  return 0;
}
