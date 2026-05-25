// error_budget.cc — Virtual InitializeParams + error budget
#include "TMinimization.h"

#include "Math/Functor.h"

#include <cmath>
#include <iostream>

static double sigma_a = 0.3, sigma_b = 0.5;

double chi2_with_syst(const double* par) {
  double mu = par[0];
  double eps_a = par[1];
  double eps_b = par[2];
  double chi2 = (mu - 5.0 + eps_a * sigma_a + eps_b * sigma_b);
  chi2 = chi2 * chi2 / (0.2 * 0.2);
  chi2 += eps_a * eps_a + eps_b * eps_b;
  return chi2;
}

// Custom minimizer with application-specific parameter initialization
class MyMinimizer : public TMini<MinimizerBackend::Minuit2> {
public:
  using TMini::TMini;

  unsigned int InitializeParams() override {
    unsigned int idx = 0;
    AddParameterInfo(idx, "mu", 4.5, 0.01, "stat");
    AddParameterInfo(idx, "eps_a", 0.0, 0.01, "syst_a");
    AddParameterInfo(idx, "eps_b", 0.0, 0.01, "syst_b");
    paras_groups["stat"].insert(0);
    return idx;
  }
};

int main() {
  MyMinimizer mini(-1);
  mini.InitializeParams();

  ROOT::Math::Functor f(&chi2_with_syst, 3);
  mini.SetFunction(f);
  mini.Fit();

  std::vector<double> best = mini.GetParametersVector();
  std::cout << "Best fit: mu=" << best[0] << std::endl;

  mini.ExportErrorBudget(best, "add", "error_budget.txt");
  return 0;
}
