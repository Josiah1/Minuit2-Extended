// state_snapshot.cc — Save/restore RAII pattern
#include "TMinimization.h"

#include "Math/Functor.h"

#include <cmath>
#include <iostream>

double simple_chi2(const double* par) {
  return (par[0] - 3.0) * (par[0] - 3.0) + (par[1] - 7.0) * (par[1] - 7.0);
}

int main() {
  TMini<> mini(-1);
  mini.SetParameter(0, "x", 1.0, 0.1);
  mini.SetParameter(1, "y", 5.0, 0.1);

  ROOT::Math::Functor f(&simple_chi2, 2);
  mini.SetFunction(f);
  mini.Fit();

  std::cout << "After fit: x=" << mini.GetParameter(0) << " y=" << mini.GetParameter(1)
            << std::endl;

  // Save state
  auto snapshot = mini.SaveState();

  // Perturb
  mini.SetParValue(0, 100.0);
  mini.FixParameter(1, 0.0);
  std::cout << "Perturbed: x=" << mini.GetParameter(0) << " y=" << mini.GetParameter(1)
            << std::endl;

  // Restore
  mini.RestoreState(snapshot);
  std::cout << "Restored:  x=" << mini.GetParameter(0) << " y=" << mini.GetParameter(1)
            << std::endl;

  return 0;
}
