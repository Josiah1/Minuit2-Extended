// grid_scan.cc — 2D grid scan -> TH2D
#include "TMinimization.h"

#include "Math/Functor.h"
#include "TFile.h"

#include <cmath>
#include <iostream>

double banana_chi2(const double* par) {
  double x = par[0], y = par[1];
  return (1.0 - x) * (1.0 - x) + 10.0 * (y - x * x) * (y - x * x);
}

int main() {
  TMini<> mini(-1);
  mini.SetParameter(0, "x", 0.5, 0.01);
  mini.SetParameter(1, "y", 0.5, 0.01);

  ROOT::Math::Functor f(&banana_chi2, 2);
  mini.SetFunction(f);
  mini.Fit();

  GridScanResult grid = mini.GridScan(0, 1, 40, 40, -0.5, 2.5, -1.0, 4.0);
  grid.SaveTSV("grid_scan.tsv");

  TFile output("grid_scan_output.root", "RECREATE");
  TH2D* h = grid.ToTH2D("chi2_grid");
  if (h) {
    h->Write();
    delete h;
  }

  std::cout << "Grid scan: best=(" << grid.best_x << ", " << grid.best_y
            << ") chi2_min=" << grid.chi2_min << std::endl;
  return 0;
}
