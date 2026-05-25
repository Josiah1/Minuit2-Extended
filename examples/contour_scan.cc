// contour_scan.cc — 2D contour + 1D chi2 scan
#include "TMinimization.h"

#include "Math/Functor.h"
#include "TFile.h"
#include "TGraph.h"

#include <cmath>
#include <iostream>

double parabolic_chi2(const double* par) {
  double x = par[0] - 3.0;
  double y = par[1] - 5.0;
  return x * x / (0.5 * 0.5) + y * y / (1.0 * 1.0) + 0.3 * x * y;
}

int main() {
  TMini<> mini(2);
  mini.SetParameter(0, "x", 2.0, 0.01);
  mini.SetParameter(1, "y", 4.0, 0.01);

  ROOT::Math::Functor f(&parabolic_chi2, 2);
  mini.SetFunction(f);
  mini.Fit();

  // 1, 2, 3-sigma contours
  for (size_t sigma = 1; sigma <= 3; ++sigma)
    mini.CreateStandardContour(sigma, 100, 0, 1);

  mini.GetContourGraphs("contour_scan_output.root", 100);

  // 1D chi2 scan
  const unsigned int nstep = 50;
  std::vector<double> sx(nstep), sy(nstep);
  mini.Scan_Chi2(0, nstep, 1.0, 5.0, sx.data(), sy.data(), "chi2scan_x.txt");

  std::cout << "Contour + scan output written." << std::endl;
  return 0;
}
