#include "Minuit2/Minuit2Minimizer.h"
#include "Math/Functor.h"

class TMini:public ROOT::Minuit2::Minuit2Minimizer
{

public:
  TMini(ROOT::Math::Functor &Chi2Functor, int level=1, double tolerance=1e-8)
      :ROOT::Minuit2::Minuit2Minimizer(ROOT::Minuit2::kMigrad)
  {
    SetPrintLevel(level);
    SetStrategy(1); // 0- cursory, 1- default, 2- thorough yet no more successful
    SetMaxFunctionCalls(500000);
    SetMaxIterations(500000);
    SetTolerance(tolerance);  // tolerance*2e-3 = edm precision
    SetPrecision(1e-18); // precision in the target function
    SetFunction(Chi2Functor);
  }
  TMini(const TMini&) {}
  ~TMini() {}

  void GetContour(int CL, unsigned int Npoints, unsigned int parI, unsigned int parJ, double arrayI[Npoints], double arrayJ[Npoints])
  {

    /// https://www.fourmilab.ch/rpkp/experiments/analysis/chiCalc.html
    switch (CL)
    {
      case 1:
        // 68.3% confidence level
        SetErrorDef(2.2977); // 1 - 0.317 with ndf = 2
        break;
      case 2:
        // 95.5% confidence level
        SetErrorDef(6.2021); // 1 - 0.045 with ndf = 2
        break;
      case 3:
        // 99.7% confidence level
        SetErrorDef(11.6182); // 1 - 0.003 with ndf = 2
        break;
      default:
        // 68.3% confidence level
        SetErrorDef(2.2977); // 1 - 0.317 with ndf = 2
    }

    Contour(parI, parJ, Npoints, arrayI, arrayJ);
  }

};
