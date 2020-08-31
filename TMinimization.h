#include "Minuit2/Minuit2Minimizer.h"
#include "Math/Functor.h"

class TMini
{

public:
  TMini(ROOT::Minuit2::Minuit2Minimizer* mini_)
  {
    mini = mini_;

    mini->SetPrintLevel(1);
    mini->SetStrategy(1); // 0- cursory, 1- default, 2- thorough yet no more successful
    mini->SetMaxFunctionCalls(500000);
    mini->SetMaxIterations(500000);
    mini->SetTolerance(1e-8);  // tolerance*2e-3 = edm precision
    mini->SetPrecision(1e-18); // precision in the target function

    //ROOT::Math::Functor Chi2Functor(fcn, npars);
    //mini->SetFunction(Chi2Functor);
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
        mini->SetErrorDef(2.2977); // 1 - 0.317 with ndf = 2
        break;
      case 2:
        // 95.5% confidence level
        mini->SetErrorDef(6.2021); // 1 - 0.045 with ndf = 2
        break;
      case 3:
        // 99.7% confidence level
        mini->SetErrorDef(11.6182); // 1 - 0.003 with ndf = 2
        break;
      default:
        // 68.3% confidence level
        mini->SetErrorDef(2.2977); // 1 - 0.317 with ndf = 2
    }

    mini->Contour(parI, parJ, Npoints, arrayI, arrayJ);
  }

private:
  ROOT::Minuit2::Minuit2Minimizer* mini;
};
