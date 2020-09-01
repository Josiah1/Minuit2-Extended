/**
 * @file TMinimization.h
 * @brief  An easy-used minuit2 package.
 * @author Jinjing Li, josiahleeoaa@outlook.com
 * @version 1.0.0
 * @date 2020-09-01
 */

#include "Minuit2/Minuit2Minimizer.h"
#include "Math/Functor.h"

class Function
{
public:
  double operator()(const double* par);
} fcn;

class TMini : public ROOT::Minuit2::Minuit2Minimizer
{

public:
  //Default constructer
  TMini(ROOT::Math::Functor& Chi2Functor, int level = 1, double tolerance = 1e-8)
      : ROOT::Minuit2::Minuit2Minimizer(ROOT::Minuit2::kMigrad)
  {
    SetPrintLevel(level);
    SetStrategy(1); // 0- cursory, 1- default, 2- thorough yet no more successful
    SetMaxFunctionCalls(500000);
    SetMaxIterations(500000);
    SetTolerance(tolerance); // tolerance*2e-3 = edm precision
    SetPrecision(1e-18);     // precision in the target function
    SetFunction(Chi2Functor);
  }
  TMini(const TMini&) {}
  ~TMini() {}

  //Get parameter list
  const double* GetParameters() const
  {
    return X();
  }

  //Get one of the parameter
  const double GetParameter(int i) const
  {
    return X()[i];
  }

  //Get error list
  const double* GetParErrors() const
  {
    return Errors();
  }

  //Get one of the error
  const double GetParError(int i) const
  {
    return Errors()[i];
  }

  //Set Parameter
  bool SetParameter(unsigned int ivar, const std::string& name, double val, double step)
  {
    return SetVariable(ivar, name, val, step);
  }

  bool SetParLimits(unsigned int ivar, double lower, double upper)
  {
    return SetVariableLimits(ivar, lower, upper);
  }

  bool SetParValue(unsigned int ivar, double val)
  {
    return SetVariableValue(ivar, val);
  }

  bool SetParValues(const double* val)
  {
    return SetVariableValues(val);
  }

  bool SetConstant(unsigned int ival, const std::string& name, double var)
  {
    return SetFixedVariable(ival, name, var);
  }

  bool FixParameter(unsigned int ivar)
  {
    return FixVariable(ivar);
  }

  bool FreeParameter(unsigned int ivar)
  {
    return ReleaseVariable(ivar);
  }

  double GetFcn() const
  {
    return MinValue();
  }

  //Reminder
  //bool GetMinosError(unsigned int i, double &errLow, double &errUp, int=0);
  //bool Minimize();
  //bool Scan (unsigned int i, unsigned int &nstep, double *x, double *y, double xmin=0, double xmax=0);
  //unsigned int NFree () const;
  //double Correlation (unsigned int i, unsigned int j) const;
  //mini->PrintResults();

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
