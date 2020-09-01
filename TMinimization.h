/**
 * @file TMinimization.h
 * @brief  An easy-used minuit2 package.
 * @author Jinjing Li, josiahleeoaa@outlook.com
 * @version 1.0.0
 * @date 2020-09-01
 */

#include "Minuit2/Minuit2Minimizer.h"
#include "Math/Functor.h"
#include "Math/QuantFuncMathCore.h"
#include "Math/PdfFuncMathCore.h"
#include "Math/ProbFuncMathCore.h"
#include "TFile.h"
#include "TGraph.h"
#include "TColor.h"

class Fcn
{
public:
  double operator()(const double* par);
};

class TMini : public ROOT::Minuit2::Minuit2Minimizer
{

public:
  // Default constructer
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

  // Get parameter list
  const double* GetParameters() const { return X(); }

  // Get one of the parameter
  const double GetParameter(int i) const { return X()[i]; }

  // Get error list
  const double* GetParErrors() const { return Errors(); }

  // Get one of the error
  const double GetParError(int i) const { return Errors()[i]; }

  // Set Parameter
  bool SetParameter(unsigned int ivar, const std::string& name, double val, double step)
  {
    return SetVariable(ivar, name, val, step);
  }

  bool SetParLimits(unsigned int ivar, double lower, double upper)
  {
    return SetVariableLimits(ivar, lower, upper);
  }

  bool SetParValue(unsigned int ivar, double val) { return SetVariableValue(ivar, val); }

  bool SetParValues(const double* val) { return SetVariableValues(val); }

  bool SetConstant(unsigned int ival, const std::string& name, double var)
  {
    return SetFixedVariable(ival, name, var);
  }

  bool FixParameter(unsigned int ivar) { return FixVariable(ivar); }

  bool FreeParameter(unsigned int ivar) { return ReleaseVariable(ivar); }

  double GetFcn() const { return MinValue(); }

  // Reminder
  // bool GetMinosError(unsigned int i, double &errLow, double &errUp, int=0);
  // bool Minimize();
  // bool Scan (unsigned int i, unsigned int &nstep, double *x, double *y, double xmin=0, double xmax=0);
  // unsigned int NFree () const; double Correlation (unsigned int i, unsigned int j) const;
  // mini->PrintResults();

  //Given the measured X² value for a set of experiments with a degree of freedom d, the probability of the result being due to chance.
  //0.683 <- 2.2977, 0.955 <- 6.2021, 0.997 <- 11.6182, with n=2
  //0.683 <- 1,      0.955 <- 4,      0.997 <- 9,       with n=1
  //1 - P value = confidence level
  static double ConfidenceLevel(double chi2, int n){
    return ROOT::Math::chisquared_cdf(chi2,n);
  }
  void CreateContour(double confi_level, unsigned int Npoints, unsigned int parI, unsigned int parJ, double *array_I, double *array_J)
  {
    //To determine the chi-square value indicating a probability Q of non-chance occurrence for an experiment with d degrees of freedom
    //0.683 -> 2.2977, 0.955 -> 6.2021, 0.997->11.6182, with n=2
    //0.683 -> 1,      0.955 -> 4,      0.997->9,       with n=1
    SetErrorDef(ROOT::Math::chisquared_quantile(confi_level,2));
    ROOT::Minuit2::Minuit2Minimizer::Contour(parI, parJ, Npoints, array_I, array_J);
  }
  void CreateStandardContour(int nsigma, unsigned int Npoints, unsigned int parI, unsigned int parJ)
  {

    /// https://www.fourmilab.ch/rpkp/experiments/analysis/chiCalc.html
    // 68.3% confidence level // 1 - 0.317 with ndf = 2
    // 95.5% confidence level // 1 - 0.045 with ndf = 2
    // 99.7% confidence level // 1 - 0.003 with ndf = 2
    double chi2_quant[3] = {2.2977, 6.2021, 11.6182};
    switch (nsigma)
    {
      case 1:
      case 2:
      case 3:
        SetErrorDef(chi2_quant[nsigma - 1]);

        contourEdge_I[nsigma - 1] = new double[Npoints];
        contourEdge_J[nsigma - 1] = new double[Npoints];

        ROOT::Minuit2::Minuit2Minimizer::Contour(parI, parJ, Npoints, contourEdge_I[nsigma - 1], contourEdge_J[nsigma - 1]);

        break;
      default:
        std::cout << "Invalid nsigma option" << std::endl;
        return;
    }
  }
  void GetContourGraphs(const char* name, unsigned int Npoints)
  {
    TFile* oufile = new TFile(name, "RECREATE");
    TGraph* contour_nsigma[3];
    int colors[3] = {
        TColor::GetColor("#ff9999"),
        TColor::GetColor("#99ff99"),
        TColor::GetColor("#9999ff"),
    };
    int color_nsigma[3] = {kRed, kBlue, kGreen + 1};

    for (int i = 0; i < 3; i++)
    {
      contour_nsigma[i] = new TGraph(Npoints, contourEdge_I[i], contourEdge_J[i]);
      contour_nsigma[i]->SetLineColor(colors[i]);
      contour_nsigma[i]->SetLineWidth(2);
      contour_nsigma[i]->SetMarkerColor(color_nsigma[i]);
      contour_nsigma[i]->SetFillColor(colors[i]);
      contour_nsigma[i]->Write(Form("contour_%dsigma", i + 1));
    }
    oufile->Close();
  }

public:
  double* contourEdge_I[3];
  double* contourEdge_J[3];
};
