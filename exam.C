/**
 * @file exam.C
 * @brief  An example to use the minuit2 package
 * @author Jinjing Li, josiahleeoaa@outlook.com
 * @version 1.0.0
 * @date 2020-09-01
 */
#include "TMinimization.h"
#include "TRandom3.h"
#include "TMath.h"
#include "TH1D.h"

#include <vector>
#include <iostream>

using namespace std;

vector<double>* data = new vector<double>;

double Function::operator()(const double* par)
{
  double res   = 0;
  double mu    = par[0];
  double sigma = par[1];

  //Gausaian likelihood function
  for (int i = 0; i < data->size(); i++)
  {
    res += 0.5 * TMath::Log(1 / pow(sigma, 2)) - pow(data->at(i) - mu, 2) / 2 / sigma / sigma;
  }
  return -2 * res;
}

//g++ -o exam exam.C `root-config --cflags --libs` -lminuit2
int main(void)
{

  //Prepare data
  TRandom3 gr(0);
  TH1D* h = new TH1D("h", "", 100, -5, 5);
  for (int i = 0; i < 1000; i++)
  {
    data->push_back(gr.Gaus(0, 1));
    h->Fill(data->back());
  }

  //Usage

  ROOT::Math::Functor Chi2Functor(fcn, 2);
  TMini* mini = new TMini(Chi2Functor);

  mini->SetParameter(0, "mu", 0, 1e-8);
  mini->SetParameter(1, "sigma", 1, 1e-8);

  mini->Minimize();

  const double* pars = mini->GetParameters();
  const double* errs = mini->GetParErrors();

  cout << Form("MinValue= %.2f", mini->GetFcn()) << endl;

  //Example of parameter scan
  unsigned int nstep = 100;
  double x[nstep];
  double y[nstep];
  mini->Scan(0, nstep, x, y, pars[0] - 3 * errs[0], pars[0] + 3 * errs[0]);
  for (int i = 0; i < 100; i++)
  {
    cout << Form("%.2f %.2f", x[i], y[i] - mini->GetFcn()) << endl;
  }

  return 0;
}
