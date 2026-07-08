/**
 * @file TMinimization.h
 * @brief Template-based Minuit2/Minuit1 minimization wrapper with modern C++ API.
 * @author Jinjing Li, josiahleeoaa@outlook.com
 * @version 2.0
 * @date 2020-09-01
 *
 * Provides a unified interface for both Minuit2 and legacy TMinuit backends via
 * compile-time template selection. Default backend is Minuit2.
 */

#ifndef TMINIMIZATION_H
#define TMINIMIZATION_H

#include <algorithm>
#include <array>
#include <cmath>
#include <ctime>
#include <fstream>
#include <functional>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <memory>
#include <numeric>
#include <set>
#include <string>
#include <unistd.h>
#include <vector>

#include "Fit/ParameterSettings.h"
#include "Math/Functor.h"
#include "Math/PdfFuncMathCore.h"
#include "Math/ProbFuncMathCore.h"
#include "Math/QuantFuncMathCore.h"
#include "Minuit2/Minuit2Minimizer.h"
#include "Minuit2/MinimumState.h"
#include "Minuit2/MnTraceObject.h"
#include "Minuit2/MnUserParameterState.h"
#include "TColor.h"
#include "TFile.h"
#include "TGraph.h"
#include "TH2.h"
#include "TMinuitMinimizer.h"
#include "TRandom3.h"
#include "TString.h"

// =============================================================================
// Enumerations
// =============================================================================

enum class MinimizerBackend { Minuit1, Minuit2 };

enum class FitStatus { Success, SuccessApprox, NotConverged, Failed, InvalidInput };

enum class Algorithm { Migrad, Simplex, Combined, Fumili2 };

// =============================================================================
// Backend Traits
// =============================================================================

namespace detail {

template <MinimizerBackend B>
struct BackendBase;

template <>
struct BackendBase<MinimizerBackend::Minuit2> {
  using type = ROOT::Minuit2::Minuit2Minimizer;
  static constexpr bool kSupportsTrace = true;
  static constexpr bool kSupportsFumili = true;
};

template <>
struct BackendBase<MinimizerBackend::Minuit1> {
  using type = TMinuitMinimizer;
  static constexpr bool kSupportsTrace = false;
  static constexpr bool kSupportsFumili = false;
};

}  // namespace detail

// =============================================================================
// Forward Declarations
// =============================================================================

template <MinimizerBackend B>
class TMini;

// =============================================================================
// FitResult
// =============================================================================

struct FitResult {
  FitStatus status = FitStatus::Failed;
  double min_value = 0.0;
  double edm = 0.0;
  unsigned int nfcn = 0;
  unsigned int niter = 0;
  unsigned int nfree = 0;
  bool has_valid_covariance = false;
  bool has_accurate_covariance = false;
  std::vector<double> parameters;
  std::vector<double> errors;
  std::vector<std::string> par_names;
  std::vector<std::vector<double>> covariance;
  std::vector<std::vector<double>> correlation;

  double Parameter(size_t i) const { return parameters.at(i); }

  double Parameter(const std::string& name) const {
    for (size_t i = 0; i < par_names.size(); ++i)
      if (par_names[i] == name)
        return parameters[i];
    return std::numeric_limits<double>::quiet_NaN();
  }

  double Error(size_t i) const { return errors.at(i); }

  double Correlation(size_t i, size_t j) const {
    if (correlation.empty())
      return (i == j) ? 1.0 : 0.0;
    return correlation.at(i).at(j);
  }

  double GoodnessOfFit(unsigned int ndf) const {
    if (ndf == 0)
      return 0.0;
    return ROOT::Math::chisquared_cdf_c(min_value, ndf);
  }

  void Print(std::ostream& os = std::cout) const {
    os << "=== FitResult ===" << std::endl;
    os << "Status: " << static_cast<int>(status) << std::endl;
    os << "MinValue: " << min_value << "  EDM: " << edm << std::endl;
    os << "NFcn: " << nfcn << "  NFree: " << nfree << std::endl;
    os << "Covariance valid: " << has_valid_covariance
       << "  accurate: " << has_accurate_covariance << std::endl;
    for (size_t i = 0; i < parameters.size(); ++i) {
      os << "  " << std::setw(25) << std::left
         << (i < par_names.size() ? par_names[i] : std::to_string(i)) << " = " << std::setw(14)
         << parameters[i] << " +/- " << errors[i] << std::endl;
    }
  }

  void SaveJSON(const std::string& filename) const {
    std::ofstream out(filename);
    if (!out.is_open())
      return;
    out << "{\n";
    out << "  \"status\": " << static_cast<int>(status) << ",\n";
    out << "  \"min_value\": " << std::setprecision(15) << min_value << ",\n";
    out << "  \"edm\": " << edm << ",\n";
    out << "  \"nfcn\": " << nfcn << ",\n";
    out << "  \"nfree\": " << nfree << ",\n";
    out << "  \"parameters\": [\n";
    for (size_t i = 0; i < parameters.size(); ++i) {
      out << "    {\"name\": \"" << (i < par_names.size() ? par_names[i] : "") << "\", "
          << "\"value\": " << parameters[i] << ", "
          << "\"error\": " << errors[i] << "}" << (i + 1 < parameters.size() ? "," : "") << "\n";
    }
    out << "  ]\n}\n";
  }
};

// =============================================================================
// ProfileResult
// =============================================================================

struct ProfileResult {
  unsigned int par_index = 0;
  std::vector<double> x_values;
  std::vector<double> chi2_values;
  double best_fit = 0.0;
  double chi2_min = 0.0;
  double lower_1sigma = 0.0;
  double upper_1sigma = 0.0;
  double lower_2sigma = 0.0;
  double upper_2sigma = 0.0;

  void SaveTSV(const std::string& filename) const {
    std::ofstream out(filename);
    if (!out.is_open())
      return;
    out << "#x\tchi2\tdelta_chi2\n";
    for (size_t i = 0; i < x_values.size(); ++i)
      out << std::setprecision(10) << x_values[i] << '\t' << chi2_values[i] << '\t'
          << chi2_values[i] - chi2_min << '\n';
  }

  TGraph* ToTGraph() const {
    if (x_values.empty())
      return nullptr;
    auto* g = new TGraph(static_cast<int>(x_values.size()), x_values.data(), chi2_values.data());
    g->SetTitle(";Parameter value;#chi^{2}");
    return g;
  }
};

// =============================================================================
// GridScanResult
// =============================================================================

struct GridScanResult {
  unsigned int par_i = 0;
  unsigned int par_j = 0;
  std::vector<double> x_values;
  std::vector<double> y_values;
  std::vector<std::vector<double>> chi2_grid;
  double best_x = 0.0;
  double best_y = 0.0;
  double chi2_min = 0.0;

  void SaveTSV(const std::string& filename) const {
    std::ofstream out(filename);
    if (!out.is_open())
      return;
    out << "#x\ty\tchi2\n";
    for (size_t ix = 0; ix < x_values.size(); ++ix)
      for (size_t iy = 0; iy < y_values.size(); ++iy)
        out << std::setprecision(10) << x_values[ix] << '\t' << y_values[iy] << '\t'
            << chi2_grid[ix][iy] << '\n';
  }

  TH2D* ToTH2D(const char* name = "grid_scan") const {
    if (x_values.empty() || y_values.empty())
      return nullptr;
    int nx = static_cast<int>(x_values.size());
    int ny = static_cast<int>(y_values.size());
    double dx = (nx > 1) ? (x_values[1] - x_values[0]) / 2.0 : 0.5;
    double dy = (ny > 1) ? (y_values[1] - y_values[0]) / 2.0 : 0.5;
    auto* h = new TH2D(name, "", nx, x_values.front() - dx, x_values.back() + dx, ny,
                        y_values.front() - dy, y_values.back() + dy);
    for (int ix = 0; ix < nx; ++ix)
      for (int iy = 0; iy < ny; ++iy)
        h->SetBinContent(ix + 1, iy + 1, chi2_grid[ix][iy]);
    return h;
  }
};

// =============================================================================
// MultiStartConfig
// =============================================================================

struct MultiStartConfig {
  unsigned int n_starts = 10;
  unsigned int seed = 0;
  std::vector<std::pair<double, double>> par_ranges;
  Algorithm algo = Algorithm::Combined;
};

// =============================================================================
// ParameterHandle (fluent builder)
// =============================================================================

template <MinimizerBackend B>
class ParameterHandle {
public:
  ParameterHandle(TMini<B>* owner, unsigned int idx) : m_owner(owner), m_index(idx) {}

  ParameterHandle& SetValue(double val) {
    m_owner->SetParValue(m_index, val);
    return *this;
  }

  ParameterHandle& SetStep(double step) {
    m_owner->SetVariableStepSize(m_index, step);
    return *this;
  }

  ParameterHandle& SetLimits(double lo, double hi) {
    m_owner->SetParLimits(m_index, lo, hi);
    return *this;
  }

  ParameterHandle& Fix(double val) {
    m_owner->FixParameter(m_index, val);
    return *this;
  }

  ParameterHandle& Fix() {
    m_owner->FixParameter(m_index);
    return *this;
  }

  ParameterHandle& Release() {
    m_owner->FreeParameter(m_index);
    return *this;
  }

  ParameterHandle& SetGroup(const std::string& group) {
    m_owner->paras_groups[group].insert(m_index);
    return *this;
  }

  double Value() const { return m_owner->GetParameter(m_index); }
  double Error() const { return m_owner->GetParError(m_index); }
  std::string Name() const { return m_owner->GetParName(m_index); }
  bool IsFixed() const { return m_owner->IsFixedVariable(m_index); }
  unsigned int Index() const { return m_index; }

private:
  TMini<B>* m_owner;
  unsigned int m_index;
};

// =============================================================================
// StateSnapshot
// =============================================================================

template <MinimizerBackend B>
class StateSnapshot {
public:
  StateSnapshot() = default;

  explicit StateSnapshot(const TMini<B>& mini) { Capture(mini); }

  void Capture(const TMini<B>& mini) {
    m_values = mini.GetParametersVector();
    m_valid = !m_values.empty();
    m_ndim = mini.NDim();
    m_fixed.clear();
    for (unsigned int i = 0; i < m_ndim; ++i) {
      if (mini.IsFixedVariable(i))
        m_fixed.insert(i);
    }
  }

  void Restore(TMini<B>& mini) const {
    if (!m_valid)
      return;
    mini.SetVariableValues(m_values.data());
    for (unsigned int i = 0; i < m_ndim; ++i) {
      if (m_fixed.count(i))
        mini.FixVariable(i);
      else
        mini.ReleaseVariable(i);
    }
  }

  bool IsValid() const { return m_valid; }

private:
  std::vector<double> m_values;
  std::set<unsigned int> m_fixed;
  unsigned int m_ndim = 0;
  bool m_valid = false;
};

// =============================================================================
// TMini Template Class
// =============================================================================

template <MinimizerBackend B = MinimizerBackend::Minuit2>
class TMini : public detail::BackendBase<B>::type {
public:
  using Base = typename detail::BackendBase<B>::type;

  TMini(int level = 2, double tolerance = 1e-6) : Base([&]() {
    if constexpr (B == MinimizerBackend::Minuit2)
      return ROOT::Minuit2::kMigrad;
    else
      return ROOT::Minuit::kMigrad;
  }()) {
    this->SetPrintLevel(level);
    this->SetStrategy(1);
    this->SetMaxFunctionCalls(5000000);
    this->SetMaxIterations(5000000);
    this->SetTolerance(tolerance);
  }

  virtual ~TMini() { DisableIterationTrace(); }

  TMini(const TMini&) = delete;
  TMini& operator=(const TMini&) = delete;

  // =========================================================================
  // Virtual hook for parameter initialization (override in derived classes)
  // =========================================================================

  virtual unsigned int InitializeParams() { return 0; }

  // =========================================================================
  // Parameter Accessors
  // =========================================================================

  const double* GetParameters() const { return this->X(); }

  std::vector<double> GetParametersVector() const {
    if (!this->X())
      return {};
    return std::vector<double>(this->X(), this->X() + this->NDim());
  }

  double GetParameter(size_t i) const { return this->X()[i]; }

  const double* GetParErrors() const { return this->Errors(); }

  std::vector<double> GetParErrorsVector() const {
    if (!this->Errors())
      return {};
    return std::vector<double>(this->Errors(), this->Errors() + this->NDim());
  }

  double GetParError(size_t i) const { return this->Errors()[i]; }

  std::string GetParName(unsigned int idx) const {
    ROOT::Fit::ParameterSettings var_settings;
    this->GetVariableSettings(idx, var_settings);
    return var_settings.Name();
  }

  unsigned int GetParIndex(const std::string& name) const {
    return static_cast<unsigned int>(this->VariableIndex(name));
  }

  // =========================================================================
  // Parameter Modifiers
  // =========================================================================

  bool SetParameter(unsigned int ivar, const std::string& name, double val, double step) {
    return this->SetVariable(ivar, name, val, step);
  }

  bool SetParLimits(unsigned int ivar, double lower, double upper) {
    return this->SetVariableLimits(ivar, lower, upper);
  }

  bool SetParValue(unsigned int ivar, double val) { return this->SetVariableValue(ivar, val); }

  bool SetParValues(const double* val) { return this->SetVariableValues(val); }

  bool SetConstant(unsigned int ival, const std::string& name, double var) {
    return this->SetFixedVariable(ival, name, var);
  }

  bool SetConstant(unsigned int ival, double var) {
    this->SetVariableValue(ival, var);
    return this->FixVariable(ival);
  }

  bool FixParameter(unsigned int idx, double val) {
    if (this->IsFixedVariable(idx))
      this->ReleaseVariable(idx);
    this->SetVariableValue(idx, val);
    return this->FixVariable(idx);
  }

  bool FixParameter(unsigned int ivar) { return this->FixVariable(ivar); }

  bool FreeParameter(unsigned int ivar) { return this->ReleaseVariable(ivar); }

  // =========================================================================
  // Algorithm Selection
  // =========================================================================

  void SetAlgorithm(Algorithm algo) { m_algorithm = algo; }
  Algorithm GetAlgorithm() const { return m_algorithm; }

  bool FitWith(Algorithm algo, int print_level = -2) {
    Algorithm prev = m_algorithm;
    m_algorithm = algo;
    bool result = Fit(print_level);
    m_algorithm = prev;
    return result;
  }

  // =========================================================================
  // Iteration Trace (Minuit2 only)
  // =========================================================================

  bool EnableIterationTrace(const std::string& output_filename, int theta13_index = 0,
                            int dm2_index = 1) {
    if constexpr (detail::BackendBase<B>::kSupportsTrace) {
      DisableIterationTrace();
      if (output_filename.empty()) {
        std::cerr << "[WARN] Empty iteration trace filename, skip trace export." << std::endl;
        return false;
      }
      m_trace_stream.open(output_filename, std::ios::out | std::ios::trunc);
      if (!m_trace_stream.is_open()) {
        std::cerr << "[WARN] Cannot open iteration trace file: " << output_filename << std::endl;
        return false;
      }
      m_trace_stream << "#iter\tnfcn\tchi2\tedm\tsin2_2theta13\tdm2_param\n";
      m_trace_stream.flush();
      m_trace_object =
          std::make_unique<IterationTraceObject>(&m_trace_stream, theta13_index, dm2_index);
      this->SetTraceObject(*m_trace_object);
      m_trace_enabled = true;
      return true;
    } else {
      (void)output_filename;
      (void)theta13_index;
      (void)dm2_index;
      std::cerr << "[WARN] Iteration trace requires Minuit2 backend." << std::endl;
      return false;
    }
  }

  void DisableIterationTrace() {
    if constexpr (detail::BackendBase<B>::kSupportsTrace) {
      if (m_trace_stream.is_open()) {
        m_trace_stream.flush();
        m_trace_stream.close();
      }
    }
    m_trace_enabled = false;
  }

  bool IterationTraceEnabled() const { return m_trace_enabled; }

  // =========================================================================
  // Config-Fixed Parameters
  // =========================================================================

  void RegisterConfigFixedParam(unsigned int idx) { m_config_fixed_params.insert(idx); }

  void ApplyConfigFixedParams(const std::vector<double>& best_vals, bool fixed_at_zeros) {
    for (unsigned int idx : m_config_fixed_params)
      FixParameter(idx, fixed_at_zeros ? 0.0 : best_vals[idx]);
  }

  // =========================================================================
  // Parameter Info Registration
  // =========================================================================

  void AddParameterInfo(unsigned int& par_index, const std::string& name, double init_val,
                        double step_size, const std::string& error_group) {
    SetParameter(par_index, name, init_val, step_size);
    paras_groups[error_group].insert(par_index);
    ParameterIndices[name] = par_index;
    par_index++;
  }

  // =========================================================================
  // Parameter Handle API (new)
  // =========================================================================

  ParameterHandle<B> AddParam(const std::string& name, double val, double step) {
    unsigned int idx = static_cast<unsigned int>(this->NDim());
    SetParameter(idx, name, val, step);
    return ParameterHandle<B>(this, idx);
  }

  ParameterHandle<B> AddFixedParam(const std::string& name, double val) {
    unsigned int idx = static_cast<unsigned int>(this->NDim());
    SetConstant(idx, name, val);
    return ParameterHandle<B>(this, idx);
  }

  ParameterHandle<B> Par(unsigned int idx) { return ParameterHandle<B>(this, idx); }

  ParameterHandle<B> Par(const std::string& name) {
    unsigned int idx = GetParIndex(name);
    return ParameterHandle<B>(this, idx);
  }

  // =========================================================================
  // State Management
  // =========================================================================

  StateSnapshot<B> SaveState() const { return StateSnapshot<B>(*this); }

  void RestoreState(const StateSnapshot<B>& snap) { snap.Restore(*this); }

  // =========================================================================
  // Fitting
  // =========================================================================

  double GetFcn() const { return this->MinValue(); }

  bool Fit(int print_level = -2) {
    int old_pl = this->PrintLevel();
    if (print_level != -2)
      this->SetPrintLevel(print_level);

    if (m_algorithm != Algorithm::Migrad)
      ApplyAlgorithm(m_algorithm);
    bool res = this->Minimize();

    if (print_level != -2)
      this->SetPrintLevel(old_pl);
    return res;
  }

  // =========================================================================
  // Structured Results
  // =========================================================================

  FitResult GetResult() const {
    FitResult r;
    r.min_value = this->MinValue();
    r.edm = this->Edm();
    r.nfcn = this->NCalls();
    r.nfree = this->NFree();
    r.has_valid_covariance = (this->CovMatrixStatus() >= 1);
    r.has_accurate_covariance = (this->CovMatrixStatus() >= 3);

    if (this->Status() == 0)
      r.status = FitStatus::Success;
    else if (this->Status() == 1)
      r.status = FitStatus::SuccessApprox;
    else if (this->Status() == 2)
      r.status = FitStatus::NotConverged;
    else
      r.status = FitStatus::Failed;

    size_t n = this->NDim();
    r.parameters.resize(n);
    r.errors.resize(n);
    r.par_names.resize(n);
    for (size_t i = 0; i < n; ++i) {
      r.parameters[i] = this->X()[i];
      r.errors[i] = this->Errors()[i];
      r.par_names[i] = GetParName(static_cast<unsigned int>(i));
    }
    return r;
  }

  FitResult FitAndGetResult(int print_level = -2) {
    Fit(print_level);
    return GetResult();
  }

  // =========================================================================
  // Advanced Statistics
  // =========================================================================

  ProfileResult ProfileLikelihood(unsigned int par_idx, unsigned int npoints, double lo,
                                  double hi) {
    ProfileResult result;
    result.par_index = par_idx;
    result.best_fit = this->X()[par_idx];
    result.chi2_min = this->MinValue();
    result.x_values.resize(npoints);
    result.chi2_values.resize(npoints);

    auto snap = SaveState();

    for (unsigned int i = 0; i < npoints; ++i) {
      double x = lo + (hi - lo) * i / std::max(1u, npoints - 1);
      FixParameter(par_idx, x);
      this->SetValidError(false);
      Fit(-1);
      result.x_values[i] = x;
      result.chi2_values[i] = this->MinValue();
    }

    FreeParameter(par_idx);
    RestoreState(snap);

    // Compute confidence intervals from delta chi2
    auto find_crossing = [&](double delta_chi2, bool lower) -> double {
      double target = result.chi2_min + delta_chi2;
      const auto& xv = result.x_values;
      const auto& yv = result.chi2_values;
      size_t best_idx = 0;
      for (size_t i = 0; i < xv.size(); ++i)
        if (yv[i] < yv[best_idx])
          best_idx = i;

      if (lower) {
        for (size_t i = best_idx; i > 0; --i) {
          if (yv[i - 1] >= target && yv[i] < target) {
            double frac = (target - yv[i]) / (yv[i - 1] - yv[i]);
            return xv[i] + frac * (xv[i - 1] - xv[i]);
          }
        }
        return xv.front();
      } else {
        for (size_t i = best_idx; i + 1 < xv.size(); ++i) {
          if (yv[i] < target && yv[i + 1] >= target) {
            double frac = (target - yv[i]) / (yv[i + 1] - yv[i]);
            return xv[i] + frac * (xv[i + 1] - xv[i]);
          }
        }
        return xv.back();
      }
    };

    result.lower_1sigma = find_crossing(1.0, true);
    result.upper_1sigma = find_crossing(1.0, false);
    result.lower_2sigma = find_crossing(4.0, true);
    result.upper_2sigma = find_crossing(4.0, false);

    return result;
  }

  GridScanResult GridScan(unsigned int par_i, unsigned int par_j, unsigned int nx, unsigned int ny,
                          double xi_lo, double xi_hi, double xj_lo, double xj_hi) {
    GridScanResult result;
    result.par_i = par_i;
    result.par_j = par_j;
    result.best_x = this->X()[par_i];
    result.best_y = this->X()[par_j];
    result.chi2_min = this->MinValue();

    result.x_values.resize(nx);
    result.y_values.resize(ny);
    result.chi2_grid.resize(nx, std::vector<double>(ny, 0.0));

    for (unsigned int ix = 0; ix < nx; ++ix)
      result.x_values[ix] = xi_lo + (xi_hi - xi_lo) * ix / std::max(1u, nx - 1);
    for (unsigned int iy = 0; iy < ny; ++iy)
      result.y_values[iy] = xj_lo + (xj_hi - xj_lo) * iy / std::max(1u, ny - 1);

    auto snap = SaveState();

    for (unsigned int ix = 0; ix < nx; ++ix) {
      for (unsigned int iy = 0; iy < ny; ++iy) {
        FixParameter(par_i, result.x_values[ix]);
        FixParameter(par_j, result.y_values[iy]);
        this->SetValidError(false);
        Fit(-1);
        result.chi2_grid[ix][iy] = this->MinValue();
      }
    }

    FreeParameter(par_i);
    FreeParameter(par_j);
    RestoreState(snap);
    return result;
  }

  std::vector<std::vector<double>> GetCorrelationMatrix() const {
    size_t n = this->NDim();
    std::vector<double> cov(n * n, 0.0);
    this->GetCovMatrix(cov.data());

    std::vector<std::vector<double>> corr(n, std::vector<double>(n, 0.0));
    for (size_t i = 0; i < n; ++i) {
      for (size_t j = 0; j < n; ++j) {
        double si = std::sqrt(std::abs(cov[i * n + i]));
        double sj = std::sqrt(std::abs(cov[j * n + j]));
        corr[i][j] = (si * sj > 0.0) ? cov[i * n + j] / (si * sj) : ((i == j) ? 1.0 : 0.0);
      }
    }
    return corr;
  }

  double Chi2PerNDF(unsigned int ndf) const {
    return (ndf > 0) ? this->MinValue() / ndf : 0.0;
  }

  double PValue(unsigned int ndf) const {
    return ROOT::Math::chisquared_cdf_c(this->MinValue(), ndf);
  }

  std::pair<double, double> ConfidenceInterval(unsigned int par_idx, double delta_chi2 = 1.0) {
    double err_low = 0.0, err_high = 0.0;
    this->SetErrorDef(delta_chi2);
    this->GetMinosError(par_idx, err_low, err_high);
    return {err_low, err_high};
  }

  // =========================================================================
  // Multi-Start Optimization
  // =========================================================================

  FitResult MultiStartFit(const MultiStartConfig& config) {
    TRandom3 rng(config.seed);
    FitResult best_result;
    best_result.min_value = std::numeric_limits<double>::max();

    auto snap = SaveState();

    for (unsigned int start = 0; start < config.n_starts; ++start) {
      // Randomize starting point
      for (size_t i = 0; i < config.par_ranges.size() && i < this->NDim(); ++i) {
        if (this->IsFixedVariable(static_cast<unsigned int>(i)))
          continue;
        double lo = config.par_ranges[i].first;
        double hi = config.par_ranges[i].second;
        this->SetVariableValue(static_cast<unsigned int>(i), rng.Uniform(lo, hi));
      }

      FitWith(config.algo, -1);
      FitResult r = GetResult();

      if (r.min_value < best_result.min_value)
        best_result = r;
    }

    RestoreState(snap);

    // Final fit from best starting point
    if (!best_result.parameters.empty())
      this->SetVariableValues(best_result.parameters.data());
    Fit(-1);
    return GetResult();
  }

  // =========================================================================
  // Error Budget
  // =========================================================================

  void ExportErrorBudget(const std::vector<double>& best_vals, const std::string& method = "add",
                         const std::string& output_filename = "", bool fixed_at_zeros = false) {
    std::cout << "Method: " << method << " is specified for error budget calculation" << std::endl;

    std::ofstream fout;
    if (!output_filename.empty())
      fout.open(output_filename);

    int original_pl = this->PrintLevel();
    this->SetPrintLevel(-1);

    const auto groups_snapshot = paras_groups;

    for (const auto& pair : groups_snapshot) {
      const std::string group_name = pair.first;
      const std::set<unsigned int> group_indices = pair.second;

      this->Clear();
      size_t params_num = InitializeParams();
      if (params_num != this->NDim()) {
        std::cerr << "Error: NDim (" << this->NDim() << ") differs from parameter count ("
                  << params_num << ")" << std::endl;
        continue;
      }

      ApplyConfigFixedParams(best_vals, fixed_at_zeros);

      std::cout << "Processing error component: " << group_name << std::endl;

      for (unsigned int j = 0; j < this->NDim(); j++) {
        if (m_config_fixed_params.count(j))
          continue;

        bool is_stat = (paras_groups.count("stat") && paras_groups.at("stat").count(j));
        bool is_in_group = group_indices.count(j);

        if (method == "add") {
          if (!is_stat && !is_in_group)
            FixParameter(j, fixed_at_zeros ? 0 : best_vals[j]);
        } else if (method == "sub") {
          if (!is_stat && is_in_group)
            FixParameter(j, fixed_at_zeros ? 0 : best_vals[j]);
        }
      }

      Fit(-1);
      double err_low_0, err_high_0, err_low_1, err_high_1;
      this->GetMinosError(0, err_low_0, err_high_0);
      this->GetMinosError(1, err_low_1, err_high_1);

      std::string result_str =
          Form("%s: %.8f %.8f +%.8f; %.6f %.6f +%.6f", group_name.c_str(), this->X()[0], err_low_0,
               err_high_0, this->X()[1], err_low_1, err_high_1);
      std::cout << result_str << std::endl;
      if (fout.is_open())
        fout << result_str << std::endl;
    }

    if (fout.is_open())
      fout.close();
    this->SetPrintLevel(original_pl);
    this->Clear();
    InitializeParams();
    ApplyConfigFixedParams(best_vals, fixed_at_zeros);
  }

  // =========================================================================
  // Covariance Matrix Export
  // =========================================================================

  void ExportCovMatrix(
      const char* filename_txt = "./Output/FullCovMatrix.txt",
      const char* filename_root = "./Output/cov_matrix.root",
      std::function<bool(const std::string&)> is_degenerate = nullptr) {
    size_t n = this->NDim();
    std::vector<double> cov_matrix(n * n, 0.0);

    // Default degenerate predicate: Eps_shape parameters
    if (!is_degenerate) {
      is_degenerate = [](const std::string& pname) {
        return pname.rfind("Eps_shape", 0) == 0;
      };
    }

    std::vector<unsigned int> temp_fixed;

    this->Hesse();
    this->GetCovMatrix(cov_matrix.data());

    double diag_sum = 0.0;
    for (size_t i = 0; i < n; ++i)
      diag_sum += std::abs(cov_matrix[i * n + i]);
    if (diag_sum > 0.0) {
      for (unsigned int i = 0; i < static_cast<unsigned int>(n); ++i) {
        if (this->IsFixedVariable(i))
          continue;
        const std::string pname = GetParName(i);
        if (!is_degenerate(pname))
          continue;
        if (cov_matrix[i * n + i] == 0.0) {
          FixParameter(i, this->X()[i]);
          temp_fixed.push_back(i);
          std::cout << "[INFO] ExportCovMatrix: temporarily fixed degenerate parameter " << pname
                    << " for Hesse retry." << std::endl;
        }
      }
      if (!temp_fixed.empty()) {
        this->Hesse();
        std::fill(cov_matrix.begin(), cov_matrix.end(), 0.0);
        this->GetCovMatrix(cov_matrix.data());
      }
    }

    for (unsigned int idx : temp_fixed)
      FreeParameter(idx);

    TH2D h_cov("h_cov", ";Parameter index;Parameter index", static_cast<int>(n), 0.5,
               static_cast<int>(n) + 0.5, static_cast<int>(n), 0.5, static_cast<int>(n) + 0.5);
    std::ofstream cov_out_file(filename_txt);

    for (size_t i = 0; i < n; i++) {
      for (size_t j = 0; j < n; j++) {
        double s1 = std::sqrt(cov_matrix[i * n + i]);
        double s2 = std::sqrt(cov_matrix[j * n + j]);
        cov_out_file << Form("%.4e  ", cov_matrix[i * n + j]);
        double corr;
        if (s1 * s2 == 0.0) {
          corr = (i == j) ? 1.0 : 0.0;
        } else {
          corr = cov_matrix[i * n + j] / (s1 * s2);
        }
        h_cov.SetBinContent(static_cast<int>(i + 1), static_cast<int>(j + 1), corr);
      }
      cov_out_file << std::endl;
    }
    h_cov.SaveAs(filename_root);
  }

  // =========================================================================
  // Contours & Scans
  // =========================================================================

  static double ConfidenceLevel(double chi2, int n) {
    return ROOT::Math::chisquared_cdf(chi2, n);
  }

  void CreateContour(double confi_level, unsigned Npoints, unsigned int parI, unsigned int parJ,
                     double* array_I, double* array_J) {
    this->SetErrorDef(ROOT::Math::chisquared_quantile(confi_level, 2));
    this->Contour(parI, parJ, Npoints, array_I, array_J);
  }

  bool CreateStandardContour(size_t nsigma, unsigned int Npoints, unsigned int parI,
                             unsigned int parJ) {
    if (nsigma < 1 || nsigma > 3) {
      std::cout << "Invalid nsigma option: " << nsigma << std::endl;
      return false;
    }
    double chi2_quant[3] = {2.2977, 6.2021, 11.6182};
    this->SetErrorDef(chi2_quant[nsigma - 1]);

    contourEdge_I[nsigma - 1].assign(Npoints, 0.0);
    contourEdge_J[nsigma - 1].assign(Npoints, 0.0);

    bool ok = this->Contour(parI, parJ, Npoints, contourEdge_I[nsigma - 1].data(),
                            contourEdge_J[nsigma - 1].data());
    if (!ok) {
      contourEdge_I[nsigma - 1].clear();
      contourEdge_J[nsigma - 1].clear();
    }
    return ok;
  }

  void CreateSpecialContour(unsigned int Npoints, unsigned int parI, unsigned int parJ) {
    double chi2_90 = 4.6051;
    this->SetErrorDef(chi2_90);

    contourEdgeS_I.assign(Npoints, 0.0);
    contourEdgeS_J.assign(Npoints, 0.0);

    this->Contour(parI, parJ, Npoints, contourEdgeS_I.data(), contourEdgeS_J.data());
  }

  void GetContourGraphs(const char* name, size_t Npoints) {
    TFile oufile(name, "RECREATE");
    Color_t colors[3] = {
        static_cast<Color_t>(TColor::GetColor("#ff9999")),
        static_cast<Color_t>(TColor::GetColor("#99ff99")),
        static_cast<Color_t>(TColor::GetColor("#9999ff")),
    };
    Color_t marker_colors[3] = {kRed, kBlue, kGreen + 1};

    for (size_t i = 0; i < 3; i++) {
      if (contourEdge_I[i].empty())
        continue;
      TGraph g(static_cast<int>(Npoints), contourEdge_I[i].data(), contourEdge_J[i].data());
      g.SetLineColor(colors[i]);
      g.SetLineWidth(2);
      g.SetMarkerColor(marker_colors[i]);
      g.SetFillColor(colors[i]);
      g.Write(Form("contour_%zusigma", i + 1));
    }
  }

  void ExportSpecialContour(const char* name, size_t Npoints) {
    TFile oufile(name, "RECREATE");
    if (!contourEdgeS_I.empty()) {
      TGraph g(static_cast<int>(Npoints), contourEdgeS_I.data(), contourEdgeS_J.data());
      g.Write("contour");
    }
  }

  bool Scan_Chi2(unsigned int idx, unsigned int nstep, double x0, double x1, double* array_x,
                 double* array_y, const std::string& output_filename = "") {
    std::ofstream ofile;
    if (!output_filename.empty())
      ofile.open(output_filename);

    // Save/restore the full parameter state around the scan (same pattern as
    // ProfileLikelihood/GridScan). Without this, every parameter is left at
    // the last grid point's re-fit values, contaminating whatever runs next
    // (e.g. a second scan over another parameter).
    auto snap = SaveState();
    bool all_ok = true;

    for (unsigned int i = 0; i < nstep; i++) {
      double this_x = x0 + (x1 - x0) * i / (std::max(1u, nstep - 1));
      FixParameter(idx, this_x);
      this->SetValidError(false);
      if (!Fit(-1) || this->Status() != 0) {
        all_ok = false;
        std::cerr << "[WARN] Scan_Chi2: re-fit at par" << idx << " = " << this_x
                  << " did not converge cleanly (status=" << this->Status() << ")" << std::endl;
      }
      array_x[i] = this->X()[idx];
      array_y[i] = this->MinValue();
      if (ofile.is_open()) {
        ofile.precision(8);
        ofile << array_x[i] << "\t" << array_y[i] << std::endl;
      }
    }
    this->ReleaseVariable(idx);
    RestoreState(snap);
    return all_ok;
  }

  // =========================================================================
  // Batch Parameter Operations
  // =========================================================================

  void FixParameterGroup(const std::string& group_name, const double* values = nullptr) {
    if (paras_groups.count(group_name)) {
      for (unsigned int idx : paras_groups.at(group_name)) {
        if (values)
          this->SetVariableValue(idx, values[idx]);
        this->FixVariable(idx);
      }
    }
  }

  void FreeParameterGroup(const std::string& group_name) {
    if (paras_groups.count(group_name)) {
      for (unsigned int idx : paras_groups.at(group_name)) {
        this->ReleaseVariable(idx);
      }
    }
  }

  void FixParList(unsigned int idx_be, unsigned int idx_en, const double* params,
                  const bool* parfix, std::set<unsigned int> skip = {2}) {
    for (unsigned int i = idx_be; i < idx_en; i++) {
      if (parfix[i] || skip.count(i))
        continue;
      this->SetVariableValue(i, params[i]);
      this->FixVariable(i);
    }
  }

  void FixParList(unsigned int idx_be, unsigned int idx_en, const bool* parfix,
                  std::set<unsigned int> skip = {2}) {
    for (unsigned int i = idx_be; i < idx_en; i++) {
      if (parfix[i] || skip.count(i))
        continue;
      this->SetVariableValue(i, 0.0);
      this->FixVariable(i);
    }
  }

  void FixShapeVars(unsigned int begin, unsigned int end, const double* params) {
    for (unsigned int i = begin; i <= end; i++) {
      this->SetVariableValue(i, params[i]);
      this->FixVariable(i);
    }
  }

  void FreeParList(unsigned int idx_be, unsigned int idx_en, const bool* parfix,
                   std::set<unsigned int> skip = {}) {
    for (unsigned int i = idx_be; i < idx_en; i++) {
      if (parfix[i] || skip.count(i))
        continue;
      this->ReleaseVariable(i);
      this->SetVariableValue(i, 0.0);
    }
  }

  // =========================================================================
  // Utility / Debug
  // =========================================================================

  void PrintFixVars(size_t N) const {
    std::cout << "\nFixed variables: ";
    for (unsigned int i = 0; i < N; i++) {
      ROOT::Fit::ParameterSettings settings;
      this->GetVariableSettings(i, settings);
      if (settings.IsFixed())
        std::cout << settings.Name() << " ";
    }
    std::cout << std::endl;
  }

  void ShowVarStatus(size_t N) const {
    std::cout << "\nVariable Status:\n";
    for (unsigned int i = 0; i < N; i++) {
      ROOT::Fit::ParameterSettings settings;
      this->GetVariableSettings(i, settings);
      std::cout << settings.Name() << " \tValue: " << settings.Value()
                << " \tFixed: " << (settings.IsFixed() ? "Yes" : "No") << std::endl;
    }
  }

  // =========================================================================
  // Public Data Members
  // =========================================================================

public:
  std::array<std::vector<double>, 3> contourEdge_I;
  std::array<std::vector<double>, 3> contourEdge_J;
  std::vector<double> contourEdgeS_I;
  std::vector<double> contourEdgeS_J;

  std::map<std::string, std::set<unsigned int>> paras_groups;
  static std::map<std::string, size_t> ParameterIndices;
  std::set<unsigned int> m_config_fixed_params;

private:
  void ApplyAlgorithm(Algorithm algo) {
    if constexpr (B == MinimizerBackend::Minuit2) {
      switch (algo) {
        case Algorithm::Migrad:
          this->SetMinimizerType(ROOT::Minuit2::kMigrad);
          break;
        case Algorithm::Simplex:
          this->SetMinimizerType(ROOT::Minuit2::kSimplex);
          break;
        case Algorithm::Combined:
          this->SetMinimizerType(ROOT::Minuit2::kCombined);
          break;
        case Algorithm::Fumili2:
          this->SetMinimizerType(ROOT::Minuit2::kFumili);
          break;
      }
    } else {
      if (algo != Algorithm::Migrad && algo != Algorithm::Simplex) {
        std::cerr << "[WARN] Minuit1 backend supports only Migrad/Simplex." << std::endl;
      }
    }
  }

  Algorithm m_algorithm = Algorithm::Migrad;
  bool m_trace_enabled = false;

  // Minuit2-only trace infrastructure
  struct IterationTraceObject;
  std::ofstream m_trace_stream;
  std::unique_ptr<IterationTraceObject> m_trace_object;
};

// =============================================================================
// IterationTraceObject definition (Minuit2 backend only)
// =============================================================================

template <MinimizerBackend B>
struct TMini<B>::IterationTraceObject : public ROOT::Minuit2::MnTraceObject {
  IterationTraceObject(std::ofstream* out, int theta13_index, int dm2_index)
      : m_out(out), m_theta13_index(theta13_index), m_dm2_index(dm2_index) {}

  void operator()(int iter, const ROOT::Minuit2::MinimumState& state) override {
    if (!m_out || !m_out->is_open())
      return;

    const auto& params = state.Vec();
    double theta13 = std::numeric_limits<double>::quiet_NaN();
    double dm2 = std::numeric_limits<double>::quiet_NaN();

    if (m_theta13_index >= 0 && static_cast<size_t>(m_theta13_index) < params.size())
      theta13 = params[static_cast<size_t>(m_theta13_index)];
    if (m_dm2_index >= 0 && static_cast<size_t>(m_dm2_index) < params.size())
      dm2 = params[static_cast<size_t>(m_dm2_index)];

    (*m_out) << iter << '\t' << state.NFcn() << '\t' << std::setprecision(15) << state.Fval()
             << '\t' << state.Edm() << '\t' << theta13 << '\t' << dm2 << '\n';
  }

  std::ofstream* m_out;
  int m_theta13_index;
  int m_dm2_index;
};

// =============================================================================
// Static Member Definitions
// =============================================================================

template <MinimizerBackend B>
std::map<std::string, size_t> TMini<B>::ParameterIndices = std::map<std::string, size_t>();

// =============================================================================
// Convenience Type Alias
// =============================================================================

using TMiniMinuit2 = TMini<MinimizerBackend::Minuit2>;
using TMiniMinuit1 = TMini<MinimizerBackend::Minuit1>;

#endif  // TMINIMIZATION_H
