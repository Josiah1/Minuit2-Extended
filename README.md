# Minuit2-Extended v2.0

A modern C++17 template-based wrapper for ROOT's Minuit2 (and legacy TMinuit) minimizers. Provides a clean, type-safe API with compile-time backend selection, structured results, and advanced statistics — while preserving full backward compatibility with existing code.

## Features

- **Template backend selection** — `TMini<MinimizerBackend::Minuit2>` (default) or `TMini<MinimizerBackend::Minuit1>` with zero-cost compile-time switching
- **Virtual `InitializeParams()`** — override in derived classes for application-specific parameter registration
- **Fluent ParameterHandle** — chainable `AddParam("name", val, step).SetLimits(lo, hi).SetGroup("grp")`
- **Structured FitResult** — status, parameters, errors, covariance, JSON export
- **Profile likelihood** — 1D profile scans with automatic confidence interval extraction
- **2D grid scan** — full grid minimization with TH2D/TSV export
- **Multi-start optimization** — random restarts for global optimization
- **State snapshots** — save/restore minimizer state (RAII-compatible)
- **Iteration trace** — per-step TSV output (Minuit2 only, via `if constexpr`)
- **Error budget** — systematic component breakdown via group fix/release
- **Contour & chi2 scans** — 1/2/3-sigma contours, special contours, 1D scans

## Quick Start

```cpp
#include "TMinimization.h"
#include "Math/Functor.h"

double my_chi2(const double* par) {
    return (par[0] - 3.0) * (par[0] - 3.0) / (0.5 * 0.5);
}

int main() {
    TMini<> mini;  // defaults to Minuit2 backend
    mini.AddParam("x", 1.0, 0.1).SetLimits(0.0, 10.0);

    ROOT::Math::Functor f(&my_chi2, 1);
    mini.SetFunction(f);

    FitResult result = mini.FitAndGetResult();
    result.Print();
    result.SaveJSON("result.json");
}
```

## Derived Class Pattern

Override `InitializeParams()` for application-specific parameter setup:

```cpp
class MyMinimizer : public TMini<MinimizerBackend::Minuit2> {
public:
    using TMini::TMini;  // inherit constructors

    unsigned int InitializeParams() override {
        unsigned int idx = 0;
        AddParameterInfo(idx, "mu", 0.0, 0.01, "physics");
        AddParameterInfo(idx, "sigma", 1.0, 0.01, "systematics");
        paras_groups["stat"].insert(0);
        return idx;
    }
};
```

## Backend Selection

```cpp
// Compile-time selection via template parameter
TMini<MinimizerBackend::Minuit2> mini2;  // Minuit2 (default)
TMini<MinimizerBackend::Minuit1> mini1;  // Legacy TMinuit

// Type aliases
TMiniMinuit2 m2;  // = TMini<MinimizerBackend::Minuit2>
TMiniMinuit1 m1;  // = TMini<MinimizerBackend::Minuit1>
```

Minuit2-only features (iteration trace, Fumili2) compile to no-ops on the Minuit1 backend via `if constexpr`.

## CMake Integration

### As a subdirectory

```cmake
add_subdirectory(external/Minuit2-Extended)
target_link_libraries(MyTarget PRIVATE minuit2_extended)
```

### Building examples

```bash
cd external/Minuit2-Extended
cmake -B build -DBUILD_EXAMPLES=ON
cmake --build build
```

## API Summary

### Core
| Method | Description |
|--------|-------------|
| `TMini(level, tolerance)` | Constructor (default: print_level=2, tol=1e-6) |
| `InitializeParams()` | Virtual, override in derived classes |
| `Fit(print_level)` | Run minimization |
| `GetFcn()` | Current minimum value |

### Algorithm Selection
| Method | Description |
|--------|-------------|
| `SetAlgorithm(algo)` | Set algorithm (Migrad/Simplex/Combined/Fumili2) |
| `FitWith(algo, print_level)` | One-shot fit with specific algorithm |

### Structured Results
| Method | Description |
|--------|-------------|
| `GetResult()` | Get FitResult struct after fit |
| `FitAndGetResult()` | Fit + return result in one call |

### Parameter Handles
| Method | Description |
|--------|-------------|
| `AddParam(name, val, step)` | Add parameter, returns ParameterHandle |
| `Par(idx)` / `Par(name)` | Get handle for existing parameter |

### State Management
| Method | Description |
|--------|-------------|
| `SaveState()` | Capture current state as StateSnapshot |
| `RestoreState(snap)` | Restore from snapshot |

### Advanced Statistics
| Method | Description |
|--------|-------------|
| `ProfileLikelihood(idx, n, lo, hi)` | 1D profile scan |
| `GridScan(i, j, nx, ny, ...)` | 2D grid scan |
| `ConfidenceInterval(idx, delta_chi2)` | MINOS errors |
| `GetCorrelationMatrix()` | Full correlation matrix |
| `PValue(ndf)` | p-value |

### Error Budget & Contours
| Method | Description |
|--------|-------------|
| `ExportErrorBudget(...)` | Systematic breakdown by parameter groups |
| `ExportCovMatrix(...)` | Hesse covariance with degenerate recovery |
| `CreateStandardContour(...)` | 1/2/3-sigma contours |
| `Scan_Chi2(...)` | 1D chi2 scan |
| `MultiStartFit(config)` | Multi-start global optimization |

## Examples

| Example | Demonstrates |
|---------|-------------|
| `basic_fit.cc` | Gaussian MLE, FitResult, ParameterHandle |
| `algorithms.cc` | Migrad vs Simplex vs Combined |
| `contour_scan.cc` | 2D contour + 1D chi2 scan |
| `profile_likelihood.cc` | Profile scan + confidence intervals |
| `grid_scan.cc` | 2D grid scan to TH2D |
| `error_budget.cc` | Virtual InitializeParams + error budget |
| `multi_start.cc` | Global optimization via random restarts |
| `state_snapshot.cc` | Save/restore RAII pattern |
| `iteration_trace.cc` | TSV trace (Minuit2 backend) |

## Dependencies

- **ROOT 6.0+** (Minuit2, MathCore, Hist)
- **C++17** compiler

## Version History

- **v2.0.0** — Template backend, virtual InitializeParams, FitResult, ParameterHandle, ProfileLikelihood, GridScan, MultiStart, StateSnapshot
- **v1.5** — Iteration trace, config-fixed params, error groups
- **v1.0** — Basic Minuit2 wrapper

## License

MIT
