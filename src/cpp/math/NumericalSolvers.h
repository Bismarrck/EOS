#ifndef NUMERICAL_SOLVERS_H
#define NUMERICAL_SOLVERS_H

#include <vector>
#include <functional>

namespace EOS_Toolkit { // Or a name like EOS_Math
namespace NumericalSolvers {

// Common result struct for root finders
struct RootResult {
  double root = 0.0;
  bool converged = false;
  int iterations = 0;
};

// Brent's root-finding method.
// Finds a root of func(x) = 0 within the bracket [x1, x2].
// func(x1) and func(x2) must have opposite signs.
RootResult brents_method(
    const std::function<double(double)>& func,
    double x1,
    double x2,
    int max_iter = 100,
    double tolerance = 1e-9
);

// Type definition for the state of an ODE system (e.g., [T, other_vars...])
using ODEState = std::vector<double>;
// Type definition for the ODE system function dy/dx = f(x, y)
using ODESystem = std::function<void(const ODEState& y, ODEState& dydx, const double x)>;
// Type definition for the observer function called at each step
using ODEObserver = std::function<void(const ODEState& y, const double x)>;

// Solves an ODE system using a constant-step Runge-Kutta 4th order method.
// Integrates from x_start to x_end with step size dx.
void integrate_ode_rk4(
    ODESystem& system,
    ODEState& start_state,
    double x_start,
    double x_end,
    double dx,
    ODEObserver& observer
);

// the trapezoidal rule for ODE integration (RK2)
void integrate_ode_heun(
    ODESystem& system,
    ODEState& start_state,
    double x_start,
    double x_end,
    double dx,
    ODEObserver& observer
);

} // namespace NumericalSolvers
} // namespace EOS_Toolkit

#endif // NUMERICAL_SOLVERS_H