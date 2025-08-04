#include "NumericalSolvers.h"
#include <iostream>  // For error messages
#include <cmath>     // For std::abs
#include <limits>    // For std::numeric_limits
#include <utility>   // For std::swap

namespace EOS_Toolkit {
namespace NumericalSolvers {

// Implementation of Brent's method (copy from tools/eos-tool.cpp)
RootResult brents_method(const std::function<double(double)>& func, double x1, double x2, int max_iter, double tolerance) {
  // ... (full implementation as before) ...
 RootResult result;
  double a = x1;
  double b = x2;
  double fa = func(a);
  double fb = func(b);

  if (fa * fb >= 0) {
    std::cerr << "Error (Brent's Method): Root not bracketed. f(a)=" << fa
              << ", f(b)=" << fb << std::endl;
    result.converged = false;
    return result;
  }

  if (std::abs(fa) < std::abs(fb)) {
    std::swap(a, b);
    std::swap(fa, fb);
  }

  double c = a;  // Contrapoint
  double fc = fa;
  double d = 0;          // Step from previous iteration
  double m_flag = true;  // Bisection flag

  for (int i = 0; i < max_iter; ++i) {
    result.iterations = i + 1;

    if (std::abs(b - a) < tolerance || std::abs(fb) < tolerance) {
      result.root = b;
      result.converged = true;
      return result;
    }

    double s = 0;
    if (fa != fc && fb != fc) {
      // Inverse quadratic interpolation
      s = (a * fb * fc / ((fa - fb) * (fa - fc))) +
          (b * fa * fc / ((fb - fa) * (fb - fc))) +
          (c * fa * fb / ((fc - fa) * (fc - fb)));
    } else {
      // Secant method
      s = b - fb * (b - a) / (fb - fa);
    }

    double half_interval = (c - b) / 2.0;
    if (!((s > b && s < b + half_interval) ||
          (s < b && s > b - half_interval)) ||  // s not in bracket
        (m_flag &&
         (std::abs(s - b) >= std::abs(b - c) / 2.0)) ||  // Bisection was slow
        (!m_flag &&
         (std::abs(s - b) >= std::abs(c - d) / 2.0)) ||  // Secant/IQI was slow
        (m_flag && (std::abs(b - c) < tolerance)) ||
        (!m_flag && (std::abs(c - d) < tolerance))) {
      // Fallback to bisection
      s = (a + b) / 2.0;
      m_flag = true;
    } else {
      m_flag = false;
    }

    double fs = func(s);
    d = c;
    c = b;
    fc = fb;

    if (fa * fs < 0) {
      b = s;
      fb = fs;
    } else {
      a = s;
      fa = fs;
    }

    if (std::abs(fa) < std::abs(fb)) {
      std::swap(a, b);
      std::swap(fa, fb);
    }
  }

  std::cerr << "Warning (Brent's Method): Failed to converge after " << max_iter
            << " iterations." << std::endl;
  result.root = b;
  result.converged = false;
  return result;
}


// Implementation of ODE integration wrapper
void integrate_ode_rk4(
    ODESystem& system,
    ODEState& start_state,
    double x_start,
    double x_end,
    double dx,
    ODEObserver& observer
) {
  if (dx == 0.0) {
    std::cerr << "Error (RK4): Step size dx cannot be zero." << std::endl;
    return;
  }
  if ((dx > 0 && x_start > x_end) || (dx < 0 && x_start < x_end)) {
    std::cerr << "Warning (RK4): Integration range [x_start, x_end] is inconsistent with sign of step size dx." << std::endl;
    // Continue anyway, the loop condition will handle it.
  }

  ODEState y = start_state;
  ODEState dydx(y.size());
  ODEState k1(y.size()), k2(y.size()), k3(y.size()), k4(y.size());
  ODEState y_temp(y.size());

  double x = x_start;
  int num_steps = static_cast<int>(std::abs((x_end - x_start) / dx));

  // Call observer for the initial state
  observer(y, x);

  for (int i = 0; i < num_steps; ++i) {
    // k1 = f(x, y)
    system(y, dydx, x);
    for (size_t j = 0; j < y.size(); ++j) k1[j] = dx * dydx[j];

    // k2 = f(x + dx/2, y + k1/2)
    for (size_t j = 0; j < y.size(); ++j) y_temp[j] = y[j] + 0.5 * k1[j];
    system(y_temp, dydx, x + 0.5 * dx);
    for (size_t j = 0; j < y.size(); ++j) k2[j] = dx * dydx[j];

    // k3 = f(x + dx/2, y + k2/2)
    for (size_t j = 0; j < y.size(); ++j) y_temp[j] = y[j] + 0.5 * k2[j];
    system(y_temp, dydx, x + 0.5 * dx);
    for (size_t j = 0; j < y.size(); ++j) k3[j] = dx * dydx[j];

    // k4 = f(x + dx, y + k3)
    for (size_t j = 0; j < y.size(); ++j) y_temp[j] = y[j] + k3[j];
    system(y_temp, dydx, x + dx);
    for (size_t j = 0; j < y.size(); ++j) k4[j] = dx * dydx[j];

    // Update y
    for (size_t j = 0; j < y.size(); ++j) {
      y[j] += (k1[j] + 2.0 * k2[j] + 2.0 * k3[j] + k4[j]) / 6.0;
    }

    x += dx;
    // Call observer with the new state
    observer(y, x);
  }
}

void integrate_ode_heun(
    ODESystem& system,
    ODEState& start_state,
    double x_start,
    double x_end,
    double dx,
    ODEObserver& observer
) {
  if (dx == 0.0) {
    std::cerr << "Error (Heun): Step size dx cannot be zero." << std::endl;
    return;
  }
  // ... (check for inconsistent x_start, x_end, dx signs as in RK4) ...

  ODEState y = start_state;
  ODEState dydx(y.size());
  ODEState y_predictor(y.size());
  ODEState dydx_predictor(y.size());

  double x = x_start;
  int num_steps = static_cast<int>(std::abs((x_end - x_start) / dx));

  // Call observer for the initial state
  observer(y, x);

  for (int i = 0; i < num_steps; ++i) {
    // 1. Calculate slope at current point: k1 = f(x, y)
    system(y, dydx, x); // dydx now holds f(x, y)

    // 2. Predictor step: y_tilde = y + h * f(x, y)
    for (size_t j = 0; j < y.size(); ++j) {
      y_predictor[j] = y[j] + dx * dydx[j];
    }

    // 3. Calculate slope at predicted point: k2 = f(x+h, y_tilde)
    system(y_predictor, dydx_predictor, x + dx); // dydx_predictor now holds f(x+h, y_tilde)

    // 4. Corrector step: Update y using average of slopes
    for (size_t j = 0; j < y.size(); ++j) {
      y[j] += (dx / 2.0) * (dydx[j] + dydx_predictor[j]);
    }

    x += dx;
    // Call observer with the new state
    observer(y, x);
  }
}

} // namespace NumericalSolvers
} // namespace EOS_Toolkit