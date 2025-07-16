#include "NumericalSolvers.h"
#include <boost/numeric/odeint.hpp> // Odeint headers
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
  // Use a Runge-Kutta 4 stepper from Boost.Odeint
  boost::numeric::odeint::runge_kutta4<ODEState> stepper;

  // Use integrate_const to take constant steps
  // This will call the observer after each successful step.
  boost::numeric::odeint::integrate_const(
      stepper,
      system,
      start_state,
      x_start,
      x_end,
      dx,
      observer
  );
}

} // namespace NumericalSolvers
} // namespace EOS_Toolkit