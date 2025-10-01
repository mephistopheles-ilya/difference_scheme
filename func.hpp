#pragma once

#include <cmath>

inline double H_0 (double x)
{
    return std::cos (3. * M_PI * x) + 1.5;
}

inline double V_0 (double x)
{
    return std::sin (4. * M_PI * x);
}

inline double H (double t, double x)
{
    return std::exp (t) * (std::cos (3. * M_PI * x) + 1.5);
}

inline double V (double t, double x)
{
    return std::cos (2 * M_PI * t) * std::sin ( 4 * M_PI * x);
}

inline double f_0 (double t, double x)
{
    double dro_dt = std::exp (t) * (std::cos (3 * M_PI * x) + 1.5); 
    double drou_dx = std::exp(t) * std::cos (2 * M_PI * t) * (4 * M_PI * std::cos (3 * M_PI * x) * std::cos (4 * M_PI * x)
             - 3 * M_PI * std::sin (3 * M_PI * x) * std::sin (4 * M_PI * x)) + 6 * M_PI * std::exp (t) * std::cos (2 * M_PI * t) * cos (4 * M_PI * x);
    return dro_dt + drou_dx;
}

inline double f_lin (double t, double x, double C, double mu)
{
    double ro = std::exp (t) * (std::cos (3 * M_PI * x) + 1.5);
    double drou_dt = (std::cos (3 * M_PI * x) * std::sin (4 * M_PI * x) + 1.5 * std::sin (4 * M_PI * x))
         * std::exp (t) * (std::cos (2 * M_PI * t) - 2 * M_PI * std::sin (2 * M_PI * t));
    double drou2_dx = std::exp (t) * std::cos (2 * M_PI * t) * std::cos (2 * M_PI * t) * (-3 * M_PI * std::sin (3 * M_PI * x) * std::sin (4 * M_PI * x)
            * std::sin (4 * M_PI * x) + 8 * M_PI * std::sin (4 * M_PI * x) * std::cos (4 * M_PI * x) * std::cos (3 * M_PI * x))
        + 12 * M_PI * std::exp (t) * std::cos (2 * M_PI * t) * std::cos (2 * M_PI * t) * std::sin (4 * M_PI * x) * std::cos (4 * M_PI * x);
    double d2u_dx2 = -1 * mu * 16 * M_PI * M_PI * std::cos (2 * M_PI * t) * std::sin (4 * M_PI * x);
    double dp_dx = C * std::exp (t) * (-3 * M_PI) * std::sin (3 * M_PI * x);
    return (drou_dt + drou2_dx + dp_dx - d2u_dx2) / ro;
}

inline double f_pow (double t, double x, double gamma, double mu)
{
    double ro = std::exp (t) * (std::cos (3 * M_PI * x) + 1.5);
    double drou_dt = (std::cos (3 * M_PI * x) * std::sin (4 * M_PI * x) + 1.5 * std::sin (4 * M_PI * x))
         * std::exp (t) * (std::cos (2 * M_PI * t) - 2 * M_PI * std::sin (2 * M_PI * t));
    double drou2_dx = std::exp (t) * std::cos (2 * M_PI * t) * std::cos (2 * M_PI * t) * (-3 * M_PI * std::sin (3 * M_PI * x) * std::sin (4 * M_PI * x)
            * std::sin (4 * M_PI * x) + 8 * M_PI * std::sin (4 * M_PI * x) * std::cos (4 * M_PI * x) * std::cos (3 * M_PI * x))
        + 12 * M_PI * std::exp (t) * std::cos (2 * M_PI * t) * std::cos (2 * M_PI * t) * std::sin (4 * M_PI * x) * std::cos (4 * M_PI * x);
    double d2u_dx2 = -1 * mu * 16 * M_PI * M_PI * std::cos (2 * M_PI * t) * std::sin (4 * M_PI * x);
    double dp_dx = gamma * std::pow ( std::exp (t) * (std::cos (3 * M_PI * x) + 1.5), gamma - 1) * std::exp (t) * (-3 * M_PI) * std::sin (3 * M_PI * x);
    return (drou_dt + drou2_dx + dp_dx - d2u_dx2) / ro;
}

inline double p_lin (double x, double C)
{
    return C * x;
}

inline double p_pow (double x, double gamma)
{
    return std::pow (x, gamma);
}

