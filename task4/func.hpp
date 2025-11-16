#pragma once

#include <cmath>

static int __K__ = 0;

inline double H_027 (double x)
{
    return 2 + std::sin (__K__ * M_PI * x);
}

inline double V_027 (double /* x */)
{
    return 0;
}

inline double H_028 (double /* x */)
{
    return 1;
}

inline double V_028 (double x)
{
    return std::sin (__K__ * M_PI * x);
}



inline double f_0 (double /* t */, double /* x */)
{
    return 0;
}

inline double f (double /* t */, double /* x */, double /* C */, double /* mu */)
{
    return 0;
}

inline double p_lin (double x, double C)
{
    return C * x;
}

inline double p_pow (double x, double gamma)
{
    return std::pow (x, gamma);
}

