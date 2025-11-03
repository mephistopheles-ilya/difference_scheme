#pragma once

#include <cmath>

inline double H_024 (double x)
{
    if (x < 4.5 || x > 5.5)
    {
        return 1;
    }
    else 
    {
        return 2;
    }
}

inline double V_024 (double /* x */)
{
    return 0;
}

inline double H_025 (double /* x */)
{
    return 1;
}

inline double V_025 (double x)
{
    if (x < 4.5 || x > 5.5)
    {
        return 0;
    }
    else
    {
        return 1;
    }
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

