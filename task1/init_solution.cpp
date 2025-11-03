#include "stddef.h"
#include "init_solution.hpp"


void fill_H_initial (std::vector<double> &H_solution, double a, double b, double h, double (*H) (double))
{
    size_t sz = H_solution.size();
    for (size_t i = 1; i < sz - 1; ++i)
    {
        H_solution[i] = H (a + h * i);
    }
    H_solution[0] = H (a);
    H_solution[sz - 1] = H (b);
}


void fill_V_initial (std::vector<double> &V_solution, double a, double b, double h, double (*V) (double), bool stiking)
{
    size_t sz = V_solution.size();
    for (size_t i = 1; i < sz - 1; ++i)
    {
        V_solution[i] = V (a + h * i);
    }
    if (stiking == true)
    {
        V_solution[0] = 0;
        V_solution[sz - 1] = 0;
    }
    else
    {
        //TODO if no stiking
        V_solution[0] = V (a);
        V_solution[sz - 1] = V (b);
    }
}
