#include <stddef.h>
#include "solve.hpp"



void solve_tree_diag (std::vector<double> &up_diag, std::vector<double> &diag, std::vector<double> &low_diag
        , std::vector<double> &rhs, std::vector<double> &H_solution)
{
    size_t n = diag.size();
    up_diag[0] = up_diag[0] / diag[0];
    rhs[0] = rhs[0] / diag[0];
    for (size_t i = 1; i < n - 1; ++i)
    {
        double den = diag[i] - low_diag[i - 1] * up_diag[i - 1];
        up_diag[i] = up_diag[i] / den;
        rhs[i] = (rhs[i] - low_diag[i - 1] * rhs[i - 1])/ den;
    }

    H_solution[n - 1] = (rhs[n - 1] - low_diag[n - 2] * rhs[n - 2]) / (diag[n - 1] - low_diag[n - 2] * up_diag[n - 2]);

    for (size_t i = n - 1; i-- > 0;)
    {
        H_solution[i] = rhs[i] - up_diag[i] * H_solution[i + 1];
    }
}


