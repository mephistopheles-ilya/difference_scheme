#include <stddef.h>
#include "fill_matrix.hpp"

void fill_H_matrix (std::vector<double> &up_diag, std::vector<double> &diag, std::vector<double> &low_diag
        , std::vector<double> &rhs, std::vector<double> &H_solution_prev, std::vector<double> &V_solution_prev
        , double a, double b, double h, double t, double tau, double (*f) (double, double), double h_left)
{
    size_t sz = diag.size();
    double CFL = tau / h;
    for (size_t m = 1; m < sz - 1; ++m)
    {
        low_diag[m - 1] = - 0.25 * CFL * (V_solution_prev[m] + V_solution_prev[m - 1]);
        diag[m] = 1;
        up_diag[m] = 0.25 * CFL * (V_solution_prev[m] + V_solution_prev[m + 1]);
        rhs[m] = H_solution_prev[m] - 0.25 * CFL * H_solution_prev[m] * (V_solution_prev[m + 1] - V_solution_prev[m - 1]) + tau * f (t,  a + m * h);
    }
    diag[0] = 1;
    up_diag[0] = 0;
    rhs[0] = h_left;
    diag[sz - 1] = 1 + 0.5 * CFL * V_solution_prev[sz - 1] ;
    low_diag[sz - 2] = -0.5 * CFL * V_solution_prev[sz - 2];
    rhs[sz - 1] = H_solution_prev[sz - 1] - 0.5 * CFL * H_solution_prev[sz - 1] * (V_solution_prev[sz - 1] - V_solution_prev[sz - 2]) - 0.5 * CFL
        * (H_solution_prev[sz - 3] * V_solution_prev[sz - 3] - 2 * H_solution_prev[sz - 2] * V_solution_prev[sz - 2] + H_solution_prev[sz - 1] * V_solution_prev[sz - 1] - 0.5
        * (H_solution_prev[sz - 4] * V_solution_prev[sz - 4] - 2 * H_solution_prev[sz - 3] * V_solution_prev[sz - 3] + H_solution_prev[sz - 2] * V_solution_prev[sz - 2])
        + H_solution_prev[sz - 1] * (V_solution_prev[sz - 3] - 2 * V_solution_prev[sz - 2] + V_solution_prev[sz - 1] - 0.5
        * (V_solution_prev[sz - 4] - 2 * V_solution_prev[sz - 3] + V_solution_prev[sz - 2]))) + tau * f(t, b);
}


void fill_V_matrix (std::vector<double> &up_diag, std::vector<double> &diag, std::vector<double> &low_diag
        , std::vector<double> &rhs, std::vector<double> &H_solution_prev, std::vector<double> &H_solution, std::vector<double> &V_solution_prev
        , double a, double /*b*/, double h, double t, double tau, double (*f) (double, double, double , double), double mu, double (*p) (double, double)
        , double pressure_param, double v_left)
{
    size_t sz = diag.size();
    double CFL = tau / h;
    for (size_t m = 1; m < sz - 1; ++m)
    {
        low_diag[m - 1] = - (1./3. * CFL * H_solution[m] * V_solution_prev[m] + 1./3. * CFL * H_solution[m - 1] * V_solution_prev[m - 1] + mu * CFL / h);
        diag[m] = H_solution[m] + 2 * mu * CFL / h;
        up_diag[m] = 1./3. * CFL * H_solution[m + 1] * V_solution_prev[m + 1] + 1./3. * CFL * H_solution[m] * V_solution_prev[m] - mu * CFL / h;
        rhs[m] = H_solution_prev[m] * V_solution_prev[m] - 1./6. * CFL * V_solution_prev[m] * V_solution_prev[m] * (H_solution[m + 1] - H_solution[m - 1])
             - 0.5 * CFL * (p(H_solution[m + 1], pressure_param) - p(H_solution[m - 1], pressure_param)) + tau * H_solution[m] * f (t, a + m * h, pressure_param, mu);
    }
    diag[0] = 1;
    rhs[0] = v_left;
    diag[sz - 1] = 1;
    rhs[sz - 1] = 0;
    up_diag[0] = 0;
    low_diag[sz - 2] = -1;
}

