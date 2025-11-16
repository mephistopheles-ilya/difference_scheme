#include "norms.hpp"

#include <cmath>
#include <stdio.h>

void print_solution (char *file_name, std::vector<double> &solution)
{
    FILE* file = fopen(file_name, "w");
    size_t sz = solution.size();
    for (size_t i = 0; i < sz; ++i)
    {
        fprintf(file, "%e\n", solution[i]); 
    }
    fclose(file);
}


double calc_C_norm(std::vector<double> &vec)
{
    size_t sz = vec.size();
    double max = 0;
    for (size_t i = 0; i < sz; ++i)
    {
        if (std::abs(vec[i]) > max)
        {
            max = std::abs(vec[i]);
        }
    }
    return max;
}

double C_norm (std::vector<double> &vec1, std::vector<double> &vec2)
{
    size_t sz = vec1.size();
    double max = 0;
    for (size_t i = 0; i < sz; ++i)
    {
        max = std::max (max, std::fabs(vec1[i] - vec2[i]));
    }
    return max;
}


double calc_L_norm(std::vector<double> &vec, double h)
{
    size_t sz = vec.size();
    double sum = 0;
    for (size_t i = 1; i < sz - 1; ++i)
    {
        sum += vec[i] * vec[i];
    }
    sum *= h;
    sum += 0.5 * h * (vec[0] * vec[0] + vec[sz - 1] * vec[sz - 1]);
    sum = std::sqrt(sum);
    return sum;
}

double calc_W_norm (std::vector<double> &vec, double h)
{
    size_t sz = vec.size();
    double sum = 0;
    for (size_t i = 1; i < sz - 1; ++i)
    {
        sum += vec[i] * vec[i];
    }
    sum *= h;
    sum += 0.5 * h * (vec[0] * vec[0] + vec[sz - 1] * vec[sz - 1]);
    double sum1 = 0;
    for (size_t i = 0; i < sz - 1; ++i)
    {
        double dv_dx = (vec[i + 1] - vec[i]) / h;
        sum1 += dv_dx * dv_dx;
    }
    sum1 *= h;
    return std::sqrt (sum + sum1);
}

void sub_real_solution (std::vector<double>& solution, double t, double a, double b, double h, double (*f) (double, double))
{
    size_t sz = solution.size();
    for (size_t i = 0; i < sz - 1; ++i)
    {
        solution[i] -= f (t, a + i * h);
    }
    solution[sz - 1] -= f (t, b);
}

double calc_discrepancy (std::vector<double> &up_diag, std::vector<double> &diag, std::vector<double> &low_diag, std::vector<double> &rhs
        , std::vector<double> &solution)
{
    size_t n = diag.size();
    rhs[0] -= solution[0] * diag[0] + solution[1] * up_diag[0];
    for (size_t i = 1; i < n - 1; ++i)
    {
        rhs[i] -= solution[i - 1] * low_diag[i - 1] + solution[i] * diag[i] + solution[i + 1] * up_diag[i];
    }
    rhs[n - 1] -= low_diag[n - 2] * solution[n - 2] + diag[n - 1] * solution[n - 1];

    double norm = 0;
    for (size_t i = 0; i < n; ++i)
    {
        norm += rhs[i] * rhs[i];
    }
    return std::sqrt(norm);
}

double stability_norm (std::vector<double> &H_solution, std::vector<double> &V_solution)
{
    double H_av = 0;
    double V_av = 0;
    double sz = H_solution.size();
    for (size_t i = 1; i < sz - 1; ++i)
    {
        H_av += H_solution[i];
    }
    H_av /= sz;
    double H_norm = 0;
    double V_norm = 0;
    for (size_t i = 0; i < sz; ++i)
    {
        H_norm = std::max(fabs(H_solution[i] - H_av), H_norm);
        V_norm = std::max(fabs(V_solution[i] - V_av), V_norm);
    }
    return std::max(H_norm, V_norm);
}

double calc_mass(std::vector<double> &H_solution)
{
    double sum = 0;
    for (auto H: H_solution)
    {
        sum += H;
    }
    return sum;
}

