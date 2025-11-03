#pragma once

#include <vector>

void fill_H_matrix (std::vector<double> &up_diag, std::vector<double> &diag, std::vector<double> &low_diag
        , std::vector<double> &rhs, std::vector<double> &H_solution_prev, std::vector<double> &V_solution_prev
        , double a, double b, double h, double t, double tau, double (*f) (double, double));

void fill_V_matrix (std::vector<double> &up_diag, std::vector<double> &diag, std::vector<double> &low_diag
        , std::vector<double> &rhs, std::vector<double> &H_solution_prev, std::vector<double> &H_solution, std::vector<double> &V_solution_prev
        , double a, double b, double h, double t, double tau, double (*f) (double, double, double, double), double mu, double (*p) (double, double)
        , double pressure_param );
