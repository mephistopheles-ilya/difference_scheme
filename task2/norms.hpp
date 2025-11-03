#pragma once

#include <vector>

double calc_C_norm (std::vector<double> &vec);
double calc_L_norm (std::vector<double> &vec, double h);
double calc_W_norm (std::vector<double> &vec, double h);


void sub_real_solution (std::vector<double>& solution, double t, double a, double b, double h, double (*f) (double, double));
double calc_discrepancy (std::vector<double> &up_diag, std::vector<double> &diag, std::vector<double> &low_diag, std::vector<double> &rhs
        , std::vector<double> &solution);
void print_solution (char *file_name, std::vector<double> &solution);
