#pragma once

#include <vector>


void fill_H_initial (std::vector<double> &H_solution, double a, double b, double h, double (*H) (double));
void fill_V_initial (std::vector<double> &V_solution, double a, double b, double h, double (*V) (double), bool stiking = true);
