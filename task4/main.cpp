#include <algorithm>
#include <cmath>
#include <stdio.h>
#include <iostream>
#include <vector>
#include <fenv.h>

#include "init_solution.hpp"
#include "fill_matrix.hpp"
#include "solve.hpp"
#include "func.hpp"
#include "norms.hpp"



static const double max_time = 6000;
static const double epsilon = 1e-4;
static const int stab_const = 50;


int main(int argc, char *argv[])
{
    //feenableexcept(FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW | FE_UNDERFLOW);

    if (argc != 8)
    {
        printf("Wrong number of arguments\n");
        printf("Usage: %s h tau mu pp is_lin_p v_left h_left \n", argv[0]);
        return 1;
    }

    double h = 0;
    double tau = 0;
    double mu = 0;
    double pp = 0;
    int  is_lin_p = 0;
    double v_left = 0;
    double h_left = 0;

    if (!(sscanf(argv[1], "%lf", &h) == 1
       && sscanf(argv[2], "%lf", &tau) == 1
       && sscanf(argv[3], "%lf", &mu) == 1
       && sscanf(argv[4], "%lf", &pp) == 1
       && sscanf(argv[5], "%d",  &is_lin_p) == 1
       && sscanf(argv[6], "%lf", &v_left) == 1
       && sscanf(argv[7], "%lf",  &h_left) == 1
       ))
    {
        printf("Wrong type of arguments\n");
        printf("Usage: %s h tau mu pp is_lin_p v_left h_left \n", argv[0]);
        return 2;
    }

    double (*P) (double, double) = nullptr;
    double (*F) (double, double, double, double) = nullptr;
    F = f;
    if (is_lin_p != 0)
    {
        P = p_lin;
    }
    else 
    {
        P = p_pow;
    }

    double x_a = 0;
    double x_b = 10;
    double t_a = 0;
    size_t x_N = std::floor ((x_b - x_a) / h);

    std::vector<double> up_diag (x_N, 0);
    std::vector<double> diag (x_N + 1, 0);
    std::vector<double> low_diag (x_N, 0);
    std::vector<double> rhs (x_N + 1, 0);

    std::vector<double> H_solution_prev (x_N + 1, 0);
    std::vector<double> H_solution (x_N + 1, 0);
    std::vector<double> V_solution (x_N + 1, 0);
    
    fill_H_initial (H_solution_prev, x_a, x_b, h, H_0);
    fill_V_initial (V_solution, x_a, x_b, h, V_0);

    std::vector<double> stab_H (x_N + 1, 0);
    std::vector<double> stab_V (x_N + 1, 0);

    int stab_count = stab_const;
    double stab_time = t_a;
    double time1 = 0, time2 = 0;
    double stab_norm = 1e64;
    int stab_index = 0;
    time1 = clock();
    for (stab_index = 0; stab_count > 0 && stab_time < max_time;) 
    {
        stab_time = t_a + stab_index * tau;

        fill_H_matrix (up_diag, diag, low_diag, rhs, H_solution_prev /* n */, V_solution /* n */, x_a, x_b, h, stab_time, tau, f_0, h_left); 
        solve_tree_diag (up_diag, diag, low_diag, rhs, H_solution /* n + 1 */);

        fill_V_matrix (up_diag, diag, low_diag, rhs, H_solution_prev /* n */, H_solution /* n + 1 */, V_solution /* n */, x_a, x_b, h, stab_time, tau, F, mu, P, pp, v_left); 
        solve_tree_diag (up_diag, diag, low_diag, rhs, V_solution /* n + 1 */);

        double norm1 = C_norm(H_solution, stab_H);
        double norm2 = C_norm(V_solution, stab_V);
        stab_norm = std::max (norm1, norm2);
        if (stab_norm < epsilon)
        {
            stab_count -= 1;
        }
        else
        {
            stab_H = H_solution;
            stab_V = V_solution;
            stab_count = stab_const;
        }
#if 0
        if (stab_index % 5000 == 0)
        {
            printf("current_stab_norm = %lf\n", stab_norm);
        }
#endif
        std::swap(H_solution_prev, H_solution);
        stab_index += 1;
    }
    std::swap(H_solution_prev, H_solution);
    stab_index -= 1;
    time2 = clock();
    double calc_time = (time2 - time1)/CLOCKS_PER_SEC;

#if 0
    std::string file_name = "H_res.txt";
    print_solution(file_name.data(), H_solution_prev);
    file_name = "V_res.txt";
    print_solution(file_name.data(), V_solution);
#endif

    printf("h = %lf tau = %lf mu = %lf pp = %lf is_lin_p = %d, v_left = %lf, h_left = %lf\n", h, tau, mu, pp, is_lin_p, v_left, h_left);
    printf("stab_norm = %e\n", stab_norm);
    printf("stab_time = %lf\n", stab_time);
    printf("stab_index = %d\n", stab_index);
    printf("time = %e \n", calc_time);

    return 0;
}
