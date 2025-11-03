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





int main(int argc, char *argv[])
{
    //feenableexcept(FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW | FE_UNDERFLOW);

    if (argc != 6)
    {
        printf("Wrong number of arguments\n");
        printf("Usage: %s h tau mu pp is_lin_p \n", argv[0]);
        return 1;
    }

    double h = 0;
    double tau = 0;
    double mu = 0;
    double pp = 0;
    int  is_lin_p = 0;

    if (!(sscanf(argv[1], "%lf", &h) == 1
       && sscanf(argv[2], "%lf", &tau) == 1
       && sscanf(argv[3], "%lf", &mu) == 1
       && sscanf(argv[4], "%lf", &pp) == 1
       && sscanf(argv[5], "%d",  &is_lin_p) == 1
       ))
    {
        printf("Wrong type of arguments\n");
        printf("Usage: %s h tau mu pp is_lin_p nested_param\n", argv[0]);
        return 2;
    }

    double (*P) (double, double) = nullptr;
    double (*F) (double, double, double, double) = nullptr;
    if (is_lin_p != 0)
    {
        P = p_lin;
        F = f_lin;
    }
    else 
    {
        P = p_pow;
        F = f_pow;
    }
    

    double x_a = 0;
    double x_b = 1;
    double t_a = 0;
    double t_b = 1;
    size_t x_N = std::floor ((x_b - x_a) / h);
    size_t t_N = std::floor ((t_b - t_a) / tau);

    std::vector<double> up_diag (x_N, 0);
    std::vector<double> diag (x_N + 1, 0);
    std::vector<double> low_diag (x_N, 0);
    std::vector<double> rhs (x_N + 1, 0);

    std::vector<double> H_solution_prev (x_N + 1, 0);
    std::vector<double> H_solution (x_N + 1, 0);
    std::vector<double> V_solution (x_N + 1, 0);
    
    
    fill_H_initial (H_solution_prev, x_a, x_b, h, H_0);
    fill_V_initial (V_solution, x_a, x_b, h, V_0);

    double t = t_a;
    double time1 = 0, time2 = 0;
    time1 = clock();
    for (size_t i = 0; i < t_N + 1; ++i) 
    {
        t = (i == t_N) ? t_b : t_a + i * tau;

        fill_H_matrix (up_diag, diag, low_diag, rhs, H_solution_prev /* n */, V_solution /* n */, x_a, x_b, h, t, tau, f_0); 
        solve_tree_diag (up_diag, diag, low_diag, rhs, H_solution /* n + 1 */);

        fill_V_matrix (up_diag, diag, low_diag, rhs, H_solution_prev /* n */, H_solution /* n + 1 */, V_solution /* n */, x_a, x_b, h, t, tau, F, mu, P, pp); 
        solve_tree_diag (up_diag, diag, low_diag, rhs, V_solution /* n + 1 */);

        std::swap(H_solution_prev, H_solution);
    }
    std::swap(H_solution_prev, H_solution);
    time2 = clock();
    t = (time2 - time1)/CLOCKS_PER_SEC;


#if 0
    std::string file_name = "H_res.txt";
    print_solution(file_name.data(), H_solution_prev);
    file_name = "V_res.txt";
    print_solution(file_name.data(), V_solution);
#endif

    sub_real_solution (H_solution, t_b, x_a, x_b, h, H); 
    sub_real_solution (V_solution, t_b, x_a, x_b, h, V);

    double C_norm_H = calc_C_norm (H_solution);
    double C_norm_V = calc_C_norm (V_solution);
    double L_norm_H = calc_L_norm (H_solution, h);
    double L_norm_V = calc_L_norm (V_solution, h);
    double W_norm_H = calc_W_norm (H_solution, h);
    double W_norm_V = calc_W_norm (V_solution, h);

    printf("CFL = %lf\n", tau / h);
    printf ("h = %lf tau = %lf mu = %lf pp = %lf is_lin_p = %d\n", h, tau, mu, pp, is_lin_p);
    printf ("C_norm_H = %e  C_norm_V = %e \nL_norm_H = %e  L_norm_V = %e \nW_norm_H = %e  W_norm_V = %e \n"
            , C_norm_H, C_norm_V, L_norm_H, L_norm_V, W_norm_H, W_norm_V);
    printf ("time = %e \n", t);

    return 0;
}
