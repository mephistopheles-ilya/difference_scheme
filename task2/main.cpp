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


static double max_time = 6000;
static double epsilon = 3 * 1e-3;


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
    
    
    fill_H_initial (H_solution_prev, x_a, x_b, h, H_024);
    fill_V_initial (V_solution, x_a, x_b, h, V_024);

    double initial_massa = calc_mass(H_solution_prev);

    double stab_time = t_a;
    double time1 = 0, time2 = 0;
    double stab_norm = 1e64;
    size_t stab_index = 0;
    time1 = clock();
    for (stab_index = 0; stab_norm > epsilon && stab_time < max_time;) 
    {
        stab_time = t_a + stab_index * tau;

        fill_H_matrix (up_diag, diag, low_diag, rhs, H_solution_prev /* n */, V_solution /* n */, x_a, x_b, h, stab_time, tau, f_0); 
        solve_tree_diag (up_diag, diag, low_diag, rhs, H_solution /* n + 1 */);

        fill_V_matrix (up_diag, diag, low_diag, rhs, H_solution_prev /* n */, H_solution /* n + 1 */, V_solution /* n */, x_a, x_b, h, stab_time, tau, F, mu, P, pp); 
        solve_tree_diag (up_diag, diag, low_diag, rhs, V_solution /* n + 1 */);

        stab_norm = stability_norm(H_solution, V_solution);
#if 0
        if (stab_index % 10 == 0)
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

    double res_mass = calc_mass(H_solution);
    double delta_mass = (res_mass - initial_massa) / initial_massa;


#if 0
    std::string file_name = "H_res.txt";
    print_solution(file_name.data(), H_solution_prev);
    file_name = "V_res.txt";
    print_solution(file_name.data(), V_solution);
#endif

    printf("CFL = %lf\n", tau / h);
    printf("h = %lf tau = %lf mu = %lf pp = %lf is_lin_p = %d\n", h, tau, mu, pp, is_lin_p);
    printf("stab_norm = %e, stab_time = %lf\n, stab_index = %lu\n", stab_norm, stab_time, stab_index);
    printf("delta_massa = %e\n", delta_mass);
    printf("time = %e \n", calc_time);

    return 0;
}
