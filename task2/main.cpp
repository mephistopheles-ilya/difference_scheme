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


#define PRINT_V 1
#define PRINT_H 0

static double max_time = 6000;
static double epsilon = 3 * 1e-3;


int main(int argc, char *argv[])
{
    //feenableexcept(FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW | FE_UNDERFLOW);

    if (argc != 9)
    {
        printf("Wrong number of arguments\n");
        printf("Usage: %s h tau mu pp is_lin_p count_to_stop task_num nested_param\n", argv[0]);
        return 1;
    }

    double h = 0;
    double tau = 0;
    double mu = 0;
    double pp = 0;
    int  is_lin_p = 0;
    int count_to_stop = 0;
    int task_num = 0;
    int nested_param = 0;

    if (!(sscanf(argv[1], "%lf", &h) == 1
       && sscanf(argv[2], "%lf", &tau) == 1
       && sscanf(argv[3], "%lf", &mu) == 1
       && sscanf(argv[4], "%lf", &pp) == 1
       && sscanf(argv[5], "%d",  &is_lin_p) == 1
       && sscanf(argv[6], "%d",  &count_to_stop) == 1
       && sscanf(argv[7], "%d",  &task_num) == 1
       && sscanf(argv[8], "%d",  &nested_param) == 1
       ))
    {
        printf("Wrong type of arguments\n");
        printf("Usage: %s h tau mu pp is_lin_p count_to_stap task_num nested_param\n", argv[0]);
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
    
    
    if (task_num == 24)
    {
        fill_H_initial (H_solution_prev, x_a, x_b, h, H_024);
        fill_V_initial (V_solution, x_a, x_b, h, V_024);
    }
    if (task_num == 25)
    {
        fill_H_initial (H_solution_prev, x_a, x_b, h, H_025);
        fill_V_initial (V_solution, x_a, x_b, h, V_025);
    }


    double initial_massa = calc_mass(H_solution_prev);

    std::string file_name = "2d.txt";
    FILE* file = fopen(file_name.data(), "w");
    double stab_time = t_a;
    double time1 = 0, time2 = 0;
    double stab_norm = 1e64;
    int stab_index = 0;
    time1 = clock();
    for (stab_index = 0; stab_norm > epsilon && stab_time < max_time && (stab_index <= count_to_stop || count_to_stop < 0);) 
    {
        stab_time = t_a + stab_index * tau;

        fill_H_matrix (up_diag, diag, low_diag, rhs, H_solution_prev /* n */, V_solution /* n */, x_a, x_b, h, stab_time, tau, f_0); 
        solve_tree_diag (up_diag, diag, low_diag, rhs, H_solution /* n + 1 */);

        fill_V_matrix (up_diag, diag, low_diag, rhs, H_solution_prev /* n */, H_solution /* n + 1 */, V_solution /* n */, x_a, x_b, h, stab_time, tau, F, mu, P, pp); 
        solve_tree_diag (up_diag, diag, low_diag, rhs, V_solution /* n + 1 */);

        stab_norm = stability_norm(H_solution, V_solution);
#if 1
        if (stab_index % 1000 == 0)
        {
            printf("current_stab_norm = %lf\n", stab_norm);
        }
#endif
        if (PRINT_V)
        {
            size_t sz = V_solution.size();
            for (size_t i = 0; i < sz; ++i)
            {
                fprintf(file, "%e ", V_solution[i]); 
            }
            fprintf(file, "\n");
        }
        if (PRINT_H)
        {
            size_t sz = H_solution.size();
            for (size_t i = 0; i < sz; ++i)
            {
                fprintf(file, "%e ", H_solution[i]); 
            }
            fprintf(file, "\n");
        }


        std::swap(H_solution_prev, H_solution);
        stab_index += 1;
    }
    fclose(file);
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
    std::vector<double> H_solution_nested;
    std::vector<double> V_solution_nested;
    if (nested_param != 0)
    {
        int scaling = std::pow (2, nested_param);
        double h_nested = h / scaling;
        double tau_nested = tau / scaling;
        double t_b = tau * count_to_stop;
        size_t x_N_nested = std::floor ((x_b - x_a) / h_nested);
        size_t t_N_nested = std::floor ((t_b - t_a) / tau_nested);

        up_diag.resize (x_N_nested);
        std::fill (up_diag.begin (), up_diag.end (), 0);
        diag.resize (x_N_nested + 1);
        std::fill(diag.begin (), diag.end (), 0);
        low_diag.resize (x_N_nested);
        std::fill (low_diag.begin (), low_diag.end (), 0);
        rhs.resize (x_N_nested + 1);
        std::fill (rhs.begin (), rhs.end (), 0);

        H_solution_prev.resize (x_N_nested + 1);
        std::fill (H_solution_prev.begin (), H_solution_prev.end (), 0);
        H_solution_nested.resize (x_N_nested + 1);
        std::fill (H_solution_nested.begin (), H_solution_nested.end (), 0);
        V_solution_nested.resize (x_N_nested + 1);
        std::fill (V_solution_nested.begin (), V_solution_nested.end (), 0);

        if (task_num == 24)
        {
            fill_H_initial (H_solution_prev, x_a, x_b, h_nested, H_024);
            fill_V_initial (V_solution_nested, x_a, x_b, h_nested, V_024);
        }
        if (task_num == 25)
        {
            fill_H_initial (H_solution_prev, x_a, x_b, h_nested, H_025);
            fill_V_initial (V_solution_nested, x_a, x_b, h_nested, V_025);
        }
        for (size_t i = 0; i < t_N_nested + 1; ++i)
        {
            double t = (i == t_N_nested) ? t_b : t_a + i * tau_nested;
            fill_H_matrix (up_diag, diag, low_diag, rhs, H_solution_prev /* n */, V_solution_nested /* n */, x_a, x_b, h_nested, t, tau_nested, f_0); 
            solve_tree_diag (up_diag, diag, low_diag, rhs, H_solution_nested /* n + 1 */);

            fill_V_matrix (up_diag, diag, low_diag, rhs, H_solution_prev /* n */, H_solution_nested /* n + 1 */, V_solution_nested /* n */, x_a, x_b, h_nested, t
                    , tau_nested, F, mu, P, pp); 
            solve_tree_diag (up_diag, diag, low_diag, rhs, V_solution_nested /* n + 1 */);

            std::swap(H_solution_prev, H_solution_nested);
        }

        std::swap(H_solution_prev, H_solution_nested);
        
        for (size_t i = 0; i <  x_N + 1; ++i)
        {
            H_solution[i] -= H_solution_nested[i*scaling];
            V_solution[i] -= V_solution_nested[i*scaling];
        }

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
    }
    else
    {
        printf("h = %lf tau = %lf mu = %lf pp = %lf is_lin_p = %d\n", h, tau, mu, pp, is_lin_p);
        printf("stab_norm = %e\n", stab_norm);
        printf("stab_time = %lf\n", stab_time);
        printf("stab_index = %d\n", stab_index);
        printf("delta_massa = %e\n", delta_mass);
        printf("time = %e \n", calc_time);
    }

    return 0;
}
