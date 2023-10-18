#ifndef economy_hpp
#define economy_hpp

#include <iostream>

class Economy{

    public:

        // Parameters for this economy;
        int B_grid_size;       // Number of points in the grid for the bond price.
        double B_grid_min;     // Minimum value of the bond price.
        double B_grid_max;     // Maximum value of the bond price.
        int Y_grid_size;       // Number of points in the grid for the income.
        double Y_default;      // Maximum income under default.
        double Beta;           // Discount factor.
        double Gamma;          // Risk aversion.
        double R;              // Interest rate.
        double Rho;            // Persistence of the income.
        double Sigma;          // Standard deviation of the income.
        double Theta;          // Probability of a re-entry.
        double Alpha;          // Recovery on defaulted debt.
        double Tol;            // Tolerance for the convergence.
        int Max_iter;          // Maximum number of iterations.
        double M;              // Number of standard deviations for the income grid.

        // Set pointers:
        double* Y_grid;
        double* Y_grid_default;
        double* B_grid;
        double* P;
        double* V;
        double* V_r;
        double* V_d;
        double* Q;
        int* B_policy;
        double* D_policy;

        // Constructor:
        Economy(int b_grid_size,  double b_grid_min,  double b_grid_max,  int y_grid_size,  double y_default,  double beta,  double gamma,  double r,  double rho,  double sigma,  double theta,  double alpha,  double tol,  int max_iter,  double m, double* ptr_y_grid, double* ptr_y_grid_default, double* ptr_b_grid, double* ptr_p_grid, double* ptr_v, double* ptr_v_r, double* ptr_v_d, double* ptr_q, int* ptr_b_policy, double* ptr_d_policy);

        // Initialize economy:
        int initialize_economy();

        // Solution Methods:
            // Guess:
            void guess_vd_vr_q();

            // Update v and default policy:
            void update_v_and_default_policy();

            // Update price:
            void update_price();

            // Update value at default:
            void update_vd();

            // Update value of repayment and bond policy:
            void update_vr_and_bond_policy();

            // Solve:
            int solve_model();
};

#endif