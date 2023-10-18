/**
 * This codes performs a simple extension on the Arellano (2008) model. 
 * We propose a model in which recovery rates on defaulted bonds can vary.
 * Author: Lucas Belmudes 10/17/2023.
*/

#include <iostream>
#include "economy.hpp"
#include "initialization.hpp"
#include "auxiliary.hpp"

int main(){

    // Set parameters:
    const int b_grid_size = 251;        // Number of points in the grid for the bond price.
    const double b_grid_min = -0.5;     // Minimum value of the bond price.
    const double b_grid_max = 0.00;     // Maximum value of the bond price.
    const int y_grid_size = 20;         // Number of points in the grid for the income.
    const double y_default = 0.969;     // Maximum income under default.
    const double beta = 0.953;          // Discount factor.
    const double gamma = 2;             // Risk aversion.
    const double r = 0.017;             // Interest rate.
    const double rho = 0.945;           // Persistence of the income.
    const double sigma = 0.025;         // Standard deviation of the income.
    const double theta = 0.282;         // Probability of a re-entry.
    const double alpha = 0.00;          // Recovery on defaulted debt.
    const int max_iter = 1000;          // Maximum number of iterations.
    const double tol = 1e-6;            // Tolerance for the convergence.
    const double m = 3;                 // Number of standard deviations for the income grid.

    // Set pointers to store the results of the model:
    double* y_grid = new double[y_grid_size];
    double* y_grid_default = new double[y_grid_size];
    double* b_grid = new double[b_grid_size];
    double* p = new double[y_grid_size*y_grid_size];
    double* v = new double[y_grid_size*b_grid_size];
    double* v_r = new double[y_grid_size*b_grid_size];
    double* v_d = new double[y_grid_size*b_grid_size];
    double* q = new double[y_grid_size*b_grid_size];
    int* b_policy = new int[y_grid_size*b_grid_size];
    double* d_policy = new double[y_grid_size*b_grid_size];

  
    // Create an instance of the Economy class:
    Economy economy(b_grid_size, b_grid_min, b_grid_max, y_grid_size, y_default, beta, gamma, r, rho, sigma, theta, alpha, tol, max_iter, m, y_grid, y_grid_default, b_grid, p, v, v_r, v_d, q, b_policy, d_policy);

    // Initialize economy:
    if (economy.initialize_economy() == EXIT_SUCCESS){
        std::cout << "Economy initialized successfully." << std::endl;
    } else {
        std::cout << "Economy initialization failed." << std::endl;
    }

    // Solve the model:
    economy.solve_model();

    // Display results:
    std::cout << "Bond price" << std::endl;
    display_matrix(economy.Q, economy.Y_grid_size, economy.B_grid_size);
    std::cout << "Value function" << std::endl;
    display_matrix(economy.V, economy.Y_grid_size, economy.B_grid_size);





    // Free memory:
    delete[] y_grid;
    delete[] y_grid_default;
    delete[] b_grid;
    delete[] p;
    delete[] v;
    delete[] v_r;
    delete[] v_d;
    delete[] q;
    delete[] b_policy;
    delete[] d_policy;

    return 0;

}