#include "economy.hpp"
#include "Initialization.hpp"
#include "auxiliary.hpp"
#include <iostream>
#include <cmath>
#include <mex.h>
#include <omp.h>

//! Note: For computations we first iterate over the low recovery bond and then iterate over the high recovery bond.
//! Vectors are then V=[(b_low1;b_high1;y_1), (b_low2;b_high1;y_1),...,(b_lown;b_high1;y_1),(b_low1;b_high2;y_1),...

Economy::Economy(int b_grid_size, double b_grid_min, double b_grid_max, int y_grid_size, double y_default, double beta, double gamma, double r, double rho, double sigma, double theta, double alpha_low, double alpha_high, double tol, int max_iter, double m, double* ptr_y_grid, double* ptr_y_grid_default, double* ptr_b_grid, double* ptr_p_grid, double* ptr_v, double* ptr_v_r, double* ptr_v_d, double* ptr_q_low, double* ptr_q_high, int* ptr_b_policy_low, int* ptr_b_policy_high, double* ptr_d_policy){
  
    // Parameters:
    B_grid_size = b_grid_size;            // Number of points in the grid for the bond price.   
    B_grid_min = b_grid_min;              // Minimum value of the bond price.
    B_grid_max = b_grid_max;              // Maximum value of the bond price.
    Y_grid_size = y_grid_size;            // Number of points in the grid for the income.
    Y_default = y_default;                // Maximum income under default.
    Beta = beta;                          // Discount factor.
    Gamma = gamma;                        // Risk aversion.
    R = r;                                // Interest rate.
    Rho = rho;                            // Persistence of the income.
    Sigma = sigma;                        // Standard deviation of the income.
    Theta = theta;                        // Probability of a re-entry.
    Alpha_low = alpha_low;                // Recovery rate for the low recovery debt.
    Alpha_high = alpha_high;              // Recovery rate for the high recovery debt.
    Tol = tol;                            // Tolerance for the convergence.
    Max_iter = max_iter;                  // Maximum number of iterations.
    M = m;                                // Number of standard deviations for the income grid.
        
    // Name the pointer to the arrays we will be working:
    Y_grid = ptr_y_grid;                  // Income grid.
    Y_grid_default = ptr_y_grid_default;  // Income grid for the default state.
    B_grid = ptr_b_grid;                  // Bond price grid.
    P = ptr_p_grid;                       // Transition matrix.
    V = ptr_v;                            // Value function.
    V_r = ptr_v_r;                        // Value function under re-entry.
    V_d = ptr_v_d;                        // Value function under default.
    Q_low = ptr_q_low;                    // Price for the low recovery debt.
    Q_high = ptr_q_high;                  // Price for the high recovery debt.
    B_policy_low = ptr_b_policy_low;      // Bond policy for the low recovery debt.
    B_policy_high = ptr_b_policy_high;    // Bond policy for the high recovery debt.
    D_policy = ptr_d_policy;              // Default policy.
    
}

// Create grids and store it in the space previously allocated:
int Economy::initialize_economy(){
    // Create the grid for the income:
    create_bond_grids(B_grid, B_grid_size, B_grid_max, B_grid_min);
    if (B_grid[B_grid_size - 1] > Tol || B_grid[B_grid_size - 1] < -Tol){
        mexPrintf("Error: the bond grid is not correctly initialized..\n");
        return EXIT_FAILURE;
    }
    // Create the grid for the income and probability matrix:
    create_income_and_prob_grids(Y_grid, P, Y_grid_size, Sigma, Rho, M);
    for (int i = 0; i < Y_grid_size; i++){
        if (Y_grid[i]<=0){
            mexPrintf("Error: the income grid is not correctly initialized.\n");
            return EXIT_FAILURE;
        }
    }
    for (int i=0; i< Y_grid_size; i++){
        double prob = 0;
        for (int j=0; j<Y_grid_size; j++){
            prob += P[i*Y_grid_size+j];
        }
        if (prob > 1+Tol || prob < 1-Tol){
            mexPrintf("Error: the probability matrix is not correctly initialized\n");
            return EXIT_FAILURE;
        }
    }
    // Create the income under default:
    create_income_under_default(Y_grid_default, Y_grid, Y_grid_size, Y_default);
    return EXIT_SUCCESS;
}


// Guess value function at default, value at reentry and price:
void Economy::guess_vd_vr_q(){
    for (int i=0; i<Y_grid_size; i++)
    {
        for (int j=0; j<B_grid_size; j++)
        {
            for (int z=0; z<B_grid_size; z++)
            {
                V_d[i*(B_grid_size*B_grid_size)+j*B_grid_size+z] = -20;
                V_r[i*(B_grid_size*B_grid_size)+j*B_grid_size+z] = -20;
                Q_low[i*(B_grid_size*B_grid_size)+j*B_grid_size+z] = 1/(1+R);
                Q_high[i*(B_grid_size*B_grid_size)+j*B_grid_size+z] = 1/(1+R);
            }
        }
    }
}


// Update value function and default policy:
void Economy::update_v_and_default_policy(){
    for (int i=0; i<Y_grid_size; i++)
    {
        for (int j=0; j<B_grid_size; j++)
        {
            for (int z=0; z<B_grid_size; z++)
            {
                double V_r_aux = V_r[i*(B_grid_size*B_grid_size)+j*B_grid_size+z];
                double V_d_aux = V_d[i*(B_grid_size*B_grid_size)+j*B_grid_size+z];
                if (V_r_aux>=V_d_aux)
                {
                    V[i*(B_grid_size*B_grid_size)+j*B_grid_size+z] = V_r_aux;
                    D_policy[i*(B_grid_size*B_grid_size)+j*B_grid_size+z] = 0;
                } else {
                    V[i*(B_grid_size*B_grid_size)+j*B_grid_size+z] = V_d_aux;
                    D_policy[i*(B_grid_size*B_grid_size)+j*B_grid_size+z] = 1;
                }    
            }
        }
    }
}
 
// Update prices given a default policies;
void Economy::update_price(){
    for (int i=0; i<Y_grid_size; i++)
    {
        for (int j=0; j<B_grid_size; j++)
        {
            for (int z=0; z<B_grid_size; z++)
            {
                double aux_low = 0;
                double aux_high = 0;
                for (int i_prime = 0; i_prime < Y_grid_size; i_prime++)
                {
                    aux_low += P[i*Y_grid_size+i_prime] * ((1-D_policy[i_prime*(B_grid_size*B_grid_size)+j*B_grid_size+z]) + D_policy[i_prime*(B_grid_size*B_grid_size)+j*B_grid_size+z] * Alpha_low) *  (1/(1+R));
                    aux_high += P[i*Y_grid_size+i_prime] * ((1-D_policy[i_prime*(B_grid_size*B_grid_size)+j*B_grid_size+z]) + D_policy[i_prime*(B_grid_size*B_grid_size)+j*B_grid_size+z] * Alpha_high) *  (1/(1+R));
                }
                Q_low[i*(B_grid_size*B_grid_size)+j*B_grid_size+z] = aux_low;
                Q_high[i*(B_grid_size*B_grid_size)+j*B_grid_size+z] = aux_high;
            }
        }
    }
}

// Update value at default:
void Economy::update_vd(){
    double* Vd0 = new double[Y_grid_size * B_grid_size* B_grid_size];      // Store initial value function at default:
    copy_vector(V_d, Vd0, Y_grid_size * B_grid_size * B_grid_size);
    for (int i=0; i<Y_grid_size; i++)
    {
        for (int j=0; j<B_grid_size; j++)
        {
            for (int z=0; z<B_grid_size; z++)
            {
                double E_V = 0;
                double E_Vd = 0;
                for (int i_prime = 0; i_prime < Y_grid_size; i_prime++)
                {
                    E_V += P[i*Y_grid_size+i_prime] * V[i_prime*(B_grid_size*B_grid_size)+(B_grid_size-1)*B_grid_size+(B_grid_size-1)];           // Expected value given exclusion and zero debt.
                    E_Vd += P[i*Y_grid_size+i_prime] * V_d[i_prime*(B_grid_size*B_grid_size)+(B_grid_size-1)*B_grid_size+(B_grid_size-1)];         // Expected value given exclusion and zero debt.      
                }
                V_d[i*(B_grid_size*B_grid_size)+j*B_grid_size+z] = utility(Y_grid_default[i] + Alpha_low * B_grid[z] + Alpha_high * B_grid[j], Gamma, Tol) + Beta * (Theta * E_V + (1-Theta) * E_Vd); // Payoff of recovery in current period.
            }
        }
    }
    delete[] Vd0;
}

// Update value of repayment and bond policy:
void Economy::update_vr_and_bond_policy(){
    for (int i=0; i<Y_grid_size; i++)
    {
        for (int j=0; j<B_grid_size; j++)
        {
            for (int z=0; z<B_grid_size; z++) // Condition on (y,b_low,b_high).
            {
                double aux_v = -1000000000;                          
                for (int x_low = 0; x_low<B_grid_size; x_low++)
                {
                    for (int x_high = 0; x_high<B_grid_size; x_high++)
                    {  
                        double E_V_rx = 0;                                      // Expected continuation value of repayment following policy x_low and x_high.
                        double c = Y_grid[i] - Q_low[i*(B_grid_size*B_grid_size)+x_high*B_grid_size+x_low] * B_grid[x_low] - Q_high[i*(B_grid_size*B_grid_size)+x_high*B_grid_size+x_low] * B_grid[x_high] + B_grid[j] + B_grid[z];
                        if (c > Tol){
                            for (int i_prime = 0; i_prime < Y_grid_size; i_prime++){
                                E_V_rx += P[i*Y_grid_size+i_prime] * V[i_prime*(B_grid_size*B_grid_size)+x_high*B_grid_size+x_low];
                            }
                            double temp = utility(c, Gamma, Tol) + Beta * E_V_rx;
                            if (temp >= aux_v){
                                aux_v = temp;
                                B_policy_low[i*(B_grid_size*B_grid_size)+j*B_grid_size+z] = x_low;
                                B_policy_high[i*(B_grid_size*B_grid_size)+j*B_grid_size+z] = x_high;
                            }
                        }
                }
                    V_r[i*(B_grid_size*B_grid_size)+j*B_grid_size+z] = aux_v;
                    if (aux_v == -1000000000)
                    {
                        B_policy_low[i*(B_grid_size*B_grid_size)+j*B_grid_size+z] = -5;
                        B_policy_high[i*(B_grid_size*B_grid_size)+j*B_grid_size+z] = -5;
                    }
                }
            }
        }
    }
}

// Solve the model:
int Economy::solve_model(){
    // Initialize economy:
    if (initialize_economy() == EXIT_SUCCESS){
        mexPrintf("Economy initialized successfully.\n");
    } else {
        mexPrintf("Economy initialization failed.\n");
        return EXIT_FAILURE;
    }
    // Guess value functions and prices and copy initial values:
    guess_vd_vr_q();    
    double* Vd0 = new double[Y_grid_size *  B_grid_size * B_grid_size];     // Store initial value function at default:
    copy_vector(V_d, Vd0, Y_grid_size * B_grid_size * B_grid_size);
    double* Vr0 = new double[Y_grid_size *  B_grid_size * B_grid_size];     // Store initial value function at reentry:
    copy_vector(V_r, Vr0, Y_grid_size * B_grid_size * B_grid_size);
    double* Q0_low = new double[Y_grid_size *  B_grid_size * B_grid_size];      // Store initial default probability:
    copy_vector(Q_low, Q0_low, Y_grid_size * B_grid_size * B_grid_size);
    double* Q0_high = new double[Y_grid_size *  B_grid_size * B_grid_size];      // Store initial default probability:
    copy_vector(Q_high, Q0_high, Y_grid_size * B_grid_size * B_grid_size);

    // Initialize difference between value functions:
    int iter = 0;          
    double diff_q_low = 1;
    double diff_q_high = 1;
    double diff_vd = 1;
    double diff_vr = 1;
    double aux_q_low, aux_q_high, aux_vd, aux_vr;

    while (iter < Max_iter){
      
        update_v_and_default_policy();                  // Update v and default policy:
        update_price();                                 // update price:
        update_vd();                                    // update value at default:
        update_vr_and_bond_policy();                    // update value of repayment and bond policy:
        
        diff_q_low = 0;                                 // Difference between prices.
        diff_q_high = 0;                                // Difference between prices.
        diff_vd = 0;                                    // Difference between value function at default.
        diff_vr = 0;                                    // Difference between value function at reentry.
        
        for (int i=0; i<Y_grid_size; i++)
        {
            for (int j=0; j<B_grid_size; j++)
            {
                for (int z=0; z<B_grid_size; z++)
                {
                    aux_q_low = fabs(Q0_low[i*B_grid_size*B_grid_size+j*B_grid_size+z] - Q_low[i*B_grid_size*B_grid_size+j*B_grid_size+z]);
                    aux_q_high = fabs(Q0_high[i*B_grid_size*B_grid_size+j*B_grid_size+z] - Q_high[i*B_grid_size*B_grid_size+j*B_grid_size+z]);
                    aux_vd = fabs(V_d[i*B_grid_size*B_grid_size+j*B_grid_size+z] - Vd0[i*B_grid_size*B_grid_size+j*B_grid_size+z]);
                    aux_vr = fabs(V_r[i*B_grid_size*B_grid_size+j*B_grid_size+z] - Vr0[i*B_grid_size*B_grid_size+j*B_grid_size+z]);
                    if (aux_q_low > diff_q_low){
                        diff_q_low = aux_q_low;
                    }
                    if (aux_q_high > diff_q_high){
                        diff_q_high = aux_q_high;
                    }
                    if (aux_vd > diff_vd){
                        diff_vd = aux_vd;
                    }
                    if (aux_vr > diff_vr){
                        diff_vr = aux_vr;
                    }
                }
            }
        }

        if (diff_q_low < Tol && diff_q_high < Tol && diff_vd < Tol && diff_vr < Tol){
            mexPrintf("Iteration: %d\n", iter);
            mexPrintf("Difference between value function at default: %f\n", diff_vd);
            mexPrintf("Difference between value function at reentry: %f\n", diff_vr);
            mexPrintf("Difference between low prices: %f\n", diff_q_low);
            mexPrintf("Difference between high prices: %f\n", diff_q_high);
            mexPrintf("Convergence achieved after: %d\n", iter);
            // Free memory:
            delete[] Vd0;
            delete[] Vr0;
            delete[] Q0_low;
            delete[] Q0_high;
            return EXIT_SUCCESS;

        } else {
            if (iter % 25 == 0){
                mexPrintf("Iteration: %d\n", iter);
                mexPrintf("Difference between value function at default: %f\n", diff_vd);
                mexPrintf("Difference between value function at reentry: %f\n", diff_vr);
                mexPrintf("Difference between low prices: %f\n", diff_q_low);
                mexPrintf("Difference between high prices: %f\n", diff_q_high);
            }
            // Update value functions and prices:
            copy_vector(V_d, Vd0, Y_grid_size * B_grid_size * B_grid_size);
            copy_vector(V_r, Vr0, Y_grid_size * B_grid_size * B_grid_size);
            copy_vector(Q_low, Q0_low, Y_grid_size * B_grid_size * B_grid_size);
            copy_vector(Q_high, Q0_high, Y_grid_size * B_grid_size * B_grid_size);
            iter += 1;
        }
    }
    // If the model does not converge:
    mexPrintf("Convergence not achieved after: %d\n", iter);
    mexPrintf("Difference between value function at default: %f\n", diff_vd);
    mexPrintf("Difference between value function at reentry: %f\n", diff_vr);
    mexPrintf("Difference between low prices: %f\n", diff_q_low);
    mexPrintf("Difference between high prices: %f\n", diff_q_high);
    // Free memory:
    delete[] Vd0;
    delete[] Vr0;
    delete[] Q0_low;
    delete[] Q0_high;
    return EXIT_FAILURE;
}
