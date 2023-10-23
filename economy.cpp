#include "economy.hpp"
#include "Initialization.hpp"
#include "auxiliary.hpp"
#include <iostream>
#include <cmath>


Economy::Economy(int b_grid_size,  double b_grid_min,  double b_grid_max,  int y_grid_size,  double y_default,  double beta,  double gamma,  double r,  double rho,  double sigma,  double theta,  double alpha,  double tol,  int max_iter,  double m, double* ptr_y_grid, double* ptr_y_grid_default, double* ptr_b_grid, double* ptr_p_grid, double* ptr_v, double* ptr_v_r, double* ptr_v_d, double* ptr_q, int* ptr_b_policy, double* ptr_d_policy){
  
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
    Alpha = alpha;                        // Recovery on defaulted debt.
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
    Q = ptr_q;                            // Default probability.
    B_policy = ptr_b_policy;              // Bond price policy.
    D_policy = ptr_d_policy;              // Default policy.
}

// Create grids and store it in the space previously allocated:
int Economy::initialize_economy(){
    // Create the grid for the income:
    create_bond_grids(B_grid, B_grid_size, B_grid_max, B_grid_min);
    if (B_grid[B_grid_size - 1] > Tol || B_grid[B_grid_size - 1] < -Tol){
        std::cout << "Error: the bond grid is not correctly initialized." << std::endl;
        return EXIT_FAILURE;
    }
    // Create the grid for the income and probability matrix:
    create_income_and_prob_grids(Y_grid, P, Y_grid_size, Sigma, Rho, M);
    for (int i = 0; i < Y_grid_size; i++){
        if (Y_grid[i]<=0){
            std::cout << "Error: the income grid is not correctly initialized." << std::endl;
            return EXIT_FAILURE;
        }
    }
    for (int i=0; i< Y_grid_size; i++){
        double prob = 0;
        for (int j=0; j<Y_grid_size; j++){
            prob += P[i*Y_grid_size+j];
        }
        if (prob > 1+Tol || prob < 1-Tol){
            std::cout << "Error: the probability matrix is not correctly initialized." << std::endl;
            return EXIT_FAILURE;
        }
    }
    // Create the income under default:
    create_income_under_default(Y_grid_default, Y_grid, Y_grid_size, Y_default);
    return EXIT_SUCCESS;
}

// Guess value function at default, value at reentry and price:
void Economy::guess_vd_vr_q(){
    for (int i=0; i<Y_grid_size; i++){
        for (int j=0; j<B_grid_size; j++){
            //V_d[i*B_grid_size+j] = utility( Y_grid_default[i] + Alpha * B_grid[j], Gamma, Tol) + (Beta/(1-Beta)) * utility(Y_grid_default[i], Gamma, Tol);
            V_d[i*B_grid_size+j] = 0 ;
            //V_r[i*B_grid_size+j] = utility( R * B_grid[j] + Y_grid[i], Gamma, Tol) + (Beta/(1-Beta)) * utility(R * B_grid[j] + Y_grid[i], Gamma, Tol);
            V_r[i*B_grid_size+j] = 0;
            //Q[i*B_grid_size+j] = 1/(1+R) * ((double)j/((double)B_grid_size-1));
            Q[i*B_grid_size+j] = 1/(1+R);
        }
    }
}

// Update value function and default policy:
void Economy::update_v_and_default_policy(){
    for (int i=0; i<Y_grid_size; i++){
        for (int j=0; j<B_grid_size; j++){
            double V_r_aux = V_r[i*B_grid_size+j];
            double V_d_aux = V_d[i*B_grid_size+j];
            if (V_r_aux >= V_d_aux){
                V[i*B_grid_size+j] = V_r_aux;
                D_policy[i*B_grid_size+j] = 0;
            } else {
                V[i*B_grid_size+j] = V_d_aux;
                D_policy[i*B_grid_size+j] = 1;
            }
        }
    }
}

// Update price given a default policy:
void Economy::update_price(){
    for (int i=0; i<Y_grid_size; i++){
        for (int j=0; j<B_grid_size; j++){
        double aux = 0;
            for (int i_prime = 0; i_prime < Y_grid_size; i_prime++){
                aux += P[i*Y_grid_size+i_prime] * ((1-D_policy[i_prime*B_grid_size+j])+D_policy[i_prime*B_grid_size+j]*Alpha)* (1/(1+R));          
            }
        Q[i*B_grid_size+j] = aux; 
        }
    }
}

// Update value at default:
void Economy::update_vd(){
    double* Vd0 = new double[Y_grid_size * B_grid_size];
    copy_vector(V_d, Vd0, Y_grid_size * B_grid_size);
    for (int i=0; i<Y_grid_size; i++){
        for (int j=0; j<B_grid_size; j++){
            double E_V = 0;
            double E_Vd = 0;
            for (int i_prime = 0; i_prime < Y_grid_size; i_prime++){
                E_V += P[i*Y_grid_size+i_prime] * V[i_prime*B_grid_size+(B_grid_size-1)];    // Expected value given reentry and zero debt.
                E_Vd += P[i*Y_grid_size+i_prime] * Vd0[i_prime*B_grid_size+(B_grid_size-1)]; // Expected value given exclusion and zero debt.        
            }
            V_d[i*B_grid_size+j] = utility(Y_grid_default[i] + Alpha * B_grid[j], Gamma, Tol) + Beta * (Theta * E_V + (1-Theta) * E_Vd); // Payoff of recovery in current period.
        }
    }
    delete[] Vd0;
}

// Update value of repayment and bond policy:
void Economy::update_vr_and_bond_policy(){
    for (int i=0; i<Y_grid_size; i++){
        for (int j=0; j<B_grid_size; j++){
            double aux_v = -1000000000; 
            // Loop over possible bond policies.                          
            for (int x = 0; x<B_grid_size; x++){    
                double E_V_rx = 0;                  // Expected continuation value of repayment following policy x.
                double c = Y_grid[i] - Q[i*B_grid_size+x] * B_grid[x] + B_grid[j];
                if (c >= Tol){
                    for (int i_prime = 0; i_prime < Y_grid_size; i_prime++){
                        E_V_rx += P[i*Y_grid_size+i_prime] * V[i_prime*B_grid_size+x];
                    }
                    double temp = utility(c, Gamma, Tol) + Beta * E_V_rx;
                    if (temp >= aux_v){
                        aux_v = temp;
                        B_policy[i*B_grid_size+j] = x;
                    }
                }
            }
            V_r[i*B_grid_size+j] = aux_v;
        }
    }
}

// Solve the model:
int Economy::solve_model(){
    
    // Initialize economy:
    if (initialize_economy() == EXIT_SUCCESS){
        std::cout << "Economy initialized successfully." << std::endl;
    } else {
        std::cout << "Economy initialization failed." << std::endl;
        return EXIT_FAILURE;
    }

    // Guess value functions and prices and copy initial values:
    guess_vd_vr_q();    
    double* Vd0 = new double[Y_grid_size * B_grid_size];      // Store initial value function at default:
    copy_vector(V_d, Vd0, Y_grid_size * B_grid_size);
    double* Vr0 = new double[Y_grid_size * B_grid_size];      // Store initial value function at reentry:
    copy_vector(V_r, Vr0, Y_grid_size * B_grid_size);
    double* Q0 = new double[Y_grid_size * B_grid_size];       // Store initial default probability:
    copy_vector(Q, Q0, Y_grid_size * B_grid_size);

    // Initialize difference between value functions:
    int iter = 0;          
    double diff_q = 1;
    double diff_vd = 1;
    double diff_vr = 1;
    double aux_q, aux_vd, aux_vr;

    while (iter < Max_iter){
      
        update_v_and_default_policy();                  // Update v and default policy:
        update_price();                                 // update price:
        update_vd();                                    // update value at default:
        update_vr_and_bond_policy();                    // update value of repayment and bond policy:
        
        diff_q = 0;                                     // Difference between prices.
        diff_vd = 0;                                    // Difference between value function at default.
        diff_vr = 0;                                    // Difference between value function at reentry.
        
        for (int i=0; i<Y_grid_size; i++){
            for (int j=0; j<B_grid_size; j++){
                aux_q = fabs(Q[i*B_grid_size+j] - Q0[i*B_grid_size+j]);
                aux_vd = fabs(V_d[i*B_grid_size+j] - Vd0[i*B_grid_size+j]);
                aux_vr = fabs(V_r[i*B_grid_size+j] - Vr0[i*B_grid_size+j]);
                if (aux_q > diff_q){
                    diff_q = aux_q;
                }
                if (aux_vd > diff_vd){
                    diff_vd = aux_vd;
                }
                if (aux_vr > diff_vr){
                    diff_vr = aux_vr;
                }
            }
        }

        if (diff_q < Tol && diff_vd < Tol && diff_vr < Tol){
            std::cout << "Convergence achieved after " << iter << " iterations." << std::endl;
            std::cout << "Difference between value function at default: " << diff_vd << std::endl;
            std::cout << "Difference between value function at reentry: " << diff_vr << std::endl;
            std::cout << "Difference between prices: " << diff_q << std::endl;
            // Free memory:
            delete[] Vd0;
            delete[] Vr0;
            delete[] Q0;
            return EXIT_SUCCESS;

        } else {
            //if (iter % 25 == 0){
            if (false){
                std::cout << "Iteration: " << iter << std::endl;
                std::cout << "Difference between value function at default: " << diff_vd << std::endl;
                std::cout << "Difference between value function at reentry: " << diff_vr << std::endl;
                std::cout << "Difference between prices: " << diff_q << std::endl;
            }
            // Update value functions and prices:
            copy_vector(V_d, Vd0, Y_grid_size * B_grid_size);
            copy_vector(V_r, Vr0, Y_grid_size * B_grid_size);
            copy_vector(Q, Q0, Y_grid_size * B_grid_size);
            iter += 1;
        }
    }
    // If the model does not converge:
    std::cout << "Convergence not achieved after " << iter << " iterations." << std::endl;
    std::cout << "Difference between value function at default: " << diff_vd << std::endl;
    std::cout << "Difference between value function at reentry: " << diff_vr << std::endl;
    std::cout << "Difference between default probability: " << diff_q << std::endl;
    // Free memory:
    delete[] Vd0;
    delete[] Vr0;
    delete[] Q0;
    return EXIT_FAILURE;
}