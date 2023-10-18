#ifndef initialization_hpp
#define intialization_hpp

#include <iostream>

// Create bond grids:
void create_bond_grids(double* prt_bond_grid, int Nb, double Bmax, double Bmin);

// Create the income grid and the transition matrix:
void create_income_and_prob_grids(double* prt_y_grid, double* prt_p_grid,  int Ny,  double Sigma,  double Rho,  double M);

// Create the income grid for the default state:
void create_income_under_default(double* prt_y_grid_default, double* prt_y_grid,  int Ny,  double y_def);
#endif