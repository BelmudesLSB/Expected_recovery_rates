#include <iostream>
#include "economy.hpp"
#include "Initialization.hpp"
#include "auxiliary.hpp"
#include <cmath>

// Create bond grids:
void create_bond_grids(double* prt_bond_grid, int Nb, double Bmax, double Bmin){
    double bstep = (Bmax - Bmin)/(Nb - 1);
    for(int i = 0; i < Nb; i++){
        prt_bond_grid[i] = Bmin + i*bstep;
    }
}

// Create the income grid and the transition matrix:
void create_income_and_prob_grids(double* prt_y_grid, double* prt_p_grid,  int Ny,  double Sigma,  double Rho,  double M){
    double sigma_y = sqrt(pow(Sigma,2)/(1-pow(Rho,2)));
    double omega = (2*M*sigma_y)/(Ny-1);
    for (int i=0; i<Ny; i++){ 
        prt_y_grid[i] = (-M*sigma_y)  + omega * i;
    }   
    for (int i=0; i<Ny; i++){
        for (int j=0; j<Ny; j++){
            if (j==0 || j==Ny-1){
                if (j==0){
                    prt_p_grid[i*Ny+j] = normalCDF((prt_y_grid[0]-Rho*prt_y_grid[i]+omega/2)/Sigma);
                }
                else {
                    prt_p_grid[i*Ny+j] = 1-normalCDF((prt_y_grid[Ny-1]-Rho*prt_y_grid[i]-omega/2)/Sigma);
                }
            } else {
                prt_p_grid[i*Ny+j] = normalCDF((prt_y_grid[j]-Rho*prt_y_grid[i]+omega/2)/Sigma)-normalCDF((prt_y_grid[j]-Rho*prt_y_grid[i]-omega/2)/Sigma);
            }
        }
    }
    for (int i=0; i<Ny; i++){
        prt_y_grid[i] = exp(prt_y_grid[i]);
    }
}

// Create the income grid for the default state:
void create_income_under_default(double* prt_y_grid_default, double* prt_y_grid,  int Ny,  double y_def){
    for (int i=0; i<Ny; i++){
        if (prt_y_grid[i]>y_def){
            prt_y_grid_default[i] = y_def;
        } else {
            prt_y_grid_default[i] = prt_y_grid[i];
        }
    }
}