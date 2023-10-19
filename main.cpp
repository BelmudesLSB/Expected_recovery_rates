/**
 * This codes performs a simple extension on the Arellano (2008) model. 
 * We propose a model in which recovery rates on defaulted bonds can vary.
 * Authors: Lucas Belmudes & Angelo Mendes 10/17/2023.
*/

#include <iostream>
#include <mex.h>
#include "economy.hpp"
#include "initialization.hpp"
#include "auxiliary.hpp"

void mexFunction(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[]){

    // Check for the proper number of arguments
    if (nrhs != 1 || nlhs != 1) {
        mexErrMsgIdAndTxt("MyToolbox:arrayProduct:nrhs", "One input required.");
    }

    // Read the input parameters from MATLAB
    const mxArray* paramsStruct = prhs[0];

    // Load the parameters from MATLAB into C++ variables:
    const int b_grid_size = static_cast<int>(mxGetScalar(mxGetField(paramsStruct, 0, "b_grid_size")));
    const double b_grid_min = static_cast<double>(mxGetScalar(mxGetField(paramsStruct, 0, "b_grid_min")));
    const double b_grid_max = static_cast<double>(mxGetScalar(mxGetField(paramsStruct, 0, "b_grid_max")));
    const int y_grid_size = static_cast<int>(mxGetScalar(mxGetField(paramsStruct, 0, "y_grid_size")));
    const double y_default = static_cast<double>(mxGetScalar(mxGetField(paramsStruct, 0, "y_default")));
    const double beta = static_cast<double>(mxGetScalar(mxGetField(paramsStruct, 0, "beta")));
    const double gamma = static_cast<double>(mxGetScalar(mxGetField(paramsStruct, 0, "gamma")));
    const double r = static_cast<double>(mxGetScalar(mxGetField(paramsStruct, 0, "r")));
    const double rho = static_cast<double>(mxGetScalar(mxGetField(paramsStruct, 0, "rho")));
    const double sigma = static_cast<double>(mxGetScalar(mxGetField(paramsStruct, 0, "sigma")));
    const double theta = static_cast<double>(mxGetScalar(mxGetField(paramsStruct, 0, "theta")));
    const double alpha = static_cast<double>(mxGetScalar(mxGetField(paramsStruct, 0, "alpha")));
    const double tol = static_cast<double>(mxGetScalar(mxGetField(paramsStruct, 0, "tol")));
    const int max_iter = static_cast<int>(mxGetScalar(mxGetField(paramsStruct, 0, "max_iter")));
    const double m = static_cast<double>(mxGetScalar(mxGetField(paramsStruct, 0, "m")));

    // Create the pointer in MATLAB to store the results:
    // Endogenous variables:
    mxArray* Q_m = mxCreateDoubleMatrix(y_grid_size * b_grid_size, 1, mxREAL);
    mxArray* V_m = mxCreateDoubleMatrix(y_grid_size * b_grid_size, 1, mxREAL);
    mxArray* V_r_m = mxCreateDoubleMatrix(y_grid_size * b_grid_size, 1, mxREAL);
    mxArray* V_d_m = mxCreateDoubleMatrix(y_grid_size * b_grid_size, 1, mxREAL);
    mxArray* B_policy_m = mxCreateDoubleMatrix(y_grid_size * b_grid_size, 1, mxREAL);
    mxArray* D_policy_m = mxCreateDoubleMatrix(y_grid_size * b_grid_size, 1, mxREAL);
    // Exogenous variables:
    mxArray* Y_grid_m = mxCreateDoubleMatrix(y_grid_size, 1, mxREAL);
    mxArray* Y_grid_default_m = mxCreateDoubleMatrix(y_grid_size, 1, mxREAL);
    mxArray* B_grid_m = mxCreateDoubleMatrix(b_grid_size, 1, mxREAL);
    mxArray* P_m = mxCreateDoubleMatrix(y_grid_size * y_grid_size, 1, mxREAL);

    
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

    mexPrintf("Initialization done.\n");
    // Solve the model:
    economy.solve_model();

    mexPrintf("Solution found .\n");
    
    // Copy the results to the matlab pointers:
    double* Q_m_aux = mxGetPr(Q_m);
    double* V_m_aux = mxGetPr(V_m);
    double* V_r_m_aux = mxGetPr(V_r_m);
    double* V_d_m_aux = mxGetPr(V_d_m);
    double* B_policy_m_aux = mxGetPr(B_policy_m);
    double* D_policy_m_aux = mxGetPr(D_policy_m);
    double* Y_grid_m_aux = mxGetPr(Y_grid_m);
    double* Y_grid_default_m_aux = mxGetPr(Y_grid_default_m);
    double* B_grid_m_aux = mxGetPr(B_grid_m);
    double* P_m_aux = mxGetPr(P_m);

    mexPrintf("Copying results to MATLAB.\n");
   
    copy_vector(economy.Q, Q_m_aux, y_grid_size*b_grid_size);
    copy_vector(economy.V, V_m_aux, y_grid_size*b_grid_size);
    copy_vector(economy.V_r, V_r_m_aux, y_grid_size*b_grid_size);
    copy_vector(economy.V_d, V_d_m_aux, y_grid_size*b_grid_size);
    copy_vector(economy.B_policy, B_policy_m_aux, y_grid_size*b_grid_size);
    copy_vector(economy.D_policy, D_policy_m_aux, y_grid_size*b_grid_size);
    copy_vector(economy.Y_grid, Y_grid_m_aux, y_grid_size);
    copy_vector(economy.Y_grid_default, Y_grid_default_m_aux, y_grid_size);
    copy_vector(economy.B_grid, B_grid_m_aux, b_grid_size);
    copy_vector(economy.P, P_m_aux, y_grid_size*y_grid_size);

    mexPrintf("Exporting results to MATLAB.\n");

    // Set the pointers to the MATLAB structure:
    const char* fieldNames[10] = {"Q", "V", "V_r", "V_d", "B_policy", "D_policy", "Y_grid", "Y_grid_default", "B_grid", "P"};
    plhs[0] = mxCreateStructMatrix(1, 1, 10, fieldNames);
    mxSetField(plhs[0], 0, "Q", Q_m);
    mxSetField(plhs[0], 0, "V", V_m);
    mxSetField(plhs[0], 0, "V_r", V_r_m);
    mxSetField(plhs[0], 0, "V_d", V_d_m);
    mxSetField(plhs[0], 0, "B_policy", B_policy_m);
    mxSetField(plhs[0], 0, "D_policy", D_policy_m);
    mxSetField(plhs[0], 0, "Y_grid", Y_grid_m);
    mxSetField(plhs[0], 0, "Y_grid_default", Y_grid_default_m);
    mxSetField(plhs[0], 0, "B_grid", B_grid_m);
    mxSetField(plhs[0], 0, "P", P_m);

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
}