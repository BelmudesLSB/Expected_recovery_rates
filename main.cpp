/**
 * We propose an extension to Arellano (2008) in which recovery rates on defaulted bonds can vary depending on the bond.
 * Sovereign can issue 2 types of bond, low and high recovery rate.
 * Authors: Lucas Belmudes & Angelo Mendes 10/23/2023.
*/

#include <iostream>
#include <mex.h>
#include <omp.h>
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
    const double alpha_low = static_cast<double>(mxGetScalar(mxGetField(paramsStruct, 0, "alpha_low")));
    const double alpha_high = static_cast<double>(mxGetScalar(mxGetField(paramsStruct, 0, "alpha_high")));
    const double tol = static_cast<double>(mxGetScalar(mxGetField(paramsStruct, 0, "tol")));
    const int max_iter = static_cast<int>(mxGetScalar(mxGetField(paramsStruct, 0, "max_iter")));
    const double m = static_cast<double>(mxGetScalar(mxGetField(paramsStruct, 0, "m")));

    // Create the pointer in MATLAB to store the results:
    // Endogenous variables:
    mxArray* Q_m_low = mxCreateDoubleMatrix(y_grid_size * b_grid_size * b_grid_size, 1, mxREAL);
    mxArray* Q_m_high = mxCreateDoubleMatrix(y_grid_size * b_grid_size * b_grid_size, 1, mxREAL);
    mxArray* V_m = mxCreateDoubleMatrix(y_grid_size * b_grid_size * b_grid_size, 1, mxREAL);
    mxArray* V_r_m = mxCreateDoubleMatrix(y_grid_size * b_grid_size * b_grid_size, 1, mxREAL);
    mxArray* V_d_m = mxCreateDoubleMatrix(y_grid_size * b_grid_size * b_grid_size, 1, mxREAL);
    mxArray* B_policy_m_low = mxCreateDoubleMatrix(y_grid_size * b_grid_size * b_grid_size, 1, mxREAL);
    mxArray* B_policy_m_high = mxCreateDoubleMatrix(y_grid_size * b_grid_size * b_grid_size, 1, mxREAL);
    mxArray* D_policy_m = mxCreateDoubleMatrix(y_grid_size * b_grid_size * b_grid_size, 1, mxREAL);
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
    double* v = new double[y_grid_size*b_grid_size*b_grid_size];
    double* v_r = new double[y_grid_size*b_grid_size*b_grid_size];
    double* v_d = new double[y_grid_size*b_grid_size*b_grid_size];
    double* q_low = new double[y_grid_size*b_grid_size*b_grid_size];
    double* q_high = new double[y_grid_size*b_grid_size*b_grid_size];
    int* b_policy_low = new int[y_grid_size*b_grid_size*b_grid_size];
    int* b_policy_high = new int[y_grid_size*b_grid_size*b_grid_size];
    double* d_policy = new double[y_grid_size*b_grid_size*b_grid_size];
  
    // Create an instance of the Economy class:
    Economy economy(b_grid_size, b_grid_min, b_grid_max, y_grid_size, y_default, beta, gamma, r, rho, sigma, theta, alpha_low, alpha_high, tol, max_iter, m, y_grid, y_grid_default, b_grid, p, v, v_r, v_d, q_low, q_high, b_policy_low, b_policy_high, d_policy);
    
    mexPrintf("Initialization done.\n");
    // Solve the model:
    //economy.initialize_economy();
    //economy.guess_vd_vr_q();
    //economy.update_v_and_default_policy();
    //economy.update_price(); 
    //economy.update_vd();
    //economy.update_vr_and_bond_policy();
    mexPrintf("Solving the model.\n"); 
    economy.solve_model();
    
    mexPrintf("Copying results to MATLAB.\n");
    copy_vector(q_low, mxGetPr(Q_m_low), y_grid_size*b_grid_size*b_grid_size);
    copy_vector(q_high, mxGetPr(Q_m_high), y_grid_size*b_grid_size*b_grid_size);
    copy_vector(v, mxGetPr(V_m), y_grid_size*b_grid_size*b_grid_size);
    copy_vector(v_r, mxGetPr(V_r_m), y_grid_size*b_grid_size*b_grid_size);
    copy_vector(v_d, mxGetPr(V_d_m), y_grid_size*b_grid_size*b_grid_size);
    copy_vector(b_policy_low, mxGetPr(B_policy_m_low), y_grid_size*b_grid_size*b_grid_size);
    copy_vector(b_policy_high, mxGetPr(B_policy_m_high), y_grid_size*b_grid_size*b_grid_size);
    copy_vector(d_policy, mxGetPr(D_policy_m), y_grid_size*b_grid_size*b_grid_size);
    copy_vector(y_grid, mxGetPr(Y_grid_m), y_grid_size);
    copy_vector(y_grid_default, mxGetPr(Y_grid_default_m), y_grid_size);
    copy_vector(b_grid, mxGetPr(B_grid_m), b_grid_size);
    copy_vector(p, mxGetPr(P_m), y_grid_size*y_grid_size);

    mexPrintf("Exporting results to MATLAB.\n");

    // Set the pointers to the MATLAB structure:
    const char* fieldNames[12] = {"Q_low", "Q_high", "V", "V_r", "V_d", "B_policy_low", "B_policy_high", "D_policy", "Y_grid", "Y_grid_default", "B_grid", "P"};
    plhs[0] = mxCreateStructMatrix(1, 1, 12, fieldNames);
    mxSetField(plhs[0], 0, "Q_low", Q_m_low);
    mxSetField(plhs[0], 0, "Q_high", Q_m_high);
    mxSetField(plhs[0], 0, "V", V_m);
    mxSetField(plhs[0], 0, "V_r", V_r_m);
    mxSetField(plhs[0], 0, "V_d", V_d_m);
    mxSetField(plhs[0], 0, "B_policy_low", B_policy_m_low);
    mxSetField(plhs[0], 0, "B_policy_high", B_policy_m_high);
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
    delete[] q_low;
    delete[] q_high;
    delete[] b_policy_low;
    delete[] b_policy_high;
    delete[] d_policy;
}