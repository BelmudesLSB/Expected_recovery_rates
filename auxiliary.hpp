#ifndef auxiliary_hpp
#define auxiliary_hpp

#include <iostream>

// Normal cumulative distribution function:
double normalCDF(double x);

// Display vector:
void display_vector(double* prt_vector, int size);

// Display matrix of doubles size (n,m):
void display_matrix(double* prt_matrix, int n, int m);

// Display matrix of integers size (n,m):
void display_matrix_int(int* prt_matrix, int n, int m);

// Utility function:
double utility(double c,  double gamma, double c_lb);

// Copy vector:
void copy_vector(double* prt_vector, double* prt_vector_copy, int size);

// Copy vector:
void copy_vector(int* prt_vector, double* prt_vector_copy, int size);
#endif