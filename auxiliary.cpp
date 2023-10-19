#include <iostream>
#include "auxiliary.hpp"
#include <cmath>
#include <mex.h>

// Normal cumulative distribution function:
double normalCDF(double x){
    return std::erfc(-x / std::sqrt(2)) / 2;
}

// Display vector:
void display_vector(double* prt_vector, int size){
    for (int i=0; i<size; i++){
        std::cout << prt_vector[i] << std::endl;
    }
}

// Display matrix of doubles of size (n,m):
void display_matrix(double* prt_matrix, int n, int m){
    for (int i=0; i<n; i++){
        for (int j=0; j<m; j++){
            std::cout << prt_matrix[i*m+j] << " ";
        }
        std::cout << std::endl;
    }
}

// Display matrix of integers of size (n,m):
void display_matrix_int(int* prt_matrix, int n, int m){
    for (int i=0; i<n; i++){
        for (int j=0; j<m; j++){
            std::cout << prt_matrix[i*m+j] << " ";
        }
        std::cout << std::endl;
    }
}

// Utility function:
double utility(double c,  double gamma, double c_lb){
    if (c>=c_lb){
        return pow(c,1-gamma)/(1-gamma);
    } else {
        std::cout << "!!!!!!!!!!!!!!!!!!!!!!!!!!!" << std::endl;
        std::cout << "!! Negative consumption  !!" << std::endl;
        std::cout << "!!!!!!!!!!!!!!!!!!!!!!!!!!!" << std::endl;
        return -1000000;
    }
}

// Copy vector:
void copy_vector(double* prt_vector, double* prt_vector_copy, int size){
    for (int i=0; i<size; i++){
        prt_vector_copy[i] = prt_vector[i];
    }
}

// Copy vector:
void copy_vector(int* prt_vector, double* prt_vector_copy, int size){
    for (int i=0; i<size; i++){
        prt_vector_copy[i] = static_cast<double>(prt_vector[i]);
    }
}

