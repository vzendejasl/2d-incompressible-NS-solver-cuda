#ifndef FUNCTIONCALLS_H
#define FUNCTIONCALLS_H

#include "functionCalls.h"
#include <complex>
#include <cuda_runtime.h>

using namespace std;

void computeRHS(double *U_in, double *V_in);

int returnIndex(int i, int j);
void time_advance_RK3(double *U_in,double *V_in, double *U_out, double *V_out, double *k1mod, double *k2mod);

void scalar_matrix_product(double *A, double *B, double c);
void add_matrices(double *A, double *B, double *C);
void subtract_matrices(double *A, double *B, double *C);
void double2complex(double *A, complex<double> *B);
void scalar_matrix_complex_product(complex<double> *A, complex<double> *B, double c);
void fft_shift(complex<double> *out);
int FFT2D(complex<double> *c,int nx,int ny,int dir);
int FFT(int dir,int m,double *x,double *y);
int Powerof2(int n,int *m,int *twopm);

void compute_curl(double *U_in, double *V_in, double *curlz);

void compute_projection_step(double *U_star,double *V_star, double *U_n1, double *V_n1, double *k1mod, double *k2mod);

void computeErrorRHS(double *U_RHS_EXACT,double *V_RHS_EXACT, 
    double *U_out,double *V_out,
    double *MaxErrorx, double *MaxErrory);

#endif