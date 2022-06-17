#ifndef CUDAHELPER_CUH
#define CUDAHELPER_CUH
 
#include "cuda_header.cuh"
#include <cuda_runtime.h>
 
#include "../include/variableCalls.h"

void call_cudacomputeRHS(double *U_in,double *V_in,double *U_out, double *V_out,
    double DX,double DY,double LY, double LX, double RE,int NX,int NY);
void cuda_time_advance_RK3(double *d_U_in, double *d_V_in, double *d_U_out, double *d_V_out, 
    double *d_k1mod, double *d_k2mod, double DT, double DX, double DY, int NX, int NY);
void cuda_compute_projection_step(double *d_U_star, double *d_V_star, double *d_U_n1, double *d_V_n1, 
    double *d_k1mod, double *d_k2mod, double DX, double DY, int NX, int NY);
void CheckKernel();

#endif
