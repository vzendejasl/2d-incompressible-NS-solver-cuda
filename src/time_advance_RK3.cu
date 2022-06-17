#include "../include/functionCalls.h" 
#include "../include/variableCalls.h"

#include <cuda_runtime.h>
#include "cuda_header.cuh"
#include <cublas_v2.h>
#include "cuda_helper.cuh"
#include <cufft.h>

void cuda_time_advance_RK3(double *d_U_in, double *d_V_in, double *d_U_out, double *d_V_out, 
    double *d_k1mod, double *d_k2mod, double DT, double DX, double DY, int NX, int NY){

    cublasHandle_t handle;
    cublasCreate(&handle);

    double one = 1;
    double zero = 0;
    double two = 2;
    double four = 4;
    double minusone = -1;
    double onehalf = 1/2;
    double onesixth = 1/6;

    // Allocate device memory
    double *d_U_k1;
    double *d_V_k1;
    double *d_U_k2;
    double *d_V_k2;
    double *d_U_k3;
    double *d_V_k3;
    double *d_U_temp;
    double *d_V_temp;

    cudaMalloc(&d_U_k1,   NX*NY*sizeof(double));
    cudaMalloc(&d_V_k1,   NX*NY*sizeof(double));
    cudaMalloc(&d_U_k2,   NX*NY*sizeof(double));
    cudaMalloc(&d_V_k2,   NX*NY*sizeof(double));
    cudaMalloc(&d_U_k3,   NX*NY*sizeof(double));
    cudaMalloc(&d_V_k3,   NX*NY*sizeof(double));
    cudaMalloc(&d_U_temp, NX*NY*sizeof(double));
    cudaMalloc(&d_V_temp, NX*NY*sizeof(double));


    ///// First RK3 step /////
    call_cudacomputeRHS(d_U_in,d_V_in,d_U_temp,d_V_temp,
        DX,DY,LY,LX,RE,NX,NY);

    // U_k1 = U_temp * dt
    cublasDgeam(handle, CUBLAS_OP_N, CUBLAS_OP_N,
        NX, NY,
        &DT, d_U_temp, NX, 
        &zero,nullptr, NX, 
        d_U_k1, NX);

    // V_k1 = V_temp * dt
    cublasDgeam(handle, CUBLAS_OP_N, CUBLAS_OP_N,
        NX, NY, 
        &DT, d_V_temp, NX, 
        &zero, nullptr,NX,
        d_V_k1, NX);
    ///////////////////////////


    ///// Second RK3 step /////
    // U_k2 = U_in + U_k1/2
    cublasDgeam(handle, CUBLAS_OP_N, CUBLAS_OP_N,
        NX, NY, 
        &one, d_U_in, NX,
        &onehalf, d_U_k1, NX,
        d_U_k2, NX);

    // V_k2 = V_in + V_k1/2
    cublasDgeam(handle, CUBLAS_OP_N, CUBLAS_OP_N,
        NX, NY, 
        &one, d_V_in, NX, 
        &onehalf, d_V_k1,NX,
        d_V_k2, NX);

    cuda_compute_projection_step(d_U_k2,d_V_k2,d_U_k2,d_V_k2,d_k1mod,d_k2mod,DX,DY,NX,NY);
    call_cudacomputeRHS(d_U_k2,d_V_k2,d_U_temp,d_V_temp,DX,DY,LY,LX,RE,NX,NY);

    // U_k2 = U_temp * dt
    cublasDgeam(handle, CUBLAS_OP_N, CUBLAS_OP_N,
        NX, NY, 
        &DT, d_U_temp, NX, 
        &zero, nullptr, NX,
        d_U_k2, NX);

    // V_k2 = V_temp * dt
    cublasDgeam(handle, CUBLAS_OP_N, CUBLAS_OP_N,
        NX, NY,
        &DT, d_V_temp, NX, 
        &zero, nullptr,NX,
        d_V_k2, NX);
    ///////////////////////////


    ///// Third RK3 step /////
    // U_k3 = U_in - U_k1 + 2*U_k2
    cublasDgeam(handle, CUBLAS_OP_N, CUBLAS_OP_N,
        NX, NY, 
        &one, d_U_in, NX, 
        &minusone,d_U_k1,NX, 
        d_U_temp, NX);
    cublasDgeam(handle, CUBLAS_OP_N, CUBLAS_OP_N,
        NX, NY, 
        &one, d_U_temp, NX, 
        &two,d_U_k2,NX,
        d_U_k3, NX);

    // V_k3 = V_in - V_k1 + 2*V_k2
    cublasDgeam(handle, CUBLAS_OP_N, CUBLAS_OP_N,
        NX, NY, 
        &one, d_V_in, NX, 
        &minusone, d_V_k1, NX,
        d_V_temp, NX);
    cublasDgeam(handle, CUBLAS_OP_N, CUBLAS_OP_N,
        NX, NY, 
        &one, d_V_temp, NX, 
        &two,d_V_k2,NX,
        d_V_k3, NX);

    cuda_compute_projection_step(d_U_k3,d_V_k3,d_U_k3,d_V_k3,d_k1mod,d_k2mod,DX,DY,NX,NY);
    call_cudacomputeRHS(d_U_k3,d_V_k3,d_U_temp,d_V_temp,DX,DY,LY,LX,RE,NX,NY);

    // U_k3 = U_temp * dt
    cublasDgeam(handle, CUBLAS_OP_N, CUBLAS_OP_N,
        NX, NY, 
        &DT, d_U_temp, NX, 
        &zero,nullptr,NX,
        d_U_k3, NX);

    // V_k3 = V_temp * dt
    cublasDgeam(handle, CUBLAS_OP_N, CUBLAS_OP_N,
        NX, NY,
        &DT, d_V_temp, NX, 
        &zero,nullptr,NX,
        d_V_k3, NX);
    ///////////////////////////


    ///// U_out, V_out /////
    // U_out = U_in + (U_k1 + 4*U_k2 + U_k3)/6
    cublasDgeam(handle, CUBLAS_OP_N, CUBLAS_OP_N,
        NX, NY, 
        &one, d_U_k1, NX, 
        &four, d_U_k2,NX,
        d_U_temp, NX);
    cublasDgeam(handle, CUBLAS_OP_N, CUBLAS_OP_N,
        NX, NY, 
        &one, d_U_temp, NX, 
        &one, d_U_k3,NX,
        d_U_temp, NX);
    cublasDgeam(handle, CUBLAS_OP_N, CUBLAS_OP_N,
        NX, NY, 
        &onesixth, d_U_temp, NX, 
        &one, d_U_in, NX,
        d_U_out, NX);

    // V_out = V_in + (V_k1 + 4*V_k2 + V_k3)/6
    cublasDgeam(handle, CUBLAS_OP_N, CUBLAS_OP_N,
        NX, NY, 
        &one, d_V_k1, NX, 
        &four, d_V_k2,NX,
        d_V_temp, NX);
    cublasDgeam(handle, CUBLAS_OP_N, CUBLAS_OP_N,
        NX, NY,
        &one, d_V_temp, NX, 
        &one, d_V_k3, NX,
        d_V_temp, NX);
    cublasDgeam(handle, CUBLAS_OP_N, CUBLAS_OP_N,
        NX, NY, 
        &onesixth, d_V_temp, NX, 
        &one, d_V_in, NX,
        d_V_out, NX);

    cuda_compute_projection_step(d_U_out,d_V_out,d_U_out,d_V_out,d_k1mod,d_k2mod,DX,DY,NX,NY);
    ///////////////////////////

    // Free memory
    cudaFree(d_U_k1);
    cudaFree(d_V_k1);
    cudaFree(d_U_k2);
    cudaFree(d_V_k2);
    cudaFree(d_U_k3);
    cudaFree(d_V_k3);
    cudaFree(d_U_temp);
    cudaFree(d_V_temp);

    cublasDestroy(handle);
}
