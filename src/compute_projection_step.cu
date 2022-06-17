#include <cstdio>
#include <cuda_runtime.h>
#include "cuda_header.cuh"
#include "cuda_helper.cuh"
#include "../include/functionCalls.h"

#include <complex>
#include <iostream>
#include <cmath>
#include <cufft.h>

#define IDX2R(i,j,N) (((i)*(N))+(j))

using namespace std;

__device__ int returnIndexCuda(int i, int j, int NX){
    return i*NX +j;
}
__global__ void cuda_compute_div(double *d_U_star, double *d_V_star, cufftComplex *d_div, 
        double DX, double DY, int NX, int NY) {

    for(int index = blockIdx.x * blockDim.x + threadIdx.x; index < NX*NY; index += blockDim.x * gridDim.x){
        int i = index/NX;
        int j = index%NX;
        int i1  = (i + 1) % NX;
        int j1  = (j + 1) % NY;

        d_div[returnIndexCuda(i,j,NX)].x = ((d_U_star[returnIndexCuda(i1,j,NX)]-d_U_star[returnIndexCuda(i,j,NX)])/DX
                                       + (d_V_star[returnIndexCuda(i,j1,NX)]-d_V_star[returnIndexCuda(i,j,NX)])/DY);
    }
}

__global__ void cuda_fftshift_2D(cufftComplex *data, int NX, int NY) {
    int i = threadIdx.y + blockDim.y * blockIdx.y;
    int j = threadIdx.x + blockDim.x * blockIdx.x;

    if (i < NX && j < NY) {
        double a = pow(-1.0, (i+j)&1);
        data[IDX2R(i,j,NY)].x *= a;
        data[IDX2R(i,j,NY)].y *= a;
    }
}

__global__ void cuda_compute_Phat(cufftComplex *d_f_hat, cufftComplex *d_P_hat, double *k1mod, double *k2mod, int NX, int NY) {

    for(int index = blockIdx.x * blockDim.x + threadIdx.x; index < NX*NY; index += blockDim.x * gridDim.x) {
        int i = index/NX;
        int j = index%NX;

        d_P_hat[returnIndexCuda(i,j,NX)].x = - d_f_hat[returnIndexCuda(i,j,NX)].x / (k1mod[returnIndexCuda(i,j,NX)]+k2mod[returnIndexCuda(i,j,NX)]);
        d_P_hat[returnIndexCuda(i,j,NX)].y = - d_f_hat[returnIndexCuda(i,j,NX)].y / (k1mod[returnIndexCuda(i,j,NX)]+k2mod[returnIndexCuda(i,j,NX)]);
    }
}

__global__ void cuda_set_Phat_middle(cufftComplex *d_P_hat, int NX, int NY) {
    d_P_hat[returnIndexCuda(NX/2,NY/2,NX)].x = 0;
    d_P_hat[returnIndexCuda(NX/2,NY/2,NX)].y = 0;
}

__global__ void cuda_P_complex2double(cufftComplex *d_P_hat, double *d_P, int NX, int NY) {
    for(int index = blockIdx.x * blockDim.x + threadIdx.x; index < NX*NY; index += blockDim.x * gridDim.x) {
        int i = index/NX;
        int j = index%NX;

        d_P[returnIndexCuda(i,j,NX)] = d_P_hat[returnIndexCuda(i,j,NX)].x;

    }
}

__global__ void cuda_update_UV(double *d_P, double *d_U_star, double *d_V_star, double *d_U_n1, double *d_V_n1, double DX, double DY, int NX, int NY) {

    for(int index = blockIdx.x * blockDim.x + threadIdx.x; index < NX*NY; index += blockDim.x * gridDim.x){
        int i = index/NX;
        int j = index%NX;
        int i1  = (i + 1) % NX;
        int j1  = (j + 1) % NY;

        d_U_n1[returnIndexCuda(i1,j,NX)] = -(d_P[returnIndexCuda(i1,j,NX)]-d_P[returnIndexCuda(i,j,NX)])/DX + d_U_star[returnIndexCuda(i1,j,NX)];
        d_V_n1[returnIndexCuda(i,j1,NX)] = -(d_P[returnIndexCuda(i,j1,NX)]-d_P[returnIndexCuda(i,j,NX)])/DY + d_V_star[returnIndexCuda(i,j1,NX)];
    }
}

void cuda_compute_projection_step(double *d_U_star, double *d_V_star, double *d_U_n1, double *d_V_n1, 
    double *d_k1mod, double *d_k2mod, double DX, double DY, int NX, int NY) {

    cufftComplex *d_div;
    cudaMalloc(&d_div, NX*NY*sizeof(cufftComplex));

    cuda_compute_div<<<1024,32>>>(d_U_star, d_V_star, d_div, DX, DY, NX, NY);

    cufftComplex *d_f_hat;
    cufftComplex *d_P_hat;
    cudaMalloc(&d_f_hat, NX*NY*sizeof(cufftComplex));
    cudaMalloc(&d_P_hat, NX*NY*sizeof(cufftComplex));

    cufftHandle plan;
    cufftPlan2d(&plan, NX, NY, CUFFT_C2C);
    cufftExecC2C(plan, d_div, d_f_hat, CUFFT_FORWARD);

    CheckKernel();
    cuda_fftshift_2D<<<1024,32>>>(d_f_hat, NX, NY);
    CheckKernel();

    CheckKernel();
    cuda_compute_Phat<<<1024,32>>>(d_f_hat, d_P_hat, d_k1mod, d_k2mod, NX, NY);
    CheckKernel();

    CheckKernel();
    cuda_set_Phat_middle<<<1,32>>>(d_P_hat,NX,NY);
    CheckKernel();

    CheckKernel();
    cuda_fftshift_2D<<<1024,32>>>(d_P_hat,NX,NY);
    CheckKernel();

    cufftExecC2C(plan, d_P_hat, d_P_hat, CUFFT_INVERSE);

    double *d_P;
    cudaMalloc(&d_P, NX*NY*sizeof(double));

    CheckKernel();
    cuda_P_complex2double<<<1024,32>>>(d_P_hat,d_P,NX,NY);
    CheckKernel();

    CheckKernel();
    cuda_update_UV<<<1024,32>>>(d_P,d_U_star,d_V_star,d_U_n1,d_V_n1,DX,DY,NX,NY);
    CheckKernel();

    cudaFree(d_div);
    cudaFree(d_f_hat);
    cudaFree(d_P_hat);
    cufftDestroy(plan);
    cudaFree(d_P);
    
}
