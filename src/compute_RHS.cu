#include <cstdio>
#include <cuda_runtime.h>
#include "cuda_header.cuh"
#include "cuda_helper.cuh"
#include "../include/functionCalls.h"

__device__ int returnIndexCudaRHS(int i, int j, int NX){
    return i*NX +j;
}

__global__ void cuda_computeRHS(double *U_in,double *V_in,double *U_out, double *V_out,
        double DX,double DY,double LY, double LX, double RE,int NX,int NY){

    for(int index = blockIdx.x * blockDim.x + threadIdx.x; index < NX*NY; index += blockDim.x * gridDim.x){
        int i = index/NX;
        int j = index%NX;

        int i1  = (i + 1) % NX;
        int i2  = (i + 2) % NX;
        int i_1 = (i - 1 + NX)% NX;

        int j1  = (j + 1) % NY;
        int j2  = (j + 2) % NY;
        int j_1 = (j - 1 + NY) % NY;

        double  duudx = 0.5 * (-U_in[returnIndexCudaRHS(i_1, j, NX)]*U_in[returnIndexCudaRHS(i_1, j, NX)] + U_in[returnIndexCudaRHS(i1, j, NX)]*U_in[returnIndexCudaRHS(i1, j, NX)])/DX;
        double  dvvdy = 0.5 * (-V_in[returnIndexCudaRHS(i, j_1, NX)]*V_in[returnIndexCudaRHS(i, j_1, NX)] + V_in[returnIndexCudaRHS(i, j1, NX)]*V_in[returnIndexCudaRHS(i, j1, NX)])/DY;

        double  du2dxdx = (U_in[returnIndexCudaRHS(i1, j, NX)]  - 2*U_in[returnIndexCudaRHS(i, j, NX)] + U_in[returnIndexCudaRHS(i_1, j, NX)])/(DX*DX);
        double  du2dydy = (U_in[returnIndexCudaRHS(i, j1, NX)]  - 2*U_in[returnIndexCudaRHS(i, j, NX)] + U_in[returnIndexCudaRHS(i, j_1, NX)])/(DY*DY);
        double  dv2dxdx = (V_in[returnIndexCudaRHS(i1, j, NX)]  - 2*V_in[returnIndexCudaRHS(i, j, NX)] + V_in[returnIndexCudaRHS(i_1, j, NX)])/(DX*DX);
        double  dv2dydy = (V_in[returnIndexCudaRHS(i, j1, NX)]  - 2*V_in[returnIndexCudaRHS(i, j, NX)] + V_in[returnIndexCudaRHS(i, j_1, NX)])/(DY*DY);

        double duvdy = 0.5 * (U_in[returnIndexCudaRHS(i, j1, NX)] * 0.25 * (V_in[returnIndexCudaRHS(i_1, j2, NX)] + V_in[returnIndexCudaRHS(i, j2, NX)] + 
            V_in[returnIndexCudaRHS(i_1, j1, NX)] + V_in[returnIndexCudaRHS(i, j1, NX)]))/DY - 0.5 * (U_in[returnIndexCudaRHS(i, j_1, NX)] * 0.25 *(V_in[returnIndexCudaRHS(i, j, NX)] + V_in[returnIndexCudaRHS(i_1, j, NX)]
            + V_in[returnIndexCudaRHS(i, j_1, NX)] + V_in[returnIndexCudaRHS(i_1, j_1, NX)]))/DY;

        double duvdx = 0.5 * (V_in[returnIndexCudaRHS(i, j, NX)] * 0.25 * (U_in[returnIndexCudaRHS(i1, j, NX)] + U_in[returnIndexCudaRHS(i2, j, NX)] + 
            U_in[returnIndexCudaRHS(i1, j_1, NX)] + U_in[returnIndexCudaRHS(i2, j_1, NX)]))/DX - 0.5 * (V_in[returnIndexCudaRHS(i_1, j, NX)] * 0.25 *(U_in[returnIndexCudaRHS(i, j, NX)] + U_in[returnIndexCudaRHS(i_1, j, NX)]
            + U_in[returnIndexCudaRHS(i_1, j_1, NX)] + U_in[returnIndexCudaRHS(i, j_1, NX)]))/DX;

        U_out[returnIndexCudaRHS(i, j, NX)] = -duvdy - duudx + 1/RE * (du2dxdx + du2dydy);
        V_out[returnIndexCudaRHS(i, j, NX)] = -duvdx - dvvdy + 1/RE * (dv2dxdx + dv2dydy);
    }
}

void call_cudacomputeRHS(double *U_in, double *V_in, double *U_out, double *V_out,
        double DX,double DY,double LY, double LX, double RE,int NX, int NY) {

    CheckKernel();
    cuda_computeRHS<<<1024,32>>>(U_in,V_in, U_out,V_out,DX, DY,LY,LX,RE,NX,NY); 
    CheckKernel();  // Copy array back to host
    cudaThreadSynchronize();
}

void CheckKernel(){
    cudaError_t err = cudaGetLastError();
    if (err != cudaSuccess)
    printf("Error: %s\n", cudaGetErrorString(err));
}