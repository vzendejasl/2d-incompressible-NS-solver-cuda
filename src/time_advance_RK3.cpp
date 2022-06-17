#include "../include/functionCalls.h" 
#include "../include/variableCalls.h"
#include "cuda_helper.cuh"
#include <cuda_runtime.h>

void time_advance_RK3(double *U_in,double *V_in, double *U_out, double *V_out, double *k1mod, double *k2mod) {

    double *U_k1   = new double [NX * NY];
    double *V_k1   = new double [NY * NY];
    double *U_k2   = new double [NX * NY];
    double *V_k2   = new double [NY * NY];
    double *U_k3   = new double [NX * NY];
    double *V_k3   = new double [NY * NY];
    double *U_temp = new double [NY * NY];
    double *V_temp = new double [NX * NY];

    for (int i = 0; i < NX; ++i) {
        for (int j = 0; j < NY; ++j) {
            U_temp[returnIndex(i,j)] = U_in[returnIndex(i,j)];
        }
    }
    for (int i = 0; i < NX; ++i) {
        for (int j = 0; j < NY; ++j) {
            V_temp[returnIndex(i,j)] = V_in[returnIndex(i,j)];
        }
    }

    /// First RK3 step ///
    computeRHS(U_temp,V_temp);
    for (int i = 0; i < NX; ++i) {
        for (int j = 0; j < NY; ++j) {
            U_k1[returnIndex(i,j)] = U_temp[returnIndex(i,j)] * DT;
            V_k1[returnIndex(i,j)] = V_temp[returnIndex(i,j)] * DT;
        }
    }
    //////////////////////

    /// Second RK3 step ///
    // U_k2 = U_in + U_k1/2
    // V_k2 = V_in + V_k1/2
    for (int i = 0; i < NX; ++i) {
        for (int j = 0; j < NY; ++j) {
            U_k2[returnIndex(i,j)] = U_in[returnIndex(i,j)] + U_k1[returnIndex(i,j)]/2;
            V_k2[returnIndex(i,j)] = V_in[returnIndex(i,j)] + V_k1[returnIndex(i,j)]/2;
        }
    }
    compute_projection_step(U_k2,V_k2,U_k2,V_k2,k1mod,k2mod);
    computeRHS(U_k2,V_k2);
    for (int i = 0; i < NX; ++i) {
        for (int j = 0; j < NY; ++j) {
            U_k2[returnIndex(i,j)] = U_k2[returnIndex(i,j)] * DT;
            V_k2[returnIndex(i,j)] = V_k2[returnIndex(i,j)] * DT;
        }
    }
    ///////////////////////

    /// Third RK3 step ///
    // U_k3 = U_in - U_k1 + 2*U_k2
    // V_k3 = V_in - V_k1 + 2*V_k2
    for (int i = 0; i < NX; ++i) {
        for (int j = 0; j < NY; ++j) {
            U_k3[returnIndex(i,j)] = U_in[returnIndex(i,j)] - U_k1[returnIndex(i,j)] + 2*U_k2[returnIndex(i,j)];
            V_k3[returnIndex(i,j)] = V_in[returnIndex(i,j)] - V_k1[returnIndex(i,j)] + 2*V_k2[returnIndex(i,j)];
        }
    }
    compute_projection_step(U_k3,V_k3,U_k3,V_k3,k1mod,k2mod);
    computeRHS(U_k3,V_k3);
    for (int i = 0; i < NX; ++i) {
        for (int j = 0; j < NY; ++j) {
            U_k3[returnIndex(i,j)] = U_k3[returnIndex(i,j)] * DT;
            V_k3[returnIndex(i,j)] = V_k3[returnIndex(i,j)] * DT;
        }
    }
    //////////////////////
    
    /// Compute U,V out ///
    // U_out = U_in + (U_k1 + 4*U_k2 + U_k3)/6
    // V_out = V_in + (V_k1 + 4*V_k2 + V_k3)/6
    for (int i = 0; i < NX; ++i) {
        for (int j = 0; j < NY; ++j) {
            U_out[returnIndex(i,j)] = U_in[returnIndex(i,j)] + (U_k1[returnIndex(i,j)] + 4*U_k2[returnIndex(i,j)] + U_k3[returnIndex(i,j)])/6;
            V_out[returnIndex(i,j)] = V_in[returnIndex(i,j)] + (V_k1[returnIndex(i,j)] + 4*V_k2[returnIndex(i,j)] + V_k3[returnIndex(i,j)])/6;
        }
    }
    compute_projection_step(U_out,V_out,U_out,V_out,k1mod,k2mod);
    ///////////////////////

    // Free memory
    delete [] U_temp;
    delete [] V_temp;
    delete [] U_k1;
    delete [] V_k1;
    delete [] U_k2;
    delete [] V_k2;
    delete [] U_k3;
    delete [] V_k3;
}
