#include <iostream>
#include <fstream>
#include <stdio.h>
#include <cmath>
#include <chrono>
#include "include/functionCalls.h"
#include "include/variableCalls.h"
#include "src/cuda_helper.cuh"

#include <cuda_runtime.h>

// Simulation parameters
int NX = 512; // Grid points in x
int NY = NX; // Grid points in y
double PI = 3.14159265358;
double LX = 2*PI; // Domain Length in x
double LY = LX; // Domain length in y
double RE = 5000.0; // Reynolds number
double DX = LX/NX; // Grid spacing in x
double DY = LY/NY; // Grid spacing in y
double DT = 0.001; // Timestep
int num_iter = 10; // Total number of iterations

int TRUE = 1;
int FALSE = 0;

// Parameters for setting Taylor Green Vortex IC
int A = 1;
int a = 1;
int B = -1;
int b = 1;

using namespace std;
using namespace std::chrono;

int main() {

    // Allocate memory for matrices
    double *x_umom = new double [NX];
    double *y_vmom = new double [NY];
    double *y_umom = new double [NY];
    double *x_vmom = new double [NX];

    double *X_umom = new double [NX * NY];
    double *Y_umom = new double [NX * NY];
    double *X_vmom = new double [NX * NY];
    double *Y_vmom = new double [NX * NY];
    double *U_umom = new double [NX * NY];
    double *V_vmom = new double [NX * NY];
    double *U_in   = new double [NX * NY];
    double *V_in   = new double [NX * NY];
    double *U_out  = new double [NX * NY];
    double *V_out  = new double [NX * NY];
    double *k1     = new double [NX * NY];
    double *k2     = new double [NX * NY];
    double *k1_mod = new double [NX * NY];
    double *k2_mod = new double [NX * NY];
    double *CURLZ  = new double [NX * NY];

    // Allocate device memory
    double *d_U_in;
    double *d_V_in;
    double *d_U_out;
    double *d_V_out;
    double *d_k1_mod;
    double *d_k2_mod;
        
    cudaMalloc(&d_U_in,  NX *NY * sizeof(double));
    cudaMalloc(&d_V_in,  NX *NY * sizeof(double));
    cudaMalloc(&d_U_out, NX *NY * sizeof(double));
    cudaMalloc(&d_V_out, NX *NY * sizeof(double));
    cudaMalloc(&d_k1_mod, NX *NY * sizeof(double));
    cudaMalloc(&d_k2_mod, NX *NY * sizeof(double));
    
    // Initialize to zero
    cudaMemset(d_U_in, 0, NX *NY * sizeof(double));
    cudaMemset(d_V_in, 0, NX *NY * sizeof(double));
    cudaMemset(d_U_out,0, NX *NY * sizeof(double));
    cudaMemset(d_V_out,0, NX *NY * sizeof(double));
    cudaMemset(d_k1_mod,0, NX *NY * sizeof(double));
    cudaMemset(d_k1_mod,0, NX *NY * sizeof(double));

    // Set up grid
    for(int j = 0; j < NY; j++){
        y_umom[j] = j*DY + DY/2;
        x_vmom[j] = j*DY + DY/2;
    }
    for(int i = 0; i < NX; i++){
        x_umom[i] = i*DX;
        y_vmom[i] = i*DX;
    }

    // "meshgrid"
    for(int i = 0; i < NX; i++){
        for(int j = 0; j < NY; j++){
           Y_umom[returnIndex(i,j)]= y_umom[j];
           X_umom[returnIndex(i,j)]= x_umom[i];
           X_vmom[returnIndex(i,j)]= x_vmom[i];
           Y_vmom[returnIndex(i,j)]= y_vmom[j];
        }
    }

    // Compute modified wave numbers for central finite-difference scheme
    for(int i = 0; i < NX; i++){
        for(int j = 0; j < NY; j++){
            k1[returnIndex(i,j)] = - NX/2 + i;
            k2[returnIndex(i,j)] = - NX/2 + j;
            k1_mod[returnIndex(i,j)] = 2 * (1 - cos(k1[returnIndex(i,j)] * DX))/(DX * DX);
            k2_mod[returnIndex(i,j)] = 2 * (1 - cos(k2[returnIndex(i,j)] * DY))/(DY * DY);
         }
     }


///// Taylor Green Vortex case /////
#if 1

    // Set up U(x,y), V(x,y) IC for Taylor Green Vortex
    for(int i = 0; i < NX; i++){
        for(int j = 0; j < NY; j++){
            U_in[returnIndex(i,j)] = A*sin(a*X_umom[returnIndex(i,j)])*cos(b*Y_umom[returnIndex(i,j)]);       
            V_in[returnIndex(i,j)] = B*cos(a*X_vmom[returnIndex(i,j)])*sin(b*Y_vmom[returnIndex(i,j)]);       
        }
    }
    
    cudaMemcpy(d_U_in, U_in, NX * NY * sizeof(double), 
        cudaMemcpyHostToDevice);
    cudaMemcpy(d_V_in, V_in, NX * NY * sizeof(double), 
        cudaMemcpyHostToDevice);
    cudaMemcpy(d_k1_mod, k1_mod, NX * NY * sizeof(double), 
        cudaMemcpyHostToDevice);
    cudaMemcpy(d_k2_mod, k2_mod, NX * NY * sizeof(double), 
        cudaMemcpyHostToDevice);

#endif


///// RUN CPU VERSION /////
#if 1
    cout << "Running CPU Code...\n";
    auto startCPU = high_resolution_clock::now();

    /// Time integration loop ///
    for (int t = 0; t < num_iter; t++) {
        time_advance_RK3(U_in,V_in,U_in,V_in,k1_mod,k2_mod);
        cout << "Iteration # " << t << "\n";
    }
    /////////////////////////////

    auto stopCPU = high_resolution_clock::now();
    auto durationCPU = duration_cast<milliseconds>(stopCPU - startCPU);
 
    cout << "Time taken by CPU: "
         << durationCPU.count() << " milliseconds" << endl;
    cout<< "---------------------------------\n";

    // Compute curl for final solution and write to CPU_output.txt
    compute_curl(U_in,V_in,CURLZ);

    ofstream CPU_out;
    CPU_out.open("CPU_output.txt");

    for (int i = 0; i < NX; ++i){
        for (int j = 0; j < NY; ++j){
            CPU_out << CURLZ[returnIndex(i,j)]<< " ";
        }
        CPU_out << "\n";
    }
    CPU_out.close();

#endif


///// RUN GPU VERSION /////
#if 1
    cout << "Running GPU Code...\n";
    auto startGPU = high_resolution_clock::now();
 
    /// Time integration loop ///
    for (int t = 0; t < num_iter; t++){
        cuda_time_advance_RK3(d_U_in, d_V_in, d_U_out, d_V_out, d_k1_mod, d_k2_mod, 
        DT, DX, DY, NX, NY);
        cout << "Iteration # " << t << "\n";
    }
    /////////////////////////////

    auto stopGPU = high_resolution_clock::now();
    auto durationGPU = duration_cast<milliseconds>(stopGPU - startGPU);
 
    cout << "Time taken by GPU: "
         << durationGPU.count() << " milliseconds" << endl;

    cudaMemcpy(U_in,d_U_out, NX * NY * sizeof(double), 
        cudaMemcpyDeviceToHost);
    cudaMemcpy(V_in, d_V_out, NX * NY * sizeof(double), 
        cudaMemcpyDeviceToHost);
    
    // Compute curl for final solution and write to GPU_output.txt
    compute_curl(U_in,V_in,CURLZ);

    ofstream GPU_out;
    GPU_out.open("GPU_output.txt");

    for (int i = 0; i < NX; ++i){
        for (int j = 0; j < NY; ++j){
            GPU_out << CURLZ[returnIndex(i,j)]<< " ";
        }
        GPU_out << "\n";
    }
    GPU_out.close();

#endif

// compute RHS Test
#if 0
    // Allocate host memory
    double *MaxErrorx = (double *) malloc(sizeof (double));
    double *MaxErrory = (double *) malloc(sizeof (double));
    *MaxErrorx = 0.0;
    *MaxErrory = 0.0;
    double *U_RHS_EXACT = new double [NX*NY];
    double *V_RHS_EXACT = new double [NX*NY];

    // compute exact solutions to compare
    // approx solutions with
    for(int i = 0; i < NX; i++){
        for(int j = 0; j < NY; j++){
           U_RHS_EXACT[returnIndex(i,j)] = 2.0*cos(X_umom[returnIndex(i,j)])*sin(X_umom[returnIndex(i,j)])
               - cos(X_umom[returnIndex(i,j)])*cos(Y_umom[returnIndex(i,j)]) - 1/RE*cos(X_umom[returnIndex(i,j)]);
           V_RHS_EXACT[returnIndex(i,j)] = -2.0*cos(Y_vmom[returnIndex(i,j)])*sin(Y_vmom[returnIndex(i,j)])
               + sin(X_vmom[returnIndex(i,j)])*sin(Y_vmom[returnIndex(i,j)]) - 1/RE*sin(Y_vmom[returnIndex(i,j)]);
           U_umom[returnIndex(i,j)] = cos(X_umom[returnIndex(i,j)]);
           V_vmom[returnIndex(i,j)] = sin(Y_vmom[returnIndex(i,j)]);
           U_in[returnIndex(i,j)] = U_umom[returnIndex(i,j)]; 
           V_in[returnIndex(i,j)] = V_vmom[returnIndex(i,j)]; 
        }
    }

    // store variables in U_in and v_in
    for(int i = 0; i < NX; i++){
        for(int j = 0; j < NY; j++){
           U_umom[returnIndex(i,j)] = cos(X_umom[returnIndex(i,j)]);
           V_vmom[returnIndex(i,j)] = sin(Y_vmom[returnIndex(i,j)]);
           U_in[returnIndex(i,j)] = U_umom[returnIndex(i,j)]; 
           V_in[returnIndex(i,j)] = V_vmom[returnIndex(i,j)]; 
        }
    }

    // Copy host to device data
    cudaMemcpy(d_U_in, U_in, NX * NY * sizeof(double), 
        cudaMemcpyHostToDevice);
    cudaMemcpy(d_V_in, V_in, NX * NY * sizeof(double), 
        cudaMemcpyHostToDevice);

    // call compure RHS function
    computeRHS(U_in, V_in);
    computeErrorRHS(U_RHS_EXACT,V_RHS_EXACT,U_in,V_in, MaxErrorx, MaxErrory);
    printf("CPP---------------------\n");
    printf("maxErrorx %4.7f\n", *MaxErrorx);
    printf("maxErrory %4.7f\n", *MaxErrory);

    //printf("---------------------\n");
    //for(int i = 0; i < NX; i++){
    //    for(int j = 0; j < NY; j++){    
    //        printf("U_out[%d][%d]: [%4.4f]\n",
    //            i,j,U_out[returnIndex(i,j)]);
    //    }
    //}
    //printf("CPP---------------------\n");

    call_cudacomputeRHS(d_U_in,d_V_in,d_U_out,d_V_out,
         DX, DY, LY, LX, RE, NX, NY);
    cudaMemcpy(V_out, d_V_out, NX * NY * sizeof(double), 
        cudaMemcpyDeviceToHost);
    cudaMemcpy(U_out, d_U_out, NX * NY * sizeof(double), 
        cudaMemcpyDeviceToHost);
    //call_cudacomputeRHS(d_U_in,d_V_in,
    //     DX, cmake --build .DY, LY, LX, RE, NX, NY);
    //cudaMemcpy(V_out, d_V_in, NX * NY * sizeof(double), 
    //    cudaMemcpyDeviceToHost);
    //cudaMemcpy(U_out, d_U_in, NX * NY * sizeof(double), 
    //    cudaMemcpyDeviceToHost);
    computeErrorRHS(U_RHS_EXACT,V_RHS_EXACT,U_out,V_out, MaxErrorx, MaxErrory);
    printf("CUDA---------------------\n");
    printf("maxErrorx %4.7f\n", *MaxErrorx);
    printf("maxErrory %4.7f\n", *MaxErrory);
    printf("---------------------\n");
    delete [] U_RHS_EXACT; 
    delete [] V_RHS_EXACT; 
    cudaFree(d_U_in);
    cudaFree(d_V_in);
    cudaFree(d_U_out);
    cudaFree(d_V_out);
#endif

    delete [] x_umom;
    delete [] y_vmom;
    delete [] y_umom;
    delete [] x_vmom;
    delete [] X_umom;
    delete [] Y_umom;
    delete [] X_vmom;     
    delete [] Y_vmom;     
    delete [] U_umom;
    delete [] V_vmom;
    delete [] U_in;
    delete [] V_in;
    delete [] U_out;
    delete [] V_out;
    delete [] CURLZ;
    cudaFree(d_U_in);
    cudaFree(d_V_in);
    cudaFree(d_U_out);
    cudaFree(d_V_out);
    cudaFree(d_k1_mod);
    cudaFree(d_k1_mod);
}