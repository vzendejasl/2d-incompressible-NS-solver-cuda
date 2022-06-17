#include <cmath>
#include <stdio.h>

#include "../include/functionCalls.h"
#include "../include/variableCalls.h"

#include "../src/cuda_helper.cuh"
void computeErrorRHS(double *U_RHS_EXACT,double *V_RHS_EXACT, 
    double *U_out,double *V_out,
    double *MaxErrorx, double *MaxErrory){

    double errorx = 0.0;
    double errory = 0.0;
    for(int i = 0; i < NX; i++){
        for(int j = 0; j < NY; j++){
           errorx = fabs(U_RHS_EXACT[returnIndex(i,j)] - U_out[returnIndex(i,j)]);
           errory = fabs(V_RHS_EXACT[returnIndex(i,j)] - V_out[returnIndex(i,j)]);

           if (errorx > *MaxErrorx){
               *MaxErrorx = errorx;
           }
           if (errory > *MaxErrory){
               *MaxErrory = errory;
           }
        }
    }
}

//void computeErrorRHS(const Matrix& U_RHS_EXACT,const Matrix& V_RHS_EXACT, 
//    const Matrix& U_out,const Matrix& V_out,
//    double *MaxErrorx, double *MaxErrory){
//
//    double errorx = 0.0;
//    double errory = 0.0;
//    for(int i = 0; i < NX; i++){
//        for(int j = 0; j < NY; j++){
//           errorx = fabs(U_RHS_EXACT(i, j) - U_out(i,j));
//           errory = fabs(V_RHS_EXACT(i, j) - V_out(i,j));
//
//           if (errorx > *MaxErrorx){
//               *MaxErrorx = errorx;
//           }
//           if (errory > *MaxErrory){
//               *MaxErrory = errory;
//           }
//        }
//    }
//}