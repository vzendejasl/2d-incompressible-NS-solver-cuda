#include <iostream>
#include <stdio.h>
#include <cmath>

#include "../include/functionCalls.h"
#include "../include/variableCalls.h"
 
int NX = 100;
int NY = NX;
double PI = 3.14159265358;
double LX = 2*PI;
double LY = LX;

double RE = 1.0;
double DX = LX/NX;
double DY = LY/NY;

// int main()
// {

//     double *MaxErrorx = (double *) malloc(sizeof (double));
//     double *MaxErrory = (double *) malloc(sizeof (double));
//     *MaxErrorx = 0.0;
//     *MaxErrory = 0.0;

//     double *x_umom = (double *) malloc(NX * sizeof (double));
//     double *y_vmom = (double *) malloc(NY * sizeof (double));
//     double *y_umom = (double *) malloc(NY * sizeof (double));
//     double *x_vmom = (double *) malloc(NX * sizeof (double));

//     double **X_umom = (double **) malloc(NX * sizeof (double *));
//     double **Y_umom = (double **) malloc(NY * sizeof (double *));
//     double **X_vmom = (double **) malloc(NX * sizeof (double *));
//     double **Y_vmom = (double **) malloc(NY * sizeof (double *));

//     double **U_RHS_EXACT = (double **) malloc(NX * sizeof (double *));
//     double **V_RHS_EXACT = (double **) malloc(NY * sizeof (double *));

//     double **U_umom = (double **) malloc(NX * sizeof (double *));
//     double **V_vmom = (double **) malloc(NY * sizeof (double *));
//     double **U_in   = (double **) malloc(NX * sizeof (double *));
//     double **V_in   = (double **) malloc(NY * sizeof (double *));
//     double **U_out   = (double **) malloc(NX * sizeof (double *));
//     double **V_out   = (double **) malloc(NY * sizeof (double *));

//     for (int i = 0; i < NX; i++){
//         X_umom[i]      = (double *)malloc(NX*sizeof (double *));
//         X_vmom[i]      = (double *)malloc(NX*sizeof (double *));
//         U_RHS_EXACT[i] = (double *)malloc(NX*sizeof (double *));
//         U_umom[i]      = (double *)malloc(NX*sizeof (double *));
//         U_in[i]        = (double *)malloc(NX*sizeof (double *));
//         U_out[i]       = (double *)malloc(NX*sizeof (double *));
//     }
//     for (int i = 0; i < NY; i++){
//         Y_umom[i]      = (double *)malloc(NY*sizeof (double *));
//         Y_vmom[i]      = (double *)malloc(NY*sizeof (double *));
//         V_RHS_EXACT[i] = (double *)malloc(NY*sizeof (double *));
//         V_vmom[i]      = (double *)malloc(NY*sizeof (double *));
//         V_in[i]        = (double *)malloc(NY*sizeof (double *));
//         V_out[i]       = (double *)malloc(NY*sizeof (double *));
//     }

//     if (!x_umom){
//         printf("Out of memory\n");
//         exit(-1);
//         }
//     if (!y_vmom){
//         printf("Out of memory\n");
//         exit(-1);
//         }
//     if (!y_umom){
//         printf("Out of memory\n");
//         exit(-1);
//         }
//     if (!x_vmom){
//         printf("Out of memory\n");
//         exit(-1);
//         }
//     if (!X_umom){
//         printf("Out of memory\n");
//         exit(-1);
//         }
//     if (!X_vmom){
//         printf("Out of memory\n");
//         exit(-1);
//         }
//     if (!Y_umom){
//         printf("Out of memory\n");
//         exit(-1);
//         }
//     if (!Y_vmom){
//         printf("Out of memory\n");
//         exit(-1);
//         }
//     if (!U_RHS_EXACT){
//         printf("Out of memory\n");
//         exit(-1);
//         }
//     if (!V_RHS_EXACT){
//         printf("Out of memory\n");
//         exit(-1);
//         }
//     if (!U_umom){
//         printf("Out of memory\n");
//         exit(-1);
//         }
//     if (!V_vmom){
//         printf("Out of memory\n");
//         exit(-1);
//         }
//     if (!U_out){
//         printf("Out of memory\n");
//         exit(-1);
//         }
//     if (!V_out){
//         printf("Out of memory\n");
//         exit(-1);
//         }


//     // set up grid space
//     for(int j = 0; j < NY; j++){
//         y_umom[j] = j*DY + DY/2;
//         x_vmom[j] = j*DY + DY/2;
//     }
//     for(int i = 0; i < NX; i++){
//         x_umom[i] = i*DX;
//         y_vmom[i] = i*DX;
//     }

//     // "meshgrid"
//     for(int i = 0; i < NX; i++){
//         for(int j = 0; j < NY; j++){
//            Y_umom[i][j] = y_umom[j];
//            X_umom[i][j] = x_umom[i];
           
//            X_vmom[i][j] = x_vmom[i];
//            Y_vmom[i][j] = y_vmom[j];
//         }
//     }

//     // compute exact solutions to compare
//     // approx solutions with
//     for(int i = 0; i < NX; i++){
//         for(int j = 0; j < NY; j++){
//            U_RHS_EXACT[i][j] = 2.0*cos(X_umom[i][j])*sin(X_umom[i][j])
//                - cos(X_umom[i][j])*cos(Y_umom[j][j]) - 1/RE*cos(X_umom[i][j]);
//            V_RHS_EXACT[i][j] = -2.0*cos(Y_vmom[i][j])*sin(Y_vmom[i][j])
//                + sin(X_vmom[i][j])*sin(Y_vmom[j][j]) - 1/RE*sin(Y_vmom[i][j]);
//            U_umom[i][j] = cos(X_umom[i][j]);
//            V_vmom[i][j] = sin(Y_vmom[i][j]);
//            U_in[i][j] = U_umom[i][j]; 
//            V_in[i][j] = V_vmom[i][j]; 
//         }
//     }

//     // store variables in U_in and v_in
//     for(int i = 0; i < NX; i++){
//         for(int j = 0; j < NY; j++){
//            U_umom[i][j] = cos(X_umom[i][j]);
//            V_vmom[i][j] = sin(Y_vmom[i][j]);
//            U_in[i][j] = U_umom[i][j]; 
//            V_in[i][j] = V_vmom[i][j]; 
//         }
//     }

//     // call compure RHS function
//     // returns U_out and Vout
//     computeRHS(U_in, V_in, U_out, V_out);

//     computeErrorRHS(U_RHS_EXACT,V_RHS_EXACT,U_out,V_out, MaxErrorx, MaxErrory);

//     printf("maxErrorx %4.7f\n", *MaxErrorx);
//     printf("maxErrory %4.7f\n", *MaxErrory);

//     // Free Variables
//     free(x_umom);
//     free(y_umom);
//     free(x_vmom);
//     free(y_vmom);

//     for (int i = 0; i < NX; i++){
//         free(U_RHS_EXACT[i]);
//         free(V_RHS_EXACT[i]);
//         free(U_in[i]);
//         free(V_in[i]);
//         free(X_umom[i]);
//         free(Y_umom[i]);
//         free(X_vmom[i]);
//         free(Y_vmom[i]);
//         free(U_out[i]);
//         free(V_out[i]);
//     }

//     free(U_RHS_EXACT);
//     free(V_RHS_EXACT);
//     free(U_in);
//     free(V_in);
//     free(X_umom);
//     free(Y_umom);
//     free(X_vmom);
//     free(Y_vmom);
    
// }

// print statements for dubugging

//for(int i = 0; i < NX; i++){
//    for(int j = 0; j < NY; j++){    
//        printf("U_out[%d][%d]: [%4.4f]\n",
//            i,j,U_out[i][j]);
//        printf("V_out[%d][%d]: [%4.4f]\n",
//           i,j,V_out[i][j]);
//    }
//}

//printf("i:%d\n", i);
//printf("x_umom[%d]: %4.4f\n",i, x_umom[i]);
//printf("x_umom[%d]: %4.4f\n",i, x_umom[i]);
//printf("U_RHS_EXACT[%d][%d]: [%4.4f]\n",
//    i,j,U_RHS_EXACT[i][j]);