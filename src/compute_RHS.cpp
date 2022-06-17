#include "../include/functionCalls.h" 
#include "../include/variableCalls.h"

void computeRHS(double *U_out, double *V_out){

    double *U_in      = new double [NY * NX];
    double *V_in      = new double [NX * NY];

    for (int i = 0; i < NX; ++i) {
        for (int j = 0; j < NY; ++j) {
            U_in[returnIndex(i,j)]= U_out[returnIndex(i,j)];
        }
    }
    for (int i = 0; i < NX; ++i) {
        for (int j = 0; j < NY; ++j) {
            V_in[returnIndex(i,j)]= V_out[returnIndex(i,j)];
        }
    }
       double duvdy = 0.0;
       double duvdx = 0.0;
       double duudx = 0.0;
       double dvvdy = 0.0;
       double du2dxdx = 0.0;
       double du2dydy = 0.0;
       double dv2dxdx = 0.0;
       double dv2dydy = 0.0;
            
       int i1  = 0;
       int i2  = 0;
       int i_1 = 0;
       int j1  = 0;
       int j2  = 0;
       int j_1 = 0;

    for(int i = 0; i < NX; i++){
        i1  = i + 1;
        i2  = i + 2;
        i_1 = i - 1;

        if (i == 0){
            i_1 = NX - 1;
        } else if( i == NX - 2){
            i2 = 0;
        } else if( i == NX - 1){
            i1 = 0;
            i2 = 1;
        }

        for(int j = 0; j < NY; j++){
            j1  = j + 1;
            j2  = j + 2;
            j_1 = j - 1;
            
            if (j == 0){
                j_1 = NY - 1;
            } else if( j == NY - 2){
                j2 = 0;
            } else if( j == NY - 1){
                j1 = 0;
                j2 = 1;
            }
            
            duudx = 0.5 * (-U_in[returnIndex(i_1, j)]*U_in[returnIndex(i_1, j)] + U_in[returnIndex(i1, j)]*U_in[returnIndex(i1, j)])/DX;
            dvvdy = 0.5 * (-V_in[returnIndex(i, j_1)]*V_in[returnIndex(i, j_1)] + V_in[returnIndex(i, j1)]*V_in[returnIndex(i, j1)])/DY;

            du2dxdx = (U_in[returnIndex(i1, j)]  - 2*U_in[returnIndex(i, j)] + U_in[returnIndex(i_1, j)])/(DX*DX);
            du2dydy = (U_in[returnIndex(i, j1)]  - 2*U_in[returnIndex(i, j)] + U_in[returnIndex(i, j_1)])/(DY*DY);
            dv2dxdx = (V_in[returnIndex(i1, j)]  - 2*V_in[returnIndex(i, j)] + V_in[returnIndex(i_1, j)])/(DX*DX);
            dv2dydy = (V_in[returnIndex(i, j1)]  - 2*V_in[returnIndex(i, j)] + V_in[returnIndex(i, j_1)])/(DY*DY);

            duvdy = 0.5 * (U_in[returnIndex(i, j1)] * 0.25 * (V_in[returnIndex(i_1, j2)] + V_in[returnIndex(i, j2)] + 
                V_in[returnIndex(i_1, j1)] + V_in[returnIndex(i, j1)]))/DY - 0.5 * (U_in[returnIndex(i, j_1)] * 0.25 *(V_in[returnIndex(i, j)] + V_in[returnIndex(i_1, j)]
                + V_in[returnIndex(i, j_1)] + V_in[returnIndex(i_1, j_1)]))/DY;
            duvdx = 0.5 * (V_in[returnIndex(i, j)] * 0.25 * (U_in[returnIndex(i1, j)] + U_in[returnIndex(i2, j)] + 
                U_in[returnIndex(i1, j_1)] + U_in[returnIndex(i2, j_1)]))/DX - 0.5 * (V_in[returnIndex(i_1, j)] * 0.25 *(U_in[returnIndex(i, j)] + U_in[returnIndex(i_1, j)]
                + U_in[returnIndex(i_1, j_1)] + U_in[returnIndex(i, j_1)]))/DX;

            U_out[returnIndex(i, j)] = -duvdy - duudx + 1/RE * (du2dxdx + du2dydy);
            V_out[returnIndex(i, j)] = -duvdx - dvvdy + 1/RE * (dv2dxdx + dv2dydy);

        }

    }

    delete [] U_in;
    delete [] V_in;
}