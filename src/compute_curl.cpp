#include "../include/functionCalls.h" 
#include "../include/variableCalls.h"

void compute_curl(double *U_in,double *V_in, double *curlz) {

    // Compute dvdx
    for (int i = 0; i < NX-1; ++i) {
        for (int j = 0; j < NY; ++j) {
            curlz[returnIndex(i,j)] = (V_in[returnIndex(i+1,j)]
                -V_in[returnIndex(i,j)])/DX;
        }
    }
    // Handle periodic BC for dvdx
    for (int j = 0; j < NY; ++j) {
        curlz[returnIndex(NX - 1,j)] = (V_in[returnIndex(0,j)]
             -V_in[returnIndex(NX - 1,j)])/DX;
    }

    // Compute dudy
    for (int i = 0; i < NX; ++i) {
        for (int j = 0; j < NY-1; ++j) {
            curlz[returnIndex(i,j)] -= (U_in[returnIndex(i,j+1)]
                -U_in[returnIndex(i,j)])/DY;
        }
    }
    // Handle periodic BC for dudy
    for (int i = 0; i < NX; ++i) {
        curlz[returnIndex(i,NY - 1)] -= (U_in[returnIndex(i,0)]
             -U_in[returnIndex(i,NY-1)])/DY;
    }

}
