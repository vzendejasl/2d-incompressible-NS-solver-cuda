#include "../include/functionCalls.h" 
#include "../include/variableCalls.h"
#include <complex>
#include <iostream>

using namespace std;

void compute_projection_step(double *U_star,double *V_star, 
    double *U_n1, double *V_n1, double *k1mod, double *k2mod) {

    double *div = new double [NX * NY];
    double *P   = new double [NX * NY];

    // Compute divergence (dudx + dvdy)
    // Compute dudx
    for (int i = 0; i < NX-1; ++i) {
        for (int j = 0; j < NY; ++j) {
            div[returnIndex(i,j)] = (U_star[returnIndex(i+1,j)]
                -U_star[returnIndex(i,j)])/DX;
        }
    }
    // Handle periodic BC for dudx
    for (int j = 0; j < NY; ++j) {
        div[returnIndex(NX-1,j)] = (U_star[returnIndex(0,j)]
            -U_star[returnIndex(NX-1,j)])/DX;
    }

    // Compute dvdy
    for (int i = 0; i < NX; ++i) {
       for (int j = 0; j < NY-1; ++j) {
           //dvdy
            div[returnIndex(i,j)] += (V_star[returnIndex(i,j+1)]
                -V_star[returnIndex(i,j)])/DY;
       }
    }
    // Handle period BC for dvdy
    for (int i = 0; i < NX; ++i) {
        div[returnIndex(i,NY - 1)] += (V_star[returnIndex(i,0)]
            -V_star[returnIndex(i,NY - 1)])/DY;
    }

    complex<double> *f_hat = new complex<double> [NX * NY];
    complex<double> *P_hat = new complex<double> [NX * NY];

    double2complex(div,f_hat);
    FFT2D(f_hat,NX,NY,1);
    fft_shift(f_hat);

    for (int i = 0; i < NX; ++i) {
        for (int j = 0; j < NY; ++j) {
            P_hat[returnIndex(i,j)].real(-real(f_hat[returnIndex(i,j)]) 
                / (k1mod[returnIndex(i,j)]+k2mod[returnIndex(i,j)]));
            P_hat[returnIndex(i,j)].imag(-imag(f_hat[returnIndex(i,j)]) 
                / (k1mod[returnIndex(i,j)]+k2mod[returnIndex(i,j)]));
        }
    }
    
    P_hat[returnIndex(NX/2,NY/2)].real(0);
    P_hat[returnIndex(NX/2,NY/2)].imag(0);

    fft_shift(P_hat);

    FFT2D(P_hat,NX,NY,-1);

    for (int i = 0; i < NX; ++i) {
        for (int j = 0; j < NY; ++j) {
            P[returnIndex(i,j)] = real(P_hat[returnIndex(i,j)]);
        }
    }

    for (int i = 0; i < NX-1; ++i) {
        for (int j = 0; j < NY; ++j) {
            U_n1[returnIndex(i+1,j)]= -(P[returnIndex(i+1,j)]-P[returnIndex(i,j)])/DX
                + U_star[returnIndex(i+1,j)];
        }
    }
    for (int j = 0; j < NY; ++j) {
        U_n1[returnIndex(0,j)]= -(P[returnIndex(0,j)]-P[returnIndex(NX-1,j)])/DX
            + U_star[returnIndex(0,j)];
    }

    for (int i = 0; i < NX; ++i) {
        for (int j = 0; j < NY-1; ++j) {
            V_n1[returnIndex(i,j+1)]= -(P[returnIndex(i,j+1)]-P[returnIndex(i,j)])/DY
                + V_star[returnIndex(i,j+1)];
        }
    }
    for (int i = 0; i < NX; ++i) {
        V_n1[returnIndex(i,0)]= -(P[returnIndex(i,0)]-P[returnIndex(i,NY-1)])/DY
            + V_star[returnIndex(i,0)];
    }
    
    delete [] div;
    delete [] P;
    delete [] f_hat;
    delete [] P_hat;
}
