#include "define_matrices.h"
#include <iostream>
#include <fstream>
#include <iomanip>
#include <stdlib.h>
#include <cmath>
#include "time.h"

void define_matrices(double rho_min,double rho_max, int n, int interacting, double** A, double** R, double omega){
    //dimension A=n-1; dimension R=n-1;
double h = (rho_max - rho_min)/n;
double *rho = new double[n+1];
rho[0] = rho_min;
rho[n] = rho_max;
for(int i=1; i<n; i++){
    rho[i] = i*h;
}

//fill matrix of eigenvectors
for(int p=0; p<n-1; p++){
        for(int q=0; q<n-1; q++){
            if(p==q){
                R[p][q] = 1;
            }
            else{
                R[p][q] = 0;
            }
        }
}


//fill matrix A our our problem depending on the flag interacting: if it's 0, it'll set the non-interacting case,
//else if it's 1, it'll set the interacting case
for(int i=0; i<n-1; i++){
    for(int j=0; j<n-1; j++){
        if(i==j&& interacting==0){
            A[i][j] = rho[i+1]*rho[i+1]+2.0/(h*h);
        }
        else if(i==j&&interacting==1){
            A[i][j] = omega*omega*rho[i+1]*rho[i+1]+2.0/(h*h) + 1.0/rho[i+1];
        }
        else if(i==j+1 || i==j-1){
            A[i][j]  =  -1.0/(h*h);
        }
        else{
            A[i][j]  = 0;
        }
    }
}
}
