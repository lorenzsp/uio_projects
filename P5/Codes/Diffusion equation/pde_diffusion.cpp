#include "pde_diffusion.h"
#include <iostream>
#include <fstream>
#include <iomanip>
#include <stdlib.h>
#include <cmath> //libreria migliore
#include <string>
#include <array>
#include "data_analysis.h"
#include "pde_diffusion.h"

using namespace std;
using std::array;


// declaration of the number of points in the interval
// and of the alpha factor = dt/(dx*dx)
pde_diffusion::pde_diffusion(int n, double a)
{
    N = n;
    alpha = a;
    r_cn = 2-2*alpha;
    r_i = 1 + 2*alpha;
    r_e = 1 - 2*alpha;
    r_2d = 1+4*alpha;
}


// solver for the tridiagonal matrix of our PDE solvers
// B is the diagonal element, A and C are the off-diagonal elements
// S is the vecotor we want to find, Q is the knwon vector
// n is the dimensionality of S and Q
void tridiagonal_solver(double A, double B, double C, int n, double* S, double* Q){

    // define the new diagonal element
    double *d = new double[n];
    d[0] = 0;
    d[1] = B;
    // new Q vector
    double *q_new = new double[n];
    q_new[0] = 0;
    q_new[1] = Q[1];

    // forward substitution
    for(int i=2; i<n-2; i++){
        d[i] = B - A*C / d[i-1];
        q_new[i] = Q[i] - A*q_new[i-1] / d[i-1];
    }


    d[n-2]=B-A*C/d[n-3];
    q_new[n-2] = Q[n-2] - A*q_new[n-3] / d[n-3];
    S[n-2] = (-C +q_new[n-2]) / d[n-2];
    // backward substitution
    for(int i=n-3; i>=1; i--){
        S[i] = (q_new[i] - C*S[i+1])/ d[i];
    }

    delete[] d;
    delete[] q_new;

}


// solve numerically the two-dimensional diffusion equation
void pde_diffusion::two_dimension(double** u, double tolerance, int cutoff){

    // variables declaration
    double** v = new double*[N+1];
    double** temp = new double*[N+1];
    for(int i=0;i<N+1;i++){
        v[i]=new double[N+1];
        temp[i]=new double[N+1];
    }

    double delta = 0.0;

    // define the actual solution
    for(int i=0; i<N+1;i++){
        for(int j=0; j<N+1;j++){
            v[i][j]=u[i][j];

        }
    }

    // declare number of iterations
    int iterations = 0;
    // declare the variable difference for the iterative solver
    double difference=tolerance+1;

    // iterative Jacobi's method
    while((iterations<cutoff)&&(difference>tolerance)){
        difference=0;

        // define a guess
        for(int i=0; i<N+1;i++){
            for(int j=0; j<N+1;j++){
                temp[i][j]=u[i][j];
            }
        }

        // alorithm to define the solution
        for(int i=1;i<N;i++){
            for(int j=1;j<N;j++){
                delta=alpha*(temp[i+1][j]+temp[i-1][j]+temp[i][j-1]+temp[i][j+1]);
                u[i][j]=(delta+v[i][j])/(r_2d);
                // definition of the difference bewteen the actual solution and the previous one
                difference+=fabs(u[i][j]-temp[i][j]);
            }
        }

        iterations++;
        // if the iterations are bigger than the chosen limit cutoff
        if(iterations>cutoff) cout<<"Failure"<<endl;
        difference/=pow(N-1,2.0);
    }

    for(int i=0; i<N+1; i++) {
        delete[] v[i];
        delete[] temp[i];
    }
    delete[] v;
    delete[] temp;
}

// solve numerically the one-dmensional PDE with fixed boundary conditions
// using Crank Nicolson method
void pde_diffusion::Crank_nicolson(double* u)
{
    double* v = new double[N+1];

    for(int i=1; i<N; i++){
        v[i] = alpha*u[i-1] + r_cn*u[i] + alpha*u[i+1];
    }

    tridiagonal_solver(-alpha, 2+2*alpha, -alpha, N+1, u, v);

    delete[] v;
}


// solve numerically the one-dmensional PDE with fixed boundary conditions
// using Implicit method
void pde_diffusion::Implicit(double* u){

    double* v = new double[N+1];
    for(int i=0; i<N+1;i++){
        v[i]=u[i];
    }

    tridiagonal_solver(-alpha,r_i,-alpha,N+1,u,v);

    delete[] v;
}


// solve numerically the one-dmensional PDE with fixed boundary conditions
// using Explicit method
void pde_diffusion::Explicit(double* u){
    double* v = new double[N+1];
    for(int i=0; i<N+1;i++){
        v[i]=u[i];
    }
    for(int i=1; i<N; i++){
        u[i] = alpha*v[i-1] + r_e*v[i] + alpha*v[i+1];
    }

    delete[] v;
}

// analytical solution of the one-dimensional diffusion equation at a given time
// with boundary conditions: solution[0]=0 and solution[N]=1
void pde_diffusion::Solution1d(double* solution, double time)
{
    // define a trunction of the fourier series
    int truncation = 1000;

    double C = M_PI/N;
    double B = M_PI*M_PI*time;
    // one dimension solution
    for(int i=0; i<N+1; i++){
        solution[i] =(double) i/N;
        for(int j=1; j<=truncation; j++){
            if(j%2==0) solution[i] += 2*sin(i*j*C)*exp(-j*j*B)/(j*M_PI);
            else solution[i] += -2*sin(i*j*C)*exp(-j*j*B)/(j*M_PI);
        }
    }

}

// analytical solution of the two-dimensional diffusion equation at a given time
// with boundary conditions: solution[0][0]=solution[N][0]=solution[0][N]=solution[N][N]=0
void pde_diffusion::Solution2d(double** solution, double time)
{
    double B = M_PI/N;
    double C = 2*M_PI*M_PI*time;
    for(int i=0; i<N+1; i++){
        for(int j=0; j<N+1; j++){
            solution[i][j] = sin(i*B)*sin(j*B)*exp(-C);
        }
    }
}

// two functions to look for the max relative error
double pde_diffusion::max_relative_error(double* u, double* s){
    double max = fabs((u[3]-s[3])/s[3]);

    for(int i=3; i<N-2; i++){
        double value = fabs((u[i]-s[i])/s[i]);
        if(max<value){
            max = value;
        }
    }
    return max;
}

double pde_diffusion::max_relative_error_2d(double** u, double** s){
    double max = fabs((u[3][3]-s[3][3])/s[3][3]);

    for(int i=3; i<N-2; i++){
        for(int j=3; j<N-2; j++){
            double value = fabs((u[i][j]-s[i][j])/s[i][j]);
            if(max<value){
                max = value;
            }
        }
    }
    return max;
}
