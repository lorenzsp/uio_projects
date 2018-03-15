#include <iostream>
#include <fstream>
#include <iomanip>
#include <stdlib.h>
#include <cmath>
#include "time.h"
#include "lib.h"

using namespace std;
ofstream ofile;

int main()
{
    int n=200; //define mesh points, size matrix

    double* d = new double[n+1]; //initialize diagonal elements of our problem matrix
    double* e = new double[n]; //initialize upper and lower diagonal elements of our problem matrix
    d[0]=e[0]=0; //set zero the first value of d and e to have easy counting

    double ** R = new double*[n]; // initialize matrix of eigenvectors
    for(int i = 0; i < n; i++){
        R[i] = new double[n];
    }
    for(int i=0; i<n; i++){
        for(int j=0; j<n; j++){
            if(i==j){
                R[i][j] = 1;
            }
            else{
                R[i][j] = 0;
            }
        }
    }

    double rho_0 = 0;
    double rho_n = 25;
    double h = (rho_n - rho_0)/n;
    double *rho = new double[n+1];
    rho[0] = rho_0;
    rho[n]= rho_n;
    for(int i=1; i<n; i++){
        rho[i] = i*h;
    }

    double *V = new double[n+1]; //initialize potential
    V[0] = rho[0];
    V[n] = 0;
    double omega=0.25; //define omega
    for(int i = 1; i<n; i++){
        V[i] = rho[i]*rho[i]; //potential for non-interacting case
    }

//    for(int i = 1; i<n; i++){
//        V[i] = rho[i]*rho[i]*omega*omega + 1/rho[i]; //potential for interacting case
//    }


    for(int i =1; i<n; i++){
        d[i]=2.0/(h*h)+V[i]; //define diagonal elements of our problem matrix
    }

    for(int i =1; i<n-1; i++){
        e[i]=-1.0/(h*h); //define upper and lower diagonal elements of our problem matrix
    }

    //start time
    clock_t start, finish;
    start = clock();

    tqli(d, e, n, R); //library function which take d,e,n,R and gives the eigenvalues saved in d

    //stop time
    finish = clock();
    double timeused = (double) (finish - start)/(CLOCKS_PER_SEC );
    cout <<"Run time: "<< timeused <<endl; //print run time

    //order eigenvalues
    double a,b;
    int it=0;
    while(it<n){
    for(int i=1; i<n; i++){
       a = d[i-1];
       b = d[i];
      if(b < a){
             d[i] = a;
             d[i-1] = b;
       }
    } it++;
    }

    cout <<"Eigenvalues:"<< endl;
    for(int i=1;i<4; i++){
        cout << d[i] << endl; //print first three eigenvalues
    }

//free memory
    for (int i = 0; i < n; i++){
        delete[] R[i];
        }

    delete[] R;
    delete[] V;
    delete[] d;
    delete[] e;
    delete[] rho;

    return 0;
}
