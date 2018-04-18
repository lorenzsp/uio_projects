#include "catch.hpp"
#include "jacobi.h"
#include "define_matrices.h"
#include <stdlib.h>
#include <cmath>

//function which calculate the Frobenius norm for a already defined matrix A
double frob(double ** A, int n){
    double off_A;
    for(int i=0; i<n-1; i++){
        for(int j=0; j<n-1; j++){
                off_A+=A[i][j]*A[i][j];
        }
    }
    return sqrt(off_A);
}


TEST_CASE("Eigenvalue:"){
    int n=200;
    int iterations=0;
    double timeused=0;
    double rho_min=0;
    double rho_max=25;
    int interacting=0;
    int eps = 1.0e-8;
    int method=0;
    double omega=0.25;

        double* eig_v = new double[n-1];
        double ** A = new double*[n-1];
        for(int i = 0; i < n-1; i++){
            A[i] = new double[n-1];
        }
        double ** R = new double*[n-1];
        for(int i = 0; i < n-1; i++){
            R[i] = new double[n-1];
        }
 define_matrices(rho_min, rho_max, n, interacting, A, R, omega);
 jacobi(n, A,R, eig_v, eps, method,iterations, timeused);

 REQUIRE(round(eig_v[0])==3);
 REQUIRE(round(eig_v[1])==7);
 REQUIRE(round(eig_v[2])==11);

 delete[] eig_v;
 for (int i = 0; i < n-1; i++){
     delete[] A[i];
     delete[] R[i];
 }

 delete[] A;
 delete[] R;

 double* eig_v2 = new double[n-1];
 double ** A2 = new double*[n-1];
 for(int i = 0; i < n-1; i++){
     A2[i] = new double[n-1];
 }
 double ** R2 = new double*[n-1];
 for(int i = 0; i < n-1; i++){
     R2[i] = new double[n-1];
 }

 method=1;
 define_matrices(rho_min, rho_max, n, interacting, A2, R2, omega);
 jacobi(n, A2,R2, eig_v2, eps, method,iterations, timeused);
 REQUIRE(round(eig_v2[0])==3);
 REQUIRE(round(eig_v2[1])==7);
 REQUIRE(round(eig_v2[2])==11);

 delete[] eig_v2;
 for (int i = 0; i < n-1; i++){
     delete[] A2[i];
     delete[] R2[i];
 }

 delete[] A2;
 delete[] R2;
}

TEST_CASE("Eigenvalue2:"){
    //define parameters to test, we chose n=365, because for the cyclic Jacobi we get an approximation to 1.25 only after n=361

    int n=365;
    int iterations=0;
    double timeused=0;
    double rho_min=0;
    double rho_max=25;
    int interacting;
    int eps = 1.0e-8;
    int method;
    double omega=0.25;

        double* eig_v = new double[n-1];
        double ** A = new double*[n-1];
        for(int i = 0; i < n-1; i++){
            A[i] = new double[n-1];
        }
        double ** R = new double*[n-1];
        for(int i = 0; i < n-1; i++){
            R[i] = new double[n-1];
        }
 interacting=1;
 method=0;
 define_matrices(rho_min, rho_max, n, interacting, A, R, omega);
 jacobi(n, A,R, eig_v, eps, method,iterations, timeused);
 REQUIRE(round(eig_v[0]*1e2)/(1e2)== Approx(1.25));

 delete[] eig_v;
 for (int i = 0; i < n-1; i++){
     delete[] A[i];
     delete[] R[i];
 }

 delete[] A;
 delete[] R;

 double* eig_v2 = new double[n-1];
 double ** A2 = new double*[n-1];
 for(int i = 0; i < n-1; i++){
     A2[i] = new double[n-1];
 }
 double ** R2 = new double*[n-1];
 for(int i = 0; i < n-1; i++){
     R2[i] = new double[n-1];
 }

 method=1;
 define_matrices(rho_min, rho_max, n, interacting, A2, R2, omega);
 jacobi(n, A2,R2, eig_v2, eps, method,iterations, timeused);
 REQUIRE(round(eig_v2[0]*1e2)/(1e2)== Approx(1.25));
 delete[] eig_v2;
 for (int i = 0; i < n-1; i++){
     delete[] A2[i];
     delete[] R2[i];
 }

 delete[] A2;
 delete[] R2;

}

TEST_CASE("Testing conservation of Frobenius norm:"){

    //define parameters to test
    int n=10;
    int iterations=0;
    double timeused=0;
    double rho_min=0;
    double rho_max=10;
    int interacting=0;
    int eps = 1.0E-8;
    int method=0;
    double omega=0.25;

        double* eig_v = new double[n-1];
        double ** A = new double*[n-1];
        for(int i = 0; i < n-1; i++){
            A[i] = new double[n-1];
        }
        double ** R = new double*[n-1];
        for(int i = 0; i < n-1; i++){
            R[i] = new double[n-1];
        }
 define_matrices(rho_min, rho_max, n, interacting, A, R, omega);
 double f1=frob(A,n);
 jacobi(n, A,R, eig_v, eps, method,iterations, timeused);
 double f2=frob(A,n);

 REQUIRE(round(f1)==round(f2));

 delete[] eig_v;
 for (int i = 0; i < n-1; i++){
     delete[] A[i];
     delete[] R[i];
 }

 delete[] A;
 delete[] R;


 double* eig_v2 = new double[n-1];
 double ** A2 = new double*[n-1];
 for(int i = 0; i < n-1; i++){
     A2[i] = new double[n-1];
 }
 double ** R2 = new double*[n-1];
 for(int i = 0; i < n-1; i++){
     R2[i] = new double[n-1];
 }
 method=1;
 define_matrices(rho_min, rho_max, n, interacting, A2, R2, omega);
 f1=frob(A2,n);
 jacobi(n, A2,R2, eig_v2, eps, method,iterations, timeused);
 f2=frob(A2,n);

 REQUIRE(round(f1)==round(f2));

 delete[] eig_v2;
 for (int i = 0; i < n-1; i++){
     delete[] A2[i];
     delete[] R2[i];
 }

 delete[] A2;
 delete[] R2;

}



TEST_CASE("Testing conservation of Frobenius norm 2:"){

    //define parameters to test
    int n=10;
    int iterations=1;
    double timeused=0;
    double rho_min=0;
    double rho_max=10;
    int interacting=0;
    int eps = 1.0E-8;
    int method=0;
    double omega=0.25;

        double* eig_v = new double[n-1];
        double ** A = new double*[n-1];
        for(int i = 0; i < n-1; i++){
            A[i] = new double[n-1];
        }
        double ** R = new double*[n-1];
        for(int i = 0; i < n-1; i++){
            R[i] = new double[n-1];
        }
 define_matrices(rho_min, rho_max, n, interacting, A, R, omega);
 double f1=frob(A,n);
 jacobi(n, A,R, eig_v, eps, method,iterations, timeused);
 double f2=frob(A,n);

 REQUIRE(round(f1)==round(f2));

 delete[] eig_v;
 for (int i = 0; i < n-1; i++){
     delete[] A[i];
     delete[] R[i];
 }

 delete[] A;
 delete[] R;


 double* eig_v2 = new double[n-1];
 double ** A2 = new double*[n-1];
 for(int i = 0; i < n-1; i++){
     A2[i] = new double[n-1];
 }
 double ** R2 = new double*[n-1];
 for(int i = 0; i < n-1; i++){
     R2[i] = new double[n-1];
 }
 method=1;
 define_matrices(rho_min, rho_max, n, interacting, A2, R2, omega);
 f1=frob(A2,n);
 jacobi(n, A2,R2, eig_v2, eps, method,iterations, timeused);
 f2=frob(A2,n);

 REQUIRE(round(f1)==round(f2));

 delete[] eig_v2;
 for (int i = 0; i < n-1; i++){
     delete[] A2[i];
     delete[] R2[i];
 }

 delete[] A2;
 delete[] R2;

}
