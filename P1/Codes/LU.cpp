#include <iostream>
#include <fstream>
#include <iomanip>
#include <stdlib.h>
#include <cmath>
#include "time.h"
#include <armadillo>

using namespace std;
using namespace arma;

ofstream ofile; /*print to file*/

//function used
double f(double x){return 100*exp(-10*x);};

int main(int argc, char* argv[])
{

    // we rewrite the poisson equation -u''(x) = f(x)
    // as a set of linear equations A v = r with A matrix (n,n)
    int n; // number of points in the interval (0,1)
    double s; //step from one n to the other
    int q; // number of general loop's repetions
                            // it will repeat the gaussian elimination q times

    if (argc<=3){
             cout <<"You forgot something: you have to n (int, number of starting points); t (double, step from one n to the other); "
                    "q (int, number of general loop's repetions)." << endl;
             exit(1);
    } else {
        n=atoi(argv[1]);
        s=atof(argv[2]);
        q=atoi(argv[3]);
    }
    if (s*n==(int) s*n && s>0 && n>0 && q>0) {
             cout <<"Your setup is good." << endl;
        } else {
             cout <<"Your setup is wrong: you must discretize the interval (0,1) and n,t, q must be positive." << endl;
             exit(1);
    }

    mat time(q,2); // time matrix: it's useful to have a file with elapsed time
    for(int ii=0;ii<=q-1;ii++){
    //Declaration and allocation of variables with Armadillo
    mat A=zeros<mat>(n+2,n+2);
    mat C=randu<mat>(n+2,n+2);
    vec x=randu<vec>(n+2);
    vec r=randu<vec>(n+2);

    double h=1.0/(n+1.0); //number of steps

    int i,j; //indexes

    //time calculation
        clock_t t;
        t=clock();

    //filling the tridiagonal matrix A
    for(i=1;i<=n+1;i++){
        for(j=1;j<=n+1;j++){
            if(i==j){
                A(i,j)=2.0;
            } else if (j==i+1||(i>=2&&j==i-1)){
                A(i,j)=-1.0;
            } else  {
                A(i,j)=0.0;
            }
        }
    }
    //A.print(); /*print to check*/

    //numerical solution
    for (i=0; i<=n+1; i++){
        x[i]=i*h;
        r[i]=h*h*f(x[i]); /*arranging the right side of equation -u''=r(x)=h^2*f(x)*/
//        cout << "x"<<i<<":"<<x[i]<< endl;
    }

    vec v=solve(A,r); /*numerical solution with LU, and I print it*/
    //v.print();

//Test
//      mat L, U;
//      lu(L,U,A);
//      (A-L*U).print("Test of LU decomposition");

//computing time
    time(ii,0)=n;
    time(ii,1)=(float) (clock()-t)/CLOCKS_PER_SEC;

      n=n*s;
     }

//    Write time matrix to file
     ofile.open("time_LU.txt");
     ofile << setiosflags(ios::showpoint | ios::uppercase);
     ofile << " n:          time: " << endl;
     for (int i=0;i<=q-1;i++) {
     ofile << setw(15) << setprecision(8) << time(i,0);
     ofile << setw(15) << setprecision(8) << time(i,1) << endl;
     }
  ofile.close();

time.print();
    return 0;
}
