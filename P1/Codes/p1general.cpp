#include <iostream>
#include <fstream>
#include <iomanip>
#include <stdlib.h>
#include <math.h>
#include "time.h"
#include <armadillo>
#include <string.h>

using namespace std;
using namespace arma;

ofstream ofile; /*print to file*/

// functions used
double solution(double x){return 1.0-(1.0-exp(-10))*x-exp(-10*x);}; // analytical solution
double f(double x){return 100*exp(-10*x);}; // right side of poisson equation -u''(x)=f(x)=100*exp(-10*x)

int main(int argc, char* argv[]){
    // we rewrite the poisson equation -u''(x) = f(x)
    // as a set of linear equations A v = r with A matrix (n,n)
    int n; // number of points in the interval (0,1)
    double s; //step from one n to the other
    int q; // number of general loop's repetions
                            // it will repeat the gaussian elimination q times


    double bb; // we assume the same value for all diagonal elements of A
    double cc; // we assume the same value for all upper-diagonal elements of A
    double aa; // we assume the same value for all lower-diagonal elements of A

    if (argc<=6){
             cout <<"You forgot something: you have to n (int, number of starting points); t (double, step from one n to the other); "
                    "q (int, number of general loop's repetions); a(int), b(int), c(int), number in the upper-mid-lower diagonals." << endl;
             exit(1);
    } else {
        n=atoi(argv[1]);
        s=atof(argv[2]);
        q=atoi(argv[3]);
        bb=atof(argv[4]);
        cc=atof(argv[5]);
        aa=atof(argv[6]);

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
    // we use n+1 to have a better reading of the index
    vec a(n+1);	// upper-diagonal elements of Matrix A
    vec b(n+1); // diagonal elements of Matrix A
    vec c(n+1); // lower-diagonal elements A
    vec x(n+2); // indipendent variable of the Poisson equation -u''(x) = f(x)
    vec u(n+2); /* analytical solution */
    vec r(n+1); //vector representing right side of Poisson equation
    vec r_t(n+1); /*vectors with _t represents tilde vectors, as our convention */
    r_t[0]=0;
    vec b_t(n+1);
    b_t[0]=0;
    vec v(n+2); /*numerical solution*/
    v[0]=0;

    double h=1.0/(n+1.0); //number of steps
    int i; //index

    //time calculation
    clock_t t;
    t=clock();

    //analytical solution
    for (i=0; i<=n+1; i++){
        x[i]=i*h;
    //cout << "x"<<i<<":"<<x[i]<< endl; /*print steps*/
    }

    for (i=0; i<n+1; i++){
        u[i]=solution(x[i]);
    //cout << "u"<<i<<":"<<u[i]<< endl; /*print analytical solution*/
    }
    //numerical solution
     for (i=1; i<=n; i++){
        r[i]=h*h*f(x[i]); /*arranging the right side of equation -u''=r(x)=h^2*f(x)*/
    }

    /*Filling a,b,c*/
     /* so I start with index 1 */
    for (i=1;i<=n;i++){
        if(i==1){
//            cout << "b[i]=" <<endl;
//            cin >> bb;
            b[i]=bb;
//            cout << "a[i]=" <<endl;
//            cin >> aa;
            a[i]=aa;
//            cout << "c[i]=" <<endl;
//            cin >> cc;
            c[i]=cc;
        } else if(i==n) {
            b[i]=b[i-1];
            a[i]=0;
            c[i]=0;
        } else {
            b[i]=b[i-1];
            a[i]=a[i-1];
            c[i]=c[i-1];
        }
//       cout << "b"<<i<<":"<<b[i]<< endl;
//       cout << "a"<<i<<":"<<a[i]<< endl;
//       cout << "c"<<i<<":"<<c[i]<< endl; //to check values
    }

    //forward substitution
    /*Declaring vectors*/

    b_t[1]=b[1];
    r_t[1]=r[1];
    for(i=2;i<=n; i++){
        b_t[i]=b[i]-(a[i-1]*a[i-1])/(b_t[i-1]);
        r_t[i]=r[i]-a[i-1]*r_t[i-1]/(b_t[i-1]);
    }
//    for(i=1;i<=n;i++){
//      cout << "b_t"<<i<<":"<<b_t[i]<< endl;
//      cout << "r"<<i<<":"<<r[i]<< endl;
//      cout << "r_t"<<i<<":"<<r_t[i]<< endl; //to check values
//    }

    //backward substitution
    v[n]=r_t[n]/b_t[n];
    for(i=n-1;i>0; i--){
        v[i]=(r_t[i] -a[i]*v[i+1])/(b_t[i]);
    }

//    for(i=0; i<n+1;i++){
//        cout << "v"<<i<<"="<<v[i]<< endl; /*print numerical solution*/
//    }

// Computing time
     time(ii,0)=n;
     time(ii,1)=(float) (clock()-t)/CLOCKS_PER_SEC;

//Incrementing value of n
     n = s*n;

// part of the code to save to file the value of the analytical and numerical solutions, only if you set one loop.
  ofile.open("dat_p1_1000.txt");
  ofile << setiosflags(ios::showpoint | ios::uppercase);
  ofile << " x:              u(x):           v(x): " << endl;
  for (int i=0;i<=n+1;i++) {
  ofile << setw(15) << setprecision(8) << x[i];
  ofile << setw(15) << setprecision(8) << u[i];
  ofile << setw(15) << setprecision(8) << v[i] << endl;
  }
  ofile.close();

}

//    Write time matrix to file
    ofile.open("time_general.txt");
    ofile << setiosflags(ios::showpoint | ios::uppercase);
    ofile << " n:           time: " << endl;
    for (int i=0;i<=q-1;i++) {
    ofile << setw(15) << setprecision(8) << time(i,0);
    ofile << setw(15) << setprecision(8) << time(i,1) << endl;
    }
 ofile.close();

time.print();

return 0;
}


