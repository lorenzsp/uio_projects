#include <iostream>
#include <fstream>
#include <iomanip>
#include <stdlib.h>
#include <cmath> //libreria migliore
#include "time.h"
#include <armadillo>

using namespace std;
using namespace arma;


ofstream ofile; /*print to file*/

// functions used
double solution(double x){return 1.0-(1.0-exp(-10))*x-exp(-10*x);};
double f(double x){return 100*exp(-10*x);};


int main(int argc, char* argv[]){
    // we rewrite the poisson equation -u''(x) = f(x)
    // as a set of linear equations A v = r with A matrix (n,n)
    int n; // number of points in the interval (0,1)
    double s; //step from one n to the other
    int q; // number of general loop's repetions
                            // it will repeat the gaussian elimination q times


    if (argc<=3){
             cout <<"You forgot something: you have to n (int, number of starting points); t (double, step from one n to the other); "
                    "q (int, number of general loop's repetions)"<< endl;
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


    //mat time(q,2); // time matrix: it's useful to have a file with elapsed time. n, time
    mat max_err(q,3); // matrix about maximum relative error. n,log10(h),max relative error

    //You have to choose a matrix before starting. You can choose no matrix, but you have to comment the last step about writing matrices to file.
    for(int ii=0;ii<=q-1;ii++){


    //Declaration and allocation of variables with Armadillo
    // we use n+1 to have a better reading of the index
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

    // time calculation
//    clock_t t;
//    t=clock();

    //analytical solution
    for (i=0; i<=n+1; i++){
        x[i]=i*h;
        //cout << "x"<<i<<":"<<x[i]<< endl; /*print steps*/
    }

    for (i=0;i<n+1; i++){
        u[i]=solution(x[i]);
        //cout << "u"<<i<<":"<<u[i]<< endl;  /*print analytical solution*/
    }

    //numerical solution
     for (i=1; i<=n; i++){
        r[i]=h*h*f(x[i]); /*arranging the right side of equation -u''=r(x)=h^2*f(x)*/
    }

    //Forward substitution
    /*Declaring vectors*/
    r_t[1]=r[1];
    for(i=2;i<=n; i++){
        b_t[i]=(i+1.0)/i;
        r_t[i]=r[i]+(i-1.0)*r_t[i-1]/(i);
    }
    //for(i=1;i<=n;i++){
//      cout << "b_t"<<i<<":"<<b_t[i]<< endl;
//      cout << "r"<<i<<":"<<r[i]<< endl;
//      cout << "r_t"<<i<<":"<<r_t[i]<< endl; //to check values
    //}

    //backward substitution
    v[n]=r_t[n]/b_t[n];
    for(i=n-1;i>0; i--){
                v[i]=(r_t[i]+v[i+1])*i/(i+1);
    }
//    for(i=0; i<=n+1;i++){
//       cout << "v"<<i<<"="<<v[i]<< endl; /*print numerical solution*/
//    }


// Computing time
//     time(ii,0)=n;
//     time(ii,1)=(float) (clock()-t)/CLOCKS_PER_SEC;

//    computing relative errors
    double *err=new double [n];
    //double *err = new double [n];
    for(int k=0;k<n;k++){
    err[k]= log10(abs((v[k+1]-u[k+1])/(u[k+1])));
    }
     //err.print();  //to check values

//    computing maximum relative error
    double MAX=err[0];
        for(int k=1;k<n;k++){
            if(err[k]>MAX)
                MAX=err[k];
        }
    max_err(ii,0)=n;
    max_err(ii,1)=log10(h);
    max_err(ii,2)=MAX;

    delete [] err;

//Incrementing value of n
  n = s*n;
  }

// part of the code to save to file the value of the analytical and numerical solutions, only if you set one loop.
//     ofile.open("dat_p1.txt");
//     ofile << setiosflags(ios::showpoint | ios::uppercase);
//     ofile << " x: u(x): v(x): " << endl;
//     for (int i=0;i<=n+1;i++) {
//     ofile << setw(15) << setprecision(8) << x[i];
//     ofile << setw(15) << setprecision(8) << u[i];
//     ofile << setw(15) << setprecision(8) << v[i] << endl;
//     }
//  ofile.close();

//    Write matrix about maximum relative error to file
     ofile.open("data_error.txt");
                ofile << setiosflags(ios::showpoint | ios::uppercase);
               ofile << " n:          log10(h):          log10(\varepsilon) " << endl;
               for (int i=0;i<=q-1;i++) {
                 ofile << setw(15) << setprecision(8) << max_err(i,0);
                 ofile << setw(15) << setprecision(8) << max_err(i,1);
                 ofile << setw(15) << setprecision(8) << max_err(i,2) << endl;
                 }
              ofile.close();
    max_err.print(); //print matrix to check

//    Write time matrix to file
    /*ofile.open("time_particular.txt");
      ofile << setiosflags(ios::showpoint | ios::uppercase);
      ofile << " n:          time: " << endl;
      for (int i=0;i<=q-1;i++) {
      ofile << setw(15) << setprecision(8) << time(i,0);
      ofile << setw(15) << setprecision(8) << time(i,1) << endl;
      }
      ofile.close();

time.print();*/ //print matrix to check

return 0;
}

