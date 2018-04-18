#include <iostream>
#include <fstream>
#include <iomanip>
#include <stdlib.h>
#include <cmath>
#include "time.h"
#include "define_matrices.h"
#include "jacobi.h"

using namespace std;
//using namespace arma;
ofstream ofile;



int main(int argc, char *argv[]) //code to find eigenvalues and eigenvectors of a matrix A in the non-interacting case or interacting case with brute force Jacobi or cyclic Jacobi.
{
    //define variables
    int n;
    //int p_p=0;
    int loop;
    double rho_min;
    double rho_max;
    //double rho_max
    int interacting; //flag to decide if consider interacting or non-interacting case:
                       //matrix A'll be defined as consequence
    double eps; //define tolerance
    int method; //0 bf // 1 cyc
    double omega; //define omega in the harmonic potential

    if (argc<=8){
             cout <<"You forgot something: you have to insert n (int, number of mesh points); loop (int, number of general loop's repetions, to perform one loop insert 2); "
                    "rho_min (double), rho_max (double), interacting (int, 0 if you want non-interacting case, 1 if you want the interacting one),"
                    "eps (double, tolerance), method (int, 0 if you want brute-force Jacobi, 1 if you want cyclic jacobi),"
                    "omega (double, pulsations for harmonic potentials) "<< endl;
             exit(1);
    } else {
        n=atoi(argv[1]);
        loop=atoi(argv[2]);
        rho_min=atof(argv[3]);
        rho_max=atof(argv[4]);
        interacting=atoi(argv[5]);
        eps=atof(argv[6]);
        method=atoi(argv[7]);
        omega=atof(argv[8]);
    }
    if ((interacting==0 || interacting==1) && (method==0 || method==1) && n>0 &&loop>1) {
             cout <<"Your setup is good." << endl;
        } else {
             cout <<"Your setup is wrong: n and loop must be positive and loop>1, while interacting and method should be 0 or 1." << endl;
             exit(1);
    }



    double* repeat = new double[loop]; //vector to save interaction performed during the loop
    double* point_numb = new double[loop]; //vector to save n during the loop
    double* t_elapsed = new double[loop];  //vector to save run-time of the code during the loop
//    double* first_eig = new double[loop]; //vector to save first eigenvalues during the loop
//    double* rho = new double[loop]; //vector to save rho during the loop

    for(int ii = 1; ii < loop; ii++){
        int t=n*ii;
//        rho_max=5*ii;
        int iterations=0;
        double timeused=0;

        //initialize vector of eigenvalues
        double* eig_v = new double[t-1];

        //initialize matrix of our problem
        double ** A = new double*[t-1];
        for(int i = 0; i < t-1; i++){
            A[i] = new double[t-1];
        }

        //initialize matrix of eigenvectors
        double ** R = new double*[t-1];
        for(int i = 0; i < t-1; i++){
            R[i] = new double[t-1];
        }

        //function to fill A depending on the choice of rho_min,rho_max,n,interacting,omega
        define_matrices(rho_min,rho_max,t,interacting, A,R, omega);

        //function which performs Jacobi method brute force if method=0, or cyclic if method=1 and it orders
        //the eigenvalues and eigenvectors
        jacobi(t,A,R, eig_v, eps, method,iterations,timeused);

//       print eigenvalues
        cout <<"Eigenvalues: "<<endl;
        for (int i = 0; i < 3; i++){
            cout << eig_v[i] << endl; //it prints just the first three eigenvalues
        }

        // part of code to answer the question: how many points t we need to have
        // in order to get the lowest three
        // eigenvalues with approximately four leading digits
        // after the decimal point? case of non-interacting

//        if(interacting==0){
//            double first;
//            first = round(eig_v[0]*1e4)/1e4;
//            double second;
//            second = round(eig_v[1]*1e4)/1e4;
//            double third;
//            third = round(eig_v[2]*1e4)/1e4;

//            if(round(first) == 3){
//                if(round(second) == 7){
//                    if(round(third) == 11){
//                        p_p = t;
//                        cout <<"appros"<<first <<endl;
//                        cout <<second <<endl;
//                        cout <<third <<endl;
//                        cout <<p_p <<endl;
//                    }
//                }
//            }
//        }


        repeat[ii-1] = iterations;
        point_numb[ii-1] = t;
        t_elapsed[ii-1] = timeused;
//        first_eig[ii-1]=eig_v[0];
//        rho[ii-1]=rho_max;

// free memory
        for (int i = 0; i < t-1; i++){
            delete[] A[i];
            delete[] R[i];
        }

        delete[] A;
        delete[] R;
        delete[] eig_v;
    }

    // write to file n,run time, iterations for any run in the loop
    ofile.open("iterations_and_runtime.txt");
    ofile << setiosflags(ios::showpoint | ios::uppercase);
    ofile << " numb:       repeat:      time:" << endl;
    for (int i=0; i < loop-1; i++) {
        ofile << setw(15) << setprecision(8) << point_numb[i] << "\t";
        ofile << setw(15) << setprecision(8) << repeat[i] << "\t";
        ofile << setw(15) << setprecision(8) << t_elapsed[i] << endl;
    }

    ofile.close();

     // write to file rho_max,first_eigenvalue for any run in the loop
//    ofile.open("first_eig_on_rho.txt");
//    ofile << setiosflags(ios::showpoint | ios::uppercase);
//    ofile << "rho:              first_eig:" << endl;
//    for (int i=0; i < loop-1; i++) {
//        ofile << setw(15) << setprecision(8) << rho[i] << "\t";
//        ofile << setw(15) << setprecision(8) << first_eig[i] << endl;
//    }

    // write to file rho_max,first_eigenvalue for any run in the loop
//    ofile.open("first_eig_on_n.txt");
//    ofile << setiosflags(ios::showpoint | ios::uppercase);
//    ofile << "rho:              first_eig:" << endl;
//    for (int i=0; i < loop-1; i++) {
//        ofile << setw(15) << setprecision(8) << point_numb[i] << "\t";
//        ofile << setw(15) << setprecision(8) << first_eig[i] << endl;
//    }

//    ofile.close();

// free memory
    delete[] point_numb;
    delete[] repeat;
    delete[] t_elapsed;

//    delete[] rho;
//    delete[] first_eig;

    return 0;

}



