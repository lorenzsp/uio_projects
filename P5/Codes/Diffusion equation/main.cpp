#include <iostream>
#include <cmath>
#include "pde_diffusion.h"
#include "data_analysis.h"
#include <string.h>
using namespace std;


int main()
{
    // define general parameters valid for both
    // one and two dimensions
    // number of points in the spatial coordinate(s)
    //  x: 0 -> 1       y: 0 -> 1
    int N = 100;
    // time step size
    double dt=2e-4;
    // tolerance for the iterative two dimensional solver
    double tolerance=1e-10;
    // maximum iteration that we want in the iterative two dimensional solver
    int cutoff=1e5;

    // class used to print data and print to screen
    data_analysis* dat = new data_analysis;

    dat->open_file("test.txt");

    // loop over all the interested spatial step size
    for(N=1000; N<=1000; N*=10){
        // step size of the spatial coordinates
        double h = 1.0/N;

        // loop over all the interested temporal step size
        for(dt = 1e-3; dt<=1e-3; dt*=5){

            // summerize the information
            double alpha = dt/(h*h);
            dat->printscreen(dt, "dt ");
            dat->printscreen(h, "h ");
            dat->printscreen(alpha, "alpha");

            //initial conditions for one dimensional
            double* I = new double[N+1];
            double* E = new double[N+1];
            double* C = new double[N+1];
            double* s = new double[N+1];
            for(int i=0;i<N;i++){
                I[i] = 0;
                E[i] = 0;
                C[i] = 0;
                s[i] = 0;
            }
            I[N] = 1;
            E[N] = 1;
            C[N] = 1;
            s[N] = 1;

            // intial conditions for two dimensional case
            double** u_2d = new double*[N+1];
            double** s_2d = new double*[N+1];
            for(int i=0;i<N+1;i++){
                u_2d[i]=new double[N+1];
                s_2d[i]=new double[N+1];
            }

            for(int i=0; i<N+1;i++){
                for(int j=0; j<N+1;j++){
                    u_2d[i][j] = sin(j*h*M_PI)*sin(i*h*M_PI);
                }
            }

            // declaration of the class used to solve PDE diffusion equation
            pde_diffusion* solver = new pde_diffusion(N, alpha);

            // final time of the simulation
            double t_final = 1;
            // time declaration
            double time=0;
            for(double time=dt; time<=t_final ; time +=dt){
                // one dimension
                solver->Implicit(I);
                solver->Explicit(E);
                solver->Solution1d(s, time);
                solver->Crank_nicolson(C);
                // two dimensions
                //solver->two_dimension(u_2d, tolerance, cutoff);
                //solver->Solution2d( s_2d, time);

            }

            for(int l=0; l<N+1; l++){
                dat->write(I[l]);
                dat->add_column();

            }
            dat->new_row();

            for(int l=0; l<N+1; l++){
                dat->write(C[l]);
                dat->add_column();

            }

            dat->new_row();

            delete[] I;
            delete[] s;
            delete[] E;
            delete[] C;

            delete solver;
            for(int i=0; i<N+1; i++) {
                delete[] u_2d[i];
                delete[] s_2d[i];
            }
            delete[] u_2d;
            delete[] s_2d;


        } // dt loop

    } // N loop

    dat->close_file();
    cout << '\a';

    delete dat;
    return 0;
}
