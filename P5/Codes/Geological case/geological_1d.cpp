#include <iostream>
#include <fstream>
#include <iomanip>
#include <stdlib.h>
#include <string.h>
#include <array>
#include <cmath>

//Declaration of constant variables
#define N 120
#define time_resol 100000

using namespace std;
using std::array;
ofstream ofile;

//Declaration of functions
void tridiagonal_solver(double, double, double, array<double, N+1> &, array<double, N+1> &);
double heat(double time, double, bool);

int main()
{
    //set variables
    double initial_y=0;
    double final_y=120;
    double step_y=(final_y-initial_y)/(N);;

    //set scaling of variable
    double from_ms_Gy=(365*24*3600*1e9);
    double k=2.5;
    double cp=1000;
    double rho=3.5e3;
    double cost_time=(rho*cp)*1e6/(k*from_ms_Gy);
    double final_time=1/cost_time;
    double initial_time=0;
    double step_time=(final_time-initial_time)/(time_resol);
    double alpha=step_time/(step_y*step_y);

    //set arrays for position and time
    array<double, N+1> y;
    array<double, time_resol+1> t;
    for(int i=0; i<N+1;i++){
        y[i]=i*step_y;
    }

    for(int i=0; i<time_resol+1;i++){
        t[i]=i*step_time;
    }

    //set array to solve Au=v;
    array<double, N+1> u{}; //set arrays to zero to have the first point to be zero
    array<double, N+1> v{};

    //set final point to be 1292, we are scaling the whole solution of -8 to have T(0)=0, T(120)=1292 as boundary conditions
    u[N]=v[N]=1300-8;

    //set variables useful for implicit scheme
    double a,b, c;
    a=c=-alpha;
    b=1+2*alpha;
    bool flag=0; //flag is set to zero to avoid radioactive decay

    // set output if you want to write to file the first case without radioactive decay, which is only helpful in our case to reach a steady state useful for the following case with radioactive decay. In this case, you have to comment the lines 98-124.
    string fileout = "geological_steady_state3.txt";
    ofile.open(fileout);

    //First loop for to reach the steady state (without radioactive decay), which will become our input for the geological case
    for(double time=initial_time;time<=final_time+step_time;time+=step_time){
        double coord_y=initial_y;
        for(int i=0; i<N+1;i++){
            v[i]+=step_time*heat(time,coord_y,flag); //head added with respect to the position in the slab
            coord_y+=step_y;
        }
        v[0] = 0;
        v[N] =1300-8;
        tridiagonal_solver(a,b,c,u,v); //recall to tridiagonal solver from Project 1 to solve the tridiagonal equations

        u[0]=0;
        u[N]=1300-8;

        for(int i=1; i<N+1;i++){
            v[i]=u[i];
        }
        ofile << setiosflags(ios::showpoint | ios::uppercase);
        // ofile << " time:              x:           u(x): " << endl;
        for (int i=0;i<N+1;i++) {
            ofile << setw(15) << setprecision(8) << y[i]; //write positions and solutions to file
            u[i]=u[i]+8; //scale up the solutions as we would want T(0)=8, T(120)=1300
            ofile << setw(15) << setprecision(8) << u[i] << endl;
        }

    }

    ofile.close();

    // set output if you want to write to file the case with radioactive decay.
    fileout = "geological_final.txt";
    ofile.open(fileout);


    flag=1; //flag is set to one to have radioactive decay
    for (int i=0;i<N+1;i++) {
        v[i]=u[i]-8;
    }
    for(double time=initial_time;time<=final_time+step_time;time+=step_time){
        double coord_y=initial_y;
        for(int i=0; i<N+1;i++){
            v[i]+=step_time*heat(time,coord_y,flag); //added heat depending on the position in the slab and also on time
            coord_y+=step_y;
        }
        v[0] = 0;
        v[N] =1300-8;
        tridiagonal_solver(a,b,c,u,v); //recall to tridiagonal solver from Project 1 to solve the tridiagonal equations

        u[0]=0;
        u[N]=1300-8;

        for(int i=1; i<N+1;i++){
            v[i]=u[i];
        }
        ofile << setiosflags(ios::showpoint | ios::uppercase);
        // ofile << " time:              x:           u(x): " << endl;
        for (int i=0;i<N+1;i++) {
            ofile << setw(15) << setprecision(8) << y[i]; //write positions and solutions to file
            u[i]=u[i]+8; //scale up the solutions as we would want T(0)=8, T(120)=1300
            ofile << setw(15) << setprecision(8) << u[i] << endl;
        }
    }

    ofile.close();
    return 0;
}

void tridiagonal_solver(double A, double B, double C, array<double, N+1> &S, array<double, N+1> &Q){
    //triadiagonal equations solver from Project 1 set for the conditions of our problem
    //double *d = new double[n];
    array<double, N+1> d;
    d[0] = 0;
    d[1] = B;
    //double *q_new = new double[n];
    array<double, N+1> q_new;
    q_new[0] = 0;
    q_new[1] = Q[1];


    for(int i=2; i<N-1; i++){
        d[i] = B - A*C / d[i-1];
        q_new[i] = Q[i] - A*q_new[i-1] / d[i-1];
    }

    d[N-1]=B-A*C/d[N-2];
    q_new[N-1] = Q[N-1] - A*q_new[N-2] / d[N-2];
    S[N-1] = (-C*1292 +q_new[N-1]) / d[N-1];
    for(int i=N-2; i>=1; i--){
        S[i] = (q_new[i] - C*S[i+1])/ d[i];
    }

}

double heat(double step_time, double y, bool flag){
    //it's a function used to give the heat with respect to time (if flag=1) and with respect to the position. If flag=0 we have the case without radioactive decay, if flag=1 we consider al the radioactive decay in the last 80 km.
    //set constant
    double from_ms_Gy=(365*24*3600*1e9);
    double k=2.5;
    double cp=1000;
    double rho=3.5e3;
    double cost_time=(rho*cp)*1e6/(k*from_ms_Gy);

    //set arrays for heat
    array<double, 3> Q{1.4, 0.35, 0.05};
    array<double, 3> B;
    for(int i=0;i<=2;i++) B[i]=Q[i]/k;

    //set decay time of the three elements
    double decay_time_U=4.47/(cost_time);
    double decay_time_Th=14.0/(cost_time);
    double decay_time_K=1.25/(cost_time);

    //set the additional heat provided by the elements
    double additional_heat=0.5/k;

    // set constant in the exponential of the radioactive decay with respect to halflives of the elements
    double lambdaU=log(2)/decay_time_U;
    double lambdaTh=log(2)/decay_time_Th;
    double lambdaK=log(2)/decay_time_K;

    //set additional heat per step time
    double additional_heat_U=(additional_heat*0.4)*exp(-step_time*lambdaU);
    double additional_heat_Th=(additional_heat*0.4)*exp(-step_time*lambdaTh);
    double additional_heat_K=(additional_heat*0.2)*exp(-step_time*lambdaK);

    //give heat with respect to the position in the slab and with respect to time for the last 80Km if flag=1 because we are considering the radioactive decay;
    // if flag=0, we don't consider the radioactive decay and we give a constant heat also for the last 80km
    if(y>=40&&y<=120){
        if(flag==1){ //case with radioactive decay
            return B[2]+additional_heat_K+additional_heat_Th+additional_heat_U;
        } else if (flag==0) return B[2];
    }
    else if (y>=20&&y<40) return B[1];
    else if (y>=0&&y<20) return B[0];
}



