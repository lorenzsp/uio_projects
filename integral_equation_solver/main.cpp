#include<iostream>  //Header file for cin & cout
#include<cmath>  //Header file for mathematical operartions
#include<iomanip> //Header file for precession
#include<math.h>
#include <stdlib.h>     /* abs */
#include"legendre_pkg.h"
#include"data_analysis.h"
using namespace std;  //calling the standard directory

//Given Function of Integration
long double f(long double x, long double mu, long double T)
{
    long double d;
    d=pow(x,0.5)/(1+exp((x-mu)*T));
    return d;
}



//Main Function
// solver of an integral equation
int main()
{
    int n,i;
    double c,d;
    n=18;
    c=-1;
    d=1;
    /*
    cout<<"Enter the value of n for Pn(x) : \n";
    cin>>n;
    cout<<"Enter the lower limit a of integration : \n";
    cin>>c;
    cout<<"Enter the upper limit b of integration : \n";
    cin>>d;
    */

    long double a[n],y[n],w[n],g=0;

    legendre_pkg integral(n, a, y, w);

    // integral from a to b --> 0 to infty
    // new mesh points
    double* X=new double[n];
    for(int i=0; i<n; i++){
        X[i] =tan(M_PI* (1+ y[i])/4);
    }
    // new weights
    double* M=new double[n];
    for(int i=0; i<n; i++){
        M[i] = (M_PI/4)*w[i]/(cos(M_PI*(1+y[i])/4)*cos(M_PI*(1+y[i])/4));
    }

    long double mu=0;

    data_analysis* dat = new data_analysis;
    dat->open_file("temperature_mu");
    // ratio T/T_F
    long double ratio=0;
    // energy
    long double E =0;
    for(ratio=0.001; ratio<1.6; ratio=ratio+0.0001){
        for( mu=-2; mu<2; mu=mu+0.001){

            g=0;
            E=0;
            //final integration
            for(i=0;i<n;i++){
                g+=M[i]*f(X[i],mu,1/ratio);
                E+=g*X[i] ;
            }

            if(abs(g-0.666666)/(0.666666)<0.001){
                dat->write(ratio);
                dat->add_column();
                dat->write(mu);
                dat->add_column();
                dat->write(E);
                dat->new_row();
                //dat->printscreen(g, "integral value ");
                //dat->printscreen(abs(g-0.666666)/(0.666666), "relative error");
                mu=30;
            }

        }

    }
    dat->close_file();


    // heat capacity

    //g=g*(d-c)/2;
    cout<<"The Value of Integration is = "<<setprecision(10)<<g<<endl;
    cout<<mu<<endl;
    return 0;
}


/*
// the final roots according to the chosen interval are in u[]
for(i=0;i<n;i++)
{
    u[i]=((d-c)*y[i]/2)+(c+d)/2;
}

cout<<"Roots\t\t\t\t"<<"Weights\n";

for(i=0;i<n;i++)
{
    cout<<setprecision(15)<<y[i]<<"\t\t"<<setprecision(15)<<w[i]<<endl;
}
*/

/*OUTPUT
Enter the value of n for Pn(x) :
6
Enter the lower limit a of integration :
1
Enter the upper limit b of integration :
10

The Legendre's Polynomial is : -0.3125 + (6.5625) X^2 + (-19.6875) X^4 + (14.4375) X^6
Roots    Weights
0.932469514203152  0.17132449237917
0.661209386466265  0.360761573048139
0.238619186083197  0.467913934572691
-0.238619186083197  0.467913934572691
-0.661209386466265  0.360761573048139
-0.932469514203152  0.17132449237917
The Value of Integration is = 2499.75
*/
