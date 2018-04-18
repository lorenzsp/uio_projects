#include "legendre_pkg.h"
#include <string>
#include <ios>
#include <iostream>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <stdlib.h>
using namespace std;



//Factorial Function
legendre_pkg::legendre_pkg(int N, long double a[N], long double y[N], long double w[N])
{
    n= N;
    m=n%2;
    poly_coeff(a);
    //Roots of Pn(x) with newton method
    long double z[n];
    long double l, v, s;

    for(int i=0; i<n; i++){
        z[i]=cos(M_PI*(i+0.75)/(n+0.5));
        l=z[i];
        do
        {
            s=l-(pn(a,l)/dn(a,l));
            v=l;
            l=s;
        }
        while(fabs(l-v)>0.0000000000000001);

        y[i]=l;
        w[i]=2/((1-pow(l,2))*(dn(a,l)*dn(a,l)));
    }
}

//For Legendre's Polynomial Pn(x)
long double legendre_pkg::pn(long double* a,long double x)
{

    int i;
    long double p=0;
    if(m==0)
    {
        for(i=0;i<= n ;i=i+2)
        {
            if(x==0)
                break;
            p+=a[i]*pow(x,i);
        }
    }
    else
    {
        for(i=1;i<=n;i=i+2)
        {
            p+=a[i]*pow(x,i);
        }
    }
    return p;
}

//Derivative of Pn(x)
long double legendre_pkg::dn(long double* a,long double x)
{
    int i;
    long double p=0;
    if(m==0)
    {
        for(i=0;i<=n;i=i+2)
        {
            if(x==0)
                break;
            p+=i*a[i]*pow(x,i-1);
        }
    }
    else
    {
        for(i=1;i<=n;i=i+2)
        {
            p+=i*a[i]*pow(x,i-1);
        }
    }
    return p;
}



long double legendre_pkg::fact(int N)
{
    int i;
    long double f=1;
    for(i=2;i<=N;i++)
    {
        f*=i;
    }
    return f;
}


// polynomial coefficients
void legendre_pkg::poly_coeff(long double *a)
{
    double N;
    if(m==0)
    {
        N=n/2;
    }
    else
    {
        N=(n-1)/2;
    }

    for(int i=0;i<=N;i++)
    {
        a[n-2*i]=(pow(-1,i)*fact(2*n-2*i))/(pow(2,n)*fact(i)*fact(n-i)*fact(n-2*i));
    }
}


