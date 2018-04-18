#ifndef LEGENDRE_PKG_H
#define LEGENDRE_PKG_H
#include <iostream>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <stdlib.h>


class legendre_pkg
{
public:
    legendre_pkg(int N, long double a[N], long double y[N], long double w[N]);
    int n, m;
    long double fact(int N);
    void poly_coeff(long double *a);
    long double pn(long double* a,long double x);
    long double dn(long double* a,long double x);
};

#endif // LEGENDRE_PKG_H
