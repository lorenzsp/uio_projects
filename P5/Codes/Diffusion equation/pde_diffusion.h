#ifndef PDE_DIFFUSION_H
#define PDE_DIFFUSION_H


class pde_diffusion
{
public:
    int N;
    double alpha;
    double r_cn, r_i, r_e, r_2d;
    void Crank_nicolson(double* u);
    pde_diffusion(int n, double a);
    void Implicit(double *u);
    void Explicit(double *u);
    void two_dimension(double **u, double tolerance, int cutoff);
    void Solution1d(double *solution, double time);
    void Solution2d(double **solution, double time);
    double max_relative_error(double *u, double *s);
    double max_relative_error_2d(double **u, double **s);
};

#endif // PDE_DIFFUSION_H
