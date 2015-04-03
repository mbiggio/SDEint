#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <math.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_cblas.h>
#include "Solver.h"
#include "Rossler.h"

int main(int argc, char* argv[])
{
    int k;

    FILE* f_out;
    f_out = fopen("output_ROSSLER.dat","w");

    FILE* f_in;
    f_in = fopen("input_ROSSLER.dat","r");

    double pars[4];
    double Tspan[2];
    Tspan[0] = 0.0;
    double dt;
    int print_step;

    if(!ParseData_ROSSLER(pars,Tspan,&dt,&print_step,f_in)) return 1;
    fclose(f_in);

    int steps;
    steps = ceil(Tspan[1]/dt);
    double x0[3] = {1.0, 0.0, 0.0};

    // Allocation of the state variables
    double** x = (double**)malloc(steps*sizeof(double*));
    for (k = 0; k < steps; k++) x[k] = malloc(3*sizeof(double));
    double* t = (double*)malloc(steps*sizeof(double));

    // Allocation and Initialization of Random Number Generator
    const gsl_rng_type * T;
    gsl_rng * r;
    gsl_rng * r_lyap;
    gsl_rng_env_setup();
    T = gsl_rng_default;
    r = gsl_rng_alloc(T);
    r_lyap = gsl_rng_alloc(T);

    // Parameters for the Calculation of the Maximum LE
    double tau = 100*dt;
    double alpha = 1e-2;
    double lyapexp = 0;
      
    SdeSolverLyap(Stepper_RK4,Update_VF_ROSSLER,Update_Diffusion_ROSSLER,Update_Diffusion_Deriv_ROSSLER,x,t,Tspan,x0,3,1,dt,pars,r,r_lyap,tau,alpha,&lyapexp,1);
    WriteResults(f_out,x,t,steps,3,print_step);

    printf("Computed value of maximum Lyapunov Exponent: %f\n",lyapexp);

    // Deallocation of the state variables
    for (k = 0; k < steps; k++) free(x[k]);
    free(x);
    fclose(f_out);
    return 0;
}
