#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <math.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_cblas.h>
#include "Solver.h"
#include "HR.h"

int main(int argc, char* argv[])
{
    FILE* f_out;
    f_out = fopen("output_HR.dat","w");

    FILE* f_in;
    f_in = fopen("input_HR.dat","r");

    double pars[6];
    double Tspan[2];
    Tspan[0] = 0.0;
    double dt;
    int print_step;

    if(!ParseData_HR(pars,Tspan,&dt,&print_step,f_in)) return 1;
    fclose(f_in);

    double x0[3] = {1.0, 0.0, 0.0};

    // Allocation and Initialization of Random Number Generator
    gsl_rng_env_setup();
    gsl_rng * r;
    gsl_rng * r_lyap;
    r = gsl_rng_alloc(gsl_rng_default);
    r_lyap = gsl_rng_alloc(gsl_rng_default);

    // Seed of the Random Number Generators
    gsl_rng_set(r,1);
    gsl_rng_set(r_lyap,1);

    // Parameters for the Calculation of the Maximum LE
    double tau = 100*dt;
    double alpha = 1e-6;
    double lyapexp = 0;
      
    SdeSolverLyapLog(Stepper_RK4,Update_VF_HR,Update_Diffusion_HR,Update_Diffusion_Deriv_HR,Tspan,x0,3,1,dt,pars,r,r_lyap,tau,alpha,&lyapexp,f_out,print_step,1);
    // SdeSolverLyapLog(Stepper_EM,Update_VF_HR,Update_Diffusion_HR,Update_Diffusion_Deriv_HR,Tspan,x0,3,1,dt,pars,r,r_lyap,tau,alpha,&lyapexp,f_out,print_step);

    printf("Computed value of maximum Lyapunov Exponent: %0.12f\n",lyapexp);

    fclose(f_out);
    return 0;
}
