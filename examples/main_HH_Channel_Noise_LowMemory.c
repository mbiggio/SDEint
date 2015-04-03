#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <math.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_cblas.h>
#include "Solver.h"
#include "HH_Channel_Noise.h"


int main(int argc, char* argv[])
{

    FILE* f_in;
    f_in = fopen("input_HH_CN.dat","r");

    double pars[10];
    double Tspan[2];
    Tspan[0] = 0.0;
    double dt;
    int print_step;
    char outfile[50];
    double th = -20.0;

    if(!ParseData_HH_Channel_Noise(pars,Tspan,outfile,&dt,&print_step,f_in)) return 1;
    fclose(f_in);

    FILE* f_out;
    f_out = fopen(outfile,"w");

    double x0[15] = {0.1, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};

    // Allocation and Initialization of Random Number Generator
    const gsl_rng_type * T;
    gsl_rng * r;
    gsl_rng_env_setup();
    T = gsl_rng_default;
    r = gsl_rng_alloc(T);
    gsl_rng_set(r,1);

    #ifdef LYAP
        // Parameters for the Calculation of the Maximum LE
        gsl_rng * r_lyap;
        r_lyap = gsl_rng_alloc(T);
        gsl_rng_set(r_lyap,1);
        double tau = 10*dt;
        double alpha = 1e-8;
        double lyapexp = 0;
        SdeSolverLyapLog(Stepper_EM,Update_VF_HH_Channel_Noise,Update_Diffusion_HH_Channel_Noise,Update_Diffusion_Deriv_HH_Channel_Noise,Tspan,x0,15,12,dt,pars,r,r_lyap,tau,alpha,&lyapexp,f_out,print_step);
        printf("Computed value of maximum Lyapunov Exponent: %f\n",lyapexp);
    #else
        SdeSolverSpikeDetection(Stepper_EM,Update_VF_HH_Channel_Noise,Update_Diffusion_HH_Channel_Noise,Update_Diffusion_Deriv_HH_Channel_Noise,Tspan,x0,15,12,dt,pars,r,f_out,th,1);
    #endif


    fclose(f_out);
    return 0;
}
