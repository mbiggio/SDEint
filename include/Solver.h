#ifndef SOLVER_H
#define SOLVER_H
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_cblas.h>
#include <math.h>

#define GREEN "\x1B[32m"
#define RED "\x1B[31m"
#define NORMAL "\x1B[0m"


// These typedefs are simply nicknames for the whole function headers
typedef void (*VF)(double* state, double* f, double* pars);
typedef void (*DIFF)(double* state, double** g, double* pars);
typedef void (*DIFF_DER)(double* state, double*** g_der, double* pars);
typedef void (*STEPPER)(VF Update_VF, DIFF Update_Diffusion, DIFF_DER Update_Diffusion_Deriv, double* new_state, double* old_state, double* pars, double dt, gsl_rng* r, int dim, int dim_noise, double* VF, double** g, double*** gd);

// Fourth-Order Stochastic Runge-Kutta Stepper (PRE 70, 017701, 2014)
void Stepper_RK4(VF Update_VF, DIFF Update_Diffusion, DIFF_DER Update_Diffusion_Deriv, double* new_state, double* old_state, double* pars, double dt, gsl_rng* r, int dim, int dim_noise, double* VF, double** g, double*** gd);

// Euler-Maruyama Scheme
void Stepper_EM(VF Update_VF, DIFF Update_Diffusion, DIFF_DER Update_Diffusion_Deriv, double* new_state, double* old_state, double* pars, double dt, gsl_rng* r, int dim, int dim_noise, double* VF, double** g, double*** gd);

// Auxiliary Function for RK4
void Update_f(VF Update_VF, DIFF Update_Diffusion, DIFF_DER Update_Diffusion_Deriv, double* state, double* f, double* pars, double dt, int dim, int dim_noise, gsl_rng* r, double* VF, double** g, double*** gd);

// Generic SDE Solver
void SdeSolver(STEPPER StepFunc, VF Update_VF, DIFF Update_Diffusion, DIFF_DER Update_Diffusion_Deriv, double** x, double* t, double* Tspan, double* x0, int dim, int dim_noise,double dt, double* pars, gsl_rng* r);

// Generic SDE Solver, Computes Maximum Lyapunov Exponent
void SdeSolverLyap(STEPPER StepFunc, VF Update_VF, DIFF Update_Diffusion, DIFF_DER Update_Diffusion_Deriv, double** x, double* t, double* Tspan, double* x0, int dim, int dim_noise, double dt, double* pars, gsl_rng* r, gsl_rng* r_lyap, double tau, double alpha, double* lyapexp, int print);

// Generic SDE Solver, directly logs the solution to file
void SdeSolverLog(STEPPER StepFunc, VF Update_VF, DIFF Update_Diffusion, DIFF_DER Update_Diffusion_Deriv, double* Tspan, double* x0, int dim, int dim_noise, double dt, double* pars, gsl_rng* r, FILE* f_out, int print_steps);

// Generic SDE Solver, directly logs the solution to file
void SdeSolverLyapLog(STEPPER StepFunc, VF Update_VF, DIFF Update_Diffusion, DIFF_DER Update_Diffusion_Deriv, double* Tspan, double* x0, int dim, int dim_noise, double dt, double* pars, gsl_rng* r, gsl_rng* r_lyap, double tau, double alpha, double* lyapexp, FILE* f_out, int print_steps, int print);

// Generic SDE Solver, only compute LE
void SdeSolverLyap_NoLog(STEPPER StepFunc, VF Update_VF, DIFF Update_Diffusion, DIFF_DER Update_Diffusion_Deriv, double* Tspan, double* x0, int dim, int dim_noise, double dt, double* pars, gsl_rng* r, gsl_rng* r_lyap, double tau, double alpha, double* lyapexp, int print);

// Generic SDE Solver, detects spikes and logs spike times to file
void SdeSolverSpikeDetection(STEPPER StepFunc, VF Update_VF, DIFF Update_Diffusion, DIFF_DER Update_Diffusion_Deriv, double* Tspan, double* x0, int dim, int dim_noise, double dt, double* pars, gsl_rng* r, FILE* f_out, double th, int print);

// Output Results to File
void WriteResults(FILE* f_out, double**x, double* t, int steps, int dim, int print_steps);

// Auxiliary Functions for Lyapunov Exponents Computation
double Distance(double* x, double* y, int dim);
double Length(double* x, int dim);

#endif
