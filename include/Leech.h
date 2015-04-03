#ifndef LEECH_H
#define LEECH_H
#include <stdio.h>
#include <math.h>

void Update_VF_LEECH(double* state, double* f, double* pars);
void Update_Diffusion_LEECH(double* state, double** g, double* pars);
void Update_Diffusion_Deriv_LEECH(double* state, double*** g, double* pars);
int ParseData_LEECH(double* pars, double* Tspan, double* dt_addr, int* print_step_addr, int* seed_addr, FILE* f_in);

#endif
