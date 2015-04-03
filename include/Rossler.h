#ifndef ROSSLER_H
#define ROSSLER_H
#include <stdio.h>

void Update_VF_ROSSLER(double* state, double* f, double* pars);
void Update_Diffusion_ROSSLER(double* state, double** g, double* pars);
void Update_Diffusion_Deriv_ROSSLER(double* state, double*** g, double* pars);
int ParseData_ROSSLER(double* pars, double* Tspan, double* dt_addr, int* print_step_addr, FILE* f_in);

#endif
