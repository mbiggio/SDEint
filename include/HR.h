#ifndef HR_H
#define HR_H
#include <stdio.h>

void Update_VF_HR(double* state, double* f, double* pars);
void Update_Diffusion_HR(double* state, double** g, double* pars);
void Update_Diffusion_Deriv_HR(double* state, double*** g, double* pars);
int ParseData_HR(double* pars, double* Tspan, double* dt_addr, int* print_step_addr, FILE* f_in);

#endif
