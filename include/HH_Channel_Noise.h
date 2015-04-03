#ifndef HH_CHANNEL_NOISE_H
#define HH_CHANNEL_NOISE_H
#include <stdio.h>
#include <math.h>
#include <stdlib.h>

void Update_VF_HH_Channel_Noise(double* state, double* f, double* pars);
void Update_Diffusion_HH_Channel_Noise(double* state, double** g, double* pars);
void Update_Diffusion_Deriv_HH_Channel_Noise(double* state, double*** g, double* pars);
int ParseData_HH_Channel_Noise(double* pars, double* Tspan, char* outfile, double* dt_addr, int* print_step_addr, FILE* f_in);
int ParseData_HH_Channel_Noise_MT(double* pars, double* Tspan, double* dt_addr, int* print_step_addr, int* seed_addr, FILE* f_in);

#endif
