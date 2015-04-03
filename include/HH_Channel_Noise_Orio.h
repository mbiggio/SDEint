#ifndef HH_CHANNEL_NOISE_ORIO_H
#define HH_CHANNEL_NOISE_ORIO_H
#include <stdio.h>
#include <math.h>
#include <stdlib.h>

void Update_VF_HH_Channel_Noise_Orio(double* state, double* f, double* pars);
void Update_Diffusion_HH_Channel_Noise_Orio(double* state, double** g, double* pars);
void Update_Diffusion_Deriv_HH_Channel_Noise_Orio(double* state, double*** g, double* pars);
int ParseData_HH_Channel_Noise_Orio(double* pars, double* Tspan, char* outfile, double* dt_addr, int* print_step_addr, FILE* f_in);

#endif
