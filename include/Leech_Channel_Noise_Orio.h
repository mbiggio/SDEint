#ifndef LEECH_CHANNEL_NOISE_ORIO_H
#define LEECH_CHANNEL_NOISE_ORIO_H
#include <stdio.h>
#include <math.h>
#include <stdlib.h>

void Update_VF_Leech_Channel_Noise_Orio(double* state, double* f, double* pars);
void Update_Diffusion_Leech_Channel_Noise_Orio(double* state, double** g, double* pars);
void Update_Diffusion_Deriv_Leech_Channel_Noise_Orio(double* state, double*** g, double* pars);
int ParseData_Leech_Channel_Noise_Orio(double* pars, double* Tspan, char* outfile, double* dt_addr, int* print_step_addr, FILE* f_in);
int ParseData_Leech_Channel_Noise_Orio_MT(double* pars, double* Tspan, double* dt_addr, int* print_step_addr, int* seed_addr, FILE* f_in);

#endif
