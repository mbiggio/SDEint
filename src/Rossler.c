#include "Rossler.h"


void Update_VF_ROSSLER(double* state, double* f, double* pars)
{
    // State Variables
    double x1 = state[0];
    double x2 = state[1];
    double x3 = state[2];

    // Parameters
    double a     = pars[0];
    double b     = pars[1];
    double c     = pars[2];
    
    f[0] = -x2 - x3;
    f[1] = x1 + a*x2;
    f[2] = b + x3*x1 - c*x3;
}

void Update_Diffusion_ROSSLER(double* state, double** g, double* pars)
{
    double D = pars[3];

    g[0][0] = D;
    g[1][0] = 0.0;
    g[2][0] = 0.0;
}

void Update_Diffusion_Deriv_ROSSLER(double* state, double*** g, double* pars)
{
    g[0][0][0] = 0.0;
    g[0][0][1] = 0.0;
    g[0][0][2] = 0.0;
    g[1][0][0] = 0.0;
    g[1][0][1] = 0.0;
    g[1][0][2] = 0.0;
    g[2][0][0] = 0.0;
    g[2][0][1] = 0.0;
    g[2][0][2] = 0.0;
}

int ParseData_ROSSLER(double* pars, double* Tspan, double* dt_addr, int* print_step_addr, FILE* f_in)
{
    fscanf(f_in, "----------------------\n");
    fscanf(f_in, "-- Model Parameters --\n");
    fscanf(f_in, "----------------------\n");
    if(!fscanf(f_in,"a\t%lf\n",pars)) 
    {
        printf("Errore nel leggere i dati!\n");
        return 0;
    }
    if(!fscanf(f_in,"b\t%lf\n",pars+1))
    {
        printf("Errore nel leggere i dati!\n");
        return 0;
    }
    if(!fscanf(f_in,"c\t%lf\n",pars+2))
    {
        printf("Errore nel leggere i dati!\n");
        return 0;
    }
    if(!fscanf(f_in,"D\t%lf\n",pars+3))
    {
        printf("Errore nel leggere i dati!\n");
        return 0;
    }
    fscanf(f_in, "----------------------------\n");
    fscanf(f_in, "-- Integration Parameters --\n");
    fscanf(f_in, "----------------------------\n");
    if(!fscanf(f_in,"T\t%lf\n",Tspan+1))
    {
        printf("Errore nel leggere i dati!\n");
        return 0;
    }
    if(!fscanf(f_in,"dt\t%lf\n",dt_addr))
    {
        printf("Errore nel leggere i dati!\n");
        return 0;
    }
    if(!fscanf(f_in,"print_step\t%d\n",print_step_addr))
    {
        printf("Errore nel leggere i dati!\n");
        return 0;
    }
    return 1;
}
