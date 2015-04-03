#include "HR.h"


void Update_VF_HR(double* state, double* f, double* pars)
{
    // State Variables
    double x = state[0];
    double y = state[1];
    double z = state[2];

    // Parameters
    double b     = pars[0];
    double I     = pars[1];
    double mu    = pars[2];
    double s     = pars[3];
    double xrest = pars[4];
    
    f[0] = y - x*x*x + b*x*x + I - z;
    f[1] = 1.0 - 5.0*x*x - y ;
    f[2] = mu*(s*(x-xrest) - z);
}

void Update_Diffusion_HR(double* state, double** g, double* pars)
{
    double D = pars[5];

    g[0][0] = D;
    g[1][0] = 0.0;
    g[2][0] = 0.0;
}

void Update_Diffusion_Deriv_HR(double* state, double*** g, double* pars)
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

int ParseData_HR(double* pars, double* Tspan, double* dt_addr, int* print_step_addr, FILE* f_in)
{
    fscanf(f_in, "----------------------\n");
    fscanf(f_in, "-- Model Parameters --\n");
    fscanf(f_in, "----------------------\n");
    if(!fscanf(f_in,"b\t%lf\n",pars)) 
    {
        printf("Errore nel leggere i dati!\n");
        return 0;
    }
    if(!fscanf(f_in,"I\t%lf\n",pars+1))
    {
        printf("Errore nel leggere i dati!\n");
        return 0;
    }
    if(!fscanf(f_in,"mu\t%lf\n",pars+2))
    {
        printf("Errore nel leggere i dati!\n");
        return 0;
    }
    if(!fscanf(f_in,"s\t%lf\n",pars+3))
    {
        printf("Errore nel leggere i dati!\n");
        return 0;
    }
    if(!fscanf(f_in,"xrest\t%lf\n",pars+4))
    {
        printf("Errore nel leggere i dati!\n");
        return 0;
    }
    if(!fscanf(f_in,"D\t%lf\n",pars+5))
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
