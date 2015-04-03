#include "Leech.h"

double sigmoid(double V, double alpha, double beta)
{
    return 1.0/(1 + exp(alpha*(V + beta)));
}

void Update_VF_LEECH(double* state, double* f, double* pars)
{
    // State Variables
    double V = state[0];
    double h = state[1];
    double m = state[2];

    // Parameters
    double C      = pars[ 0];
    double Iapp   = pars[ 1];
    double gNa    = pars[ 2];
    double gK2    = pars[ 3];
    double gL     = pars[ 4];
    double ENa    = pars[ 5];
    double EK     = pars[ 6];
    double EL     = pars[ 7];
    double tauNa  = pars[ 8];
    double tauK2  = pars[ 9];
    double Vshift = pars[10];
    
    double ninf = sigmoid(V,-150.0,0.0305);
    double hinf = sigmoid(V, 500.0,0.0333);
    double minf = sigmoid(V,-83.0,0.018 + Vshift);

    double INa = gNa *ninf*ninf*ninf *h* (V - ENa);
    double IK2 = gK2 *m*m*               (V - EK );
    double IL  = gL  *                   (V - EL );

    f[0] = -(INa + IK2 + IL - Iapp)/C;
    f[1] =  (hinf - h)/tauNa;
    f[2] =  (minf - m)/tauK2;
}

void Update_Diffusion_LEECH(double* state, double** g, double* pars)
{
    double D = pars[11];

    g[0][0] = D;
    g[1][0] = 0.0;
    g[2][0] = 0.0;
}

void Update_Diffusion_Deriv_LEECH(double* state, double*** g, double* pars)
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

int ParseData_LEECH(double* pars, double* Tspan, double* dt_addr, int* print_step_addr, int* seed_addr, FILE* f_in)
{
    fscanf(f_in, "----------------------\n");
    fscanf(f_in, "-- Model Parameters --\n");
    fscanf(f_in, "----------------------\n");
    if(!fscanf(f_in,"C\t%lf\n",pars)) 
    {
        printf("Errore nel leggere i dati!\n");
        return 0;
    }
    if(!fscanf(f_in,"Iapp\t%lf\n",pars+1))
    {
        printf("Errore nel leggere i dati!\n");
        return 0;
    }
    if(!fscanf(f_in,"gNa\t%lf\n",pars+2))
    {
        printf("Errore nel leggere i dati!\n");
        return 0;
    }
    if(!fscanf(f_in,"gK2\t%lf\n",pars+3))
    {
        printf("Errore nel leggere i dati!\n");
        return 0;
    }
    if(!fscanf(f_in,"gL\t%lf\n",pars+4))
    {
        printf("Errore nel leggere i dati!\n");
        return 0;
    }
    if(!fscanf(f_in,"ENa\t%lf\n",pars+5))
    {
        printf("Errore nel leggere i dati!\n");
        return 0;
    }
    if(!fscanf(f_in,"EK\t%lf\n",pars+6))
    {
        printf("Errore nel leggere i dati!\n");
        return 0;
    }
    if(!fscanf(f_in,"EL\t%lf\n",pars+7))
    {
        printf("Errore nel leggere i dati!\n");
        return 0;
    }
    if(!fscanf(f_in,"tauNa\t%lf\n",pars+8))
    {
        printf("Errore nel leggere i dati!\n");
        return 0;
    }
    if(!fscanf(f_in,"tauK2\t%lf\n",pars+9))
    {
        printf("Errore nel leggere i dati!\n");
        return 0;
    }
    if(!fscanf(f_in,"Vshift\t%lf\n",pars+10))
    {
        printf("Errore nel leggere i dati!\n");
        return 0;
    }
    if(!fscanf(f_in,"D\t%lf\n",pars+11))
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
    if(!fscanf(f_in,"seed\t%d\n",seed_addr))
    {
        printf("Errore nel leggere i dati!\n");
        return 0;
    }
    return 1;
}
