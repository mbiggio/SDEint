#include "HH_Channel_Noise.h"

void Update_VF_HH_Channel_Noise(double* state, double* f, double* pars)
{
    // State Variables
    double V = state[0];
    double m = state[1];     
    double h = state[2];
    double n = state[3];

    // Sodium Stochastic Processes
    double xi_Na_1 = state[4];
    double xi_Na_2 = state[5];
    double xi_Na_3 = state[6];
    double xi_Na_4 = state[7];
    double xi_Na_5 = state[8];
    double xi_Na_6 = state[9];
    double xi_Na_7 = state[10];

    // Potassium Stochastic Processes
    double xi_K_1 = state[11];
    double xi_K_2 = state[12];
    double xi_K_3 = state[13];
    double xi_K_4 = state[14];

    // Parameters
    double C      = pars[ 0];
    double Iapp   = pars[ 1];
    double gNa    = pars[ 2];
    double gK2    = pars[ 3];
    double gL     = pars[ 4];
    double ENa    = pars[ 5];
    double EK     = pars[ 6];
    double EL     = pars[ 7];
    
    double alpham = -0.1*(V + 40.0)/(exp(-0.1*(V + 40.0)) - 1.0);
    double betam  = 4.0*exp(-(V + 65.0)/18.0);
    double alphah = 0.07*exp(-(V + 65.0)/20.0);
    double betah  = 1.0/(exp(-0.1*(V + 35.0)) + 1.0);
    double alphan = -0.01*(V + 55.0)/(exp(-0.1*(V + 55.0)) - 1.0);
    double betan  = 0.125*exp(-(V + 65.0)/80.0);

    
    double taum = 1.0/(alpham + betam);
    double tauh = 1.0/(alphah + betah);
    double taun = 1.0/(alphan + betan);

    // Sodium Time Constants
    double tau_Na_1 =     tauh;
    double tau_Na_2 =     taum;
    double tau_Na_3 =     taum/2.0;
    double tau_Na_4 =     taum/3.0;
    double tau_Na_5 =     taum*tauh/(taum + tauh);
    double tau_Na_6 =     taum*tauh/(taum + 2.0*tauh);
    double tau_Na_7 =     taum*tauh/(taum + 3.0*tauh);

    // Potassium Time Constants
    double tau_K_1 =     taun;
    double tau_K_2 =     taun/2.0;
    double tau_K_3 =     taun/3.0;
    double tau_K_4 =     taun/4.0;

    double INa = gNa *fabs(m*m*m*h + xi_Na_1 + xi_Na_2 + xi_Na_3 + xi_Na_4 + xi_Na_5 + xi_Na_6 + xi_Na_7) * (ENa - V);
    double IK2 = gK2 *fabs(n*n*n*n + xi_K_1  + xi_K_2  + xi_K_3  + xi_K_4)                                * (EK  - V);
    //double INa = gNa *(m*m*m*h + xi_Na_1 + xi_Na_2 + xi_Na_3 + xi_Na_4 + xi_Na_5 + xi_Na_6 + xi_Na_7) * (ENa - V);
    //double IK2 = gK2 *(n*n*n*n + xi_K_1  + xi_K_2  + xi_K_3  + xi_K_4)                                * (EK  - V);
    double IL  = gL  *         (EL  - V);

    f[0]  = (INa + IK2 + IL + Iapp)/C;
    f[1]  = alpham*(1.0 - m) - betam*m;
    f[2]  = alphah*(1.0 - h) - betah*h;
    f[3]  = alphan*(1.0 - n) - betan*n;
    f[4]  = -xi_Na_1/tau_Na_1;
    f[5]  = -xi_Na_2/tau_Na_2;
    f[6]  = -xi_Na_3/tau_Na_3;
    f[7]  = -xi_Na_4/tau_Na_4;
    f[8]  = -xi_Na_5/tau_Na_5;
    f[9]  = -xi_Na_6/tau_Na_6;
    f[10] = -xi_Na_7/tau_Na_7;
    f[11] = -xi_K_1/tau_K_1;
    f[12] = -xi_K_2/tau_K_2;
    f[13] = -xi_K_3/tau_K_3;
    f[14] = -xi_K_4/tau_K_4;
}

void Update_Diffusion_HH_Channel_Noise(double* state, double** g, double* pars)
{
    // State Variables
    double V = state[0];
    double m = state[1];
    double h = state[2];
    double n = state[3];

    double gNa         = pars[ 2];
    double gK2         = pars[ 3];
    double length      = pars[ 8];
    double diameter    = pars[ 9];
    double gs          = pars[10];
    //double N_Na        = (M_PI*diameter*length)*gNa/gs * 1e-2;
    //double N_K         = (M_PI*diameter*length)*gK2/gs * 1e-2;
    double N_Na        = (M_PI*diameter*length + 2.0*M_PI*diameter/2.0*diameter/2.0)*gNa/gs * 1e-1;
    double N_K         = (M_PI*diameter*length + 2.0*M_PI*diameter/2.0*diameter/2.0)*gK2/gs * 1e-1;

    double alpham = -0.1*(V + 40.0)/(exp(-0.1*(V + 40.0)) - 1.0);
    double betam  = 4.0*exp(-(V + 65.0)/18.0);
    double alphah = 0.07*exp(-(V + 65.0)/20.0);
    double betah  = 1.0/(exp(-0.1*(V + 35.0)) + 1.0);
    double alphan = -0.01*(V + 55.0)/(exp(-0.1*(V + 55.0)) - 1.0);
    double betan  = 0.125*exp(-(V + 65.0)/80.0);

    
    double taum = 1.0/(alpham + betam);
    double tauh = 1.0/(alphah + betah);
    double taun = 1.0/(alphan + betan);

    // Sodium Variances
    double sigma2_Na_1 =     m*m*m*m*m*m*h*(1.0-h)/N_Na;
    double sigma2_Na_2 = 3.0*m*m*m*m*m*h*h*(1.0-m)/N_Na;
    double sigma2_Na_3 = 3.0*m*m*m*m*h*h*(1.0-m)*(1.0-m)/N_Na;
    double sigma2_Na_4 =     m*m*m*h*h*(1.0-m)*(1.0-m)*(1.0-m)/N_Na;
    double sigma2_Na_5 = 3.0*m*m*m*m*m*h*(1.0-m)*(1.0-h)/N_Na;
    double sigma2_Na_6 = 3.0*m*m*m*m*h*(1.0-m)*(1.0-m)*(1.0-h)/N_Na;
    double sigma2_Na_7 =     m*m*m*h*(1.0-m)*(1.0-m)*(1.0-m)*(1.0-h)/N_Na;

    // Sodium Time Constants
    double tau_Na_1 =     tauh;
    double tau_Na_2 =     taum;
    double tau_Na_3 =     taum/2.0;
    double tau_Na_4 =     taum/3.0;
    double tau_Na_5 =     taum*tauh/(taum + tauh);
    double tau_Na_6 =     taum*tauh/(taum + 2.0*tauh);
    double tau_Na_7 =     taum*tauh/(taum + 3.0*tauh);

    // Potassium Variances
    double sigma2_K_1 = 4.0*n*n*n*n*n*n*n*(1.0-n)/N_K;
    double sigma2_K_2 = 6.0*n*n*n*n*n*n*(1.0-n)*(1.0-n)/N_K;
    double sigma2_K_3 = 4.0*n*n*n*n*n*(1.0-n)*(1.0-n)*(1.0-n)/N_K;
    double sigma2_K_4 =     n*n*n*n*(1.0-n)*(1.0-n)*(1.0-n)*(1.0-n)/N_K;

    // Potassium Time Constants
    double tau_K_1 =     taun;
    double tau_K_2 =     taun/2.0;
    double tau_K_3 =     taun/3.0;
    double tau_K_4 =     taun/4.0;

    // Synaptic Noise
    double D = pars[11];

    // Initialization
    int k,s;
    for(k = 0; k < 14; k++)
        for(s = 0; s < 11; s++) g[k][s] = 0.0;

    g[ 0][ 0] = D;
    g[ 4][ 1] = sqrt(2.0*sigma2_Na_1/tau_Na_1);
    g[ 5][ 2] = sqrt(2.0*sigma2_Na_2/tau_Na_2);
    g[ 6][ 3] = sqrt(2.0*sigma2_Na_3/tau_Na_3);
    g[ 7][ 4] = sqrt(2.0*sigma2_Na_4/tau_Na_4);
    g[ 8][ 5] = sqrt(2.0*sigma2_Na_5/tau_Na_5);
    g[ 9][ 6] = sqrt(2.0*sigma2_Na_6/tau_Na_6);
    g[10][ 7] = sqrt(2.0*sigma2_Na_7/tau_Na_7);
    g[11][ 8] = sqrt(2.0*sigma2_K_1/tau_K_1);
    g[12][ 9] = sqrt(2.0*sigma2_K_2/tau_K_2);
    g[13][10] = sqrt(2.0*sigma2_K_3/tau_K_3);
    g[14][11] = sqrt(2.0*sigma2_K_4/tau_K_4);

}

void Update_Diffusion_Deriv_HH_Channel_Noise(double* state, double*** g, double* pars)
{
    int r,s,t;
    for (r = 0; r < 14; r++)
    {
        for (s = 0; s < 11; s++)
        {
            for (t = 0; t < 14; t++)
                g[r][s][t] = 0.0;
        }
    }
}

int ParseData_HH_Channel_Noise(double* pars, double* Tspan, char* outfile, double* dt_addr, int* print_step_addr, FILE* f_in)
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
    if(!fscanf(f_in,"length\t%lf\n",pars+8))
    {
        printf("Errore nel leggere i dati!\n");
        return 0;
    }
    if(!fscanf(f_in,"diameter\t%lf\n",pars+9))
    {
        printf("Errore nel leggere i dati!\n");
        return 0;
    }
    if(!fscanf(f_in,"gs\t%lf\n",pars+10))
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
    if(!fscanf(f_in,"output_file\t%s\n",outfile))
    {
        printf("Errore nel leggere i dati!\n");
        return 0;
    }
    return 1;
}


int ParseData_HH_Channel_Noise_MT(double* pars, double* Tspan, double* dt_addr, int* print_step_addr, int* seed_addr, FILE* f_in)
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
    if(!fscanf(f_in,"length\t%lf\n",pars+8))
    {
        printf("Errore nel leggere i dati!\n");
        return 0;
    }
    if(!fscanf(f_in,"gs\t%lf\n",pars+10))
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
