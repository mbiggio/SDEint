#include "Leech_Channel_Noise_Orio.h"

void Update_VF_Leech_Channel_Noise_Orio(double* state, double* f, double* pars)
{
    //****** State Variables ******//
    // Membrane Potential //
    double V    = state[ 0];
    // Potassium Channel //
    double n0   = state[ 1];     
    double n1   = state[ 2];
    double n2   = state[ 3];
    // Sodium Channel //
    double m0h0 = state[ 4];     
    double m1h0 = state[ 5];
    double m2h0 = state[ 6];
    double m3h0 = state[ 7];    
    double m0h1 = state[ 8];    
    double m1h1 = state[ 9];        
    double m2h1 = state[10];            
    double m3h1 = state[11];            
    
    //****** Parameters ******//
    double C      = pars[ 0];
    double Iapp   = pars[ 1];
    double gNa    = pars[ 2];
    double gK2    = pars[ 3];
    double gL     = pars[ 4];
    double ENa    = pars[ 5];
    double EK2    = pars[ 6];
    double EL     = pars[ 7];
    double Vshift = pars[ 8];    

    //****** Time Constants   ******//
    //double tau_nK2 = 0.057 + 0.043/(1.0 + exp(200.0*(V + 0.035)));
    //double tau_mNa = 0.0001;            
    //double tau_hNa = 0.004 + 0.006/(1.0 + exp(500.0*(V + 0.028))) + 0.01/(cosh(300.0*(V + 0.027)));
    double tau_nK2 = 0.25;
    double tau_mNa = 0.0001;            
    double tau_hNa = 0.0405;            
    //****** Asymptotic Values  ******//
    //double nK2_inf = 1.0/(1.0 + exp(- 83.0*(V + 0.02 )));
    //double mNa_inf = 1.0/(1.0 + exp(-150.0*(V + 0.029)));
    //double hNa_inf = 1.0/(1.0 + 2.0*exp(180.0*(V + 0.047)) + exp(500.0*(V + 0.047)));    
    double nK2_inf = 1.0/(1.0 + exp(- 83.0*(V + 0.018 + Vshift)));
    double mNa_inf = 1.0/(1.0 + exp(-150.0*(V + 0.0305        )));
    double hNa_inf = 1.0/(1.0 + exp( 500.0*(V + 0.0333        )));
    //****** Transition Rates ******//
    double alpham = mNa_inf/tau_mNa;
    double betam  = (1.0 - mNa_inf)/tau_mNa;
    double alphah = hNa_inf/tau_hNa;
    double betah  = (1.0 - hNa_inf)/tau_hNa;
    double alphan = nK2_inf/tau_nK2;
    double betan  = (1.0 - nK2_inf)/tau_nK2;

    //****** Membrane Currents ******//
    double INa = gNa * m3h1 * (ENa - V);
    double IK2 = gK2 * n2   * (EK2 - V);
    double IL  = gL  *        (EL  - V);

    //****** Vector Field ******//
    // Membrane Potential //
    f[ 0] = (INa + IK2 + IL + Iapp)/C                                                                        ;
    // Potassium Channel //    
    f[ 1] = -2.0*alphan*n0   +     betan*n1                                                                  ;
    f[ 2] =  2.0*alphan*n0   + 2.0*betan*n2   -     alphan*n1   -     betan*n1                               ;
    f[ 3] =      alphan*n1   - 2.0*betan*n2                                                                  ;
    // Sodium Channel //        
    f[ 4] = -3.0*alpham*m0h0 +     betam*m1h0 -     alphah*m0h0 +     betah*m0h1                             ;
    f[ 5] =  3.0*alpham*m0h0 -     betam*m1h0 - 2.0*alpham*m1h0 + 2.0*betam*m2h0 - alphah*m1h0 + betah*m1h1  ;
    f[ 6] =  2.0*alpham*m1h0 - 2.0*betam*m2h0 -     alpham*m2h0 + 3.0*betam*m3h0 - alphah*m2h0 + betah*m2h1  ;    
    f[ 7] =      alpham*m2h0 - 3.0*betam*m3h0 -     alphah*m3h0 +     betah*m3h1                             ;
    f[ 8] = -3.0*alpham*m0h1 +     betam*m1h1 +     alphah*m0h0 -     betah*m0h1                             ;
    f[ 9] =  3.0*alpham*m0h1 -     betam*m1h1 - 2.0*alpham*m1h1 + 2.0*betam*m2h1 + alphah*m1h0 - betah*m1h1  ;
    f[10] =  2.0*alpham*m1h1 - 2.0*betam*m2h1 -     alpham*m2h1 + 3.0*betam*m3h1 + alphah*m2h0 - betah*m2h1  ;
    f[11] =      alpham*m2h1 - 3.0*betam*m3h1 +     alphah*m3h0 -     betah*m3h1                             ;
}

void Update_Diffusion_Leech_Channel_Noise_Orio(double* state, double** g, double* pars)
{
    //****** State Variables ******//
    // Membrane Potential //
    double V    = state[ 0];
    // Potassium Channel //
    double n0   = state[ 1];     
    double n1   = state[ 2];
    double n2   = state[ 3];
    // Sodium Channel //
    double m0h0 = state[ 4];     
    double m1h0 = state[ 5];
    double m2h0 = state[ 6];
    double m3h0 = state[ 7];    
    double m0h1 = state[ 8];    
    double m1h1 = state[ 9];        
    double m2h1 = state[10];            
    double m3h1 = state[11];            

    //****** Parameters ******//
    double Vshift = pars[ 8];        
    double NNa    = pars[ 9];
    double NK     = pars[10];    
    
    //****** Time Constants   ******//
    //double tau_nK2 = 0.057 + 0.043/(1.0 + exp(200.0*(V + 0.035)));
    //double tau_mNa = 0.0001;            
    //double tau_hNa = 0.004 + 0.006/(1.0 + exp(500.0*(V + 0.028))) + 0.01/(cosh(300.0*(V + 0.027)));
    double tau_nK2 = 0.25;
    double tau_mNa = 0.0001;            
    double tau_hNa = 0.0405;            
    //****** Asymptotic Values  ******//
    //double nK2_inf = 1.0/(1.0 + exp(- 83.0*(V + 0.02 )));
    //double mNa_inf = 1.0/(1.0 + exp(-150.0*(V + 0.029)));
    //double hNa_inf = 1.0/(1.0 + 2.0*exp(180.0*(V + 0.047)) + exp(500.0*(V + 0.047)));    
    double nK2_inf = 1.0/(1.0 + exp(- 83.0*(V + 0.018 + Vshift)));
    double mNa_inf = 1.0/(1.0 + exp(-150.0*(V + 0.0305        )));
    double hNa_inf = 1.0/(1.0 + exp( 500.0*(V + 0.0333        )));
    //****** Transition Rates ******//
    double alpham = mNa_inf/tau_mNa;
    double betam  = (1.0 - mNa_inf)/tau_mNa;
    double alphah = hNa_inf/tau_hNa;
    double betah  = (1.0 - hNa_inf)/tau_hNa;
    double alphan = nK2_inf/tau_nK2;
    double betan  = (1.0 - nK2_inf)/tau_nK2;

    double D = pars[11];

    //****** Initialization of Diffusion Matrix ******//
    int k,s;
    for(k = 0; k < 12; k++)
        for(s = 0; s < 13; s++) g[k][s] = 0.0;

    //****** Membrane Potential ******//
    g[ 0][ 0] = D;
    //****** Potassium Channel ******//    
    g[ 1][ 1] =   sqrt(fabs(2.0*alphan*n0   +     betan*n1  )) / sqrt(NK );
    g[ 2][ 1] = - sqrt(fabs(2.0*alphan*n0   +     betan*n1  )) / sqrt(NK );
    g[ 2][ 2] =   sqrt(fabs(    alphan*n1   + 2.0*betan*n2  )) / sqrt(NK );
    g[ 3][ 2] = - sqrt(fabs(    alphan*n1   + 2.0*betan*n2  )) / sqrt(NK );
    //****** Sodium Channel ******//        
    g[ 4][ 3] =   sqrt(fabs(3.0*alpham*m0h0 +     betam*m1h0)) / sqrt(NNa);
    g[ 4][ 6] =   sqrt(fabs(    alphah*m0h0 +     betah*m0h1)) / sqrt(NNa);
    g[ 5][ 3] = - sqrt(fabs(3.0*alpham*m0h0 +     betam*m1h0)) / sqrt(NNa);
    g[ 5][ 4] =   sqrt(fabs(2.0*alpham*m1h0 + 2.0*betam*m2h0)) / sqrt(NNa);
    g[ 5][ 7] =   sqrt(fabs(    alphah*m1h0 +     betah*m1h1)) / sqrt(NNa);
    g[ 6][ 4] = - sqrt(fabs(2.0*alpham*m1h0 + 2.0*betam*m2h0)) / sqrt(NNa);
    g[ 6][ 5] =   sqrt(fabs(    alpham*m2h0 + 3.0*betam*m3h0)) / sqrt(NNa);
    g[ 6][ 8] =   sqrt(fabs(    alphah*m2h0 +     betah*m2h1)) / sqrt(NNa);
    g[ 7][ 5] = - sqrt(fabs(    alpham*m2h0 + 3.0*betam*m3h0)) / sqrt(NNa);
    g[ 7][ 9] =   sqrt(fabs(    alphah*m3h0 +     betah*m3h1)) / sqrt(NNa);
    g[ 8][ 6] = - sqrt(fabs(    alphah*m0h0 +     betah*m0h1)) / sqrt(NNa);
    g[ 8][10] =   sqrt(fabs(3.0*alpham*m0h1 +     betam*m1h1)) / sqrt(NNa);
    g[ 9][10] = - sqrt(fabs(3.0*alpham*m0h1 +     betam*m1h1)) / sqrt(NNa);
    g[ 9][11] =   sqrt(fabs(2.0*alpham*m1h1 + 2.0*betam*m2h1)) / sqrt(NNa);
    g[ 9][ 7] = - sqrt(fabs(    alphah*m1h0 +     betah*m1h1)) / sqrt(NNa);
    g[10][11] = - sqrt(fabs(2.0*alpham*m1h1 + 2.0*betam*m2h1)) / sqrt(NNa);
    g[10][12] =   sqrt(fabs(    alpham*m2h1 + 3.0*betam*m3h1)) / sqrt(NNa);
    g[10][ 8] = - sqrt(fabs(    alphah*m2h0 +     betah*m2h1)) / sqrt(NNa);
    g[11][12] = - sqrt(fabs(    alpham*m2h1 + 3.0*betam*m3h1)) / sqrt(NNa);
    g[11][ 9] = - sqrt(fabs(    alphah*m3h0 +     betah*m3h1)) / sqrt(NNa);    
}

void Update_Diffusion_Deriv_Leech_Channel_Noise_Orio(double* state, double*** g, double* pars)
{
    int r,s,t;
    for (r = 0; r < 12; r++)
    {
        for (s = 0; s < 13; s++)
        {
            for (t = 0; t < 12; t++)
                g[r][s][t] = 0.0;
        }
    }
}

int ParseData_Leech_Channel_Noise_Orio(double* pars, double* Tspan, char* outfile, double* dt_addr, int* print_step_addr, FILE* f_in)
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
    if(!fscanf(f_in,"EK2\t%lf\n",pars+6))
    {
        printf("Errore nel leggere i dati!\n");
        return 0;
    }
    if(!fscanf(f_in,"EL\t%lf\n",pars+7))
    {
        printf("Errore nel leggere i dati!\n");
        return 0;
    }
    if(!fscanf(f_in,"Vshift\t%lf\n",pars+8))
    {
        printf("Errore nel leggere i dati!\n");
        return 0;
    }    
    if(!fscanf(f_in,"NNa\t%lf\n",pars+9))
    {
        printf("Errore nel leggere i dati!\n");
        return 0;
    }    
    if(!fscanf(f_in,"NK\t%lf\n",pars+10))
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

int ParseData_Leech_Channel_Noise_Orio_MT(double* pars, double* Tspan, double* dt_addr, int* print_step_addr, int* seed_addr, FILE* f_in)
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
    if(!fscanf(f_in,"EK2\t%lf\n",pars+6))
    {
        printf("Errore nel leggere i dati!\n");
        return 0;
    }
    if(!fscanf(f_in,"EL\t%lf\n",pars+7))
    {
        printf("Errore nel leggere i dati!\n");
        return 0;
    }
    if(!fscanf(f_in,"Vshift\t%lf\n",pars+8))
    {
        printf("Errore nel leggere i dati!\n");
        return 0;
    }    
    if(!fscanf(f_in,"NNa\t%lf\n",pars+9))
    {
        printf("Errore nel leggere i dati!\n");
        return 0;
    }    
    if(!fscanf(f_in,"NK\t%lf\n",pars+10))
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
