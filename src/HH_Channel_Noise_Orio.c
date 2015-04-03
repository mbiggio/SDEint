#include "HH_Channel_Noise_Orio.h"

void Update_VF_HH_Channel_Noise_Orio(double* state, double* f, double* pars)
{
    //****** State Variables ******//
    // Membrane Potential //
    double V    = state[ 0];
    // Potassium Channel //
    double n0   = state[ 1];     
    double n1   = state[ 2];
    double n2   = state[ 3];
    double n3   = state[ 4];    
    double n4   = state[ 5];    
    // Sodium Channel //
    double m0h0 = state[ 6];     
    double m1h0 = state[ 7];
    double m2h0 = state[ 8];
    double m3h0 = state[ 9];    
    double m0h1 = state[10];    
    double m1h1 = state[11];        
    double m2h1 = state[12];            
    double m3h1 = state[13];            
    
    //****** Parameters ******//
    double C      = pars[ 0];
    double Iapp   = pars[ 1];
    double gNa    = pars[ 2];
    double gK     = pars[ 3];
    double gL     = pars[ 4];
    double ENa    = pars[ 5];
    double EK     = pars[ 6];
    double EL     = pars[ 7];
    
    //****** Transition Rates ******//
    double alpham = 0.1*(V + 40.0)/(1.0 - exp(-(V + 40.0)/10.0));
    double betam  = 4.0*exp(-(V + 65.0)/18.0);
    double alphah = 0.07*exp(-(V + 65.0)/20.0);
    double betah  = 1.0/(1.0 + exp(-(V + 35.0)/10.0));
    double alphan = 0.01*(V + 55.0)/(1.0 - exp(-(V + 55.0)/10.0));
    double betan  = 0.125*exp(-(V + 65.0)/80.0);

    //****** Membrane Currents ******//
    double INa = gNa * m3h1 * (ENa - V);
    double IK2 = gK  * n4   * (EK  - V);
    double IL  = gL  *        (EL  - V);

    //****** Vector Field ******//
    // Membrane Potential //
    f[ 0] = (INa + IK2 + IL + Iapp)/C                                                                        ;
    // Potassium Channel //    
    f[ 1] = -4.0*alphan*n0   +     betan*n1;
    f[ 2] =  4.0*alphan*n0   -     betan*n1   - 3.0*alphan*n1   + 2.0*betan*n2                               ;
    f[ 3] =  3.0*alphan*n1   - 2.0*betan*n2   - 2.0*alphan*n2   + 3.0*betan*n3                               ;    
    f[ 4] =  2.0*alphan*n2   - 3.0*betan*n3   -     alphan*n3   + 4.0*betan*n4                               ;
    f[ 5] =      alphan*n3   - 4.0*betan*n4                                                                  ;
    // Sodium Channel //        
    f[ 6] = -3.0*alpham*m0h0 +     betam*m1h0 -     alphah*m0h0 +     betah*m0h1                             ;
    f[ 7] =  3.0*alpham*m0h0 -     betam*m1h0 - 2.0*alpham*m1h0 + 2.0*betam*m2h0 - alphah*m1h0 + betah*m1h1  ;
    f[ 8] =  2.0*alpham*m1h0 - 2.0*betam*m2h0 -     alpham*m2h0 + 3.0*betam*m3h0 - alphah*m2h0 + betah*m2h1  ;    
    f[ 9] =      alpham*m2h0 - 3.0*betam*m3h0 -     alphah*m3h0 +     betah*m3h1                             ;
    f[10] = -3.0*alpham*m0h1 +     betam*m1h1 +     alphah*m0h0 -     betah*m0h1                             ;
    f[11] =  3.0*alpham*m0h1 -     betam*m1h1 - 2.0*alpham*m1h1 + 2.0*betam*m2h1 + alphah*m1h0 - betah*m1h1  ;
    f[12] =  2.0*alpham*m1h1 - 2.0*betam*m2h1 -     alpham*m2h1 + 3.0*betam*m3h1 + alphah*m2h0 - betah*m2h1  ;
    f[13] =      alpham*m2h1 - 3.0*betam*m3h1 +     alphah*m3h0 -     betah*m3h1                             ;
}

void Update_Diffusion_HH_Channel_Noise_Orio(double* state, double** g, double* pars)
{
    //****** State Variables ******//
    // Membrane Potential //
    double V    = state[ 0];
    // Potassium Channel //
    double n0   = state[ 1];     
    double n1   = state[ 2];
    double n2   = state[ 3];
    double n3   = state[ 4];    
    double n4   = state[ 5];    
    // Sodium Channel //
    double m0h0 = state[ 6];     
    double m1h0 = state[ 7];
    double m2h0 = state[ 8];
    double m3h0 = state[ 9];    
    double m0h1 = state[10];    
    double m1h1 = state[11];        
    double m2h1 = state[12];            
    double m3h1 = state[13];            

    //****** Parameters ******//
    double NNa    = pars[ 8];
    double NK     = pars[ 9];    

    
    //****** Transition Rates ******//
    double alpham = 0.1*(V + 40.0)/(1.0 - exp(-(V + 40.0)/10.0));
    double betam  = 4.0*exp(-(V + 65.0)/18.0);
    double alphah = 0.07*exp(-(V + 65.0)/20.0);
    double betah  = 1.0/(1.0 + exp(-(V + 35.0)/10.0));
    double alphan = 0.01*(V + 55.0)/(1.0 - exp(-(V + 55.0)/10.0));
    double betan  = 0.125*exp(-(V + 65.0)/80.0);

    double D = pars[10];

    //****** Initialization of Diffusion Matrix ******//
    int k,s;
    for(k = 0; k < 14; k++)
        for(s = 0; s < 15; s++) g[k][s] = 0.0;

    //****** Membrane Potential ******//
    g[ 0][ 0] = D;
    //****** Potassium Channel ******//    
    g[ 1][ 1] =   sqrt(fabs(4.0*alphan*n0   +     betan*n1  )) / sqrt(NK );
    g[ 2][ 1] = - sqrt(fabs(4.0*alphan*n0   +     betan*n1  )) / sqrt(NK );
    g[ 2][ 2] =   sqrt(fabs(3.0*alphan*n1   + 2.0*betan*n2  )) / sqrt(NK );
    g[ 3][ 2] = - sqrt(fabs(3.0*alphan*n1   + 2.0*betan*n2  )) / sqrt(NK );
    g[ 3][ 3] =   sqrt(fabs(2.0*alphan*n2   + 3.0*betan*n3  )) / sqrt(NK );
    g[ 4][ 3] = - sqrt(fabs(2.0*alphan*n2   + 3.0*betan*n3  )) / sqrt(NK );
    g[ 4][ 4] =   sqrt(fabs(    alphan*n3   + 4.0*betan*n4  )) / sqrt(NK );
    g[ 5][ 4] = - sqrt(fabs(    alphan*n3   + 4.0*betan*n4  )) / sqrt(NK );
    //****** Sodium Channel ******//        
    g[ 6][ 5] =   sqrt(fabs(3.0*alpham*m0h0 +     betam*m1h0)) / sqrt(NNa);
    g[ 6][ 8] =   sqrt(fabs(    alphah*m0h0 +     betah*m0h1)) / sqrt(NNa);
    g[ 7][ 5] = - sqrt(fabs(3.0*alpham*m0h0 +     betam*m1h0)) / sqrt(NNa);
    g[ 7][ 6] =   sqrt(fabs(2.0*alpham*m1h0 + 2.0*betam*m2h0)) / sqrt(NNa);
    g[ 7][ 9] =   sqrt(fabs(    alphah*m1h0 +     betah*m1h1)) / sqrt(NNa);
    g[ 8][ 6] = - sqrt(fabs(2.0*alpham*m1h0 + 2.0*betam*m2h0)) / sqrt(NNa);
    g[ 8][ 7] =   sqrt(fabs(    alpham*m2h0 + 3.0*betam*m3h0)) / sqrt(NNa);
    g[ 8][10] =   sqrt(fabs(    alphah*m2h0 +     betah*m2h1)) / sqrt(NNa);
    g[ 9][ 7] = - sqrt(fabs(    alpham*m2h0 + 3.0*betam*m3h0)) / sqrt(NNa);
    g[ 9][11] =   sqrt(fabs(    alphah*m3h0 +     betah*m3h1)) / sqrt(NNa);
    g[10][ 8] = - sqrt(fabs(    alphah*m0h0 +     betah*m0h1)) / sqrt(NNa);
    g[10][12] =   sqrt(fabs(3.0*alpham*m0h1 +     betam*m1h1)) / sqrt(NNa);
    g[11][12] = - sqrt(fabs(3.0*alpham*m0h1 +     betam*m1h1)) / sqrt(NNa);
    g[11][13] =   sqrt(fabs(2.0*alpham*m1h1 + 2.0*betam*m2h1)) / sqrt(NNa);
    g[11][ 9] = - sqrt(fabs(    alphah*m1h0 +     betah*m1h1)) / sqrt(NNa);
    g[12][13] = - sqrt(fabs(2.0*alpham*m1h1 + 2.0*betam*m2h1)) / sqrt(NNa);
    g[12][14] =   sqrt(fabs(    alpham*m2h1 + 3.0*betam*m3h1)) / sqrt(NNa);
    g[12][10] = - sqrt(fabs(    alphah*m2h0 +     betah*m2h1)) / sqrt(NNa);
    g[13][14] = - sqrt(fabs(    alpham*m2h1 + 3.0*betam*m3h1)) / sqrt(NNa);
    g[13][11] = - sqrt(fabs(    alphah*m3h0 +     betah*m3h1)) / sqrt(NNa);    
}

void Update_Diffusion_Deriv_HH_Channel_Noise_Orio(double* state, double*** g, double* pars)
{
    int r,s,t;
    for (r = 0; r < 14; r++)
    {
        for (s = 0; s < 15; s++)
        {
            for (t = 0; t < 14; t++)
                g[r][s][t] = 0.0;
        }
    }
}

int ParseData_HH_Channel_Noise_Orio(double* pars, double* Tspan, char* outfile, double* dt_addr, int* print_step_addr, FILE* f_in)
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
    if(!fscanf(f_in,"gK\t%lf\n",pars+3))
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
    if(!fscanf(f_in,"NNa\t%lf\n",pars+8))
    {
        printf("Errore nel leggere i dati!\n");
        return 0;
    }    
    if(!fscanf(f_in,"NK\t%lf\n",pars+9))
    {
        printf("Errore nel leggere i dati!\n");
        return 0;
    }        
    if(!fscanf(f_in,"D\t%lf\n",pars+10))
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

