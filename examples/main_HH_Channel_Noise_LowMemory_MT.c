#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <math.h>
#include <string.h>
#include <pthread.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_cblas.h>
#include <unistd.h>
#include "Solver.h"
#include "Stack.h"
#include "HH_Channel_Noise.h"


item* stack = NULL;
int Niters;
pthread_mutex_t mutexdata; 
int parnum = 9;

void* ThreadFunction(void* temp)
{
    printf("Started...\n");
    while(1)
    {
        thread_data p;
        char outfile[50];
        char buffer[50];
        double pars[10];

        pthread_mutex_lock (&mutexdata);
        if(Niters <= 0)
        {
            pthread_mutex_unlock (&mutexdata);
            break;
        }
        else
        {
            pop(&stack,&p);
            Niters--;
            pthread_mutex_unlock (&mutexdata);
        }

        pars[parnum] = p.val;
        strcpy(outfile,"./data_MT_prova/output_HH_diam_");
        sprintf(buffer,"%1.1f",pars[parnum]);
        strcat(outfile,buffer);
        strcat(outfile,".dat");

        FILE* f_in;
        f_in = fopen("input_HH_CN_MT.dat","r");
    
        double Tspan[2];
        Tspan[0] = 0.0;
        double dt;
        int seed;
        int print_step;
        double th = -20.0;

        if(!ParseData_HH_Channel_Noise_MT(pars,Tspan,&dt,&print_step,&seed,f_in)) 
            pthread_exit((void*) 1);
        fclose(f_in);
    
        FILE* f_out;
        f_out = fopen(outfile,"w");
  
        double x0[15] = {0.1, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};

        // Allocation and Initialization of Random Number Generator
        const gsl_rng_type * T;
        gsl_rng * r;
        gsl_rng_env_setup();
        T = gsl_rng_default;
        r = gsl_rng_alloc(T);
        if (seed < 0)
        {
            // Random Seed
            gsl_rng_set(r,time(NULL)*getpid());
        }
        else
        {
            // Deterministic Seed
            gsl_rng_set(r,seed);
        }

        SdeSolverSpikeDetection(Stepper_EM,Update_VF_HH_Channel_Noise,Update_Diffusion_HH_Channel_Noise,Update_Diffusion_Deriv_HH_Channel_Noise,Tspan,x0,15,12,dt,pars,r,f_out,th,0);

        fclose(f_out);
    }
    printf("Ended...\n");
    pthread_exit((void*) 0);
}


int main(int argc, char* argv[])
{
    int k;
    Niters = 6;
    int Nthreads = 6;
    double start = 1.0;
    double inc = 1.0;
    double diameter;
    for(k = 0; k < Niters; k++) 
    {
        diameter = start + k*inc;
        push(&stack,diameter,k); 
    }

    pthread_t* threads=(pthread_t *)malloc(Nthreads*sizeof(pthread_t));
    pthread_attr_t attr;
    pthread_attr_init(&attr);
    pthread_attr_setdetachstate(&attr,PTHREAD_CREATE_JOINABLE);

    for(k = 0; k < Nthreads; k++)
    {
        if(pthread_create(&threads[k], &attr, ThreadFunction, NULL)) 
        {
            printf("Fatal error forking thread # %d.\n", k);
            exit(1);
        }
    }

    for(k = 0; k < Nthreads; k++)
    {
        if(pthread_join(threads[k],(void **)NULL))
        {
            printf("Fatal error joining thread # %d.\n",k);
            exit(1);
        }
    }


    return 0;
}
