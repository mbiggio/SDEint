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
#include "Leech_Channel_Noise_Orio.h"

//#define LYAP

item* stack = NULL;
int Niters;
pthread_mutex_t mutexdata; 
int parnum = 8;

void* ThreadFunction(void* temp)
{
    thread_data p;
    char outfile[50];
    char buffer[50];
    double pars[12];

    pthread_t tid = pthread_self();

    printf("Started thread number %ld ...\n", (long)tid);
    while(1)
    {
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


        FILE* f_in;
        f_in = fopen("input_Leech_CN_Orio_MT.dat","r");
    
        double Tspan[2];
        Tspan[0] = 0.0;
        double dt;
        int seed;
        int print_step;
        double th = -0.020;

        if(!ParseData_Leech_Channel_Noise_Orio_MT(pars,Tspan,&dt,&print_step,&seed,f_in)) 
            pthread_exit((void*) 1);
        fclose(f_in);

        pars[parnum] = p.val;
    
        strcpy(outfile,"./data_MT_Leech/output_Leech_Vshift_");
        sprintf(buffer,"%1.6f",pars[parnum]);
        strcat(outfile,buffer);
        strcat(outfile,".dat");

        FILE* f_out;
        f_out = fopen(outfile,"w");
  
        double x0[12] = {-0.046772,0.791195,0.196222,0.012583,0.000496,0.000115,-0.000008,0.000000,0.774144,0.206001,0.018730,0.000523};

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

        SdeSolverSpikeDetection
(Stepper_EM,Update_VF_Leech_Channel_Noise_Orio,Update_Diffusion_Leech_Channel_Noise_Orio,Update_Diffusion_Deriv_Leech_Channel_Noise_Orio,Tspan,x0,12,13,dt,pars,r,f_out,th,0);        

        fclose(f_out);
    }
    printf("Ended thread number %d ...\n", (int)tid);
    pthread_exit((void*) 0);
}


int main(int argc, char* argv[])
{
    int k;
    Niters = 100;
    int Nthreads = 6;
    double start = -0.0248;
    double end = -0.0238;
    double inc = (end-start)/Niters;
    double Vshift;
    for(k = 0; k < Niters; k++) 
    {
        Vshift = start + k*inc;
        push(&stack,Vshift,k); 
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
