#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <math.h>
#include <string.h>
#include <pthread.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_cblas.h>
#include "Solver.h"
#include "Stack.h"
#include "Leech.h"
#include <unistd.h>

#define NITERS 200 

int Niters = NITERS;
item* stack = NULL;
pthread_mutex_t mutexdata; 
int parnum = 10;
double lyap[NITERS];

void* ThreadFunction(void* temp)
{
    printf("Started...\n");
    int iteration = 1;
    while(1)
    {


        printf("Iterazione Numero: %d\n",iteration);

        thread_data p;
        // char outfile[50];
        // char buffer[50];
        double pars[12];

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
        f_in = fopen("input_LEECH.dat","r");
    
        double Tspan[2];
        Tspan[0] = 0.0;
        double dt;
        int seed;
        int print_step;

        if(!ParseData_LEECH(pars,Tspan,&dt,&print_step,&seed,f_in)) 
            pthread_exit((void*) 1);
        fclose(f_in);

        pars[parnum] = p.val;

        /*
        strcpy(outfile,"./data_MT_prova/output_Leech_Vt_");
        sprintf(buffer,"%.8f",pars[parnum]);
        strcat(outfile,buffer);
        strcat(outfile,".dat");
 
        FILE* f_out;
        f_out = fopen(outfile,"w");
        */
  
        double x0[3] = {-0.0042, 0.0330, 0.3353};

        // Allocation and Initialization of Random Number Generator
        const gsl_rng_type * T;
        gsl_rng* r;
        gsl_rng* r_lyap;
        gsl_rng_env_setup();
        T = gsl_rng_default;
        r = gsl_rng_alloc(T);
        r_lyap = gsl_rng_alloc(T);
        if (seed < 0)
        {
            // Random Seed
            gsl_rng_set(r,     time(NULL)*getpid());
            gsl_rng_set(r_lyap,time(NULL)*getpid());
        }
        else
        {
            // Deterministic Seed
            gsl_rng_set(r,     seed);
            gsl_rng_set(r_lyap,seed);
        }

        double tau = 20*dt;
        double alpha = 1e-8;
        double lyapexp = 0;
      
        // SdeSolverLog(Stepper_RK4,Update_VF_HR,Update_Diffusion_HR,Tspan,x0,3,dt,pars,r,f_out,print_step);
        //SdeSolverLyapLog(Stepper_EM,Update_VF_LEECH,Update_Diffusion_LEECH,Update_Diffusion_Deriv_LEECH,Tspan,x0,3,1,dt,pars,r,r_lyap,tau,alpha,&lyapexp,f_out,print_step,0);
        SdeSolverLyap_NoLog(Stepper_EM,Update_VF_LEECH,Update_Diffusion_LEECH,Update_Diffusion_Deriv_LEECH,Tspan,x0,3,1,dt,pars,r,r_lyap,tau,alpha,&lyapexp,0);

        lyap[p.ndx] = lyapexp;
        // fclose(f_out);
        
        iteration++;
    }
    printf("Ended...\n");
    pthread_exit((void*) 0);
}


int main(int argc, char* argv[])
{
    int k;
    int Nthreads = 4;
    double start = -0.024;
    double end = -0.0223;
    double inc = (end - start)/NITERS;
    for(k = 0; k < NITERS; k++) 
        push(&stack,start + k*inc,k); 

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

    FILE* f_out;
    f_out = fopen("./data_MT_prova/lyapunov.dat","w");
    for (k = 0; k < NITERS; k++)
        fprintf(f_out,"%.15f\t%.15f\n",start+k*inc,lyap[k]);
    fclose(f_out);
        


    return 0;
}
