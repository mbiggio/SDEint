#include "Solver.h"


void Update_f(VF Update_VF, DIFF Update_Diffusion, DIFF_DER Update_Diffusion_Deriv, double* state, double* f, double* pars, double dt, int dim, int dim_noise, gsl_rng* r, double* VF, double** g, double*** gd)
{
    int k,h,l;


    double dW[dim_noise];
    Update_VF(state, VF, pars);
    Update_Diffusion(state, g, pars);
    Update_Diffusion_Deriv(state, gd, pars);
    for (h = 0; h < dim_noise; h++)
        dW[h] = gsl_ran_gaussian(r,1)*sqrt(dt);

    for (k = 0; k < dim; k++)
    {
        f[k] = VF[k]*dt;
        for (l = 0; l < dim; l++)
            for (h = 0; h < dim_noise; h++) f[k] += g[l][h]*gd[k][h][l]*dt;
        for (h = 0; h < dim_noise; h++)
            f[k] += g[k][h]*dW[h];
    }

}

void Stepper_EM(VF Update_VF, DIFF Update_Diffusion, DIFF_DER Update_Diffusion_Deriv, double* new_state, double* old_state, double* pars, double dt, gsl_rng* r, int dim, int dim_noise, double* VF, double** g, double*** gd)
{
    int k,h;
    Update_VF(old_state, VF, pars);
    Update_Diffusion(old_state, g, pars);
    double dW[dim_noise];

    for (h = 0; h < dim_noise; h++)
        dW[h] = gsl_ran_gaussian(r,1)*sqrt(dt);

    for (k = 0; k < dim; k++)
    {
        new_state[k] = old_state[k] + VF[k]*dt;
        for (h = 0; h < dim_noise; h++)
            new_state[k] += g[k][h]*dW[h];
    }
}

void Stepper_RK4(VF Update_VF, DIFF Update_Diffusion, DIFF_DER Update_Diffusion_Deriv, double* new_state, double* old_state, double* pars, double dt, gsl_rng* r, int dim, int dim_noise, double* VF, double** g, double*** gd)
{
    double K1[dim];
    double K2[dim];
    double K3[dim];
    double K4[dim];

    double tmp[dim];
    int k;

    // Compute Coefficient K1
    Update_f(Update_VF,Update_Diffusion,Update_Diffusion_Deriv,old_state, K1, pars, dt, dim, dim_noise, r, VF, g, gd);

    // Compute Coefficient K2
    for (k = 0; k < dim; k++) tmp[k] = old_state[k] + 0.5*K1[k];
    Update_f(Update_VF,Update_Diffusion,Update_Diffusion_Deriv, tmp, K2, pars, dt, dim, dim_noise, r, VF, g, gd);

    // Compute Coefficient K3
    for (k = 0; k < dim; k++) tmp[k] = old_state[k] + 0.5*K2[k];
    Update_f(Update_VF,Update_Diffusion,Update_Diffusion_Deriv, tmp, K3, pars, dt, dim, dim_noise, r, VF, g, gd);
    
    // Compute Coefficient K4
    for (k = 0; k < dim; k++) tmp[k] = old_state[k] + K3[k];
    Update_f(Update_VF,Update_Diffusion,Update_Diffusion_Deriv, tmp, K4, pars, dt, dim, dim_noise, r, VF, g, gd);

    // Compute the new state
    for (k = 0; k < dim; k++) new_state[k] = old_state[k] + (K1[k] + 2*K2[k] + 2*K3[k] + K4[k])/6.0;


}



void SdeSolver(STEPPER StepFunc, VF Update_VF, DIFF Update_Diffusion, DIFF_DER Update_Diffusion_Deriv, double** x, double* t, double* Tspan, double* x0, int dim, int dim_noise, double dt, double* pars, gsl_rng* r)
{    

    int k,h;
    // Initialization of vector field
    double* VF = (double*)malloc(dim*sizeof(double));
    // Initialization of diffusion matrix
    double**  g = (double**)malloc(dim*sizeof(double*));
    for (k = 0; k < dim; k++) g[k] = (double*)malloc(dim_noise*sizeof(double));
    // Initialization of derivatives of diffusion matrix
    double***  gd = (double***)malloc(dim*sizeof(double**));
    for (k = 0; k < dim; k++) gd[k] = (double**)malloc(dim_noise*sizeof(double*));
    for (k = 0; k < dim; k++)
        for (h = 0; h < dim_noise; h++) gd[k][h] = (double*)malloc(dim*sizeof(double));


    int steps = ceil(Tspan[1]/dt);
    
    t[0] = Tspan[0];
    for (h = 0; h < dim; h++) x[0][h] = x0[h];

    int print_std = floor(steps/100);
    printf("%s", GREEN);
    for(k = 1; k < steps; k++)
    {
        if (k % print_std == 0)
        {
            printf("\r");
            printf("[%03d/%03d]", k/print_std+1, 100);
            fflush(stdout);
        }
        StepFunc(Update_VF,Update_Diffusion,Update_Diffusion_Deriv,x[k],x[k-1],pars,dt,r,dim,dim_noise,VF,g,gd);
        t[k] = t[k-1] + dt;
    }
    printf("\n");
    printf("%s", NORMAL);


    // Deallocation of Vector Field
    free(VF);
    // Deallocation of diffusion matrix
    for (k = 0; k < dim; k ++) free(g[k]);
    free(g);
    // Deallocation of derivatives of diffusion matrix
    for (k = 0; k < dim; k ++)
        for (h = 0; h < dim_noise; h++) free(gd[k][h]);
    for (k = 0; k < dim; k++) free(gd[k]);
    free(gd);
}

double Distance(double* x, double* y, int dim)
{
    double length2 = 0;
    int k;
    for (k = 0; k < dim; k++) length2 += (x[k]-y[k])*(x[k]-y[k]);
    return sqrt(length2);
}

double Length(double* x, int dim)
{
    double length2 = 0;
    int k;
    for (k = 0; k < dim; k++) length2 += x[k]*x[k];
    return sqrt(length2);
}


void SdeSolverLyap(STEPPER StepFunc, VF Update_VF, DIFF Update_Diffusion, DIFF_DER Update_Diffusion_Deriv, double** x, double* t, double* Tspan, double* x0, int dim, int dim_noise, double dt, double* pars, gsl_rng* r, gsl_rng* r_lyap, double tau, double alpha, double* lyapexp,int print)
{   
    int k,h;
    // Initialization of vector field
    double* VF = (double*)malloc(dim*sizeof(double));
    // Initialization of diffusion matrix
    double**  g = (double**)malloc(dim*sizeof(double*));
    for (k = 0; k < dim; k++) g[k] = (double*)malloc(dim_noise*sizeof(double));
    // Initialization of derivatives of diffusion matrix
    double***  gd = (double***)malloc(dim*sizeof(double**));
    for (k = 0; k < dim; k++) gd[k] = (double**)malloc(dim_noise*sizeof(double*));
    for (k = 0; k < dim; k++)
        for (h = 0; h < dim_noise; h++) gd[k][h] = (double*)malloc(dim*sizeof(double));


    int steps = ceil(Tspan[1]/dt);
    int steps_lyap = ceil(tau/dt);
    int iterations = 0;
    double x_lyap[dim];
    double x_lyap_start[dim];
    double d0,dnew;
    
    
    *lyapexp = 0.0;
    t[0] = Tspan[0];
    for (h = 0; h < dim; h++) 
    {
        x[0][h] = x0[h];
	x_lyap_start[h] = (1+alpha)*x0[h];
	x_lyap[h] = x_lyap_start[h];
    }
    d0 = Distance(x_lyap_start,x0,dim);

    int print_std = floor(steps/100);

    if (print)
        printf("%s", GREEN);
    for(k = 1; k < steps; k++)
    {
        if (print && k % print_std == 0)
        {
            printf("\r");
            printf("[%03d/%03d]", k/print_std+1, 100);
            fflush(stdout);
        }
        StepFunc(Update_VF,Update_Diffusion,Update_Diffusion_Deriv,x[k],x[k-1],pars,dt,r,dim,dim_noise,VF,g,gd);
	StepFunc(Update_VF,Update_Diffusion,Update_Diffusion_Deriv,x_lyap,x_lyap,pars,dt,r_lyap,dim,dim_noise,VF,g,gd);
	if (k % steps_lyap == 0)
	{
	   dnew = Distance(x_lyap,x[k],dim);
	   for(h = 0; h < dim; h++) x_lyap[h] = x[k][h] + (x_lyap[h]-x[k][h])/dnew*d0;
	   *lyapexp += log(dnew/d0);
	   iterations++;
	}
        t[k] = t[k-1] + dt;
    }

    if (print)
    {
        printf("\n");
        printf("%s", NORMAL);
    }
    *lyapexp = *lyapexp/tau/iterations;


    // Deallocation of Vector Field
    free(VF);
    // Deallocation of diffusion matrix
    for (k = 0; k < dim; k ++) free(g[k]);
    free(g);
    // Deallocation of derivatives of diffusion matrix
    for (k = 0; k < dim; k ++)
        for (h = 0; h < dim_noise; h++) free(gd[k][h]);
    for (k = 0; k < dim; k++) free(gd[k]);
    free(gd);
}

void SdeSolverLog(STEPPER StepFunc, VF Update_VF, DIFF Update_Diffusion, DIFF_DER Update_Diffusion_Deriv, double* Tspan, double* x0, int dim, int dim_noise, double dt, double* pars, gsl_rng* r, FILE* f_out, int print_steps)
{    

    int k,h;
    // Initialization of vector field
    double* VF = (double*)malloc(dim*sizeof(double));
    // Initialization of diffusion matrix
    double**  g = (double**)malloc(dim*sizeof(double*));
    for (k = 0; k < dim; k++) g[k] = (double*)malloc(dim_noise*sizeof(double));
    // Initialization of derivatives of diffusion matrix
    double***  gd = (double***)malloc(dim*sizeof(double**));
    for (k = 0; k < dim; k++) gd[k] = (double**)malloc(dim_noise*sizeof(double*));
    for (k = 0; k < dim; k++)
        for (h = 0; h < dim_noise; h++) gd[k][h] = (double*)malloc(dim*sizeof(double));


    double x[dim];
    double t;
    int steps = ceil(Tspan[1]/dt);
    
    t = Tspan[0];
    for (h = 0; h < dim; h++) x[h] = x0[h];

    int print_std = floor(steps/100);
    printf("%s", GREEN);

    fprintf(f_out,"%f\t",t);
    #ifdef PRINTALL
    for (h = 0; h < dim; h++) fprintf(f_out,"%f\t",x[h]);
    #else
    fprintf(f_out,"%f\t",x[0]);
    #endif

    fprintf(f_out,"\n");

    for(k = 1; k < steps; k++)
    {
        if (k % print_std == 0)
        {
            printf("\r");
            printf("[%03d/%03d]", k/print_std+1, 100);
            fflush(stdout);
        }
        StepFunc(Update_VF,Update_Diffusion,Update_Diffusion_Deriv,x,x,pars,dt,r,dim,dim_noise,VF,g,gd);
        t = t + dt;
        if (k % print_steps == 0)
        {
            fprintf(f_out,"%f\t",t);
            #ifdef PRINTALL
            for (h = 0; h < dim; h++) fprintf(f_out,"%f\t",x[h]);
            #else
            fprintf(f_out,"%f\t",x[0]);
            #endif
            fprintf(f_out,"\n");
        }
    }
    printf("\n");
    printf("%s", NORMAL);


    // Deallocation of Vector Field
    free(VF);
    // Deallocation of diffusion matrix
    for (k = 0; k < dim; k ++) free(g[k]);
    free(g);
    // Deallocation of derivatives of diffusion matrix
    for (k = 0; k < dim; k ++)
        for (h = 0; h < dim_noise; h++) free(gd[k][h]);
    for (k = 0; k < dim; k++) free(gd[k]);
    free(gd);
}

void SdeSolverSpikeDetection(STEPPER StepFunc, VF Update_VF, DIFF Update_Diffusion, DIFF_DER Update_Diffusion_Deriv, double* Tspan, double* x0, int dim, int dim_noise, double dt, double* pars, gsl_rng* r, FILE* f_out, double th, int print)
{    

    int k,h,less_than_th = 0;
    // Initialization of vector field
    double* VF = (double*)malloc(dim*sizeof(double));
    // Initialization of diffusion matrix
    double**  g = (double**)malloc(dim*sizeof(double*));
    for (k = 0; k < dim; k++) g[k] = (double*)malloc(dim_noise*sizeof(double));
    // Initialization of derivatives of diffusion matrix
    double***  gd = (double***)malloc(dim*sizeof(double**));
    for (k = 0; k < dim; k++) gd[k] = (double**)malloc(dim_noise*sizeof(double*));
    for (k = 0; k < dim; k++)
        for (h = 0; h < dim_noise; h++) gd[k][h] = (double*)malloc(dim*sizeof(double));


    double x[dim];
    double t;
    int steps = ceil(Tspan[1]/dt);
    
    t = Tspan[0];
    for (h = 0; h < dim; h++) x[h] = x0[h];

    // Initialize less_than_th
    if (x[0] <= th) less_than_th = 1;

    int print_std = floor(steps/100);
    if(print)
    {
        printf("%s", GREEN);
    }

    for(k = 1; k < steps; k++)
    {
        if (print && k % print_std == 0)
        {
            printf("\r");
            printf("[%03d/%03d]", k/print_std+1, 100);
            fflush(stdout);
        }
        StepFunc(Update_VF,Update_Diffusion,Update_Diffusion_Deriv,x,x,pars,dt,r,dim,dim_noise,VF,g,gd);
        t = t + dt;
        if (x[0] >= th && less_than_th == 1)
        {
             less_than_th = 0;
             fprintf(f_out,"%.12f\n",t);
        }
        else if (x[0] < th && less_than_th == 0)
        {
             less_than_th = 1;
        }
    }
    
    if(print)
    {
        printf("\n");
        printf("%s", NORMAL);
    }


    // Deallocation of Vector Field
    free(VF);
    // Deallocation of diffusion matrix
    for (k = 0; k < dim; k ++) free(g[k]);
    free(g);
    // Deallocation of derivatives of diffusion matrix
    for (k = 0; k < dim; k ++)
        for (h = 0; h < dim_noise; h++) free(gd[k][h]);
    for (k = 0; k < dim; k++) free(gd[k]);
    free(gd);
}


void SdeSolverLyapLog(STEPPER StepFunc, VF Update_VF, DIFF Update_Diffusion, DIFF_DER Update_Diffusion_Deriv, double* Tspan, double* x0, int dim, int dim_noise, double dt, double* pars, gsl_rng* r, gsl_rng* r_lyap, double tau, double alpha, double* lyapexp, FILE* f_out, int print_steps, int print)
{   
    int k,h;
    // Initialization of vector field
    double* VF = (double*)malloc(dim*sizeof(double));
    // Initialization of diffusion matrix
    double**  g = (double**)malloc(dim*sizeof(double*));
    for (k = 0; k < dim; k++) g[k] = (double*)malloc(dim_noise*sizeof(double));
    // Initialization of derivatives of diffusion matrix
    double***  gd = (double***)malloc(dim*sizeof(double**));
    for (k = 0; k < dim; k++) gd[k] = (double**)malloc(dim_noise*sizeof(double*));
    for (k = 0; k < dim; k++)
        for (h = 0; h < dim_noise; h++) gd[k][h] = (double*)malloc(dim*sizeof(double));

    int steps = ceil(Tspan[1]/dt);
    int steps_lyap = ceil(tau/dt);
    int iterations = 0;
    double x_lyap[dim];
    double x_lyap_start[dim];
    double d0,dnew;
    
    double x[dim];
    double t; 
    *lyapexp = 0.0;
    t = Tspan[0];
    for (h = 0; h < dim; h++) 
    {
        x[h] = x0[h];
	x_lyap_start[h] = (1+alpha)*x0[h];
	x_lyap[h] = x_lyap_start[h];
    }
    
    int print_std = floor(steps/100);

    if (print)
    {
        printf("%s", GREEN);
    }

    fprintf(f_out,"%f\t",t);
    #ifdef PRINTALL
    for (h = 0; h < dim; h++) fprintf(f_out,"%f\t",x[h]);
    #else
    fprintf(f_out,"%f\t",x[0]);
    #endif

    fprintf(f_out,"\n");

    d0 = Distance(x_lyap_start,x0,dim);

    for(k = 1; k < steps; k++)
    {
        if (print && k % print_std == 0)
        {
            printf("\r");
            printf("[%03d/%03d]", k/print_std+1, 100);
            fflush(stdout);
        }
        StepFunc(Update_VF,Update_Diffusion,Update_Diffusion_Deriv,x,x,pars,dt,r,dim,dim_noise,VF,g,gd);
	StepFunc(Update_VF,Update_Diffusion,Update_Diffusion_Deriv,x_lyap,x_lyap,pars,dt,r_lyap,dim,dim_noise,VF,g,gd);
	if (k % steps_lyap == 0)
	{
	   dnew = Distance(x_lyap,x,dim);
	   for(h = 0; h < dim; h++) x_lyap[h] = x[h] + (x_lyap[h]-x[h])/dnew*d0;
	   *lyapexp += log(dnew/d0);
	   iterations++;
	}
        if (k % print_steps == 0)
        {
            fprintf(f_out,"%f\t",t);
            #ifdef PRINTALL
            for (h = 0; h < dim; h++) fprintf(f_out,"%f\t",x[h]);
            #else
            fprintf(f_out,"%f\t",x[0]);
            #endif
            fprintf(f_out,"\n");
        }
        t = t + dt;
    }

    if (print)
    {
        printf("\n");
        printf("%s", NORMAL);
    }
    *lyapexp = *lyapexp/tau/iterations;


    // Deallocation of Vector Field
    free(VF);
    // Deallocation of diffusion matrix
    for (k = 0; k < dim; k ++) free(g[k]);
    free(g);
    // Deallocation of derivatives of diffusion matrix
    for (k = 0; k < dim; k ++)
        for (h = 0; h < dim_noise; h++) free(gd[k][h]);
    for (k = 0; k < dim; k++) free(gd[k]);
    free(gd);
}



void SdeSolverLyap_NoLog(STEPPER StepFunc, VF Update_VF, DIFF Update_Diffusion, DIFF_DER Update_Diffusion_Deriv, double* Tspan, double* x0, int dim, int dim_noise, double dt, double* pars, gsl_rng* r, gsl_rng* r_lyap, double tau, double alpha, double* lyapexp, int print)
{   
    int k,h;
    // Initialization of vector field
    double* VF = (double*)malloc(dim*sizeof(double));
    // Initialization of diffusion matrix
    double**  g = (double**)malloc(dim*sizeof(double*));
    for (k = 0; k < dim; k++) g[k] = (double*)malloc(dim_noise*sizeof(double));
    // Initialization of derivatives of diffusion matrix
    double***  gd = (double***)malloc(dim*sizeof(double**));
    for (k = 0; k < dim; k++) gd[k] = (double**)malloc(dim_noise*sizeof(double*));
    for (k = 0; k < dim; k++)
        for (h = 0; h < dim_noise; h++) gd[k][h] = (double*)malloc(dim*sizeof(double));

    int steps = ceil(Tspan[1]/dt);
    int steps_lyap = ceil(tau/dt);
    int iterations = 0;
    double x_lyap[dim];
    double x_lyap_start[dim];
    double d0,dnew;
    
    double x[dim];
    double t; 
    *lyapexp = 0.0;
    t = Tspan[0];
    for (h = 0; h < dim; h++) 
    {
        x[h] = x0[h];
	x_lyap_start[h] = (1+alpha)*x0[h];
	x_lyap[h] = x_lyap_start[h];
    }
    
    int print_std = floor(steps/100);

    if (print)
    {
        printf("%s", GREEN);
    }

    d0 = Distance(x_lyap_start,x0,dim);

    for(k = 1; k < steps; k++)
    {
        if (print && k % print_std == 0)
        {
            printf("\r");
            printf("[%03d/%03d]", k/print_std+1, 100);
            fflush(stdout);
        }
        StepFunc(Update_VF,Update_Diffusion,Update_Diffusion_Deriv,x,x,pars,dt,r,dim,dim_noise,VF,g,gd);
	StepFunc(Update_VF,Update_Diffusion,Update_Diffusion_Deriv,x_lyap,x_lyap,pars,dt,r_lyap,dim,dim_noise,VF,g,gd);
	if (k % steps_lyap == 0)
	{
	   dnew = Distance(x_lyap,x,dim);
	   for(h = 0; h < dim; h++) x_lyap[h] = x[h] + (x_lyap[h]-x[h])/dnew*d0;
	   *lyapexp += log(dnew/d0);
	   iterations++;
	}
        t = t + dt;
    }

    if (print)
    {
        printf("\n");
        printf("%s", NORMAL);
    }
    *lyapexp = *lyapexp/tau/iterations;


    // Deallocation of Vector Field
    free(VF);
    // Deallocation of diffusion matrix
    for (k = 0; k < dim; k ++) free(g[k]);
    free(g);
    // Deallocation of derivatives of diffusion matrix
    for (k = 0; k < dim; k ++)
        for (h = 0; h < dim_noise; h++) free(gd[k][h]);
    for (k = 0; k < dim; k++) free(gd[k]);
    free(gd);
}

void WriteResults(FILE* f_out, double**x, double* t, int steps, int dim, int print_steps)
{
    int k,h;
    for (k = 0; k < steps; k++)
    {
        if (k % print_steps == 0)
        {
            fprintf(f_out,"%f\t",t[k]);
            for (h = 0; h < dim; h++) fprintf(f_out,"%f\t",x[k][h]);
            fprintf(f_out,"\n");
        }
    }
}

