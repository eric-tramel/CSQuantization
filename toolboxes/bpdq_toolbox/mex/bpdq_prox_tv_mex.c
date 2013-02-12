/* 
 * TV-regularized minimization algorithm:
 *  min 0.5||y-x||_2 + gamma ||x||_TV
 *
 * The problem is solved thanks to the forward-backward method described 
 * in: 
 *  Multiplicative Noise Removal Using L1 Fidelity on Frame Coefficients,
 *  Durand S., Fadili J. and Nikolova M.
 *  arXiv:0812.1697v1
 *
 * 
 * This file is part of BPDQ Toolbox (Basis Pursuit DeQuantizer)
 * 
 * Copyright (C) 2009, the BPDQ Team (see the file AUTHORS distributed with
 * this library) (See the notice at the end of the file.) 
*/

#include "mex.h"
#include <string.h>
#include <math.h>

/* 
 * Discrete gradient.
 * Input:
 *	- u: image of size N-by-M
 * Ouputs:
 *  - nabla_u_1[i, j] = u[i+1][j]-u[i][j]
 *  - nabla_u_2[i, j] = u[i][j+1]-u[i][j] 
 *
 */
void dis_gradient(double **nabla_u_1, double **nabla_u_2, double **u, int N, int M)
{
    int i,j;
    
    for(i=0;i<N;i++) for(j=0;j<M;j++)
	{
		if(i<N-1) nabla_u_1[i][j] = u[i+1][j]-u[i][j];
		else nabla_u_1[N-1][j] = 0;
            
        if(j<M-1) nabla_u_2[i][j] = u[i][j+1]-u[i][j];
		else nabla_u_2[i][M-1] = 0;
	}
}

/*
 * Discrete divergence (adjoint of the discrete gradient).
 * Inputs:
 *  - nabla_u_1 (see dis_gradient) (size N-by-M)
 *  - nabla_u_2 (see dis_gradient) (size N-by-M)
 * Output:
 *  - div_u
 */
void dis_divergence(double **div_u, double **nabla_u_1, double **nabla_u_2, int N, int M)
{
	int i,j;
    
    for(i=0;i<N;i++)
        for(j=0;j<M;j++)
        {            
            if (i==0) div_u[i][j] = nabla_u_1[i][j];
            else if (i==N-1) div_u[i][j] = -nabla_u_1[i-1][j];
            else div_u[i][j] = nabla_u_1[i][j]-nabla_u_1[i-1][j];
            
            if (j==0) div_u[i][j] += nabla_u_2[i][j];
            else if (j==M-1) div_u[i][j] += -nabla_u_2[i][j-1];
            else div_u[i][j] += nabla_u_2[i][j]-nabla_u_2[i][j-1];
        }    
}

/* 
 * Total variation of an image u (size N-by-M)
 */
double tv_norm(double **u, int N, int M) 
{
    double norm;
    double **nabla_u_1, **nabla_u_2;
    int i,j;
    
    /* Allocate memory */
    nabla_u_1 = (double **) mxCalloc(N, sizeof(double *));
    nabla_u_2 = (double **) mxCalloc(N, sizeof(double *));
    for (i=0;i<N;i++)
    {
        nabla_u_1[i]=(double *) mxCalloc(M, sizeof(double));
        nabla_u_2[i]=(double *) mxCalloc(M, sizeof(double));
    }
    
    /* Gradient of u */
    dis_gradient(nabla_u_1, nabla_u_2, u, N, M);
    
    /* TV norm */
    norm = 0;
    for(i=0;i<N;i++) for(j=0;j<M;j++)
		norm += sqrt(pow(nabla_u_1[i][j],2) + pow(nabla_u_2[i][j],2));    
    
    /* Free memory */
    for (i=0;i<N;i++) 
    {
        mxFree(nabla_u_1[i]);
        mxFree(nabla_u_2[i]);
    }
    mxFree(nabla_u_1);
    mxFree(nabla_u_2);
    
    return norm;
}


/*
 * TV-regularized minimization algorithm.
 * Inputs:
 *  - y: noisy image (size N-by-M)
 *  - gamma: regularization parameter
 *  - tol: minimun relative change of the objective value (stopping criterion)
 *  - it_max, it_min: minimun and maximun number of iteration
 *  - verbose: print log or not
 * Output:
 *  - x: image (size N-by-M)
 *
 */
void tv_min(double **x, double **y, double gamma, double beta, int N, int M, 
	double tol, int it_max, int it_min, int verbose)
{
    double **z_1, **z_2;
    double **z_tilde_1, **z_tilde_2;
    double **div_z, **res;
    double obj, prev_obj, rel_obj;
    double norm;    
    int i, j, iter;
    
    /* 
     * Allocate memory
     */    
    z_1 = (double **) mxCalloc(N, sizeof(double *));
    z_2 = (double **) mxCalloc(N, sizeof(double *));
    res = (double **) mxCalloc(N, sizeof(double *));
    z_tilde_1 = (double **) mxCalloc(N, sizeof(double *));
    z_tilde_2 = (double **) mxCalloc(N, sizeof(double *));
    div_z = (double **) mxCalloc(N, sizeof(double *));    
    for (i=0;i<N;i++)
    {
        z_1[i]=(double *) mxCalloc(M, sizeof(double));
        z_2[i]=(double *) mxCalloc(M, sizeof(double));
        res[i]=(double *) mxCalloc(M, sizeof(double));
        z_tilde_1[i]=(double *) mxCalloc(M, sizeof(double));
        z_tilde_2[i]=(double *) mxCalloc(M, sizeof(double));
        div_z[i]=(double *) mxCalloc(M, sizeof(double));
    }
    
    /*
     * Initialization
     */
    for (i=0;i<N;i++) for (j=0;j<M;j++)
	{
		z_1[i][j] = 0; z_2[i][j] = 0;
	}
    
    /* 
     * Main loop
     */
    iter = 1; 
    rel_obj = 0;
    prev_obj = 0;
    while (1)
    {
        /* 
		 * Current estimate of the solution and
		 * computation of the objective value 
		 */
        dis_divergence(div_z, z_1, z_2, N, M);
		obj = 0; 
		for (i=0;i<N;i++) for (j=0;j<M;j++)
		{		
			x[i][j] = y[i][j] - div_z[i][j];
			obj += pow(fabs(div_z[i][j]), 2);
		}
		obj = obj/2 + tv_norm(x, N, M);
		
		/* Relative change of the objective value between two consecutive iteration */
        if (iter>1)
			rel_obj = fabs(obj-prev_obj)/prev_obj;
        
		/* LOG */
		if (verbose)
			mexPrintf("Iter: %i, obj = %e, rel_obj = %e\n", iter, obj, rel_obj);
		
        /* Stopping criterions */        
        if (iter > it_min && (rel_obj < tol || iter > it_max))                        
            break;		
        
        /* 
		 * Compute gradient of the current estimate and
		 * projection onto the L_infinite ball of radius gamma
		 */
        dis_gradient(z_tilde_1, z_tilde_2, x, N, M);
        for (i=0;i<N;i++) for (j=0;j<M;j++)
		{
			z_1[i][j] = (z_1[i][j] - beta*z_tilde_1[i][j]);
			z_2[i][j] = (z_2[i][j] - beta*z_tilde_2[i][j]);
			norm = sqrt(pow(z_1[i][j],2) + pow(z_2[i][j],2))/gamma;
			if (norm>1)
			{
				z_1[i][j] = z_1[i][j]/norm;
                z_2[i][j] = z_2[i][j]/norm;
			}
		}        

        /* Update variables */
        prev_obj = obj;
        iter = iter + 1;
    }
    
    /*
     * FREE MEMORY
     *
     */    
    for (i=0;i<N;i++) 
    {
        mxFree(z_1[i]); mxFree(z_2[i]); 
        mxFree(res[i]); mxFree(div_z[i]); 
        mxFree(z_tilde_1[i]); mxFree(z_tilde_2[i]);
    }
    mxFree(z_1); mxFree(z_2); 
    mxFree(res); mxFree(div_z);
    mxFree(z_tilde_1); mxFree(z_tilde_2);
	
	return; 
}



void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) 
{ 
    /* Main variables */
    double *y_1D; /*Input image (1D) */
    double *x_1D; /*Output image (1D) */
    double gamma; /* Regularization parameter */
	int nb_l, nb_c; /* Size of the input image */	
	
    /* Other variables */
    double min_rel_obj; /* Minimun relative change of the objective function */
	double beta;
	int it_max, it_min; /* Minimun and maximun number of iterations */
	int verbose;    
	double **y_2D; /* Input image (2D) */
	double **x_2D; /* Output image (2D) */	
	int i, j, len;
	double *value;
	char opt[11];
    
		    
    /*
     * CHECK NUMBER OF ARGUMENTS
	 */       
    if (nrhs < 2) 
        mexErrMsgTxt("At least 2 input argument are required");
	if (nlhs < 1)
		mexErrMsgTxt("One output argument is required.");	
	
	
	
    /*
     * READ INPUT ARGUMENTS     
     */
	
	/* First argument: the image */     
		/* Check if the image is real and if the image is in double precision */
	if (mxIsComplex(prhs[0]) || !mxIsDouble(prhs[0]))
		mexErrMsgTxt("Input image must be a noncomplex array in double precision.");	
		/* Pointer to the image.*/
	y_1D = mxGetPr(prhs[0]);  
		/* Dimensions of the image. */
	nb_l = mxGetM(prhs[0]);	nb_c = mxGetN(prhs[0]);    
    
	/* Second argument: gamma */
		/* Check that gamma is real */
    if (mxIsComplex(prhs[1]) || !mxIsDouble(prhs[1]))
        mexErrMsgTxt("Gamma must be real value in double precision.");    
    gamma = *mxGetPr(prhs[1]);    
    
    /* Optional parameters */
    beta = 0.249; min_rel_obj = 0.001; it_max=500; it_min=1; verbose=0;
	if (nrhs%2 == 1)
		mexErrMsgTxt("Optional parameters should always go by pairs.");	
	for (i=2;i<=nrhs-2;i+=2) 
    {
		/* Check parameters */
        if (!mxIsChar(prhs[i]))
            mexErrMsgTxt("Optional parameter definition must be a string.");	
		if (!mxIsDouble(prhs[i+1]))
			mexErrMsgTxt("Optional parameter value must be a value in double precision.");
			
		/* Read optional parameter definition */
			/* Length of the input string. */
		len = (mxGetM(prhs[i]) * mxGetN(prhs[i])) + 1;
			/* Copy it into opt*/
		mxGetString(prhs[i], opt, len);
		/* Read optional parameter value */
		value = mxGetPr(prhs[i+1]);
				
		if(strcmp(opt, "beta") == 0) 
        {
			beta = *value;            
        }
		else if(strcmp(opt, "min_rel_obj")==0) 
        {
            min_rel_obj = *value;
        }
        else if(strcmp(opt, "it_max")==0) 
        {
            it_max = (int) *value;
        }
		else if(strcmp(opt, "it_min")==0) 
        {
            if ((*value)>1) it_min = (int) *value;
        }
		else if(strcmp(opt, "verbose")==0) 
        {
            verbose = (int) *value;
        }
		else
            mexErrMsgTxt("Unrecognized parameter.");
    }
    
    
    /* 
     * ALLOCATE MEMORY
	 */
    y_2D = (double **) mxCalloc(nb_l, sizeof(double *));
    for (i=0;i<nb_l;i++) y_2D[i]=(double *) mxCalloc(nb_c, sizeof(double));
    x_2D = (double **) mxCalloc(nb_l, sizeof(double *));    
    for (i=0;i<nb_l;i++) x_2D[i]=(double *) mxCalloc(nb_c, sizeof(double));        
    
	
    /*
     * INITIALIZATION
	 */        
		/* Transform 1D input image in a 2D image */
    for (i=0;i<nb_l;i++) for (j=0;j<nb_c;j++)
		y_2D[i][j] = y_1D[i+nb_l*j];
    
	
    /*
     * TV MINIMIZATION    
     */
    tv_min(x_2D, y_2D, gamma, beta, nb_l, nb_c, min_rel_obj, it_max, 
        it_min, verbose);
	
	/*
     * COPY SOLUTION     
     */
    plhs[0] = mxCreateDoubleMatrix(nb_l, nb_c, mxREAL);
	x_1D = mxGetPr(plhs[0]);    
    for (i=0;i<nb_l;i++) for (j=0;j<nb_c;j++)
		x_1D[i+nb_l*j] = (double) x_2D[i][j];
		
	/*
     * FREE MEMORY     
     */
    for (i=0;i<nb_l;i++) 
    {
        mxFree(y_2D[i]); mxFree(x_2D[i]);        
    }
    mxFree(y_2D); mxFree(x_2D);    
 
    return;
}

/* The BPDQ Toolbox is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *  
 * The BPDQ Toolbox is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *  
 * You should have received a copy of the GNU General Public License
 * along with The BPDQ Toolbox. If not, see
 *  <http://www.gnu.org/licenses/>. */
