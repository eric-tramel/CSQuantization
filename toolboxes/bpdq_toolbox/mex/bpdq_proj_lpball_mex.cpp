/*
  bpdq_proj_lpball_mex.cpp
  Matlab gateway for lp ball projection. 
  Call as [x_out,n]=proj_lpball_newton(y_in,c,r,p)

  where y_in, c are N-length vectors
  r is a scalar radius
  p is a positive exponent >2

  Computes the projection of y_in onto the Lp ball
  x: ||x-c||_p < r
*/

/* This file is part of BPDQ Toolbox (Basis Pursuit DeQuantizer)
 * 
 * Copyright (C) 2009, the BPDQ Team (see the file AUTHORS distributed with
 * this library) (See the notice at the end of the file.) */

#include "mex.h"
#include "LpProjector.h"
#include "LpNormalizer.h"
#include <cmath>

#include <limits>
const double Inf = std::numeric_limits<double>::infinity();

const bool verbose=false;
//const bool verbose=true;
// matlab gateway
char usage_string[]=
"call as [x_out,n]=proj_lpball_newton(y_in,c,r,p)\n";
#define y_in_pa prhs[0]
#define c_pa prhs[1]
#define r_pa prhs[2]
#define p_pa prhs[3]

static LpNormalizer *myNormalizer=NULL;
static LpProjector *myProjector=NULL;

static void my_exit(){
  myNormalizer->delete_all();
  myProjector->delete_all();
  delete myNormalizer;
  delete myProjector;
}

void create_output(mxArray *plhs[],int nlhs,const double *out1,int N,int out2){
  if (nlhs>0){
    plhs[0]=mxCreateDoubleMatrix(N,1,mxREAL);
    double *xout_pr = mxGetPr(plhs[0]);
      for(int k=0;k<N;k++)
	xout_pr[k]=out1[k];
  }
  if (nlhs>1){
    plhs[1]=mxCreateDoubleMatrix(1,1,mxREAL);
    *(mxGetPr(plhs[1]))=out2;
  }
}

extern "C" void mexFunction(int nlhs,mxArray *plhs[],int nrhs,
			    const mxArray *prhs[]){
  // check inputs
  if (nrhs!=4){
    mexPrintf("%s",usage_string);
    return;
  }
  for(int k=0;k<4;k++){
    if ( !(mxIsDouble(prhs[k]) & !mxIsComplex(prhs[k]))){
      mexErrMsgTxt("All inputs must be real, double arrays\n");
    }
  }
  if (mxGetNumberOfElements(y_in_pa) !=mxGetNumberOfElements(c_pa)){
    mexErrMsgTxt("y_in and c must have same number of elements\n");
  }
  if ( (mxGetNumberOfElements(r_pa)!=1) || (mxGetNumberOfElements(p_pa)!=1))
    mexErrMsgTxt("r and p must be scalars\n");

  mexAtExit(my_exit);
  //////////////////////////////////////////////////
  // transfer input arguments to c++ variables
  int N;
  double r,p;
  p=*mxGetPr(p_pa);
  r=*mxGetPr(r_pa);
  N=mxGetNumberOfElements(y_in_pa);
  if (p<2){
    mexErrMsgTxt("do not use this code for p<2\n");
  }  

  double *y_in, *c;
  y_in=mxGetPr(y_in_pa);
  c=mxGetPr(c_pa);
  if (verbose)
    mexPrintf("p %f, r %f, N %i\n",p,r,N);
  //////////////////////////////////////////////////
  // initialize LpNormalizer, LpProjector objects
  if (!myNormalizer){
    if (verbose)
      mexPrintf("reallocating LpNormalizer\n");    
    myNormalizer=new LpNormalizer;
  }
  myNormalizer->init(N);
  if (!myProjector){
    if (verbose)
      mexPrintf("reallocating LpProjector\n");
    myProjector=new LpProjector();
  }
  myProjector->init(N); 
  
  // center, remove sign from input
  myNormalizer->center(y_in,c,r);
  myNormalizer->remove_sign();
  // input is within lp ball, just return
  if (pnorm(myNormalizer->y_n,N,p) < 1.0){
    if (verbose)
      mexPrintf("input within lp ball\n");
    create_output(plhs,nlhs,y_in,N,0);
    return;
  }  
  // actually do calculation
  if (verbose){
    mexPrintf("input not within lp ball\n");
    double dp=pnorm(myNormalizer->y_n,N,p);
    mexPrintf("pnorm %f\n",dp);
  }
  int numiter;
  myProjector->proj_lpball_newton_normalized(myNormalizer->y_n,
					     myNormalizer->x_out,p,numiter);
  // undo sign removal, normalization
  myNormalizer->replace_sign();
  myNormalizer->uncenter(c,r);
  create_output(plhs,nlhs,myNormalizer->x_out,N,numiter);
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
