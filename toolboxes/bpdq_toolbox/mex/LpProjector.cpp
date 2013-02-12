/*
  LpProjector.cpp
  
  code definitions for class LpProjector
  and some utility functions
*/
/* This file is part of BPDQ Toolbox (Basis Pursuit DeQuantizer)
 * 
 * Copyright (C) 2009, the BPDQ Team (see the file AUTHORS distributed with
 * this library) (See the notice at the end of the file.) */

#include "LpProjector.h"
#include <cmath>
#include <limits>

const double Inf = std::numeric_limits<double>::infinity();

LpProjector::LpProjector(){
  N=0;
  initialized=false;
}

LpProjector::LpProjector(int newN){
  init(newN);
}

LpProjector::~LpProjector(){
  delete_all();
}
void LpProjector::init(int newN){
  if (N==newN && initialized)
    return;
  delete_all(); 
  N=newN;
  allocate_all();
}
      
void LpProjector::delete_all(){
  if (initialized){
    delete [] F;
    delete [] z;
    delete [] dz;
    delete [] zpm1;
    delete [] xn1;
    delete [] xn2;
    delete [] btwid;
    delete [] d;
  }
  initialized=false;
}

void LpProjector::allocate_all(){
  //mexPrintf("allocating memory in LpProjector\n");
  F=new double[N+1];
  z=new double[N+1];
  dz=new double[N+1];
  zpm1=new double[N+1];
  xn1=new double[N];
  xn2=new double[N];
  btwid=new double[N];
  d=new double[N];  
  initialized=true;
}


void LpProjector::proj_lpball_newton_normalized(const double *y, double *xout,
					      double p,int &numiter){
  double normF,normz0,mu,chi;
  double tol=1e-15;
  numiter=0;
  // special case p=2, p=Inf
  if(p==2.0){
    radial_lp_project(y,xout,N,p);
    return;
  }
  if (p==Inf){
    l_infinity_project(y,xout,N);
    return;
  }
  
  //////////////////////////////////////////////////
  // Initialization
  // xn1 : radial projection onto Lp ball
  // xn2 : L\infty projection followed by radial projection
  // pick the one closest to y
  //////////////////////////////////////////////////
  radial_lp_project(y,xn1,N,p);
  l_infinity_project(y,xn2,N);
  radial_lp_project(xn2,xn2,N,p);
  if (dpnorm(xn1,y,N,2.0) < dpnorm(xn2,y,N,2.0)) {
    vcopy(xn1,z,N);  // initialize with xn1
  } 
  else {
    vcopy(xn2,z,N); // initialize with xn2
  }
  // initialize lagrange multiplier coordinate with least squares fit
  z[N]=lsq_lambda_init(z,y,p);

  // we are initialized!
  normz0=pnorm(z,N+1,2.0);
  normF=tol*normz0+1; 

  while ( (normF/normz0) > tol ) {
    // build residual F
    for (int k=0;k<N;k++){
      zpm1[k]=pow(z[k],p-1);
      F[k]=z[k]+z[N]*zpm1[k]-y[k];
    }
    double szp=0;
    for (int k=0;k<N;k++)
      szp+=pow(z[k],p);
    F[N]=(szp-1)/p;

    normF=pnorm(F,N+1,2.0);

    // build Jacobian matrix J
    for (int k=0;k<N;k++)
      d[k]=1+z[N]*(p-1)*pow(z[k],p-2);

    vdiv(zpm1,d,btwid,N);
    mu=dotp(zpm1,btwid,N);
    chi=-dotp(btwid,F,N); // uses only first N entries of F
    
    for (int k=0;k<N;k++)
      dz[k]=-F[k]/d[k] + btwid[k]*(-F[N]-chi)/mu;
    dz[N]=(chi+F[N])/mu;
    
    for (int k=0;k<N+1;k++)
      z[k]=z[k]+dz[k];

    if (verbose){
      vprint(F,N+1);
      mexPrintf("nF %e\n",normF); }

    numiter++;
    if (numiter>max_numiter)
      mexErrMsgTxt("maximum # of iterations exceeded in proj_lpball_newton\n");
  }
  // answer is first N coordinates of z.
  vcopy(z,xout,N);
}


// little utility functions

void vprint(double *a,int N){
  for (int k=0;k<N;k++)
    mexPrintf("%e\n",a[k]);
}

// compute radial projection onto unit lp_ball
inline void radial_lp_project(const double *xin,double *xout,int N,double p){
  double xnorm=pnorm(xin,N,p);
  for (int k=0;k<N;k++)
    xout[k]=xin[k]/xnorm;
}

// projection onto l_infinity ball ... assumes xin all positive
void l_infinity_project(const double *xin,double *xout,int N){
  for (int k=0;k<N;k++)
    xout[k]=min(xin[k],1.0);
}

inline double dotp(const double *a,const double *b,int N){
  double r=0;
  for (int k=0;k<N;k++)
    r+=a[k]*b[k];
  return r;
}

// c=a/b, element wise
inline void vdiv(const double *a,const double *b,double *c,int N){
  for (int k=0;k<N;k++)
    c[k]=a[k]/b[k];
}
inline void vcopy(const double *xin,double *xout,int N){
  for (int k=0;k<N;k++)
    xout[k]=xin[k];
}

// calculates sign of x and stores in xsign, and makes x positive
void remove_sign(double *x,bool *xsign,int N){
  bool tmp;
  for (int k=0;k<N;k++){
    tmp=(x[k]>0);
    x[k]=(tmp ? x[k] : -x[k]);
    xsign[k]=tmp;
  }
}

// negates each element of x when xsign is false
void replace_sign(double *x, const bool *xsign,int N){
  for (int k=0;k<N;k++)
    x[k]=(xsign[k]? x[k]:-x[k]);
}

// pick lambda best solving lambda*z_i^(p-1) = y_i-z_i
// lambda = sum (b_i*a_i) / sum (a_i^2)
// with a_i = z_i^(p-1) and b_i = y_i-z_i
double LpProjector::lsq_lambda_init(const double *z,const double *y, double p){
  if (!initialized)
    return 0;

  double a,b,absum,a2sum;
  absum=a2sum=0;
  for (int k=0;k<N;k++){
    a=pow(z[k],p-1.0);
    b=y[k]-z[k];
    absum+=a*b;
    a2sum+=a*a;
  }
  return absum/a2sum;
}

// ||a-b||_p
double dpnorm(const double *a,const double *b,int N,double p){
  if (p==Inf){
    return dnorm_inf(a,b,N);
  }
  else{
    double r=0;
    for (int k=0;k<N;k++){
      r+=pow( fabs(a[k]-b[k]),p);
    }
    r=pow(r,1.0/p);
    return r;
  }
}

// ||a||_p
double pnorm(const double *a, int N, double p){
  if (p==Inf){
    return norm_inf(a,N);
  }
  else{
    double r=0;
    for(int k=0;k<N;k++)
      r+=pow(fabs(a[k]),p);
    return pow(r,1.0/p);
  }
}

// ||a||_inf = max(abs(a_i))
double norm_inf(const double *a, int N){
  double currmax=0.0;
  double abs_a;
  for (int k=0;k<N;k++){
    abs_a=fabs(a[k]);
    if (abs_a>currmax)
      currmax=abs_a;
  }
  return currmax;
}

// ||a-b||_inf
double dnorm_inf(const double *a,const double *b, int N){
  double currmax=0.0;
  double abs_d; // absolute value of difference
  for (int k=0;k<N;k++){
    abs_d=fabs(a[k]-b[k]);
    if (abs_d>currmax)
      currmax=abs_d;
  }
  return currmax;
}



inline double min(double a,double b){
  return (a<b)?a:b;
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
