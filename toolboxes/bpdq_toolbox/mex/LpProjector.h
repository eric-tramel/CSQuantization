/*
  LpProjector.h
  
  Defines class  LpProjector
  and assorted utility functions for working with N-element vectors
*/

/* This file is part of BPDQ Toolbox (Basis Pursuit DeQuantizer)
 * 
 * Copyright (C) 2009, the BPDQ Team (see the file AUTHORS distributed with
 * this library) (See the notice at the end of the file.) */

#ifndef LPPROJECTOR_H_
#define LPPROJECTOR_H_
#include "mex.h"


/*
  LpProjector

  LpProjector solves the problem of projecting y onto the 
  unit Lp ball, i.e. x solving
  argmin ||xout-y||_2 s.t. ||xout||_p = 1

  Intended use : 
  LpProjector->proj_lpball_newton_normalized(y,xout,p,n)
  
  y and xout are N length arrays of doubles. Before use, the class
  must be initialized with the length of the vectors

  myLpProjection->init(N)

  Repeated initialization with the same N will not reallocate memory, while
  initialization with a subsequently different N will reallocate memory.
  This is so the class can be a persistant object and incur no performance
  penalty when repeatedly initialized (as done in the matlab gateway)

  This class implements an iterative Newton-Raphson method for solving
  the lagrange multiplier equations corresponding to the constrained 
  optimization problem. 

  Iteration continues until the relative norm of the residual
  is below specified tolerance tol, or maximum number of iterations reached
  // TODO : tol is hard coded in proj_lpball_newton_normalized, it should
  // TODO : made class member
*/
class LpProjector{
 public:
  static const bool verbose=0;
  LpProjector(int newN);
  LpProjector();
  ~LpProjector();
  void proj_lpball_newton_normalized(const double *y, double *xout, 
				     double p,int &numiter);
  void init(int newN);
  void delete_all();
 private:
  int N;
  double *F,*z,*dz,*zpm1; // size N+1
  double *xn1,*xn2,*btwid,*d; //size N
  bool initialized;
  static const int max_numiter=50;
  double lsq_lambda_init(const double *z,const double *y, double p);
  void allocate_all();
};


/* 
   Utility functions for l2 and l_infinity projection. These are used
   for initializing the Lp projection.
*/
inline void radial_lp_project(const double *xin,double *xout,int N,double p);
void l_infinity_project(const double *xin,double *xout,int N);

/* Utility functions for N length vectors as arrays of doubles.
   Commented in .cpp file.
 */
double dotp(const double *a,const double *b,int N);
double dpnorm(const double *a,const double *b,int N,double p);
double dnorm_inf(const double *a,const double *b,int N);
double pnorm(const double *a, int N, double p);
double norm_inf(const double *a, int N);
inline void vcopy(const double *xin,double *xout,int N);
inline void vdiv(const double *a,const double *b,double *c,int N);
void vprint(double *a,int N);
void vnorm(const double *xin,double *xout,int N);
inline double min(double a,double b);

#endif // LPPROJECTOR_H_

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
