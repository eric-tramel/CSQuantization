/*
  LpNormalizer.cpp
  
  code definitions for class LpNormalizer
*/

/* This file is part of BPDQ Toolbox (Basis Pursuit DeQuantizer)
 * 
 * Copyright (C) 2009, the BPDQ Team (see the file AUTHORS distributed with
 * this library) (See the notice at the end of the file.) */

#include "LpNormalizer.h"
#include "mex.h"

// removes sign from y_n, stores sign in ysign
void LpNormalizer::remove_sign(){
  bool tmp;
  for(int k=0;k<N;k++){
    tmp=(y_n[k]>0);
    y_n[k]=(tmp ? y_n[k] : -y_n[k]);
    ysign[k]=tmp;
  }
}

// modifies x_out according to stored ysign
void LpNormalizer::replace_sign(){
  for (int k=0;k<N;k++){
    // x_out[k]=ysign[k] ? x_out[k]:-x_out[k];    
    if (!ysign[k])
      x_out[k]=-x_out[k];
  }
}

// stores y_n = (y_in-c)/r
void LpNormalizer::center(const double *y_in, const double *c, const double r){
  double ir=1.0/r; // possibly useless optimization 
  for (int k=0;k<N;k++){
    y_n[k]=(y_in[k]-c[k])*ir;
  }
}

// sets x_out = x_out*r + c
void LpNormalizer::uncenter(const double *c, const double r){
  for (int k=0;k<N;k++){
    x_out[k]=x_out[k]*r + c[k];
  }
}

// (re)inititializes for calculating with newN length vectors
void LpNormalizer::init(int newN){
  if (N==newN && initialized)
    return;
  delete_all();
  N=newN;
  allocate_all();
}

// release allocated memory
void LpNormalizer::delete_all(){
  if (initialized){
    delete [] y_n;
    delete [] ysign;
    delete [] x_out;
  }
  initialized=false;
}

void LpNormalizer::allocate_all(){
  y_n=new double[N];
  x_out=new double[N];
  ysign=new bool[N];
  initialized=true;
}

LpNormalizer::LpNormalizer(int newN){
  init(newN);
}

LpNormalizer::LpNormalizer(){
  N=0;
  initialized=false;
}

LpNormalizer::~LpNormalizer(){
  if (initialized)
    delete_all();
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
