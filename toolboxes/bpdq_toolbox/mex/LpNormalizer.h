/*
  LpNormalizer.h
  
  Defines class  LpNormalizer
*/

/* This file is part of BPDQ Toolbox (Basis Pursuit DeQuantizer)
 * 
 * Copyright (C) 2009, the BPDQ Team (see the file AUTHORS distributed with
 * this library) (See the notice at the end of the file.) */

#ifndef LPNORMALIZER_H_
#define LPNORMALIZER_H_

/*
  LpNormalizer

  This class provides functionality to center, normalize and 
  remove the signs (i.e. reflect into positive orthant) from
  given input. 

  It stores y_n, which should be normalized input, and
  x_out, which is array that should be passed to LpProjector
  to store output of Lp Projection.

  center computes y_n = (y_in - c)/r
  remove_sign computes y_n = |y_n|
  replace_sign, uncenter undoes these operations

  Before use, the class must be initialized with the length of 
  the vectors
  myLpNormalizer->init(N)
  Repeated initialization with the same N will not reallocate memory, while
  initialization with a subsequently different N will reallocate memory.
  This is so the class can be a persistant object and incur no performance
  penalty when repeatedly initialized (as done in the matlab gateway)

  Intended use is 
  myLpNormalizer->center(y_in,c,r)
  myLpNormalizer->remove_sign();
  
  call other functions taking myLpNormalizer->y_n as input
  and myLpNormalizer->x_out as storage for output

  myLpNormalizer->replace_sign()
  myLpNormalizer->uncenter(c,r)
  
  then myLpNormalizer->x_out contains output that has been 
  uncentered, with the signs properly replaced.

  It is an error to call replace_sign before remove_sign
*/
  
class LpNormalizer{
public:
  void remove_sign();
  void replace_sign();
  void center(const double *y_in,const double *c, const double r);
  void uncenter(const double *c, const double r);
 
  void init(int newN);  
  double *y_n;
  double *x_out;

  void delete_all();

  LpNormalizer();
  LpNormalizer(int newN);
  ~LpNormalizer();

 private:
  // reallocate memory 
  void allocate_all();
  // indicate whether memory has been allocated
  bool initialized;
  int N;
  // ysign used to store sign for remove_sign, replace_sign
  bool *ysign;
};

#endif // LPNORMALIZER_H_

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
