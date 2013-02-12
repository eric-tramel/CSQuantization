Matlab code implementing and demonstrating 
Basis Pursuit Dequantization, as described in

"Dequantizing Compressed Sensing : When Oversampling and Non-Gaussian
Constraints Combine" by Laurent Jacques, David Hammond, M. Jalal
Fadili ( submitted to IEEE Transactions on Signal Processing)

To use, first compile the cpp code implementing the mex version of
lp-ball projection by running 
  >> bpdq_make_mex
	
Then, run 
  >> bpdq_setpath
to add all subpaths to the matlab path.

Next, try running the demo's

  >> bpdq_demo_1d  
  >> bpdq_demo_2d 

which each run (a single case) of the numerical experiments described in 
Jacques et al.

Enjoy!

Directory to hold cleaned up dequantization code, for public release

Structure:

	bpdq_toolbox/ - root directory, contains demo files 
	bpdq_toolbox/one_d/ - code for one dimensional simulations
	bpdq_toolbox/two_d/ - code for two dimensional simulations
	bpdq_toolbox/common/ - code common between one and two dimensions
	bpdq_toolbox/mex/ -  mex source code  

	bpdq_toolbox/doc/ - html documentation, view
	bpdq_toolbox/doc/index.html with your favourite web browser


License:

The BPDQ Toolbox is a GPL Matlab/C/C++ Library. 

The BPDQ Toolbox is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

The BPDQ Toolbox is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with The BPDQ Toolbox.  If not, see <http://www.gnu.org/licenses/>.
