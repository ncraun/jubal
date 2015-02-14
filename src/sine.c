/*
    Copyright (C) 2014-2015 Nate Craun

    This file is part of Jubal.

    Jubal is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    Jubal is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with Jubal.  If not, see <http://www.gnu.org/licenses/>.
*/

#include <stdlib.h>
#include <math.h>
#include <string.h>

/*
    Fast/efficient sin generator.
    The key to this improvement is taking advantage of the fact that the
    arguments to sin() increase linearly by a fixed amount.

    In general, for frequency f and time t, the generator should output
    sin(2*pi*f*t), by definition of sin()

    However, in our computer system time is not continuous, it is discrete.
    The smallest unit of our discrete time is the sample. For a system with
    sampling rate 1/S, there are S samples per second, and each sample has a
    duration of 1/S seconds. 

    The generator must produce an output for each sample.

    let sin(Ai) be the value produced by the generator for sample i.

    Ai = 2*pi*f*t by definition of our sin generator.
    t = i/S, as there are S samples in one second.
    so, Ai = 2*pi*f*i/S

    let B be the difference between consecutive values of A.
    So, B = Ai - A_{i-1} = 2*pi*f*i/S - 2*pi*f*(i-1)/S = 2*pi*f*(i-(i-1))/S
    B = 2*pi*f*1/S

    Now to use this information to efficiently compute sin(Ai) given sin(A{i-1}) 
    Ai = A{i-1} + B, so sin(Ai) = sin(A{i-1} + B)

    Using the trig identity: sin(C + D) = sin(C)cos(D) + sin(D)cos(C)
    we can rewrite:

    sin(A{i-1} + B) = sin(A{i-1})cos(B) + sin(B)cos(A{i-1})
    B only depends on f and S, so B is constant for constant f and S.
    So, cos(B) and sin(B) are also constant for constant f and S, and
    will only need to be recomputed when f or S change.

    sin(A{i-1}) will have already been computed for the previous sample,
    so there is no need to recompute it with an expensive call to sin().

    But now, we need to compute cos(Ai) from cos(A{i-1})
    Using the trig identity: cos(C + D) = cos(C)cos(D) - sin(C)sin(D) 
    we can rewrite:
    cos(A{i-1} + B) = cos(A{i-1})cos(B) - sin(A{i-1})sin(B)

    In order to compute sin(Ai), we must store:
        sin(B)
        cos(B)
        sin(A{i-1})
        cos(A{i-1})

    So, we can compute successive sin values with 2 add, 4 multiplies and
    4 state values.

    Initial values:
    at sample number 0: 
        sin(0) = 0
        cos(0) = 1

    Constant values:
        sin(B) = sin(2*pi*f/S)
        cos(B) = cos(2*pi*f/S)
*/

/*
    Linear interpolation of a wavetable.
    
    given Frequency f, table length N, table W, sample rate R.

    delta t = f*N/R
    the amount of time between each sample in the table

    For a given output at a time
    Map that value to a table position.

    for output sample number n,

    x = n * f*(N/R) mod N
    s = |_ x _|

    Linear Interpolation:

    y[n] = (W[s+1] - W[s])*(x - s) - W[s]
*/

