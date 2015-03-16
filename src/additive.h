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

#ifndef JUBAL_ADDITIVE_H
#define JUBAL_ADDITIVE_H

/*
    Efficient Alias-Free Implementation of an Additive Synthesizer

    An additive synthesizer produces sound by summing together differently
    weighted sine waves of various harmonics of the fundamental. The weight
    of each sine wave determines the intesnity of that harmonic in the
    sounds spectrum.

    To avoid aliasing, sine waves with a frequency above fs/2 should not be
    mixed into the output. This gives us a max harmonic, km, of
    km = int(fs/(2*f0)). Additional aliasing concerns are possible interpolation
    errors, if implemented with a sine table.

    For efficiency reasons, there are two major concerns. We need an efficient
    way to compute sin(kx) given sin(x) and sin(x+n) given sin(x). sin(kx) will
    give us the k-th harmonic, and as the additive synthesizer is built on
    summing many harmonics it is very important to be able to do that efficiently,
    especially if many harmonics are involved.

    Direct Formula Approach:
    sin(kx) from sin(x)

    From: Weisstein, Eric W. "Multiple-Angle Formulas." From MathWorld--A Wolfram Web Resource. http://mathworld.wolfram.com/Multiple-AngleFormulas.html

    The most appropriate formula is the recursion:
    sin(nx) = 2sin[(n-1)x]cos(x) - sin[(n-2)x]
    cos(nx) = 2cos[(n-1)x]cos(x) - cos[(n-2)x]

    Base case:
        2 base cases are required:
        sin[(n-1)x] for n=2 is sin(x)
        sin[(n-2)x] for n=2 is sin(0x) = sin(0) = 0
        for n = 1, the value of sin(1x) = sin(x), which is given
        for n = 0, the value of sin(0x) = sin(0) = 0

    This allows us to compute sin(kx) from sin(x) with:
        1 variable: cos(x)
            possibly 2*cosx
        2 recursion variables: sin[(n-1)x], sin[(n-2)x] 
        2 multiplies
            can be reduced to 1 by saving 2*cosx in a constant
        1 subtraction

    sin(x+m) from sin(x)

    in this case m is the time taken by 1 sample, m = 1/fs

    Weisstein, Eric W. "Trigonometric Addition Formulas." From MathWorld--A Wolfram Web Resource. http://mathworld.wolfram.com/TrigonometricAdditionFormulas.html

    sin(x+m) = sin(x)cos(m) + sin(m)cos(x) 
    cos(x+m) = cos(x)cos(m) - sin(x)sin(m)
    base case: 
        sin(0) = 0
        cos(0) = 1
    
    This allows us to compute sin(x+m) from sin(x) with:
        2 constant: sin(m) = sin(1/fs), cos(m) = cos(1/fs)
        2 state variables: sin(x), cos(x)
        4 multiplies
        1 add
        1 subtract

    
    We now have separate horizontal (sin(x+m)) and vertical (sin(kx))
    efficiencies. The next question is how we can combine these to make an
    overall more efficient synthesizer.

    First implementation:
    Evaluate each harmonic on each sample, then recompute sin(x) and reeval.
    for each sample:
        for each harmonic:
            compute

    Horizontal
        Vertical

        K - number of harmonics
        M - number of samples
        Cost per sample: 
            K multiplies
            K subtractions
        Cost per buffer of harmonics:
            M*K multiplies
            M*K subtractions
        Cost per buffer of computing sin(x+m)
            M*4 multiplies
            M adds
            M subtractions
        Total Cost:
            M*(K+4) multiplies
            M adds
            M*(K+1) subtractions
            Cycles: (wrong)  3*M*(K+4) + M + M*(K+1)
            = 3MK + 12M + M + MK + M
            = 14M + 4MK
            = M(14+4K)

            O(MK)

        Worst Case:
            20Hz, K = 48000/(2*20) = 1200
            M = 1024 = 1024/48000 s = 0.0213 seconds
            1024(14+4*1200) = 4929536 cycles
            1 GHz machine: 4929536/10^9 = 4929536/10^9 = .00493
            int(0.0213/.00493) = 4
        
        In the worst case we could only evaluate 4 additive synths in the period!
    
    Implementation Idea:  
    Vertical
        Horizontal
    for each harmonic:
        for each sample:
            compute

    K - number of harmonics
    M - number of samples
    Cost per harmonic:
        M*4 multiplies
        M adds
        M subtractions
        

    Better Idea:
        We need to use DP to take advantage of the redundancies

    Essentially we can model the problem as filling a 2d array with the following
    pattern:


    At time t:

    m = 1/fs

    sin(t)  sin(t+m)    sin(t+2m)   sin(t+3m)   ... sin(t+Mm)
    sin(2t) sin(2t+m)   sin(2t+2m)  sin(2t+3m)  ... sin(2t+Mm)
    ...     ...         ...         ...         ... ...
    sin(Kt) sin(Kt+m)   sin(Kt+2m)  sin(Kt+3m)  ... sin(Kt+Mm)

    Let T(k, n) = sin(k*t + n*m), the entry at row k and column n in the table.
    the trig identites give us a few recurrences

    As a reminder, here are the relevant trig identities:
    sin(nx) = 2sin[(n-1)x]cos(x) - sin[(n-2)x]
    cos(nx) = 2cos[(n-1)x]cos(x) - cos[(n-2)x]
    sin(x+m) = sin(x)cos(m) + sin(m)cos(x) 
    cos(x+m) = cos(x)cos(m) - sin(x)sin(m)

    
    T(k, 0) = 2*T(k-1, 0)*cos(x) - T(k-2, 0)
    T(0, n) = T(0)*cos
*/

#endif
