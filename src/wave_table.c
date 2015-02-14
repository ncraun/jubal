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

#include <math.h>
#include <string.h>
#include "const.h"

void make_sine_table(float *tbl, size_t size, int sampleRate) {
    // Fill table with 1 period of the waveform
    // 0 to 2pi
    
    for (size_t i = 0; i < size; i++) {
        tbl[i] = sinf(2.0*PI*(i+1.0)/(float)size);
    }
}

// Band limited
/*
void make_square_table(float *tbl, size_t size, int sampleRate) {

    memset(tbl, 0, sizeof(float)*size);

    for (size_t i = 0; i < size; i++) {
        for (int k = 1; (2*k-1) < sampleRate/2; k++) {
            float h = (2.0*k-1.0);
            float t = (i+1.0)/(float)size;
            tbl[i] += sin(2.0*PI*t*h)/h;
        }
    }
}
*/

/* 
    Create square waves without aliasing by using bandwidth limited
    wavetables. Evaluate only the terms of the Fourier Series that will
    not be adding a frequency that would introduce aliasing, 1 table per
    octave.

    See the digital-alias page in jubal-doc for more information on how
    this works. 

    Slowest implementation.
    Can be made faster by:
        Reusing computed values from lower km for higher km octaves
        Using the fast incrementing sin generator.

    Table is 2d array of size octaves*buffers.
*/
void
make_square_tables(unsigned int sampleRate, float *table, int octaves, int size) {
    memset(table, 0, octaves*size*sizeof(float));
    for (int octave = 0; octave < octaves; octave++) {
        double fm = 27.5 * (1<<(octave+1));
        int km = (int)(sampleRate/(2*fm));
        for (int n = 0; n < size; n++) {
            //for (int k = 1; k <= km; k += 2) {
            for (int k = 1; k <= km; k++) {
                float t = (n + 1.0)/(float)size;
                table[octave*size+n] += sinf(k*2.0*PI*t)/(float)k; 
            }
            table[octave*size+n] *= 4.0/PI;
        }
    }
}

