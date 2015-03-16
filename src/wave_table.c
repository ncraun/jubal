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

#include <assert.h>
#include <math.h>
#include <string.h>
#include "const.h"

#include "wave_table.h"

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

struct wave_table G_SAW_TABLE;

static float
sawf_8ve(float p, float f, int sample_rate)
{
	// What octave is the freq in?
	float fm = 27.5;
	while (fm < f) {
		fm *= 2;
	}

	float v = 0.0f;
	int km = sample_rate/(2.0*fm);

	float t = p/f;
	float a = f*2.0*PI*t;
	for (int k = 1; k <= km; k++) {
		v += sinf(k*a)/k;
	}
	return 2.0f/PI * v;
}

/*
	Table size:

	I is The Increment is the number of samples skipped when reading from the table
	N is the number of entries in the table
	fs is the sample rate
	fm is the maximum fundamental frequency for the table
	km is the maximum harmonic number in the table
	ft is the frequency of the table

	ft = fs/N
	
	The increment should be less than one-half the period of the highest harmonic in
	the table

	I = (N*f)/(fs)

	The frequency of the highest harmonic is fm*km
	One half that period is 1/(2*km*fm)	

	I < 1/(2*ft*km)

	N*(fm/fs) < 1/(2*fm*km)
	N < fs/fm * 1/(2*fm*km)
	N < fs/(2*km*fm^2)
	
*/
void
fill_saw_wave_table(struct wave_table *saw, int sample_rate)
{
	// memset(saw->data, 0, sizeof(float)*WAVE_TABLES_LENGTH);
	saw->over_sample = 2;
	int curr_offset = 0;
	int N = 2048*saw->over_sample;//512;//2048;
	for (int t = 0; t < NUM_WAVE_TABLES; t++) {
		struct wave_table_entry *entry = &(saw->tables[t]);

		float fm = powf(2.0, t+1)*27.5;
		int km_from_fm = (float)sample_rate/(2*fm);
		int km_from_N = N/2;
		//int km_from_increment = ;

		int km = km_from_fm < km_from_N ? km_from_fm : km_from_N;


		entry->max_freq = fm;
		entry->len = N;
		
		float Im = N*(fm/sample_rate);
		// printf("Octave: %d km: %d fm: %f Im: %f Len: %d Off: %d\n", t, km, entry->max_freq, Im, entry->len, entry->offset);

		for (int n = 0; n < N; n++) {
			float p = (float)n/(float)N;
			float f = 1/(2.0*PI);
			entry->data[n] = sawf_8ve(p, f, sample_rate);
		}

		curr_offset += N;

		//N /= 2;
	}

	for (int t = 0; t < NUM_WAVE_TABLES; t++) {
		struct wave_table_entry *entry = &(saw->tables[t]);
		for (int n = 0; n < entry->len; n++) {
			printf("%f\n", entry->data[n]);
		}
	}

	printf("%zu\n", sizeof(*saw));
}
