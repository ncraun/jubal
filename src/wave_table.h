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

#ifndef JUBAL_WAVE_TABLE_H
#define JUBAL_WAVE_TABLE_H

#include <stdbool.h>

/*
	Wave Tables

	I have evaluated the DPW, BLIT and BLEP methods, but I have decided the wavetable
	method is the one to use. The main disadvantage of the wavetable method is its
	memory use, but that can be curbed in a few ways listed below. The wavetable is
	easier to implement, easier to understand, more extendable and higher quality
	than the other methods.
*/

struct wave_table_entry {
	float max_freq;
	// int offset;
	float data[4096];
	int len;
};

#define NUM_WAVE_TABLES 20
#define WAVE_TABLES_LENGTH 2048*NUM_WAVE_TABLES 

struct wave_table {
	int over_sample;
	struct wave_table_entry tables[NUM_WAVE_TABLES];
	//float data[WAVE_TABLES_LENGTH];
};

void fill_saw_wave_table(struct wave_table *saw, int sample_rate);

struct wave_table G_SAW_TABLE;

#endif
