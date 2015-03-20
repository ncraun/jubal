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

#ifndef JUBAL_NODES_H
#define JUBAL_NODES_H

#include "autil.h"

#define NUM_VOICES 128
#define DPW_ORDER 6
#define MAX_HARM 64

enum {SUB_CUTOFF=0, SUB_RES=1};

struct sine_voice {
    float freq;
    float time;
    int note;
    bool on;
};

struct sine {
    /* Voice stealing using array and index as a */
    struct sine_voice voices[NUM_VOICES];
    /* Round robin voice stealing */
    int voice_assign;
    float vib_rate; /* Hz */
    /* Remember to convert vib_depth from cents to Hz (as pitch is logarithmic) */
    float vib_depth; /* Cents */

    float trem_rate;
    float trem_depth;
    float atk_time;
    float sus_level;
    float rls_time;
    float port_time;
    float port_depth;

    /* Additive Synth */
    float a_freq;
    bool a_on;
    float a_amp[MAX_HARM];
    float a_time;
};

struct fm_operator {
	float m; // freq multiple
	ADSR adsr;
};

struct fm_voice {
    uint8_t on;
	uint8_t key_off;
	uint8_t note;

    float time;
    float c_freq; 
	float a_off;
	float t_off;
};

struct fm {

	int last_assign;

    float c;
    float m;

	/* Vibrato */
	float v_f; /* Hz */
	float v_d; /* Cents */

	/* Tremolo */
	float tr_f; /* Hz */
	float tr_d; /* dB */

	ADSR adsr;

	struct fm_voice voices[NUM_VOICES];
	struct fm_operator ops[6];
};

struct dpw {
	
	//uint8_t initialized;
	/* TODO: Feedback param, should it be double?*/
	float prev;
};

struct svf {
	double hp_delay;
	double bp_delay;
	double lp_delay;
};

struct ladder {
	/* TODO: These are feedback params should they be double instead of float? */
	float s1;
	float s2;
	float s3;
	float s4;
};

struct sub_voice {
	uint8_t on;			// 1
	uint8_t key_off;	// 1
	uint8_t note;		// 1
						// 1 byte wasted
	float freq;			// 4
	float a_off;		// 4
	float pulse_width;	// 4
	float time;		// 4
	// Normalized phase in range [0 to 1)
	// Should phasor be a double?
	float phasor;		// 4

	struct dpw dpw_d1;	// 4?
	struct dpw dpw_d2;	// 4?
	struct ladder ladder;	// 4*4 = 16
};

struct sub {
	uint8_t last_assign;	// 1
							// 3 bytes wasted
	//float cutoff;		// 4
	// = tan(cutoff/(2.0*sample_rate)), Remember to update when cutoff or sample_rate changes
	float prewarp_cutoff_g;	// 4
	float res;		// 4
	
	float pwm_rate;	// 8
	float pwm_amount;	// 8
	struct sub_voice voices[NUM_VOICES]; // sub_voice*8
};

union NodeData {
	struct sub sub;
	struct fm fm;
};

#endif
