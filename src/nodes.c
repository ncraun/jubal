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
#include <string.h>
#include <stdio.h>
#include <math.h>
#include <assert.h>
#include "const.h"
#include "engine.h"
#include "autil.h"
#include "wave_table.h"
#include "blep.h"

#define NUM_VOICES 8

#define DPW_ORDER 6

/*
    Voice Stealing:

    Each node has a maximum polyphony, which is the number of different notes that can be
    played simulatenously. The maximum polyphony is determined by the number of voices the
    synth has. In a hardware synth, each voice would consist of dedicated circuitry, but
    in Jubal, a voice is just a set of parameters unique to the sound generating process.
    Certain parameters are shared across all voices (such as ADSR, waveform, etc). A voice
    consists of a structure of parameters specific to that voice, such as current note,
    current time in the cycle, etc.

    When a note on (or note off) message is sent to a node there must be a policy
    for which voice will be allocated to handle the new note.

    Cometimes the node is asked to play more notes at the same time than it has voices
    to play with. In this case, a common aproach is "voice stealing", where a previous
    note will cut out, and that voice will be used to play the new note. There must be
    a policy to govern this as well.

    In Jubal, it is up to each node to implement this, but to help promote cohesiveness
    the following governance policies are suggested:

    TODO: determine policies
*/

struct sine_voice {
    float freq;
    float time;
    int note;
    bool on;
};

#define MAX_HARM 64

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

static void*
def_sine_make(int sample_rate)
{
    struct sine *data = malloc(sizeof(struct sine));
    for (int i = 0; i < NUM_VOICES; i++) {
        data->voices[i].on = false;
    }
    data->voice_assign = 0;
    data->vib_rate = 3.0f;
    data->vib_depth = 10.0f;

    data->trem_rate = 440.0f;
    data->trem_depth = .08;

    data->atk_time = 0.1f;
    data->sus_level = 1.0;
    data->rls_time = .5;

    data->port_time = .3f;
    data->port_depth = 20.0f;

    data->a_on = false;
    data->a_freq = 0.0f;
    for (int k = 1; k < MAX_HARM; k++) {
        data->a_amp[k] = (4.0/PI) * (1.0/(float)k);
    }

//    data->a_amp[1] = 1.0f;
    data->a_time = 0;

    return data;
}

static void
def_sine_free(void *data)
{
    free(data);
}

/*
    TODO: Move this function somewhere it makes more sense to be.
*/
float midi_freq(int note)
{
    return powf(2.0, (note-69.0)/12.0)*440.0;
}

static float
sawf(float a)
{
    float v = 0.0f;
    int km = 100;
    for (int k = 1; k <= km; k++) {
       v += sinf(k*a)/k;
    }
    return 4.0f/PI*v;
}

static float
sawf_bl(float p, float f, int sample_rate)
{
	float v = 0.0f;
	int km = sample_rate/(2.0*f);
	float t = p/f;
	float a = f*2.0*PI*t;
	for (int k = 1; k <= km; k++) {
		v += sinf(k*a)/k;
	}
	return 4.0f/PI * v;
}


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

	printf("km = %d fm = %f\n", km, fm);

	float t = p/f;
	float a = f*2.0*PI*t;
	for (int k = 1; k <= km; k++) {
		v += sinf(k*a)/k;
	}
	return 4.0f/PI * v;
}


static void
def_sine_eval(void *instance, float *left_input, float *right_input, float *left_output, float *right_output, int len, int sample_rate, int num_msg, RtNodeMsg *msgs)
{
 //   printf("Call\n");

    struct sine *data = instance;

    // m5
   double p, t, f, w, swp, cwp, st, ct, ckx, ckx_1, skx, skx_1, skx_2, cx, cx_1, sx, sx_1, sum, cx2;

    p = 1.0f/sample_rate;
    t = data->a_time;
    f = data->a_freq;
    w = 2 * PI * f;
    
    swp = sinf(w*p);
    cwp = cosf(w*p);
    st = sinf(t);
    ct = cosf(t);

    // TODO: change this to work when frequency changes
    // To change frequency, simply update cwp and swp in the loop
    // cwp = cosf(w2), swp = sinf(w2)
    // TODO: suppot changing amp[k] on a sample-by-sample basis
    // Common case:
    // Changing f often
    // Changing amp[k] not so much.

    /*
        sin(-x) = -sin(x)

        k
        -2  sin(-2x) = -sin(2x) = -(2sin(x)cos(x) - sin(0)) = -(2sin(x)cos(x))
        -1  sin(-x) = -sin(x)
         0  sin(0x) = 0.0f
         1  sin(x) = sx 

        computing 2*cos(x)    
    */
    int curr_msg = 0;

    cx = cosf(w*t);
    sx = sinf(w*t);
    for (int n = 0; n < len; n++) {

        while (curr_msg < num_msg && msgs[curr_msg].time == n) {
            if (msgs[curr_msg].type == NOTE_ON) {
                data->a_on = true;
                data->a_freq = midi_freq(msgs[curr_msg].note);
                f = data->a_freq;
                w = 2.0f*PI*f;
                swp = sinf(w*p);
                cwp = cosf(w*p);
                printf("anote_on\n");
            }
            else if (msgs[curr_msg].type == NOTE_OFF) {
                data->a_on = false;
                printf("anote_off\n");
            }
            curr_msg++;
        }

        if (!data->a_on) {
            left_output[n] = 0.0f;
            right_output[n] = 0.0f;  
            continue;
        }

        skx = 0.0f;
        skx_1 = -sx;
        skx_2 = -2.0f*cx*sx;
        sum = 0.0f;

        float x = 2.0f*data->a_freq*PI*(data->a_time + n*p);

        cx2 = 2.0f*cx;

        /* Max k avoiding aliasing. */
        int km = (int)(float)sample_rate/(2*data->a_freq);
        /* Max k for iteration */
        int K = km < MAX_HARM ? km : MAX_HARM;
        //printf("K: %d Freq: %f\n", K, data->a_freq*K);

        for (int k = 0; k < MAX_HARM; k++) {
            sum += data->a_amp[k]*skx; 

            skx_2 = skx_1;
            skx_1 = skx;
            skx = cx2*skx_1 - skx_2;
        }

        left_output[n] = sum;
        right_output[n] = sum;

        cx_1 = cx;
        sx_1 = sx;
        cx = cx_1*cwp - sx_1*swp;
        sx = sx_1*cwp + cx_1*swp;

    }

    data->a_time += len*p;

/*
    for (int n = 0; n < len; n++) {
        printf("%f ", left_output[n]);
    }
    printf("\n");
*/
}

static void
def_sine2_eval(void *instance, float *left_input, float *right_input, float *left_output, float *right_output, int len, int sample_rate, int num_msg, RtNodeMsg *msgs)
{
    //printf("sine eval\n");
    struct sine *data = instance;
    int curr_msg = 0;
    float sample_time = 1.0/(float)sample_rate;
    for (int f = 0; f < len; f++) {
        while (curr_msg < num_msg && msgs[curr_msg].time == f) {
            printf("msg\n");
            if (msgs[curr_msg].type == NOTE_ON) {

                bool steal = true;
                for (int i = 0; i < NUM_VOICES; i++) {
                    if (data->voices[i].on == false) {
                        data->voices[i].on = true;
                        data->voices[i].freq = midi_freq(msgs[curr_msg].note);
                        data->voices[i].time = 0;
                        data->voices[i].note = msgs[curr_msg].note;
                        steal = false;
                        break;
                    }
                }

                if (steal) {
                    int curr = data->voice_assign;

                    data->voices[curr].on = true;
                    data->voices[curr].freq = midi_freq(msgs[curr_msg].note);
                    data->voices[curr].note = msgs[curr_msg].note;
                    data->voices[curr].time = 0;

                    printf("Note On: %d\n", msgs[curr_msg].note);
                    printf("freq -> %f\n", data->voices[curr].freq);

                    data->voice_assign = (data->voice_assign + 1) % NUM_VOICES;
                }
            }
            else if (msgs[curr_msg].type == NOTE_OFF) {
                for (int i = 0; i < NUM_VOICES; i++) {
                    if (data->voices[i].note == msgs[curr_msg].note) {
                        data->voices[i].on = false;
                    }
                }
            }
            else if (msgs[curr_msg].type == CC) {
            }
            curr_msg++;
        }

        float mix = 1.0f/NUM_VOICES;
        for (int i = 0; i < NUM_VOICES; i++) {
            float val = 0.0;
            if (data->voices[i].on) {

#if 0
                float t = data->voices[i].time;

                float amp = data->sus_level;
                // Attack
                if (t < data->atk_time) {
                    float rise = (data->sus_level-0.0f);
                    float run = (data->atk_time-0.0f);
                    float slope = rise/run;
                    amp = slope*t;
                }
                // Release
                else {
                }

                float vib = data->vib_depth*cosf(2.0f*PI*data->vib_rate*t);
                /* The Vib sinewave is in cents, convert vib to Hz and add to f. */
                /* Or... should we convert beforehand, and then make the Vib sinewave be in Hz? */
                float vib_hz = data->voices[i].freq*powf(2.0f, vib/1200.0f);
                float f = data->voices[i].freq + vib_hz; 
                float f_c = 1200*log2f(f/data->voices[i].freq); // Let's see how many cents the output is...
                static float f_c_min = 100000.0f;
                static float f_c_max = -1000000.0f; 
                if (f_c < f_c_min) f_c_min = f_c;
                if (f_c > f_c_max) f_c_max = f_c;
                printf("Out pitch: %0.8f Cents Min: %0.8f Max: %0.8f\n", f_c, f_c_min, f_c_max);

                //printf("%0.8f\n", f);
                if (t < data->port_time) {
/*
                    float rise = data->port_depth;
                    float run = data->port_time;
                    float slope = rise/run;

                    f = slope*t + (data->voices[i].freq-data->port_depth);
*/
                }

                //float trem = data->trem_depth*sinf(2.0f*PI*data->trem_rate*t);
                //val = amp*mix*cosf(2*PI*f*t);

            #endif

/* Why isn't this working? */
                float a_o = mix;
                float f_o = data->voices[i].freq;
                float v_d = data->vib_depth;
                float v_r = data->vib_rate;
                float t = data->voices[i].time; 


                float vs_cents = v_d*sinf(2.0*PI*v_r*t);
                float out_freq = f_o * powf(2.0f, vs_cents/1200.0f);

                //val = a_o*cosf(2.0f*PI*(out_freq)*t);

                float w_c = 2.0f*PI*data->voices[i].freq; 
                float A_v = f_o - out_freq;
                float w_v = 2.0f*PI*data->vib_rate;
            
                //val = .5*sinf(w_c*t-A_v/w_v*cos(w_v*t));

                val = .5*sawf(w_c*t+A_v*sinf(w_v*t));
            
                // Scale
                float smin = -0.99710417;
                float smax =  -0.75658751;
                float peak_amp = smax - smin;

                //val = 2*(val/peak_amp) + 7.3;

                //val *= .5;

                //printf("%0.8f\n", data->voices[i].freq); 
            
                data->voices[i].time += sample_time;
                //printf("%0.8f %0.8f %0.8f\n", t, vib, val);
            }

            left_output[f] += val;
            right_output[f] += val;
        }
    }
}


static void*
def_disabled_make(int sample_rate)
{
    return NULL;
}

static void
def_disabled_init()
{
}

static void
def_disabled_free(void *instance)
{
}

static void
def_disabled_eval(void *def_data, void *instance, float *left_input, float *right_input, float *left_output, float *right_output, int len, int sample_rate, int num_msg, RtNodeMsg *msgs)
{
    // Noise
    //memset(left_output, 0, sizeof(float)*len);
    //memset(right_output, 0, sizeof(float)*len);
/*
    for (int i = 0; i < len; i++) {
        float r = rand();
        r /= RAND_MAX;
        r *= 2.0;
        r -= 1.0;
        r *= .05;
        left_output[i] = r;
        right_output[i] = r;
    }
*/
}

static void
def_disabled_read(int size, uint8_t *data)
{
}

static void
def_disabled_write(int size, uint8_t *data)
{
}

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
	bool initialized;
	double prev;
};

struct svf {
	double hp_delay;
	double bp_delay;
	double lp_delay;
};

struct sub_voice {
	uint8_t on;
	uint8_t key_off;
	uint8_t note;

	double phasor; // Normalized phase in range [0 to 1)

	float freq;
	float a_off;
	float t_off;

	float poly_blep_fir[4];
	double prev_dpw[DPW_ORDER];
	bool dpw_initialized;
	int dpw_delayed;
	struct dpw dpw_d1;
	struct dpw dpw_d2;
	struct svf svf;
};

struct sub {
	int last_assign;
	double cutoff;
	double prewarp_cutoff_g; // = tan(cutoff/(2.0*sample_rate)), Remember to update when cutoff or sample_rate changes
	double res;
	struct sub_voice voices[NUM_VOICES];
};

static void*
def_fm_make(int sample_rate)
{
	printf("%d\n", sizeof(struct fm));
    struct fm *data = malloc(sizeof(struct fm));

	data->c = 1.0f;
	data->m = 3.0f/5.0f;
	data->adsr.atk_t = .1;
	data->adsr.atk_l = 1.0;
	data->adsr.dec_t = .5;
	data->adsr.sus_l = 0.75;
	data->adsr.rel_t = .01;//.8;
	data->v_f = 6.0f;
	data->v_d = 1.0f;
	data->tr_f = 4.0f;
	data->tr_d = 0.0f;//0.1f;

	for (int i = 0; i < NUM_VOICES; i++) {
		data->voices[i].time = 0.0f;
		data->voices[i].on = false;
		data->voices[i].c_freq = 0.0f;
		data->voices[i].a_off = 0.0f;
		data->voices[i].key_off = 0.0f;
		data->voices[i].t_off = 0.0f;
	}
	data->last_assign = 0;

    return data;
}

static void def_fm_eval(void *def_data, void *instance, float *l_in, float *r_in, float *l_out, float *r_out, int len, int sample_rate, int num_msg, RtNodeMsg *msgs)
{
    struct fm *data = instance;
    float p = 1.0f/(float)sample_rate;

    int curr_msg = 0;
    for (int n = 0; n < len; n++) {

        while (curr_msg < num_msg && msgs[curr_msg].time == n) {
            if (msgs[curr_msg].type == NOTE_ON) {

				struct fm_voice *voice = &data->voices[data->last_assign];

				bool found = false;
				/* What to do when we receive a NOTE_ON for a note that is already on? */
				/* Should we retrigger, or just layer? right now retrigger is causing click... */
				/* First see if there's a voice already playing that note. */
				for (int i = 0; i < NUM_VOICES; i++) {
					if (data->voices[i].note == msgs[curr_msg].note) {
						voice = &data->voices[i];
						found = true;
						break;
					}
				}

				/* Then search for an free voice */
				if (!found) {
					for (int i = 0; i < NUM_VOICES; i++) {
						if (!(data->voices[i].on)) {
							voice = &data->voices[i];
							found = true;
							break;
						}
					}	
				}

				/* Otherwise, steal one */
				if (!found) {
					voice = &data->voices[data->last_assign];
					data->last_assign = (data->last_assign + 1) % NUM_VOICES;
				}

                voice->on = true;
                voice->c_freq = midi_freq(msgs[curr_msg].note);
				voice->key_off = false;
				voice->time = 0.0f;
				voice->note = msgs[curr_msg].note;
            }
            else if (msgs[curr_msg].type == NOTE_OFF) {
				for (int i = 0; i < NUM_VOICES; i++) {
					if (data->voices[i].note == msgs[curr_msg].note) {
						data->voices[i].key_off = true;
						data->voices[i].t_off = data->voices[i].time + n*p;

						// Should NOTE_OFF turn off all notes playing that voice, or just one?
						// break;
					}
				}
            }
            curr_msg++;
        }

		for (int i = 0; i < NUM_VOICES; i++) {
			struct fm_voice *voice = data->voices + i;

			if (!voice->on) {
				continue;
			}

			// 2op fm

			float t = voice->time + n*p;

			float a = eg_lin_adsr(t, &(data->adsr), voice->key_off, voice->t_off, &voice->a_off);

			if (a < 0.0f) {
				voice->on = false;
				continue;
			}

			float w_c = 2.0f*PI*voice->c_freq;
			float w_f = w_c*(data->m/data->c);

			// Vibrato
			float v_d_hz = voice->c_freq*powf(2.0f, data->v_d/1200.0f) - voice->c_freq;
			float w_v = 2.0f*PI*data->v_f;
			float w_tr = 2.0f*PI*data->tr_f;
			float norm = 1.0f/NUM_VOICES;

			// Evaluate the operator
			// Evaluate the algorithm

			float v = norm*(data->tr_d*sinf(w_tr*t) + a)*sinf(w_c*t + sinf(w_f*t) + v_d_hz*sinf(w_v*t));

			l_out[n] += v;
			r_out[n] += v;
		}
    }

	for (int i = 0; i < NUM_VOICES; i++) {
		data->voices[i].time += len*p;
	}
	
}

static void def_fm_free(void *instance)
{
}

static void def_fm_read(int size, uint8_t *data)
{
}

static void def_fm_write(int size, uint8_t *data)
{
}

static void*
def_sub_make(int sample_rate)
{
    struct sub *data = malloc(sizeof(struct sub));

	for (int i = 0; i < NUM_VOICES; i++) {
		data->voices[i].phasor = 0.0f;
		data->voices[i].on = false;
		data->voices[i].freq = 0.0f;
		data->voices[i].a_off = 0.0f;
		data->voices[i].key_off = 0.0f;
		data->voices[i].t_off = 0.0f;
		memset(data->voices[i].poly_blep_fir, 0, 4*sizeof(float));
		data->voices[i].dpw_initialized = false;
		data->voices[i].svf.hp_delay = 0.0;
		data->voices[i].svf.bp_delay = 0.0;
		data->voices[i].svf.lp_delay = 0.0;
		for (int o = 0; o < DPW_ORDER; o++) {
			data->voices[i].prev_dpw[o] = 0.0;
		}
	}
	data->last_assign = 0;
	data->cutoff = 2.0*PI*500; // Hz
	data->prewarp_cutoff_g = tan(data->cutoff/(2.0*sample_rate));
	data->res = 0.2; // 1.0 is no resonance, 0.0 is unstable

    return data;
}

/*
	Virtual Analog Oscillator Algorithms:

	The key is to avoid aliasing (frequencies > fs/2)

	Additive
		Use the fourier series for saw, square, etc
		Don't include any harmonics kf > fs/2

	Wavetable
		For each interval (octave, for example) create a table
		the table is the fourier series for the waveform evaluated up to km
		km depends on the frequency of the table, km*f < fs/2	
		table length is at least km*f*2, maybe more to oversample
		create square from the sawtooth table

	BLIT
		The derivative of a classic waveform is an impulse (train?)
		Integrating a bandlimited impulse train creates a (kind of) bandlimited classic wave
		The bandlimited impulse train (BLIT) is stored in a table
		The BLIT is created by a summation of windowed sinc functions
		Then you integrate the BLIT and mix it with the waveform at transitions.

	BLEP
		Like a BLIT, but do the integration ahead of time
		A bandwidth limited step is the integral of a bandwidth limited impulse
		Create bandwidth limited step functions into a table
		insert it whenever there is a discontinuity

	PolyBLEP
		Like a BLEP but instead of storing a table, use a polynomial to approximate the BLEP
			during the step
		Gives approximation across the 2 sample discontinuity
		1. Make polynomial to approximate sinc in the interval
		2. Integrate it
		3. Subtract the unit step (this is because the original waveform has a unit step)
			To elaborate: W = unit_step at the time point(s)
			W + polyBLEP_added = integrated polynomial approx for sinc
			unit_step + polyBLEP_added = inte. poly. app. for sinc
			polyBLEP_added = i. p. a. for sinc - unit_step
		4. Add value to the geometric/naive waveform at the time

	DPW
		What I chose, it works
		What order should I use?


	Virtual Analog Filter Algorithms:
*/

/*
	Poly BLEP approximation functions:
1:
	3 point spline, given in "Computationally Efficient Music Synthesis", Pekonen, 2007

	t e [-N, N]
	t e [-1/fs, 1/fs]
	We use N = 1/fs = 1 sample's amount of time, before and after the step

	r(t) = (-3/14)t^4 - (4/7)t^3 + (6/7)t + 1/2,  -1 <= t <= 0 
	r(t) = (3/14)t^4 - (4/7)t^3 + (6/7)t - 1/2,    0 < t <= 1
	
*/
static float
poly_blep_spline3(float t)
{

	static float c4 = 3.0f/14.0f;
	static float c3 = 4.0f/7.0f;
	static float c1 = 6.0f/7.0f;
	static float c0 = 0.5f;

	float t3 = t*t*t;
	float t4 = t3*t;

	// 0 <= t < 1
	if (t < 1.0f && t >= 0.0f) {

		return c4*t4 - c3*t3 + c1*t + -c0;
	}

	// -1 < t < 0
	else if (t > -1.0f && t < 0.0f) {
		return -c4*t4 - c3*t3 + c1*t + c0;
	}

	return 0.0f;
}

/*
	Optimal PolyBLEP from
	Filter-Based Oscillator Algorithms for Virtual Analog Synthesis	

	d is the fractional delay of the discontinuity
	it can be computed from a normalized unipolar modulo phase counter
	(from 0 to 1) using similar triangles

	p(n)   f_0*T
	---  = -----
	T*d       T

	d = (p(n))/(f_0*T)

	where p(n) is the modulo counter value after the discontinuity
	d is the fractional delay
	T is the sampling interval
	f_0T is the step size of the counter (the phase increment)


	For a rectangular pulse wave, fractional delay w/ downward discontinuity is
	computed from the value of the modulo counter after it exceeds P, the duty
	cycle. In this case P needs to be subtracted from the modulo counter value.
	When generating a triangle wave, the fractional delays of turning points
	are obtained similarly with P = 0.5
	
	D is the delay in sampling intervals
	N = order of the polynomial	

	odd-order interpolation polynomials:
		d = D - (N - 1)/2

	If polynomial is even:
		if discontinuity is closer to 
		

	Optimized PolyBLEP function (from Pekonen)
	alias-free (minimal aliasing?) range of fundamental freqs above range used in music!

	3rd order FIR filter Pub III section II

	Polynomial spline g(t) from L consecutive polynomial segments
	g(t) = Sum from k=0 to K_l of (g_{l,k})t^k,   t e [t_{l-1}, t_l],  l = 1, 2, ..., L
	t is the time variable
	K_l is the order
	g_{l,k} are the coefficients of the polynomial in the l-th subrange

	outside of definition range t e [t_0, t_L] g(t) is def to be 0

	Create an FIR filter of order N = ceil(tL - t0 - 1)
	filter coefficients bn for n = 0, 1, ..., N are given by bn(d) = g(t0 + n + d)
	where d e [0, 1) is the fractional delay which defines the sample point within
	the sampling interval.

	Optimized BLEP residual:
		b0(d) = 0.00029*d + 0.00737*d^2 + 0.00617*d^3 + 0.03871*d^4
		b1(d) = 0.05254 + 0.18783*d + 0.22663*d^2 + 0.14955*d^3 - 0.11656*d^4
		b2(d) = -0.5 + 0.62351*d + 0.02409*d^2 - 0.31670*d^3 + 0.11656*d^4
		b3(d) = -0.05254 + 0.18839*d - 0.25816*d^2 + 0.16102*d^3 - 0.03871*d^4

	N = 3
*/

/*
	Based off of "Perceptually informed synthesis of bandlimited classic waveforms"
*/

static float
poly_blep_optimal_fir(float d, float* prev) 
{
	float d2 = d*d;
	float d3 = d2*d;
	float d4 = d2*d2;

	double b0 = 0.00029*d + 0.00737*d2 + 0.00617*d3 + 0.03871*d4;
	double b1 = 0.05254 + 0.18783*d + 0.22663*d2 + 0.14955*d3 - 0.11656*d4;
	double b2 = -0.5 + 0.62351*d + 0.02409*d2 - 0.31670*d3 + 0.11656*d4;
	double b3 = -0.05254 + 0.18839*d - 0.25816*d2 + 0.16102*d3 - 0.03871*d4;

	return b0*prev[0] + b1*prev[1] + b2*prev[2] + b3*prev[3];
}

static float
poly_blep(float t)
{
	return poly_blep_spline3(t);
}

long
factorial(long a)
{
	long f = 1;
	while (a >= 1) {
		f *= a;
		a--;
	}
	return f;
}


double
dpw_saw(struct dpw *data, double phasor, double freq, int sample_rate)
{
	double dpw_P = sample_rate/freq;
	double phase_inc = freq/sample_rate;

	if (!data->initialized) {
		// Time travel
		double prev_phasor = phasor - phase_inc;
		if (prev_phasor < 0) {
			prev_phasor += 1.0;
		}
		
		double prev_saw = 2.0*prev_phasor - 1.0;	
		double prev_saw2 = prev_saw*prev_saw;

		data->prev = prev_saw2;
		data->initialized = true;
	}

	double saw = 2.0*phasor - 1.0;	
	double saw2 = saw*saw;

	/* The original response for DPW was 1 - z^-1
	In "NEW APPROACHES TO DIGITAL SUBTRACTIVE SYNTHESIS" it is shown that
	the response (1-z^-1)*(1+z^-1) has a better rolloff */
	double dpw = saw2 - data->prev;
	data->prev = saw2;

	/*
		// This is the better scaling verion
		double dpw_scale = PI/(4.0*sin(PI/dpw_P));
	*/
	// This is the faster scaling version
	double dpw_scale = dpw_P/4.0;

	double dpw_val = dpw*dpw_scale;

	return dpw_val;
}

double
dpw_tri(struct dpw *data, double phasor, double freq, int sample_rate)
{
	double dpw_P = sample_rate/freq;
	double phase_inc = freq/sample_rate;

	if (!data->initialized) {
		// Time travel
		double prev_phasor = phasor - phase_inc;
		if (prev_phasor < 0) {
			prev_phasor += 1.0;
		}
		
		double prev_saw = 2.0*prev_phasor - 1.0;
		double prev_tri_i = prev_saw*fabs(prev_saw) - prev_saw;

		data->prev = prev_tri_i;
		data->initialized = true;
	}

	double saw = 2.0*phasor - 1.0;	
	double tri_i = saw*fabs(saw) - saw;

	/* The original response for DPW was 1 - z^-1
	In "NEW APPROACHES TO DIGITAL SUBTRACTIVE SYNTHESIS" it is shown that
	the response (1-z^-1)*(1+z^-1) has a better rolloff */

	double dpw = tri_i - data->prev;
	data->prev = tri_i;

	/* triangle scaling factor is twice the saw scaling factor */
	/* This is the better scaling function */
	/* double dpw_scale = 2.0*(PI/(4.0*sin(PI/dpw_P))); */
	// This is the faster scaling function
	// 2*P/4 = P/2
	double dpw_scale = dpw_P/2.0;
	double dpw_val = dpw*dpw_scale;

	return dpw_val;
}

double
dpw_rect(struct dpw *saw_1, struct dpw *saw_2, double phasor, double freq, int sample_rate, double pw)
{

	// A pw (pulse width) rect can be formed by summing two saw waves
	// The other saw wave must be inverted and phase shifted by the pulse width
	double phasor_2 = phasor - pw;
	if (phasor_2 < 0) {
		phasor_2 += 1.0;
	}

	double saw_1_val = dpw_saw(saw_1, phasor, freq, sample_rate);
	double saw_2_val = dpw_saw(saw_2, phasor_2, freq, sample_rate);

	// Correct DC offset
	double dc_correction = 1.0/pw;
	if (pw < .5) {
		dc_correction = 1.0/(1.0-pw);
	}

	double rect_val = dc_correction*(.5*saw_1_val - .5*saw_2_val);
	return rect_val;
}

void
dpw_order_n()
{
}

void
dpw_oversample()
{
}

/* From: https://varietyofsound.wordpress.com/2011/02/14/efficient-tanh-computation-using-lamberts-continued-fraction/ */
double
fast_tanh(double x)
{
	double x2 = x*x;
	double a = (((x2 + 378.0)*x2 + 17325.0)*x2 + 135135.0)*x;
	double b = ((28.0*x2 + 3150.0)*x2 + 62370.0)*x2 + 135135.0;
	double f = a/b;
	return f;
}

static void def_sub_eval(void *def_data, void *instance, float *l_in, float *r_in, float *l_out, float *r_out, int len, int sample_rate, int num_msg, RtNodeMsg *msgs)
{
    struct sub *data = instance;
    float p = 1.0f/(float)sample_rate;

    int curr_msg = 0;
    for (int n = 0; n < len; n++) {

        while (curr_msg < num_msg && msgs[curr_msg].time == n) {
            if (msgs[curr_msg].type == NOTE_ON) {

				struct sub_voice *voice = &data->voices[data->last_assign];

				bool found = false;
				/* What to do when we receive a NOTE_ON for a note that is already on? */
				/* Should we retrigger, or just layer? right now retrigger is causing click... */
				/* First see if there's a voice already playing that note. */
				for (int i = 0; i < NUM_VOICES; i++) {
					if (data->voices[i].note == msgs[curr_msg].note) {
						voice = &data->voices[i];
						found = true;
						break;
					}
				}

				/* Then search for an free voice */
				if (!found) {
					for (int i = 0; i < NUM_VOICES; i++) {
						if (!(data->voices[i].on)) {
							voice = &data->voices[i];
							found = true;
							break;
						}
					}	
				}

				/* Otherwise, steal one */
				if (!found) {
					voice = &data->voices[data->last_assign];
					data->last_assign = (data->last_assign + 1) % NUM_VOICES;
				}

                voice->on = true;
                voice->freq = midi_freq(msgs[curr_msg].note);
				printf("Freq: %f\n", voice->freq);
				voice->key_off = false;
				voice->phasor = 0.0f;
				voice->note = msgs[curr_msg].note;
				memset(voice->poly_blep_fir, 0, 4*sizeof(float));
				for (int o = 0; o < DPW_ORDER; o++) {
					voice->prev_dpw[o] = 0;
				}
				voice->dpw_initialized = false;

				voice->dpw_d1.initialized = false;
				voice->dpw_d2.initialized = false;
				voice->dpw_d1.prev = 0.0;
				voice->dpw_d2.prev = 0.0;
				voice->svf.hp_delay = 0.0;
				voice->svf.bp_delay = 0.0;
				voice->svf.lp_delay = 0.0;
            }
            else if (msgs[curr_msg].type == NOTE_OFF) {
				for (int i = 0; i < NUM_VOICES; i++) {
					if (data->voices[i].note == msgs[curr_msg].note) {
						data->voices[i].key_off = true;
						data->voices[i].on = false;
						data->voices[i].dpw_initialized = false;
						data->voices[i].dpw_d1.initialized = false;
						data->voices[i].dpw_d2.initialized = false;
						data->voices[i].dpw_d1.prev = 0.0;
						data->voices[i].dpw_d2.prev = 0.0;
						// data->voices[i].t_off = data->voices[i].phase + n*p;
						// TODO: get amp envelope working w/ phase

						// Should NOTE_OFF turn off all notes playing that voice, or just one?
						// break;
					}
				}
            }
            curr_msg++;
        }

// Oversampling DPW

#if 1
// DPW func-call
	for (int i = 0; i < NUM_VOICES; i++) {
		struct sub_voice *voice = data->voices + i;
		if (!voice->on) {
			continue;
		}
		
		double phase_inc = voice->freq/sample_rate;

		double amp = 1.0/NUM_VOICES;
		// Saw wave
		 double dpw = dpw_saw(&(voice->dpw_d1), voice->phasor, voice->freq, sample_rate);
		// Square wave (rect w/ pw=50%)
		// double dpw = dpw_rect(&(voice->dpw_d1), &(voice->dpw_d2), voice->phasor, voice->freq, sample_rate, .15);
		// Triangle Wave
		// double dpw = dpw_tri(&(voice->dpw_d1), voice->phasor, voice->freq, sample_rate);

		// Filter
		#if 0
		// Moog Ladder
		#endif

		#if 1
		// 2-pole TPT-SVF

		// Frequency warping
		//

		double g = data->prewarp_cutoff_g;

		double R = data->res;

		// Zero-Delay loop, filter
		// TODO: can calculate all these independently if need be
		// hp = (x - 2RS1 - GS1 - S2)/(1 + 2RG + G^2)
		// bp = G*hp + S1
		// lp = G*bp + S2

		/* Linear Model */
/*
		double high_pass = (dpw - 2.0*R*voice->svf.hp_delay - g*voice->svf.hp_delay - voice->svf.bp_delay)/(1.0 + 2.0*R*g + g*g);
		double band_pass = (g*high_pass + voice->svf.hp_delay); ///(1.0+g);
		double low_pass = (g*band_pass + voice->svf.bp_delay); ///(1.0+g);
*/
		/* NonLinear Model w/ Saturation on bandpass feedback path */

		/* NonLinear Model w/ Saturation on integrators */
		double high_pass = (dpw - 2.0*R*voice->svf.hp_delay - g*voice->svf.hp_delay - voice->svf.bp_delay)/(1.0 + 2.0*R*g + g*g);
		/* TODO: replace calls to tanh w/ a Pade approximant */
		double band_pass = fast_tanh(g*high_pass + voice->svf.hp_delay); ///(1.0+g);
		double low_pass = fast_tanh(g*band_pass + voice->svf.bp_delay); ///(1.0+g);
		voice->svf.hp_delay = band_pass;
		voice->svf.bp_delay = low_pass;

		// Non-linear saturation for self-oscillation?
		// tanh()

		// printf("DPW: %f G: %f S1: %f S2: %f cutoff_warp: %f HP: %f BP: %f LP: %f\n", dpw, G, S1, S2, cutoff_warp, high_pass, band_pass, low_pass);
		// assert(low_pass <= 10.0);


		double v = amp*low_pass;
		// printf("%f\n",low_pass);
		#endif

		#if 0
		// Naive lowpass filter
		double w_c = data->cutoff/(2.0*sample_rate);
		voice->svf.lp_delay = voice->svf.lp_delay + (dpw - voice->svf.lp_delay)*w_c;
		printf("NaiveLP: %f w_c: %f\n", voice->svf.lp_delay, w_c);
		assert(!isnan(voice->svf.lp_delay));
		double v = amp*voice->svf.lp_delay;
		#endif 


		#if 0
		// 1 pole filter RC
		double g = data->cutoff/(2.0*sample_rate);

		// double u = (dpw - voice->svf.lp_delay)*g/(1.0+g);
		// double lp = u + voice->svf.lp_delay;

		double lp = (dpw*g + voice->svf.lp_delay)/(1.0 + g);
		voice->svf.lp_delay = lp;
		printf("1PoleLP: %f g: %f \n", lp, g );

		double v = amp*lp;
		#endif		

		l_out[n] += v;
		r_out[n] += v;

		voice->phasor += phase_inc;
		if (voice->phasor >= 1.0) {
			voice->phasor -= 1.0;
		}
	}
#endif

#if 0
// 2nd order DPW (no delay)
// Create a PWM rectangle by subtracting two DPW sawtooths
// Create a DPW triangle by ... (it's possible)

	for (int i = 0; i < NUM_VOICES; i++) {
		struct sub_voice *voice = data->voices + i;
		if (!voice->on) {
			continue;
		}

		double dpw_P = sample_rate/voice->freq;
		double phase_inc = voice->freq/sample_rate;

		if (!voice->dpw_data.initialized) {
			// Time travel
			double prev_phasor = voice->phasor - phase_inc;
			if (prev_phasor < 0) {
				prev_phasor += 1.0;
			}

			double saw_prev = 2.0*prev_phasor - 1.0;
			double saw2_prev = saw_prev*saw_prev;
			voice->dpw_data.prev = saw2_prev;
			voice->dpw_data.initialized = true;
		}

		double saw = 2.0*voice->phasor - 1.0;	
		double saw2 = saw*saw;

		double dpw = saw2 - voice->dpw_data.prev;
		voice->dpw_data.prev = saw2;

		double dpw_scale = PI/(4.0*sin(PI/dpw_P));

		double amp = 1.0/NUM_VOICES;
		double v = amp*dpw*dpw_scale;

		l_out[n] += v;
		r_out[n] += v;

		voice->phasor += phase_inc;
		if (voice->phasor >= 1.0) {
			voice->phasor -= 1.0;
		}
	}
#endif

#if 0
// Order N DPW

	for (int i = 0; i < NUM_VOICES; i++) {
		struct sub_voice *voice = data->voices + i;
		if (!voice->on) {
			continue;
		}
		double phase_inc = voice->freq/sample_rate;
		double dpw_P = sample_rate/voice->freq;


		// DPW introduces DPW_ORDER delays when starting
		// So if DPW needs to be initialized, "Go back in time" DPW_ORDER samples
		// and start there

#if 0
		if (!voice->dpw_initialized) {
			// Initialize by going back in time DPW_ORDER samples
			double phasor = voice->phasor - DPW_ORDER*phase_inc;
			while (phasor < 0) {
				phasor += 1.0;
			}

			double curr_dpw[DPW_ORDER];

			for (int o = 0; o < DPW_ORDER; o++) {
				voice->prev_dpw[o] = 0.0;
			}

			for (int s = 0; s < DPW_ORDER; s++) {

				double saw = 2.0*voice->phasor - 1.0;
				double saw2 = saw*saw;
				double saw4 = saw2*saw2;
				double saw6 = saw4*saw2;


				curr_dpw[DPW_ORDER-1] = saw6 - 5.0*saw4 + 7.0*saw2;
				for (int o = DPW_ORDER-2; o >= 0; o--) {
					curr_dpw[o] = curr_dpw[o+1] - voice->prev_dpw[o+1];
				}

				for (int o = 0; o < DPW_ORDER; o++) {
					voice->prev_dpw[o] = curr_dpw[o];
				}

				phasor += phase_inc;
				if (phasor >= 1.0) {
					phasor -= 1.0;
				}
			}

			voice->dpw_initialized = true;
		}
#endif

#if 1
		// Delay (w/ out the go back in time thing.)
		if (!voice->dpw_initialized) {
			voice->dpw_delayed = 0;
			voice->dpw_initialized = true;
		}

#endif

		// Main DPW	
		double curr_dpw[DPW_ORDER];

		double saw = 2.0*voice->phasor - 1.0;
		double saw2 = saw*saw;
		double saw4 = saw2*saw2;
		double saw6 = saw4*saw2;

		curr_dpw[DPW_ORDER-1] = saw6 - 5.0*saw4 + 7.0*saw2;
		for (int o = DPW_ORDER-2; o >= 0; o--) {
			curr_dpw[o] = curr_dpw[o+1] - voice->prev_dpw[o+1];
		}

		for (int o = 0; o < DPW_ORDER; o++) {
			voice->prev_dpw[o] = curr_dpw[o];
		}

		double dpw_scale = pow(PI, DPW_ORDER-1)/(factorial(DPW_ORDER)*pow(2.0*sin(PI/dpw_P), DPW_ORDER-1));
		double v = curr_dpw[0]*dpw_scale;

		if (voice->dpw_delayed <= DPW_ORDER) {
			v = 0.0;
			voice->dpw_delayed += 1;
			printf("%d\n", voice->dpw_delayed);
		}

		double amp = 1.0/NUM_VOICES;
		v *= amp;
		
		l_out[n] += v;
		r_out[n] += v;

		voice->phasor += phase_inc;
		if (voice->phasor >= 1.0) {
			voice->phasor -= 1.0;
		}
	}

#endif

#if 0
// BLEP

		struct blep *blep = def_data;

		for (int i = 0; i < NUM_VOICES; i++) {
			struct sub_voice *voice = data->voices + i;
			double phase_inc = voice->freq/sample_rate;
			

			if (!voice->on) {
				continue;
			}

			// Saw wave

			assert(voice->phasor <= 1.0);

			float mod = 0.0f;
			if (voice->phasor > 1.0 - phase_inc) {
				// On the left side of the discontinuity
				// Find distance from discontinuity
				float dist = (voice->phasor - 1.0);
				float idx = (1.0 + dist)*(blep->len/2 - 1);
				int index = idx;
				printf("L Dist: %f idx: %f Index: %d\n", index, dist, idx);
				mod = blep->table[index];
			}
			else if (voice->phasor < phase_inc) {
				// On the right side of the discontinuity
				float dist = voice->phasor; 
				float idx = (1.0 + dist)*(blep->len/2 - 1);
				int index = idx;
				printf("R Dist: %f idx: %f Index: %d\n", index, dist, idx);
				mod = blep->table[index];
			}

			// Create the naive sawtooth and modify it w/ the polyblep
			float saw = 2.0f*voice->phasor - 1.0f;
			float amp = 1.0f/(NUM_VOICES);
			float v = amp*(saw + mod);

			l_out[n] += v;
			r_out[n] += v;

			voice->phasor += phase_inc;
			if (voice->phasor >= 1.0) {
				voice->phasor -= 1.0;
			}
		}
#endif

#if 0
// PolyBLEP
		for (int i = 0; i < NUM_VOICES; i++) {
			struct sub_voice *voice = data->voices + i;
			double phase_inc = voice->freq/sample_rate;

			if (!voice->on) {
				continue;
			}

			// Saw wave
			float mod = 0.0f;
			if (voice->phasor > 1.0 - phase_inc) {
				// On the left side of the discontinuity
				// Find distance from discontinuity
				float dist = (1.0 - voice->phasor); 
			// printf("Phasor: %f Phase Inc: %f Dist: %f\n", voice->phasor, phase_inc, dist);
				mod = poly_blep(dist);
			}
			else if (voice->phasor < phase_inc) {
				// On the right side of the discontinuity
				float dist = voice->phasor; 
				mod = poly_blep(dist);
			//	printf("right\n");
			}

			// Create the naive sawtooth and modify it w/ the polyblep
			float saw = 2.0f*voice->phasor - 1.0f;
			float amp = 1.0f/(NUM_VOICES);
			float v = amp*(saw + mod);

			l_out[n] += v;
			r_out[n] += v;

			voice->phasor += phase_inc;
			if (voice->phasor >= 1.0) {
				voice->phasor -= 1.0;
				// printf("Discontinuity: Phasor: %f Inc: %f\n", voice->phasor, phase_inc);
			}
		}
		
#endif

#if 0
// FIR optimal PolyBLEP

		for (int i = 0; i < NUM_VOICES; i++) {
			struct sub_voice *voice = data->voices + i;
			double phase_inc = voice->freq/sample_rate;

			if (!voice->on) {
				continue;
			}

			// Saw wave
			float mod = 0.0f;
			if (voice->phasor >= 1.0 - phase_inc) {
				// On the left side of the discontinuity
				float d = (voice->phasor - 1.0);
				mod = poly_blep_optimal_fir(d, voice->poly_blep_fir);
			}
			else if (voice->phasor <= phase_inc) {
				// On the right side of the discontinuity
				float d = voice->phasor;
				mod = poly_blep_optimal_fir(d, voice->poly_blep_fir);
			}


			// Create the naive sawtooth and modify it w/ the polyblep
			float saw = 2.0f*voice->phasor - 1.0f;
			float corrected = saw + mod;
			float amp = 1.0f/(NUM_VOICES);
			float v = amp*corrected;

			l_out[n] += v;
			r_out[n] += v;

			// Update fir
			voice->poly_blep_fir[3] = voice->poly_blep_fir[2];
			voice->poly_blep_fir[2] = voice->poly_blep_fir[1];
			voice->poly_blep_fir[1] = voice->poly_blep_fir[0];
			voice->poly_blep_fir[0] = corrected;

			voice->phasor += phase_inc;
			if (voice->phasor >= 1.0) {
				voice->phasor -= 1.0;
				// printf("Discontinuity: Phasor: %f Inc: %f\n", voice->phasor, phase_inc);
			}
		}
		
#endif

#if 0
// Wavetable
	struct wave_table *saw = def_data;
		for (int i = 0; i < NUM_VOICES; i++) {
			struct sub_voice *voice = data->voices + i;
			double phase_inc = voice->freq/sample_rate;

			if (!voice->on) {
				continue;
			}

			// Read value from the appropriate wave table
			// Linear scan, only 10 tables so should be ok.
			int t = 0;
			while (t < NUM_WAVE_TABLES && saw->tables[t].max_freq < voice->freq) {
				t++;
			}

			// TODO: if t is out of range of the tables, just use a sine wave
			
			struct wave_table_entry *table = saw->tables + t;

			float incr = table->len*(voice->freq/(float)sample_rate);
			int km = sample_rate/(2.0*table->max_freq);

			// Linear Interpolation (1st order spline)
			double pos = voice->phasor*(table->len/saw->over_sample);
			int s = pos;
			double frac = pos - s;
			int s1 = (s+saw->over_sample) % (table->len);

			int so = s;// + table->offset;
			int s1o = s1;// + table->offset;

			printf("Table: %d Max Freq: %f Req Freq: %f Increment: %f Km: %d Tlen: %d Pos: %f S0o: %d S1o: %d\n", t, table->max_freq, voice->freq, incr, km, table->len, pos, so, s1o);
			double linear_interp = (table->data[s1o]-table->data[so])*frac + table->data[so];

			float amp = 1.0/NUM_VOICES;
			float v = amp * linear_interp;

			float time = voice->phasor/voice->freq;
			//v = amp*sawf_8ve(voice->phasor, voice->freq, sample_rate);

			// v = sawf(w_c*t+A_v*sinf(w_v*t));
	
			// printf("%f\n", v);

			l_out[n] += v;
			r_out[n] += v;

			voice->phasor += phase_inc;	
			if (voice->phasor >= 1.0) {
				voice->phasor -= 1.0;
			}
		}
#endif

	}
}

static void def_sub_free(void *instance)
{
}

static void def_sub_write(int size, uint8_t *data)
{
}

/* IDEA: Config audio params on the sysout node editor? */
/*
AudioDef def_sysout = {
    .name = "SysOut",
    .num_controls = 0,
    .make = def_disabled_make,
    .init = def_disabled_init,
    .eval = def_copy_io_eval,
    .free = def_disabled_free,
    .read = def_disabled_read,
    .write = def_disabled_write,
};
*/

AudioDef def_disabled = {
    .name = "Disabled",
    .num_controls = 0,
    .make = def_disabled_make,
    .init = def_disabled_init,
    .eval = def_disabled_eval,
    .free = def_disabled_free,
    .read = def_disabled_read,
    .write = def_disabled_write,
};

AudioDef def_fm = {
    .name = "FM",
    .num_controls = 0,
    .make = def_fm_make,
    .init = def_disabled_init,
    .eval = def_fm_eval,
    .free = def_fm_free,
    .read = def_fm_read,
    .write = def_fm_write,
};

AudioDef def_sub = {
	.name = "Subtractive",
	.data = &G_BLEP,
	.num_controls = 0,
	.make = def_sub_make,
	.init = def_disabled_init,
	.eval = def_sub_eval,
	.free = def_sub_free,
	.read = def_disabled_read,
	.write = def_sub_write
};

/* IDEA: Config audio params on the sysout node editor? */
/*
AudioDef def_sysout = {
    .name = "SysOut",
    .num_controls = 0,
    .make = def_disabled_make,
    .init = def_disabled_init,
    .eval = def_copy_io_eval,
    .free = def_disabled_free,
    .read = def_disabled_read,
    .write = def_disabled_write,
};
*/
