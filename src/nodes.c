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
#include "const.h"
#include "engine.h"

#define NUM_VOICES 5

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

struct sine {
    /* Voice stealing using array and index as a */
    struct sine_voice voices[NUM_VOICES];
    /* Round robin voice stealing */
    int voice_assign;
    float vib_rate;
    float vib_depth;
    float trem_rate;
    float trem_depth;
    float atk_time;
    float sus_level;
    float rls_time;
    float port_time;
    float port_depth;
};

static void*
def_sine_make()
{
    struct sine *data = malloc(sizeof(struct sine));
    for (int i = 0; i < NUM_VOICES; i++) {
        data->voices[i].on = false;
    }
    data->voice_assign = 0;
    data->vib_rate = 440.0f;
    data->vib_depth = .1;
    data->trem_rate = 440.0f;
    data->trem_depth = .08;

    data->atk_time = 0.1f;
    data->sus_level = 1.0;
    data->rls_time = .5;

    data->port_time = .3f;
    data->port_depth = 20.0f;

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

static void
def_sine_eval(void *instance, float *left_input, float *right_input, float *left_output, float *right_output, int len, int sample_rate, int num_msg, RtNodeMsg *msgs)
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

                float f = data->voices[i].freq;// + vib; 
                if (t < data->port_time) {
                    float rise = data->port_depth;
                    float run = data->port_time;
                    float slope = rise/run;

                    f = slope*t + (data->voices[i].freq-data->port_depth);
                }

                //float vib = data->vib_depth*sinf(2.0f*PI*data->vib_rate*t);
                //float trem = data->trem_depth*sinf(2.0f*PI*data->trem_rate*t);
                val = amp*mix*sawf(2*PI*f*t);
                data->voices[i].time += sample_time;
                //printf("%0.8f %0.8f %0.8f\n", t, vib, val);
            }

            left_output[f] += val;
            right_output[f] += val;
        }
    }
}


static void*
def_disabled_make()
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
def_disabled_eval(void *instance, float *left_input, float *right_input, float *left_output, float *right_output, int len, int sample_rate, int num_msg, RtNodeMsg *msgs)
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

AudioDef def_sine = {
    .name = "Sine",
    .num_controls = 0,
    .make = def_sine_make,
    .init = def_disabled_init,
    .eval = def_sine_eval,
    .free = def_sine_free,
    .read = def_disabled_read,
    .write = def_disabled_write,
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
