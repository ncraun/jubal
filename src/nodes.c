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

struct sine {
    float freq;
    float time;
    bool on;
};

static void*
def_sine_make()
{
    struct sine *data = malloc(sizeof(struct sine));
    data->freq = 0;
    data->time = 0;
    data->on = false;
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
                data->on = true;
                printf("Note On: %d\n", msgs[curr_msg].note);
                data->freq = midi_freq(msgs[curr_msg].note);
                printf("freq -> %f\n", data->freq);
            }
            else if (msgs[curr_msg].type == NOTE_OFF) {
                data->on = false;
            }
            else if (msgs[curr_msg].type == CC) {
            }
            curr_msg++;
        }

        float v = 0.0;
        if (data->on) {
            v = .1*sinf(2*PI*data->freq*data->time);
            data->time += sample_time;
        }

        left_output[f] = v;
        right_output[f] = v;
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
