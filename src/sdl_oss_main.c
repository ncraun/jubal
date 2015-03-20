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
b

    You should have received a copy of the GNU General Public License
    along with Jubal.  If not, see <http://www.gnu.org/licenses/>.
*/

#define _GNU_SOURCE
#include <time.h>
#include <math.h>
#include <termios.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <sys/soundcard.h>
#include <stdint.h>
#include <limits.h>
#include <unistd.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/ioctl.h>
#include <fcntl.h>
#include <pthread.h>
#include "engine.h"
#include "SDL.h"
#include "nodes.h"

#define NON_RT_POLL_NS 500000000L

#define BILLION 10000000000LU

int mygetch( ) {
  struct termios oldt,
                 newt;
  int            ch;
  tcgetattr( STDIN_FILENO, &oldt );
  newt = oldt;
  newt.c_lflag &= ~( ICANON | ECHO );
  tcsetattr( STDIN_FILENO, TCSANOW, &newt );
  ch = getchar();
  tcsetattr( STDIN_FILENO, TCSANOW, &oldt );
  return ch;
}

/*
    Jack Backend ?
*/

struct timespec start;

void
sdl_callback(void *userdata, Uint8* stream, int len)
{
    Engine *engine = userdata;
    int samples = len/sizeof(float)/2; //2 channel stereo
    float *fstr = (float*)stream;
    float l[samples];
    float r[samples];
    AuBuff buff;
    buff.left = l;
    buff.right = r;
    struct timespec now;
    clock_gettime(CLOCK_MONOTONIC, &now);

    // Find difference between start and now
    while (now.tv_sec >= 0 && start.tv_nsec > now.tv_nsec) { 
        now.tv_sec--;
        now.tv_nsec += BILLION;
    }
    uint64_t diff_sec = now.tv_sec - start.tv_sec;
    uint64_t diff_nsec = now.tv_nsec - start.tv_nsec;

    uint64_t samp_diff = (diff_sec*engine->rt.sample_rate) + (diff_nsec*engine->rt.sample_rate)/BILLION;

    jubal_callback(buff, samples, engine, samp_diff);

    for (int i = 0; i < samples; i++) {
        fstr[i*2] = l[i];
    }
    for (int i = 0; i < samples; i++) {
        fstr[i*2+1] = r[i];
    }
}

static int
sdl_key_to_note(int ch)
{
    switch(ch) {
        case SDLK_z:
            return 24;
        case SDLK_s:
            return 25;
        case SDLK_x:
            return 26;
        case SDLK_d:
            return 27;
        case SDLK_c:
            return 28;
        case SDLK_v:
            return 29;
        case SDLK_g:
            return 30;
        case SDLK_b:
            return 31;
        case SDLK_h:
            return 32;
        case SDLK_n:
            return 33;
        case SDLK_j:
            return 34;
        case SDLK_m:
            return 35;
        case SDLK_COMMA:
            return 36;
        case SDLK_q:
            return 36;
        case SDLK_2:
            return 37;
        case SDLK_w:
            return 38;
        case SDLK_3:
            return 39;
        case SDLK_e:
            return 40;
        case SDLK_r:
            return 41;
        case SDLK_5:
            return 42;
        case SDLK_t:
            return 43;
        case SDLK_6:
            return 44;
        case SDLK_y:
            return 45;
        case SDLK_7:
            return 46;
        case SDLK_u:
            return 47;
        case SDLK_i:
            return 48;
        case SDLK_9:
            return 49;
        case SDLK_o:
            return 50;
        case SDLK_0:
            return 51;
        case SDLK_p:
            return 52;
    }

    return -1;
}

static uint16_t
f32_to_u16_sat(float f)
{
	float sc = (f+1.0f)*0.5f;
	if (sc > 1.0f) { sc = 1.0f; }
	else if (sc < 0.0f) { sc = 0.0f; }

	return sc*65535;
}

int
oss_init()
{
	int snd = open("/dev/dsp", O_WRONLY);//|O_NONBLOCK);
	if (snd == -1) {
		perror("Err:");
		exit(1);
	}

	/* Set up Device */

	/* Format */
	int snd_fmt = AFMT_S16_LE;
	if (ioctl(snd, SNDCTL_DSP_SETFMT, &snd_fmt) == -1) {
		perror("Err:");
		exit(1);
	}
	if (snd_fmt != AFMT_S16_LE) {
		fprintf(stderr, "Wrong format.\n");
		exit(1);
	}

	/* Channels */
	int snd_chls = 2;
	if (ioctl(snd, SNDCTL_DSP_CHANNELS, &snd_chls) == -1) {
		perror("Err:");
		exit(1);
	}
	if (snd_chls != 2) {
		fprintf(stderr, "Not stereo.\n");
		exit(1);
	}

	/* Sampling Rate */
	int snd_rate = 48000;
	if (ioctl(snd, SNDCTL_DSP_SPEED, &snd_rate) == -1) {
		perror("sampling rate:");
		exit(1);
	}
	if (snd_rate != 48000) {
		fprintf(stderr, "Sampling rate is %d instead of 48000.\n", snd_rate);
		exit(1);
	}

	/* Buffer Size */
	int snd_buff_size = 128*sizeof(uint16_t)*2; // Blocksize in BYTES
	if (ioctl(snd, SNDCTL_DSP_SETBLKSIZE, &snd_buff_size) == -1) {
		perror("setblksize:");
		exit(1);
	}
	printf("Buffer Size: %d bytes\n", snd_buff_size);

	return snd;
}

struct sound_data {
	Engine *engine;
	int snd;
	int snd_buff_size;
	float time;
};

void*
oss_thread(void *arg) {

	struct sound_data *data = arg;

	/* Write Loop (Rand()) */
	for (;;) {
		Engine *engine = data->engine;
		int samples = data->snd_buff_size/4; //2 channel stereo
		float fstr[samples*2];
		float l[samples];
		float r[samples];
		AuBuff buff;
		buff.left = l;
		buff.right = r;
		struct timespec now;
		clock_gettime(CLOCK_MONOTONIC, &now);

		// Find difference between start and now
		while (now.tv_sec >= 0 && start.tv_nsec > now.tv_nsec) { 
			now.tv_sec--;
			now.tv_nsec += BILLION;
		}
		uint64_t diff_sec = now.tv_sec - start.tv_sec;
		uint64_t diff_nsec = now.tv_nsec - start.tv_nsec;

		uint64_t samp_diff = (diff_sec*engine->rt.sample_rate) + (diff_nsec*engine->rt.sample_rate)/BILLION;

		jubal_callback(buff, samples, engine, samp_diff);

		for (int i = 0; i < samples; i++) {
			fstr[i*2] = l[i];
		}
		for (int i = 0; i < samples; i++) {
			fstr[i*2+1] = r[i];
		}

		int16_t out[samples*2];
		for (int i = 0; i < samples*2; i++) {
			out[i] = fstr[i]*INT16_MAX;//f32_to_u16_sat(fstr[i]);
		}
		write(data->snd, out, sizeof(int16_t)*samples*2);
		// printf("output buffer.\n");
	}

	return NULL;
}

int
main()
{

    SDL_Init(SDL_INIT_VIDEO);

    SDL_Window *window;
    SDL_Renderer *renderer;
    SDL_CreateWindowAndRenderer(320, 240, 0, &window, &renderer);

    SDL_AudioSpec want, have;
    SDL_AudioDeviceID dev;

    int num_nodes = 256;
    int queue_size = 256;

	int snd = oss_init();
	struct sound_data sd;
	sd.snd = snd;
	sd.snd_buff_size = 128;

    Engine *engine = make_engine(num_nodes, 128, queue_size, 48000);
	sd.engine = engine;
	sd.time = 0.0f;

    clock_gettime(CLOCK_MONOTONIC, &start);
    engine_play(engine);
	pthread_t oss_thr;
	pthread_create(&oss_thr, NULL, oss_thread, &sd);

    /* Initiliaze sine node and connect to output */
    int sine_node = 5;
    engine_change_node_def(engine, sine_node, 3);
    engine_connect_output(engine, sine_node);
    /*
        Live Play!
    */

    bool running = true;
    int oct = 2;
    SDL_Event e;
    SDL_PauseAudioDevice(dev, 0);
	double f = 10000;
	double r = 1.0;
    while (running) {
        while (SDL_WaitEvent(&e)) {
            if (e.type == SDL_QUIT) {
                running = false;
                break;
            }
            else if (e.type == SDL_KEYDOWN && !e.key.repeat) {
				if (e.key.keysym.sym == SDLK_EQUALS) {
					oct++;
				}
				else if (e.key.keysym.sym == SDLK_MINUS) {
					oct--;
				}
				/* Filter */
				else if (e.key.keysym.sym == SDLK_LEFTBRACKET) {
					union param_data p;
					f -= 250.0;
					p.f32 = f;
					engine_node_cc(engine, sine_node, SUB_CUTOFF, p, 0);
				}
				else if (e.key.keysym.sym == SDLK_RIGHTBRACKET) {
					union param_data p;
					f += 250.0;
					p.f32 = f;
					engine_node_cc(engine, sine_node, SUB_CUTOFF, p, 0);
				}
				else if (e.key.keysym.sym == SDLK_DOWN) {
					union param_data p;
					r -= .25;
					p.f32 = r;
					engine_node_cc(engine, sine_node, SUB_RES, p, 0);
				}
				else if (e.key.keysym.sym == SDLK_UP) {
					union param_data p;
					r += .25;
					p.f32 = r;
					engine_node_cc(engine, sine_node, SUB_RES, p, 0);
				}

                int note = sdl_key_to_note(e.key.keysym.sym);
                if (note != -1) {
                    engine_note_on(engine, sine_node, note+12*oct, 0);
                    printf("Note On: %d\n", note);
                }
            }
            else if (e.type == SDL_KEYUP && !e.key.repeat) {
                int note = sdl_key_to_note(e.key.keysym.sym);
                if (note != -1) {
                    engine_note_off(engine, sine_node, note+12*oct, 0);
                    printf("Note On: %d\n", note);
                }
            }
        }
    }

	
	close(snd);
    return 0;
}
