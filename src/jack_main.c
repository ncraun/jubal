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

#define _GNU_SOURCE
#include <time.h>
#include <math.h>
#include <termios.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <jack/jack.h>
#include <pthread.h>
#include "engine.h"

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

jack_port_t *l_playback_port;
jack_port_t *r_playback_port;
struct timespec start;

int
backend_jack_buffer_size_change(jack_nframes_t buffer_size, void *arg)
{
    Engine *engine = arg;
    fprintf(stderr, "Changing buffer size from %d to %d\n", engine->rt.buffer_size, buffer_size);
    engine_change_buffer_size(engine, buffer_size);

    return 0;
}

int
backend_jack_sample_rate_change(jack_nframes_t sample_rate, void *arg)
{
    /*fprintf(stderr, "Changing sample rate not implemeted yet!\n");
    exit(1);*/
    Engine *engine = arg;
    engine_change_sample_rate(engine, sample_rate);
    return 0;
}

int
backend_jack_callback(jack_nframes_t nframes, void *arg)
{
    float *l_out = jack_port_get_buffer(l_playback_port, nframes);
    float *r_out = jack_port_get_buffer(r_playback_port, nframes);
    Engine *engine  = arg;
    AuBuff buff;
    buff.left = l_out;
    buff.right = r_out;
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
    
    jubal_callback(buff, nframes, engine, samp_diff);
    return 0;
}
int done = 0;
void*
non_rt_func(void* param)
{
    Engine *engine = param;
    struct timespec poll_time;
    poll_time.tv_sec = 0;
    poll_time.tv_nsec = NON_RT_POLL_NS;
    while (!__sync_fetch_and_add(&done, 0)) {
        // Poll the free queue when needed
        // experiment with wait time

        printf("poll\n");
        FreeMsg tmp;
        while (lfqueue_size(engine->free_q) > 0) {
            lfqueue_poll(engine->free_q, &tmp);
            printf("NONRT Free: %p!\n", tmp.data);
            //tmp.free(tmp.data);
        }

        nanosleep(&poll_time, NULL);
    }
    printf("done!\n");
    // pthread_exit vs return ?
    // pthread_exit(NULL);
    return NULL;
}

static unsigned char
key_to_note_base(int ch)
{
    int r;
    switch(ch) {
        case 'z':
            return 24;
        case 's':
            return 25;
        case 'x':
            return 26;
        case 'd':
            return 27;
        case 'c':
            return 28;
        case 'v':
            return 29;
        case 'g':
            return 30;
        case 'b':
            return 31;
        case 'h':
            return 32;
        case 'n':
            return 33;
        case 'j':
            return 34;
        case 'm':
            return 35;
        case ',':
            return 36;
        case 'q':
            return 36;
        case '2':
            return 37;
        case 'w':
            return 38;
        case '3':
            return 39;
        case 'e':
            return 40;
        case 'r':
            return 41;
        case '5':
            return 42;
        case 't':
            return 43;
        case '6':
            return 44;
        case 'y':
            return 45;
        case '7':
            return 46;
        case 'u':
            return 47;
        case 'i':
            return 48;
        case '9':
            return 49;
        case 'o':
            return 50;
        case '0':
            return 51;
        case 'p':
            return 52;
    }

    return 0;
}
int
main()
{

    jack_client_t *client = jack_client_open("Jubal 0.0.1", JackNullOption, NULL);
    uint32_t sample_rate = jack_get_sample_rate(client);
    uint32_t buffer_size = jack_get_buffer_size(client);
    int num_nodes = 256;
    int queue_size = 256;

    l_playback_port = jack_port_register(client, "left", JACK_DEFAULT_AUDIO_TYPE, JackPortIsOutput|JackPortIsTerminal, 0);
    r_playback_port = jack_port_register(client, "right", JACK_DEFAULT_AUDIO_TYPE, JackPortIsOutput|JackPortIsTerminal, 0);

    Engine *engine = make_engine(num_nodes, buffer_size, queue_size, sample_rate);

    /*Start non-rt Free(ing) thread */
    pthread_t non_rt;
    pthread_attr_t non_rt_attr;

    pthread_attr_init(&non_rt_attr);
    pthread_attr_setdetachstate(&non_rt_attr, PTHREAD_CREATE_DETACHED);
    // How big a stack do we need?
    pthread_attr_setstacksize(&non_rt_attr, 256*1024);

    pthread_create(&non_rt, &non_rt_attr, non_rt_func, engine);


    jack_set_buffer_size_callback(client, backend_jack_buffer_size_change, engine);
    jack_set_process_callback(client, backend_jack_callback, engine);
    jack_set_sample_rate_callback(client, backend_jack_sample_rate_change, engine);


    /*
        Can't connect ports before client is activated.
        NOTE: Jack's port types are confusing:
        Output means output from jack -> application, this is a capture port
        Input means input jack <- application, this is a playback port
    */
    clock_gettime(CLOCK_MONOTONIC, &start);
    engine_play(engine);
    jack_activate(client);
    const char **ports = jack_get_ports(client, NULL, NULL, JackPortIsPhysical|JackPortIsInput);
    jack_connect(client, jack_port_name(l_playback_port), ports[0]);
    jack_connect(client, jack_port_name(r_playback_port), ports[1]);


    /* Initiliaze sine node and connect to output */
    int sine_node = 5;
    engine_change_node_def(engine, sine_node, 1);
    engine_connect_output(engine, sine_node);
    /*
        Live Play!
    */
/*
    engine_note_on(engine, sine_node, 60, 0);
    sleep(1);
    engine_note_on(engine, sine_node, 65, 0);
    sleep(1);
    engine_note_on(engine, sine_node, 72, 0);
    sleep(1);
*/
    char ch = '\0';

    struct timespec slptime;
    slptime.tv_sec = 0;
    slptime.tv_nsec =  125000000L;

    int scale[] = {21, 23, 24, 26, 28, 29, 32};
    int scale_len = sizeof(scale)/sizeof(int);

    while (1) {
        ch = mygetch();
        int note = 36+key_to_note_base(ch);
        // Make ch pentatonic
        engine_note_on(engine, sine_node, note, 0);
        nanosleep(&slptime, NULL);
        engine_note_off(engine, sine_node, note, 0);
    }

    jack_deactivate(client);
    jack_client_close(client);
    __sync_add_and_fetch(&done, 1);
    pthread_exit(NULL);
    return 0;
}
