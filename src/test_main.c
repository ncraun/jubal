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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <pthread.h>
#include "engine.h"

#define NON_RT_POLL_NS 500000000L

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

int
main()
{

    uint32_t sample_rate = 48000;
    uint32_t buffer_size = 256;
    int num_nodes = 256;
    int queue_size = 256;

    Engine *engine = make_engine(num_nodes, buffer_size, queue_size, sample_rate);

    /*Start non-rt Free(ing) thread */
    pthread_t non_rt;
    pthread_attr_t non_rt_attr;

    pthread_attr_init(&non_rt_attr);
    pthread_attr_setdetachstate(&non_rt_attr, PTHREAD_CREATE_DETACHED);
    // How big a stack do we need?
    pthread_attr_setstacksize(&non_rt_attr, 256*1024);


    engine_play(engine);

    /* Initiliaze sine node and connect to output */
    int sine_node = 5;
    engine_change_node_def(engine, sine_node, 1);
    if (!engine_connect_output(engine, sine_node)) {
        printf("Error connecting output.\n");
    }
    printf("Connected!\n");
    AuBuff buff;
    buff.left = malloc(256*sizeof(float));
    buff.right = malloc(256*sizeof(float));
    for (int i = 0; i < 100; i++) {
        jubal_callback(buff, 256, engine, 0);
    }

    __sync_add_and_fetch(&done, 1);

    for (int i = 0; i < 256; i++) {
        printf("%0.2f ", buff.left[i]);
    }
    puts("\n");

    engine_free(engine);
    free(buff.left);
    free(buff.right);

    pthread_exit(NULL);
    return 0;
}

