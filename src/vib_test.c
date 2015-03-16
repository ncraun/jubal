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
#include "engine.h"

int
main()
{

    uint32_t sample_rate = 48000;
    uint32_t buffer_size = 256;
    int num_nodes = 256;
    int queue_size = 256;

    Engine *engine = make_engine(num_nodes, buffer_size, queue_size, sample_rate);

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
    engine_note_on(engine, sine_node, 48, 0);
    int num_sec = 120;
    int num_callbacks = (48000.0/256.0)*num_sec;
    float min = 10000;
    float max = -10000;
    for (int i = 0; i < num_callbacks; i++) {
        jubal_callback(buff, 256, engine, 0);
        write(2, buff.left, 256*sizeof(float));
    }

    engine_free(engine);
    free(buff.left);
    free(buff.right);

    pthread_exit(NULL);
    return 0;
}

