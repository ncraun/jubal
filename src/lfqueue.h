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

#ifndef JUBAL_LFQUEUE_H
#define JUBAL_LFQUEUE_H

/*
    Lock free message queue.

    Only use this in single reader single writer scenario.

    It's probably not as correct as it should be, but it seems to work in
    practice.

    This data structure isn't really that safe. Be careful!
*/

#include <stdint.h>
#include <stdbool.h>

#define CACHE_LINE_SIZE 64

struct lfqueue {
    /*
        TODO: Make sure head and tail are in different cache lines.
        Head and tail are written to by different threads so if they are in
        the same cache line this will cause contention, because the threads
        will invalidate each others caches.

        Assuming 64-byte cacheline
    */
    int head;
    uint8_t padding_1[CACHE_LINE_SIZE];
    int tail;
    uint8_t padding_2[CACHE_LINE_SIZE];

    int el_size;
    int capacity;

    /* Char instead of uint8_t due to strict aliasing rules */
    char* data;
};

typedef struct lfqueue Lfqueue;

/*
    capacity: number of objects in the queue.
    el_size: size of a single object in bytes
*/
Lfqueue *make_lfqueue(int capacity, int el_size);

bool lfqueue_is_empty(Lfqueue *q);

/* Call from producer thread ONLY */
void lfqueue_add(Lfqueue *q, const void *data);

/* Call from consumer thread ONLY */
void lfqueue_poll(Lfqueue *q, void *out);

/* Call from consumer thread ONLY */
/* out must be el_size big */
void lfqueue_peek(const Lfqueue *q, void *out);

void lfqueue_advance(Lfqueue *q);

int lfqueue_size(Lfqueue *q);

void free_lfqueue(Lfqueue *q);

#endif
