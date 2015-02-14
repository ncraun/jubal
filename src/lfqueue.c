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

/*
    Implementation of lfqueue using gcc sync builtins
*/

#include <stdlib.h>
#include <string.h>
#include "lfqueue.h"

Lfqueue*
make_lfqueue(int capacity, int el_size)
{
    Lfqueue *q = malloc(sizeof(Lfqueue));
    q->head = 0;
    q->tail = 0;
    q->capacity = capacity;
    q->el_size = el_size;
    q->data = malloc(el_size*capacity);

    return q;
}

bool
lfqueue_is_empty(Lfqueue *q)
{
    __sync_synchronize();
    return q->head == q->tail;
}

int
lfqueue_size(Lfqueue *q)
{
    __sync_synchronize();
    int sz = q->head - q->tail;
    return sz >= 0 ? sz : -sz;
}

// Adding to a full queue?
void
lfqueue_add(Lfqueue *q, const void *data)
{
    memcpy(q->data+(q->tail*q->el_size), data, q->el_size);
    // increase size after writing element
    __sync_synchronize();
    q->tail = (q->tail + 1) % q->capacity;
}

// What to do when polling empty queue?
void
lfqueue_poll(Lfqueue *q, void *out)
{
    memcpy(out, q->data+(q->head*q->el_size), q->el_size);
    // Write data before updating index
    __sync_synchronize();
    q->head = (q->head + 1) % q->capacity;
}

void
lfqueue_peek(const Lfqueue *q, void *out)
{
    __sync_synchronize();
    memcpy(out, q->data+(q->head*q->el_size), q->el_size);
}

void
lfqueue_advance(Lfqueue *q)
{
    __sync_synchronize();
    q->head++;
}

void
free_lfqueue(Lfqueue *q)
{
    free(q->data);
    free(q);
}
