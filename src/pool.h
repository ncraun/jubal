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
	Node Pool
	
	Pool based memory manager for Jubal's Node Store.
	
	The pool is realtime safe, but it is not thread-safe. The pool should only
	be accessed from 1 thread at a time. The pool requires a single memory
	allocation at construction time, all other operations are real-time safe.
	
	If the pool is full, no new allocations can be performed.
	
	Originally Based on "Fast Efficient Fixed-Size Memory Pool: No Loops and
	No Overhead". However, it was changed to move all the work to a loop in
	initializtion time, instead of amatorizing it across all calls as
	suggested in the paper. So it's really not using the same technique at all.
	
	Free list is maintained implicitly by the pool structure. Unused blocks are
	by definition not in use, so we use them to store the links in the free list
	without any additional memory overhead.
	
	Addresses are 16bit integers, change if you need more entries.
	
	Usage Caveats:
	
	Valgrind can't detect double frees, etc.
*/

#ifndef JUBAL_POOL_H
#define JUBAL_POOL_H

/*
	Memory Pool allocator
*/

#include "node.h"

union pool_block {
	struct fm_data;
	struct sub_data;
	
	uint16_t next;
};

struct pool {
	int capacity;
	uint16_t head;
	union pool_block *pool;
};

typedef struct pool pool;

pool* make_pool(uint16_t capacity);
uint16_t pool_alloc(pool *p);
void pool_free(pool *p, uint16_t addr);
void* pool_get(pool *p. uint16_t addr);

#endif
