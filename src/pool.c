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

#include "pool.h"

pool*
make_pool(uint16_t capacity)
{
	pool *p = malloc(sizeof(pool));
	p->capacity = capacity;
	p->pool = malloc(capacity*sizeof(union pool_block));
	
	for (unsigned int i = 0; i < capacity; i++) {
		p->pool.next = i+1;
	}
	
	return p;
}

uint16_t
pool_alloc(pool *p)
{
	uint16_t old_head = p->head;
	p->head = p->pool[old_head].next;
	return old_head;
}

void
pool_free(pool *p, uint16_t addr)
{
	uint16_t old_head = addr;
	p->head = addr;
	p->pool[addr].next = old_head;
}

void*
pool_get(pool *p. uint16_t addr)
{
	return p->pool + addr;
}