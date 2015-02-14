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
    Leak detection wrapper around malloc and free.
*/
#define _GNU_SOURCE
#include <dlfcn.h>
#include <stdio.h>

static size_t alloc = 0;

void* malloc(size_t sz)
{
    void *(*libc_malloc)(size_t) = dlsym(RTLD_NEXT, "malloc");
    alloc += sz;
    void *addr = libc_malloc(sz);
    printf("malloc %p\n", addr);
    return addr;
}

void free(void *p)
{
    void (*libc_free)(void*) = dlsym(RTLD_NEXT, "free");
    printf("free %p\n", p);
    libc_free(p);
}
