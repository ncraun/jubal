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

#ifndef JUBAL_BITSET_H
#define JUBAL_BITSET_H

/*  
    Based on the C-FAQ Bitset implementation: 

    http://c-faq.com/misc/bitsets.html

    Assumes 8 bits per byte.

    Usage:

    Allocation:
    Single:
    Static: uint32_t bs[BITSET_SIZE(n)];
    Dynamic: uint32_t *bs = malloc(BITSET_BYTES(n)); 

    Contiguous Multiple:
    Static: uint32_t bs[BITSET_SIZE(n)*m];
    Dynamic:  *bs = malloc(BITSET_BYTES(n)*m); 
    Accessing bitset b ptr uint32_t *bs: (bs + b*BITSET_SIZE(n))

    Initialization:
    All Bits False: memset(bs, 0, BITSET_BYTES(n));
    All Bits True:  memset(bs, 0xFF, BITSET_BYTES(n));

    TODO:

    Set Functions:
    union:
    intersection:
    complement:

    Iteration:
    next_set_bit
    next_clear_bit
*/

#include <stdint.h>

/* Private */
#define BITSET_SLOT(i) ((i)/32)
#define BITSET_MASK(i) (1 << ((i)%32))

/* Public */
#define BITSET_SIZE(n) ( ((n)+32-1)/32 )
#define BITSET_BYTES(n) (4*BITSET_SIZE(n))
#define BITSET_TEST(b, i) ( (b)[BITSET_SLOT(i)] & BITSET_MASK(i) )
#define BITSET_SET(b, i) ( (b)[BITSET_SLOT(i)] |= BITSET_MASK(i) )
#define BITSET_CLEAR(b, i) ( (b)[BITSET_SLOT(i)] &= ~(BITSET_MASK(i)) )
#define BITSET_FLIP(b, i) ( (b)[BITSET_SLOT(i)] ^= BITSET_MASK(i) )

int bitset_next_set_bit(uint32_t *bitset, int from, int num_bits);

#endif
