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

#include "bitset.h"

/* 
    Fallback C implementation.
    From Hacker's Delight. 
*/
static int
ctz_c(uint32_t x)
{
   int n;

   if (x == 0) return(32);
   n = 1;
   if ((x & 0x0000FFFF) == 0) {n = n +16; x = x >>16;}
   if ((x & 0x000000FF) == 0) {n = n + 8; x = x >> 8;}
   if ((x & 0x0000000F) == 0) {n = n + 4; x = x >> 4;}
   if ((x & 0x00000003) == 0) {n = n + 2; x = x >> 2;}
   return n - (x & 1);
}

int 
bitset_next_set_bit(uint32_t *bitset, int from, int num_bits)
{
    int i = BITSET_SLOT(from);
    /* UINT32_MAX is a number with all the bits in the long set to 1 */
    uint32_t slot = bitset[i] & (UINT32_MAX << from);
    while (i < BITSET_SIZE(num_bits)) {
        if (slot != 0) {
            return i*32 + __builtin_ctz(slot);
        }

        i++;
        if (i < BITSET_SIZE(num_bits)) {
            slot = bitset[i];
        }
    }

    return -1; // Out of bits
}

