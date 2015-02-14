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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "bitset.h"
#include "test.h"

int
main(int argc, char **argv)
{

    // Slot
    {
        for (int i = 0; i < 32; i++) {
            int slot = BITSET_SLOT(i);
            TEST("Bitset Slot 0", slot == 0);
        }
        for (int i = 32; i < 64; i++) {
            int slot = BITSET_SLOT(i);
            TEST("Bitset Slot 1", slot == 1);
        }
        for (int i = 64; i < 96; i++) {
            int slot = BITSET_SLOT(i);
            TEST("Bitset Slot 2", slot == 2);
        }
        for (int i = 96; i < 128; i++) {
            int slot = BITSET_SLOT(i);
            TEST("Bitset Slot 3", slot == 3);
        }
    }

    // Mask
    {
        uint32_t masks[] = {
            1, 2, 4, 8, 16, 32, 64, 128,
            256, 512, 1024, 2048, 4096, 8192, 16384, 32768, 65536,
            131072, 262144, 524288, 1048576, 2097152, 4194304, 8388608, 16777216,
            33554432, 67108864, 134217728, 268435456, 536870912, 1073741824, 2147483648};

        for (int i = 0; i < 4*32; i++) {
            TEST("Bitset Mask", BITSET_MASK(i) == masks[i%32]);
        }
    }

    // Size
    {
        TEST("Bitset Size 0", BITSET_SIZE(0)==0);
        for (int i = 1; i <= 32; i++) {
            int sz = BITSET_SIZE(i);
            TEST("Bitset Size 1", sz == 1);
        }
        for (int i = 33; i <= 64; i++) {
            int sz = BITSET_SIZE(i);
            TEST("Bitset Size 2", sz == 2);
        }
        for (int i = 65; i <= 96; i++) {
            int sz = BITSET_SIZE(i);
            TEST("Bitset Size 3", sz == 3);
        }
        for (int i = 97; i <= 128; i++) {
            int sz = BITSET_SIZE(i);
            TEST("Bitset Size 4", sz == 4);
        }
    }
    // Bytes
    {
        TEST("Bitset Bytes", BITSET_BYTES(1) == 4);
        TEST("Bitset Bytes", BITSET_BYTES(33) == 8);
        TEST("Bitset Bytes", BITSET_BYTES(65) == 12);
        TEST("Bitset Bytes", BITSET_BYTES(97) == 16);
        
    }

    // Single word bitset tests
    // Test
    {
        uint32_t bs[1] = {0};
        for (int i = 0; i < 32; i++) {
            TEST("Bitset Test", !BITSET_TEST(bs, i));
        }

        uint32_t bs2[1] = {0xFFFFFFFF};
        for (int i = 0; i < 32; i++) {
            TEST("Bitset Test", BITSET_TEST(bs2, i));
        }
    }
    // Set
    {
        uint32_t bs[1] = {0};
        for (int i = 0; i < 32; i++) {
            TEST("Bitset Test", !BITSET_TEST(bs,i));
            BITSET_SET(bs, i);
            TEST("Bitset Test", BITSET_TEST(bs,i));
        }
        TEST("Bitset Set", bs[0] == 0xFFFFFFFF);
    }
    // Clear
    {
        uint32_t bs[1] = {0xFFFFFFFF};
        for (int i = 0; i < 32; i++) {
            TEST("Bitset Test", BITSET_TEST(bs,i));
            BITSET_CLEAR(bs, i);
            TEST("Bitset Test", !BITSET_TEST(bs,i));
        }
        TEST("Bitset Clear", bs[0] == 0);
    }
    // Flip
    {
        uint32_t bs[1] = {0};
        for (int i = 0; i < 32; i++) {
            TEST("Bitset Test", !BITSET_TEST(bs,i));
            BITSET_FLIP(bs, i);
            TEST("Bitset Test", BITSET_TEST(bs,i));
        }
        TEST("Bitset Flip", bs[0] == 0xFFFFFFFF);

        uint32_t bs2[1] = {0xFFFFFFFF};
        for (int i = 0; i < 32; i++) {
            TEST("Bitset Test", BITSET_TEST(bs2,i));
            BITSET_CLEAR(bs2, i);
            TEST("Bitset Test", !BITSET_TEST(bs2,i));
        }
        TEST("Bitset Flip", bs2[0] == 0);
    }

    // Multi word bitset tests
    // Test
    {
        uint32_t bs[BITSET_SIZE(128)];
        memset(bs, 0, BITSET_BYTES(128));
        for (int i = 0; i < 128; i++) {
            TEST("Bitset Test", !BITSET_TEST(bs, i));
        }

        uint32_t bs2[BITSET_SIZE(128)];
        memset(bs, 0xFF, BITSET_BYTES(128));
        for (int i = 0; i < 32; i++) {
            TEST("Bitset Test", BITSET_TEST(bs, i));
        }
    }
    // Set
    {
        uint32_t bs[BITSET_SIZE(128)];
        memset(bs, 0, BITSET_BYTES(128));
        for (int i = 0; i < 128; i++) {
            TEST("Bitset Test", !BITSET_TEST(bs,i));
            BITSET_SET(bs, i);
            TEST("Bitset Test", BITSET_TEST(bs,i));
        }

        uint32_t bs2[BITSET_SIZE(128)];
        memset(bs2, 0xFF, BITSET_BYTES(128));
        TEST("Bitset Set", memcmp(bs, bs2, BITSET_BYTES(128))==0);
    }
    // Clear
    {
        uint32_t bs[BITSET_SIZE(128)];
        memset(bs, 0xFF, BITSET_BYTES(128));
        for (int i = 0; i < 128; i++) {
            TEST("Bitset Test", BITSET_TEST(bs,i));
            BITSET_CLEAR(bs, i);
            TEST("Bitset Test", !BITSET_TEST(bs,i));
        }
        uint32_t bs2[BITSET_SIZE(128)];
        memset(bs2, 0, BITSET_BYTES(128));
        TEST("Bitset Clear", memcmp(bs, bs2, BITSET_BYTES(128))==0);
    }
    // Flip
    {
        uint32_t bs[BITSET_SIZE(128)];
        memset(bs, 0, BITSET_BYTES(128));
        for (int i = 0; i < 128; i++) {
            TEST("Bitset Test", !BITSET_TEST(bs,i));
            BITSET_FLIP(bs, i);
            TEST("Bitset Test", BITSET_TEST(bs,i));
        }
        uint32_t bs2[BITSET_SIZE(128)];
        memset(bs2, 0xFF, BITSET_BYTES(128));
        TEST("Bitset Flip", memcmp(bs, bs2, BITSET_BYTES(128))==0);

        for (int i = 0; i < 128; i++) {
            TEST("Bitset Test", BITSET_TEST(bs2,i));
            BITSET_CLEAR(bs2, i);
            TEST("Bitset Test", !BITSET_TEST(bs2,i));
        }
        uint32_t bs3[BITSET_SIZE(128)];
        memset(bs3, 0, BITSET_BYTES(128));
        TEST("Bitset Flip", memcmp(bs2, bs3, BITSET_BYTES(128))==0);
    }
    // Next Set Bit
    {
        uint32_t bs[BITSET_SIZE(128)];
        memset(bs, 0, BITSET_BYTES(128));

        TEST("Bitset next set bit", bitset_next_set_bit(bs, 0, 128)==-1);

        BITSET_SET(bs, 0);
        BITSET_SET(bs, 3);
        BITSET_SET(bs, 40);
        BITSET_SET(bs, 101);


        TEST("Bitset next set bit", bitset_next_set_bit(bs, 0, 128)==0);
        TEST("Bitset next set bit", bitset_next_set_bit(bs, 1, 128)==3);
        TEST("Bitset next set bit", bitset_next_set_bit(bs, 2, 128)==3);
        TEST("Bitset next set bit", bitset_next_set_bit(bs, 3, 128)==3);
        TEST("Bitset next set bit", bitset_next_set_bit(bs, 4, 128)==40);
        TEST("Bitset next set bit", bitset_next_set_bit(bs, 40, 128)==40);
        TEST("Bitset next set bit", bitset_next_set_bit(bs, 41, 128)==101);
        TEST("Bitset next set bit", bitset_next_set_bit(bs, 101, 128)==101);
        TEST("Bitset next set bit", bitset_next_set_bit(bs, 102, 128)==-1);
        
    }
    return 0;
}
