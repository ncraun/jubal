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

#ifndef JUBAL_TEST_H
#define JUBAL_TEST_H

#include <stdio.h>
#include <stdlib.h>

#define TEST(name, cond) do { if (!(cond)) { fprintf(stderr, "FAIL: %s, %s at %s:%u\n", name, #cond, __FILE__, __LINE__); exit(1); } } while (0)

#endif
