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

#include "autil.h"

int
main(int argc, char *argv[])
{
   float a1[] = {1.0f, 1.0f, 1.0f}; 
   float a2[] = {2.0f, 2.0f, 2.0f}; 
   float a3[] = {3.0f, 3.0f, 3.0f}; 
   float a4[] = {4.0f, 4.0f, 4.0f}; 

    mix_scale_stereo(a1, a2, 1.0f, a3, a4, 3);

    return 0;
}
