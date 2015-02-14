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

#ifndef JUBAL_AUTIL_H
#define JUBAL_AUTIL_H

/*
    Interface for audio utility functions. The idea behind this is we can
    have implementations for these basic functions that can take advantage
    of cpu features. For example, we might have a plain C function, or
    SSE2, AVX, NEON, whatever.
*/

/* NOTE: Maybe replace w/ a mix(float*,float* that works on mono streams, and
    let the user call mix(l1, l2) and mix(r1, r2) */
void mix_scale_stereo(const float *l1, const float *r1, float scale, float *l2, float *r2, int len);

#endif
