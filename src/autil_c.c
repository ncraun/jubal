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

void
mix_scale_stereo(const float *l1, const float *r1, float scale, float *l2, float *r2, int len)
{
    /* TODO: check for aliasing issues. */
    for (int i = 0; i < len; i++) {
        l2[i] += scale*l1[i];
    }
    for (int i = 0; i < len; i++) {
        r2[i] += scale*r1[i];
    }
}

/*
	What to do if there's a key off during the attack or decay stage?
*/
float 
eg_lin_adsr(float t, ADSR *adsr, bool key_off, float t_off, float *a_off)
{
	float r;
	if (key_off) {
		if (t > (t_off + adsr->rel_t)) {
			// Off
			r =  -1.0f;
		}
		else {
			// Release
			// r = adsr->sus_l + (adsr->sus_l/adsr->rel_t)*(adsr->t_off - t);
			float m = -*a_off/adsr->rel_t;
			float b = *a_off - m*t_off;
			r = m*t + b;
		}
	}
	else {
		if (t < adsr->atk_t) {
			// Attack
			// Rise/Run
			// t on line from (0, 0) to (atk_t, atk_l)
			r = t*(adsr->atk_l)/(adsr->atk_t);
		}

		else if (t >= adsr->atk_t && t < (adsr->atk_t + adsr->dec_t)) {
			// Decay
			float m = (adsr->sus_l - adsr->atk_l)/adsr->dec_t;
			float b = adsr->atk_l - m*adsr->atk_t;

			r = m*t + b;
		}
		else {
			// Sustain
			r = adsr->sus_l;
		}

		*a_off = r;
	}

	//printf("Amp: r=%f a_off=%f t_off=%f t=%f\n", r, *a_off, t_off, t);
	return r;
}
