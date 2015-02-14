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
#include <string.h>
#include <stdlib.h>

/*
    reallocations:

    add_tracks():
        realloc(old_size + tracks_added*pattern_length)
        shift_zerofill_tracks()
    del_shift_tracks():
        new = malloc(old_size - tracks_added*pattern_length)
        shift_copy(old, new)
        free(old)

    inc_pattern_length():
        realloc(old_size + num_tracks*lines_added)
        zero_fill(start_addr+old_size,)
    dec_pattern_length()
        realloc(old_size - num_tracks*lines_removed)
*/

int*
add_tracks(int *old_pattern, int pattern_length, int old_num_tracks, int new_num_tracks)
{
    // Shift
    // Before: T1T2 T1T2 T1T2
    // After:  T1T2T3 T1T2T3 T1T2T3
    // So Every old_num_tracks, insert a gap
    // Every old_num_tracks, shift all other elements +1
    
    // This is a slower way, there will be more efficient ways
    // IT would be better to do it inplace.
    int *new_pattern = malloc(new_num_tracks*pattern_length*sizeof(int));
    for (int p = 0; p < pattern_length; p++) {
        // Copy N-1 tracks
        memcpy(new_pattern+p*new_num_tracks, old_pattern+p*old_num_tracks, sizeof(int)*old_num_tracks);
        // Set Nth track to 0
        memset(new_pattern+p*new_num_tracks+old_num_tracks, 0, sizeof(int)*(new_num_tracks-old_num_tracks));
    }
    free(old_pattern);
    return new_pattern;
}

int*
del_track(int *old_pattern, int pattern_length, int old_num_tracks, int track_del_idx)
{
    int *new_pattern = malloc((old_num_tracks-1)*pattern_length*sizeof(int)); 
    int new_num_tracks = old_num_tracks-1;
    // Copy rest of tracks over
    for (int p = 0; p < pattern_length; p++) {
        // Copy tracks [0..D)
        for (int t = 0; t < track_del_idx; t++) {
            new_pattern[p*new_num_tracks+t] = old_pattern[p*old_num_tracks+t];
        }
        // Copy tracks [D+1..End]
        for (int t = track_del_idx+1; t < old_num_tracks; t++) {
            int nt = t-1;
            new_pattern[p*new_num_tracks+nt] = old_pattern[p*old_num_tracks+t];
        }
    }
    free(old_pattern);
    return new_pattern;
}
// should this return ptr, or in param?
int*
change_pattern_length(int *old_pattern, int num_tracks, int old_length, int new_length)
{
    int old_sz = num_tracks*old_length;
    int new_sz = num_tracks*new_length;

    int *new_pattern = realloc(old_pattern, sizeof(int)*new_sz);
    if (new_pattern == NULL) { /* fail */ }
    
    if (new_length > old_length) {
        // zero empty space
        memset(new_pattern+old_sz, 0, sizeof(int)*(new_sz-old_sz));
    }

    return new_pattern;
}

/*
    1 min      60s    1 beat    10**9 ns
    ------- * ----- * ------- * ------
    B beats   1 min   L lines     1s

    60*10**9/(lpb*bpm) seconds per line
*/

/* 
	Transpose 
*/

/*
	Time Stretch/Shrink
*/

/*
	Shift	
*/

/*
	Insert
*/

/*
	Delete
*/
