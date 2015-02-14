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

#ifndef JUBAL_SONG_H
#define JUBAL_SONG_H

struct PatternCell {
    uint8_t note;
    uint8_t vel;
    Effect effect;
};

typedef struct PatternCell PatternCell;

struct Pattern {
    int num_tracks;
    int length;
    PatternCell *data;
    /* Map track # to node target */
    int *targets;
};

typedef struct Pattern Pattern;

struct Effect {
    uint8_t param;
    uint8_t value;
};

typedef struct Effect Effect;

struct SongMeta {
    char *title;
    char *author;
    char *description;
    int date;
};

typedef struct SongMeta SongMeta;

struct Song {
    SongMeta meta;
    char *title;
    char *author;
    uint8_t bpm;
    uint8_t lpb;
    uint8_t num_patterns;
    uint8_t song_len;
    uint8_t *song;
    Pattern *patterns;
};

typedef struct Song Song;

#endif
