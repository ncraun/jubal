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

#ifndef JUBAL_ENGINE_H
#define JUBAL_ENGINE_H
#include "lfqueue.h"

enum Au_BufferType {L_IN=0, R_IN=1, L_OUT=2, R_OUT=3};

/* Should be inline */
// float *au_buffer(float *buffer, int node, int buff_size, enum Au_BufferType type);
/* TODO: reduce multiplaction distributed */
#define au_buffer(BUFFER, NODE, BUFF_SIZE, TYPE) ((BUFFER)+(NODE)*(BUFF_SIZE)*4 + (TYPE)*(BUFF_SIZE))

struct RtMsg;
typedef struct RtMsg RtMsg;

/* Messages for the non RT Free thread */

struct FreeMsg {
    void *data;
    void (*free)(void *data);
};

typedef struct FreeMsg FreeMsg;

struct AudioNode {
    uint8_t def;
    //void *pdef;
    void *data;
};

typedef struct AudioNode AudioNode;


/* Messages for the Realtime thread */

enum RtNodeMsgType {NOTE_ON, NOTE_OFF, CC, CHANGE_DEF};

struct RtNodeMsg {
    enum RtNodeMsgType type;
    /* TODO: Replace w/ union */
    int node;
    int note;
    int time;
    int param;
    int value;
};

typedef struct RtNodeMsg RtNodeMsg;


enum RtSysMsgType {CHANGE_SAMPLE_RATE, CHANGE_BUFFER_SIZE, CHANGE_NUM_NODES};

struct RtSysMsg {
    enum RtSysMsgType type;
    /* Change Buffer Size Data */
    int new_buffer_size;
    /* This data depends on Buffer Size and Num Nodes */
    float *new_audio_buffer;

    /* Change Sample Rate Data */
    int new_sample_rate;

    /* Change Num Nodes Data */
    AudioNode *new_nodes;
    uint32_t *new_route;
    int new_num_nodes;
    uint32_t *new_seen;
    RtNodeMsg *new_mailbox;
    int *new_mailbox_size;
};

typedef struct RtSysMsg RtSysMsg;

enum RtGraphMsgType {CONNECT_NODE, DISCONNECT_NODE, CONNECT_OUT, DISCONNECT_OUT, MUTE_NODE, UNMUTE_NODE, UNMUTE_ALL, SOLO_NODE, UNSOLO_NODE, UNSOLO_ALL};

struct RtGraphMsg {
    enum RtGraphMsgType type;
    int node_out;
    int node_in;
    int node;
    bool swap;
};

typedef struct RtGraphMsg RtGraphMsg;

struct RtDefMsg {
    int node;
    int new_def;
    void *new_instance;
};

typedef struct RtDefMsg RtDefMsg;

struct ParamHints {
    char *name;
#if 0
    /* TODO: Find a good way to indicate an envelope should be drawn. */
    max_value;
    min_value;
    default_value;
    scale; // Logarithmic, Linear, Exponential
    control_type;
#endif
};

typedef struct ParamHints ParamHints;

struct AudioDef {
    const char *name;

    int num_controls;
    ParamHints *controls;

    /* Interface */
    void* (*make)(void);
    void (*init)(void);
    void (*eval)(void *instance, float* left_input, float *right_input, float *left_output, float *right_output, int len, int sample_rate, int num_msg, RtNodeMsg* msgs);
    void (*free)(void *instance);

    void (*read)(int size, uint8_t *data);
    void (*write)(int size, uint8_t *data);
};

typedef struct AudioDef AudioDef;

struct AuBuff {
    float *left;
    float *right;
};

typedef struct AuBuff AuBuff;

struct EngineRT {

    int *node_order;
    uint32_t *out_conn;
	uint32_t *mute;
	uint32_t *solo;
	int num_solo;

    AudioDef **defs;
    int num_defs;
    AudioNode *nodes;
    int num_nodes;
    int queue_size;
    int buffer_size;
    int sample_rate;

    // AuBuff inputs;
    // AuBuff outputs;

    float *audio_buffer;

    /* Temporary Memory used in functions */
    /* temp mem used in cycle checking and traversal */
    uint32_t *seen;
    uint32_t *visiting;
    uint32_t *done;

    RtNodeMsg *mailbox;
    int *mailbox_size;
};
typedef struct EngineRT EngineRT;

struct EngineNonRT {
    int *node_order_swap;
    /* Maybe replace this with adjacency list? */
    uint32_t *route;
    uint32_t *connected_outputs;
};
typedef struct EngineNonRT EngineNonRT;

struct Engine {

    /* Should be read/written to by RT thread only */
    EngineRT rt;
    /* Should be read/written to by non-RT thread only */
    EngineNonRT non_rt;

    /* Shared Communication? */
    bool paused;
    /* System */
    Lfqueue *sys_q;
    /* Graph */
    Lfqueue *graph_q;

    /* Params and Notes */
    Lfqueue *imm_q;
    Lfqueue *seq_q;

    /* Def */
    Lfqueue *def_q;

    Lfqueue *free_q;
};

typedef struct Engine Engine;

/*
    Engine interface.
*/
/*
    Engine Message Types

    Node:
    Note On
    Note Off
    Node CC

    Connections:
    Change Node Type
    Connect Node
    Disconnect Node

    System:

    Pause
    Resume
    *Change Buffer Size
    *Change Sample Rate
    *Change Number of Nodes

    * Engine must be in paused state to apply these actions
*/

/**/
void engine_note_on(Engine *engine, int node, int note, int time);
void engine_note_off(Engine *engine, int node, int note, int time);
void engine_node_cc(Engine *engine, int node, int param, int value, int time);
void engine_change_node_def(Engine *engine, int node, int def);

/* TODO: How can we reply back about cycle detection? */
/* Connections */
bool engine_connect_nodes(Engine *engine, int node_out, int node_in);
bool engine_disconnect_nodes(Engine *engine, int node_out, int node_in);
bool engine_connect_output(Engine *engine, int node);
bool engine_disconnect_output(Engine *engine, int node);
void engine_mute_node(Engine *engine, int node);
void engine_unmute_node(Engine *engine, int node);
void engine_solo_node(Engine *engine, int node);
void engine_unsolo_node(Engine *engine, int node);
void engine_unmute_all_nodes(Engine *engine);
void engine_unsolo_all_nodes(Engine *engine);
bool engine_remove_node(Engine *engine);
bool engine_add_node(Engine *engine, int def);

/* System */
void engine_pause(Engine *engine);
void engine_play(Engine *engine);
void engine_change_buffer_size(Engine *engine, int buffer_size);
void engine_change_num_nodes(Engine *engine, int num_nodes);
void engine_change_sample_rate(Engine *engine, int sample_rate);

/* Functions */
bool topo_sort_helper(Engine *engine, int curr, uint32_t *visiting, uint32_t *done, int *node_order, int *pos);
bool topo_sort(Engine *engine, int *node_order);
int jubal_callback(AuBuff out, unsigned long frames, Engine *engine, unsigned long start_time);

Engine* make_engine(int num_nodes, int buffer_size, int queue_size, int sample_rate);
void engine_free(Engine *engine);
#endif
