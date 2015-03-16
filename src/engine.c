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

#include <stdbool.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>

#include "engine.h"
#include "lfqueue.h"
#include "bitset.h"
#include "nodes.h"
#include "autil.h"
#include "wave_table.h"

#define AUDIO_OUT_NODE 0
#define DISABLED_DEF 0

static AudioDef *audio_defs[] = {&def_disabled, &def_disabled, &def_fm, &def_sub};

Engine*
make_engine(int num_nodes, int buffer_size, int queue_size, int sample_rate)
{
    Engine *engine = malloc(sizeof(Engine));

    engine->paused = true;

    engine->rt.defs = audio_defs;
    engine->rt.num_defs = sizeof(audio_defs);

    engine->rt.nodes = malloc(sizeof(AudioNode)*num_nodes);
    for (int i = 0; i < num_nodes; i++) {
        engine->rt.nodes[i].def = 0;
        engine->rt.nodes[i].data = NULL;
    }

    engine->non_rt.route = malloc(num_nodes*BITSET_BYTES(num_nodes));
    memset(engine->non_rt.route, 0, num_nodes*BITSET_BYTES(num_nodes));

    engine->non_rt.connected_outputs = malloc(num_nodes*BITSET_BYTES(num_nodes));
    memset(engine->non_rt.connected_outputs, 0, num_nodes*BITSET_BYTES(num_nodes));

    engine->rt.audio_buffer = malloc(num_nodes*sizeof(float)*buffer_size*4);
    memset(engine->rt.audio_buffer, 0, num_nodes*sizeof(float)*buffer_size*4);

    engine->rt.node_order = malloc(num_nodes*sizeof(int));
    memset(engine->rt.node_order, 0, num_nodes*sizeof(int));

    engine->non_rt.node_order_swap = malloc(num_nodes*sizeof(int));
    memset(engine->non_rt.node_order_swap, 0, num_nodes*sizeof(int));

    engine->rt.num_nodes = num_nodes;
    engine->rt.queue_size = queue_size;
    engine->rt.buffer_size = buffer_size;
    engine->rt.sample_rate = sample_rate;

    engine->imm_q = make_lfqueue(queue_size, sizeof(RtNodeMsg));
    engine->seq_q = make_lfqueue(queue_size, sizeof(RtNodeMsg));
    engine->free_q = make_lfqueue(queue_size, sizeof(FreeMsg));
    engine->graph_q = make_lfqueue(queue_size, sizeof(RtGraphMsg));
    engine->sys_q = make_lfqueue(queue_size, sizeof(RtSysMsg));
    engine->def_q = make_lfqueue(queue_size, sizeof(RtDefMsg));

    engine->rt.seen = malloc(BITSET_BYTES(num_nodes));
    memset(engine->rt.seen, 0, BITSET_BYTES(num_nodes));
    engine->rt.visiting = malloc(BITSET_BYTES(num_nodes));
    memset(engine->rt.visiting, 0, BITSET_BYTES(num_nodes));
    engine->rt.done = malloc(BITSET_BYTES(num_nodes));
    memset(engine->rt.done, 0, BITSET_BYTES(num_nodes));

    engine->rt.mailbox = malloc(sizeof(RtNodeMsg)*queue_size*num_nodes);

    engine->rt.mailbox_size = malloc(sizeof(int)*num_nodes);
    memset(engine->rt.mailbox_size, 0, sizeof(int)*num_nodes);

    engine->rt.out_conn = malloc(BITSET_BYTES(num_nodes));
    memset(engine->rt.out_conn, 0, BITSET_BYTES(num_nodes));

    engine->rt.mute = malloc(BITSET_BYTES(num_nodes));
    memset(engine->rt.mute, 0, BITSET_BYTES(num_nodes));

    engine->rt.solo = malloc(BITSET_BYTES(num_nodes));
    memset(engine->rt.solo, 0, BITSET_BYTES(num_nodes));

	engine->rt.num_solo = 0;

	// TODO: Handle the wavetable initialization better...
	// fill_saw_wave_table(&G_SAW_TABLE, sample_rate);
	
    return engine;
}

void
engine_free(Engine *engine)
{
    for (int i = 0; i < engine->rt.num_nodes; i++) {
        AudioNode *curr = &(engine->rt.nodes[i]);
        AudioDef *curr_def = engine->rt.defs[curr->def];
        curr_def->free(curr->data);
    }

    free_lfqueue(engine->imm_q);
    free_lfqueue(engine->seq_q);
    free_lfqueue(engine->free_q);
    free_lfqueue(engine->sys_q);
    free_lfqueue(engine->graph_q);
    free_lfqueue(engine->def_q);

    free(engine->rt.nodes);
    free(engine->non_rt.route);
    free(engine->non_rt.connected_outputs);
    free(engine->rt.audio_buffer);
    free(engine->rt.node_order);
    free(engine->non_rt.node_order_swap);
    free(engine->rt.seen);
    free(engine->rt.visiting);
    free(engine->rt.done);
    free(engine->rt.mailbox);
    free(engine->rt.mailbox_size);
    free(engine->rt.out_conn);
    free(engine);
}

/*
    Topological Sort
    Returns false if there was a cycle, true if the topo-sort worked correctly.
*/
bool
topo_sort_helper(Engine *engine, int curr, uint32_t *visiting, uint32_t *done, int *node_order, int *pos)
{
    //printf("Traverse: %d\n", node);
    uint32_t *connections = engine->non_rt.route + curr*BITSET_SIZE(engine->rt.num_nodes);
    int next = bitset_next_set_bit(connections, 0, engine->rt.num_nodes);

    BITSET_SET(visiting, curr);

    while (next >= 0) {
        if (BITSET_TEST(visiting, next)) {
            return false;
        }

        if (!BITSET_TEST(done, next)) {
            if (!topo_sort_helper(engine, next, visiting, done, node_order, pos)) {
                return false;
            }

            next = bitset_next_set_bit(connections, next+1, engine->rt.num_nodes);
        }
    }

    /*
        All children of curr have been processed.
    */

    BITSET_SET(done, curr);
    BITSET_CLEAR(visiting, curr);

    node_order[*pos] = curr;
    (*pos)++;

    return true;
}

bool
topo_sort(Engine *engine, int *node_order)
{
    // int next = bitset_next_set_bit(engine->rt.out_conn, 0, engine->rt.num_nodes);
    int pos = 0;

    uint32_t visiting[BITSET_SIZE(engine->rt.num_nodes)];
    uint32_t done[BITSET_SIZE(engine->rt.num_nodes)];

    /* Can also be done with = {0}; */
    memset(visiting, 0, BITSET_BYTES(engine->rt.num_nodes));
    memset(done, 0, BITSET_BYTES(engine->rt.num_nodes));

    /*
        while (next >= 0) {
        if (!topo_sort_helper(engine, next, visiting, done, node_order, &pos)) {
            return false;
        }
        next = bitset_next_set_bit(engine->rt.out_conn, next+1, engine->rt.num_nodes);
    }
    */

    for (int i = 0; i < engine->rt.num_nodes; i++) {
        /* TODO: replace with an implementation of bitset_next_clear_bit */
        if (!BITSET_TEST(done, i)) {
            if (!topo_sort_helper(engine, i, visiting, done, node_order, &pos)) {
                return false;
            }
        }
    }

    return true;
}

void
handle_msg(Engine *engine)
{
}

int
jubal_callback(AuBuff out, unsigned long frames, Engine *engine, unsigned long start_time)
{
    if (engine->paused) {
        memset(out.left, 0, sizeof(float)*frames);
        memset(out.right, 0, sizeof(float)*frames);
		// printf("Paused\n");
        return 0;
    }

    int buffer_end_time = frames/engine->rt.sample_rate;

    memset(engine->rt.mailbox_size, 0, engine->rt.num_nodes*sizeof(int));
    /* TODO: When polling the event queue, adjust event times to samples within the frame using the start_time. */

    while (lfqueue_size(engine->sys_q) > 0) {
        RtSysMsg tmp;
        lfqueue_poll(engine->sys_q, &tmp);
        if (tmp.type == CHANGE_SAMPLE_RATE) {
            engine->rt.sample_rate = tmp.new_sample_rate;
            /* TODO: Update time stamps on sample rate change. */
        }
        else if (tmp.type == CHANGE_BUFFER_SIZE) {
            FreeMsg free_old_buff;
            free_old_buff.data = engine->rt.audio_buffer;
            free_old_buff.free = free;

            engine->rt.audio_buffer = tmp.new_audio_buffer;
            engine->rt.buffer_size = tmp.new_buffer_size;

            lfqueue_add(engine->free_q, &free_old_buff);
        }
        else if (tmp.type == CHANGE_NUM_NODES) {
        /* TODO: not supported yet */
        /*
            engine->rt.inputs.left = tmp.data.sys.new_left_inputs;
            engine->rt.inputs.right = tmp.data.sys.new_right_inputs;
            engine->rt.outputs.left = tmp.data.sys.new_left_outputs;
            engine->rt.outputs.right = tmp.data.sys.new_right_outputs;
                engine->rt.nodes = tmp.data.sys.new_nodes;
                engine->rt.seen = tmp.data.sys.new_seen;
                engine->rt.mailbox = tmp.data.sys.new_mailbox;
                engine->rt.mailbox_size = tmp.data.sys.new_mailbox_size;
                engine->rt.num_nodes = tmp.data.sys.new_num_nodes;
            */
        }
    }

    while (lfqueue_size(engine->graph_q) > 0) {
        RtGraphMsg tmp_msg;
        lfqueue_poll(engine->graph_q, &tmp_msg);

        if (tmp_msg.type == CONNECT_NODE) {
            BITSET_SET(engine->non_rt.route+tmp_msg.node_in*BITSET_SIZE(engine->rt.num_nodes), tmp_msg.node_out);
            BITSET_SET(engine->non_rt.connected_outputs+tmp_msg.node_out*BITSET_SIZE(engine->rt.num_nodes), tmp_msg.node_in);
            if (tmp_msg.swap) {
                int *tmp = engine->rt.node_order;
                engine->rt.node_order = engine->non_rt.node_order_swap;
                engine->non_rt.node_order_swap = tmp;
            }
        }
        else if (tmp_msg.type == DISCONNECT_NODE) {
            BITSET_CLEAR(engine->non_rt.route+tmp_msg.node_in*BITSET_SIZE(engine->rt.num_nodes), tmp_msg.node_out);
            BITSET_CLEAR(engine->non_rt.connected_outputs+tmp_msg.node_out*BITSET_SIZE(engine->rt.num_nodes), tmp_msg.node_in);
            if (tmp_msg.swap) {
                int *tmp = engine->rt.node_order;
                engine->rt.node_order = engine->non_rt.node_order_swap;
                engine->non_rt.node_order_swap = tmp;
            }
        }
        else if (tmp_msg.type == CONNECT_OUT) {
            BITSET_SET(engine->rt.out_conn, tmp_msg.node);
            printf("Connect %d to out.\n", tmp_msg.node);
            if (tmp_msg.swap) {
                int *tmp = engine->rt.node_order;
                engine->rt.node_order = engine->non_rt.node_order_swap;
                engine->non_rt.node_order_swap = tmp;
            }
        }
        else if (tmp_msg.type == DISCONNECT_OUT) {
            BITSET_CLEAR(engine->rt.out_conn, tmp_msg.node);
            if (tmp_msg.swap) {
                int *tmp = engine->rt.node_order;
                engine->rt.node_order = engine->non_rt.node_order_swap;
                engine->non_rt.node_order_swap = tmp;
            }
        }
		else if (tmp_msg.type == MUTE_NODE) {
			BITSET_CLEAR(engine->rt.solo, tmp_msg.node);
			BITSET_SET(engine->rt.mute, tmp_msg.node);
		}
		else if (tmp_msg.type == UNMUTE_NODE) {
			BITSET_CLEAR(engine->rt.mute, tmp_msg.node);
		}
		else if (tmp_msg.type == UNMUTE_ALL) {
			memset(engine->rt.mute, 0, BITSET_BYTES(engine->rt.num_nodes));
		}
		else if (tmp_msg.type == SOLO_NODE) {
			BITSET_CLEAR(engine->rt.mute, tmp_msg.node);
			BITSET_SET(engine->rt.solo, tmp_msg.node);
			if (engine->rt.num_solo < engine->rt.num_nodes) {
				engine->rt.num_solo++;
			}
		}
		else if (tmp_msg.type == UNSOLO_NODE) {
			BITSET_CLEAR(engine->rt.solo, tmp_msg.node);
			if (engine->rt.num_solo > 0) {
				engine->rt.num_solo--;
			}
		}
		else if (tmp_msg.type == UNSOLO_ALL) {
			memset(engine->rt.solo, 0, BITSET_BYTES(engine->rt.num_nodes));
			engine->rt.num_solo = 0;
		}
    }

    /* NOTE: Changing definition of the node will be done before node messaging. */

    while (lfqueue_size(engine->def_q) > 0) {
        RtDefMsg tmp_msg;
        lfqueue_poll(engine->def_q, &tmp_msg);

        int new_def = tmp_msg.new_def;
        int node = tmp_msg.node;
        void *new_instance = tmp_msg.new_instance;
        void *old_instance = engine->rt.nodes[node].data;
        int old_def = engine->rt.nodes[node].def;

        engine->rt.nodes[node].def = new_def;
        engine->rt.nodes[node].data = new_instance;

        FreeMsg fi;
        fi.data = old_instance;
        fi.free = engine->rt.defs[old_def]->free;
        lfqueue_add(engine->free_q, &fi);
    }

    while (lfqueue_size(engine->imm_q) > 0) {
        RtNodeMsg tmp_msg;
        lfqueue_poll(engine->imm_q, &tmp_msg);

        /* Adjust time to 0 so it is executed at the start of the buffer. */
        tmp_msg.time = 0;
        if (engine->rt.mailbox_size[tmp_msg.node] < engine->rt.queue_size) {
            engine->rt.mailbox[tmp_msg.node*(engine->rt.queue_size) + engine->rt.mailbox_size[tmp_msg.node]] = tmp_msg;
            engine->rt.mailbox_size[tmp_msg.node]++;
        }
    }

    while (lfqueue_size(engine->seq_q) > 0) {
        RtNodeMsg tmp_msg;
        lfqueue_peek(engine->seq_q, &tmp_msg);
        if (tmp_msg.time > buffer_end_time) {
            break;
        }

        lfqueue_advance(engine->seq_q);

        if (tmp_msg.time < start_time) {
            continue;
        }

        /* Adjust time to be sample offset into the period. */
        tmp_msg.time -= start_time;

        if (engine->rt.mailbox_size[tmp_msg.node] < engine->rt.queue_size) {
            engine->rt.mailbox[tmp_msg.node*(engine->rt.queue_size) + engine->rt.mailbox_size[tmp_msg.node]] = tmp_msg;
            engine->rt.mailbox_size[tmp_msg.node]++;
        }
    }


    /* Zero ALL Buffers */
	/* Actually this is kind of slow, we should only zero the buffers we need to */
    memset(out.left, 0, sizeof(float)*frames);
    memset(out.right, 0, sizeof(float)*frames);
/*
    memset(engine->rt.audio_buffer, 0, engine->rt.buffer_size*4*engine->rt.num_nodes*sizeof(float));
*/

    /* Evaluate nodes */
	/* TODO: In cases of mute, solo, disabled def, have a way to skip to 
		the next nonmuted, not disabled, or solo node, instead of iterating through */
    for (int i = 0; i < engine->rt.num_nodes; i++) {
        int node = engine->rt.node_order[i];
        AudioNode *curr = &(engine->rt.nodes[node]);
        if (curr->def == DISABLED_DEF) {
            continue; // Don't eval empty nodes
        }

		if (BITSET_TEST(engine->rt.mute, node)) {
			// printf("Node: %d muted\n", node);
			continue; // Don't eval muted nodes
		}

		if (engine->rt.num_solo != 0 && !BITSET_TEST(engine->rt.solo, node)) {
			// printf("Node: %d not solo\n", node);
			continue; // Don't eval if solo mode is active, and curr node not solo.
		}

        AudioDef *curr_def = engine->rt.defs[curr->def];
        float *curr_left_input = au_buffer(engine->rt.audio_buffer, node, engine->rt.buffer_size, L_IN);
        float *curr_right_input = au_buffer(engine->rt.audio_buffer, node, engine->rt.buffer_size, R_IN);
        float *curr_left_output = au_buffer(engine->rt.audio_buffer, node, engine->rt.buffer_size, L_OUT);
        float *curr_right_output = au_buffer(engine->rt.audio_buffer, node, engine->rt.buffer_size, R_OUT);
		// Zero buffers
		memset(curr_left_input, 0, engine->rt.buffer_size*sizeof(float));
		memset(curr_right_input, 0, engine->rt.buffer_size*sizeof(float));
		memset(curr_left_output, 0, engine->rt.buffer_size*sizeof(float));
		memset(curr_right_output, 0, engine->rt.buffer_size*sizeof(float));
		// printf("Calling eval for node %d\n", node);
        curr_def->eval(curr_def->data, curr->data, curr_left_input, curr_right_input, curr_left_output, curr_right_output, frames, engine->rt.sample_rate, engine->rt.mailbox_size[node], engine->rt.mailbox+node*engine->rt.queue_size);

        /* Mix Evaulated node into its outputs */

        int next = bitset_next_set_bit(engine->non_rt.connected_outputs+i*BITSET_SIZE(engine->rt.num_nodes), 0, engine->rt.num_nodes);
        while (next >= 0) {
            float *next_left_input = au_buffer(engine->rt.audio_buffer, next, engine->rt.buffer_size, L_IN);
            float *next_right_input = au_buffer(engine->rt.audio_buffer, next, engine->rt.buffer_size, R_IN);
            /* Mix node's output to all other nodes that have the current node as an input. */
            // 2 channels for stereo
            // Could benefit from sse2
            // Mix curr_l/r output into next_l/r input
            mix_scale_stereo(curr_left_output, curr_right_output, 1.0f, next_left_input, next_right_input, frames);

            next = bitset_next_set_bit(engine->non_rt.connected_outputs+i*BITSET_SIZE(engine->rt.num_nodes), next, engine->rt.num_nodes);
        }
		// printf("Evaluated Node: %d\n", node);
    }
    /* All Nodes have been evaluated. Mix the connected ones into the system output. */
    int node_out = bitset_next_set_bit(engine->rt.out_conn, 0, engine->rt.num_nodes);
    while (node_out >= 0) {
        float *curr_left_output = au_buffer(engine->rt.audio_buffer, node_out, engine->rt.buffer_size, L_OUT);
        float *curr_right_output = au_buffer(engine->rt.audio_buffer, node_out, engine->rt.buffer_size, R_OUT);

        mix_scale_stereo(curr_left_output, curr_right_output, 1.0f, out.left, out.right, frames);

        node_out = bitset_next_set_bit(engine->rt.out_conn, node_out+1, engine->rt.num_nodes);
    }

    return 0;
}

/*
    While the engine is paused the queue isn't being read, so we can't just
    send a message on the queue?

    Non-RT functions
*/
void
engine_pause(Engine *engine)
{
    __sync_bool_compare_and_swap(&engine->paused, false, true);
}

void
engine_play(Engine *engine)
{
    __sync_bool_compare_and_swap(&engine->paused, true, false);
}

void
engine_change_buffer_size(Engine *engine, int buffer_size)
{
    if (buffer_size > engine->rt.buffer_size) {
        RtSysMsg msg;
        msg.type = CHANGE_BUFFER_SIZE;
        msg.new_buffer_size = buffer_size;
        msg.new_audio_buffer = malloc(sizeof(float)*buffer_size*engine->rt.num_nodes*4);
        lfqueue_add(engine->sys_q, &msg);
    }
}

void
engine_change_sample_rate(Engine *engine, int sample_rate)
{
    printf("Change sample rate from %d to %d.\n", engine->rt.sample_rate, sample_rate);
    RtSysMsg msg;
    msg.type = CHANGE_SAMPLE_RATE;
    msg.new_sample_rate = sample_rate;

/* TODO: A better way of refilling the wavetables when samplerate changes. */
	fill_saw_wave_table(&G_SAW_TABLE, sample_rate);
    lfqueue_add(engine->sys_q, &msg);
}

void
engine_change_num_nodes(Engine *engine, int num_nodes)
{
    /*
    RtMsg msg;
    msg.type = SYS;
    msg.data.sys.type = CHANGE_NUM_NODES;
    msg.data.sys.new_num_nodes = num_nodes;
    msg.data.sys.new_left_inputs = malloc(sizeof(float)*buffer_size*engine->num_nodes);
    msg.data.sys.new_right_inputs = malloc(sizeof(float)*buffer_size*engine->num_nodes);
    msg.data.sys.new_left_outputs = malloc(sizeof(float)*buffer_size*engine->num_nodes);
    msg.data.sys.new_right_outputs = malloc(sizeof(float)*buffer_size*engine->num_nodes);
    msg.data.sys.new_nodes = malloc(sizeof());
    memcpy(msg.data.sys.new_nodes,
    msg.data.sys.new_route = ;
    msg.data.sys.new_seen = ;
    msg.data.sys.new_mailbox = ;
    msg.data.sys.new_mailbox_size = ;

    // Also reallocate the node_order and node_order_swap buffers.
    */
    fprintf(stderr, "Changing Num Nodes not supported (yet)!\n");
}

void
engine_note_on(Engine *engine, int node, int note, int time)
{
    RtNodeMsg msg;

    msg.type = NOTE_ON;
    msg.node = node;
    msg.time = time;
    msg.note = note;

    if (time == 0) { lfqueue_add(engine->imm_q, &msg); }
    else { lfqueue_add(engine->seq_q, &msg); }
}

void
engine_note_off(Engine *engine, int node, int note, int time)
{
    RtNodeMsg msg;
    msg.type = NOTE_OFF;
    msg.node = node;
    msg.time = time;
    msg.note = note;

    if (time == 0) { lfqueue_add(engine->imm_q, &msg); }
    else { lfqueue_add(engine->seq_q, &msg); }
}

void
engine_node_cc(Engine *engine, int node, int param, int value, int time)
{
    RtNodeMsg msg;
    msg.type = CC;
    msg.node = node;
    msg.time = time;
    msg.param = param;
    msg.value = value;

    if (time == 0) { lfqueue_add(engine->imm_q, &msg); }
    else { lfqueue_add(engine->seq_q, &msg); }
}

void
engine_change_node_def(Engine *engine, int node, int def)
{
    RtDefMsg msg;
    msg.node = node;
    msg.new_def = def;
    msg.new_instance = engine->rt.defs[def]->make(engine->rt.sample_rate);

    lfqueue_add(engine->def_q, &msg);
}

bool
engine_connect_nodes(Engine *engine, int node_out, int node_in)
{
    /* TODO: Cycle Check */
    RtGraphMsg msg;
    msg.type = CONNECT_NODE;
    msg.node_out = node_out;
    msg.node_in = node_in;
    msg.swap = true;
    /* TODO: Might not be safe. */
    memcpy(engine->non_rt.node_order_swap, engine->rt.node_order, engine->rt.num_nodes*sizeof(int));
    if (topo_sort(engine, engine->non_rt.node_order_swap)) {
        lfqueue_add(engine->graph_q, &msg);
        return true;
    }
    return false;
}

bool
engine_disconnect_nodes(Engine *engine, int node_out, int node_in)
{
    RtGraphMsg msg;
    msg.type = DISCONNECT_NODE;
    msg.node_out = node_out;
    msg.node_in = node_in;
    msg.swap = true;

    /* TODO: Might not be safe. */
    memcpy(engine->non_rt.node_order_swap, engine->rt.node_order, engine->rt.num_nodes*sizeof(int));
    if (topo_sort(engine, engine->non_rt.node_order_swap)) {
        lfqueue_add(engine->graph_q, &msg);
        return true;
    }
    return false;
}

bool
engine_connect_output(Engine *engine, int node)
{
    RtGraphMsg msg;
    msg.type = CONNECT_OUT;
    msg.node = node;
    msg.swap = true;

    memcpy(engine->non_rt.node_order_swap, engine->rt.node_order, engine->rt.num_nodes*sizeof(int));
    if (topo_sort(engine, engine->non_rt.node_order_swap)) {
        lfqueue_add(engine->graph_q, &msg);
        return true;
    }
    return false;
}

bool
engine_disconnect_output(Engine *engine, int node)
{
    RtGraphMsg msg;
    msg.type = DISCONNECT_OUT;
    msg.node = node;
    msg.swap = true;

    memcpy(engine->non_rt.node_order_swap, engine->rt.node_order, engine->rt.num_nodes*sizeof(int));
    if (topo_sort(engine, engine->non_rt.node_order_swap)) {
        lfqueue_add(engine->graph_q, &msg);
        return true;
    }
    return false;
}

void
engine_mute_node(Engine *engine, int node)
{
	RtGraphMsg msg;
	msg.type = MUTE_NODE;
	msg.node = node;
	lfqueue_add(engine->graph_q, &msg);
}

void
engine_unmute_node(Engine *engine, int node)
{
	RtGraphMsg msg;
	msg.type = UNMUTE_NODE;
	msg.node = node;
	lfqueue_add(engine->graph_q, &msg);
}

void
engine_unmute_all(Engine *engine)
{
	RtGraphMsg msg;
	msg.type = UNMUTE_ALL;
	lfqueue_add(engine->graph_q, &msg);
}

void
engine_solo_node(Engine *engine, int node)
{
	RtGraphMsg msg;
	msg.type = SOLO_NODE;
	msg.node = node;
	lfqueue_add(engine->graph_q, &msg);
}

void
engine_unsolo_node(Engine *engine, int node)
{
	RtGraphMsg msg;
	msg.type = UNSOLO_NODE;
	msg.node = node;
	lfqueue_add(engine->graph_q, &msg);
}

void
engine_unsolo_all(Engine *engine)
{
	RtGraphMsg msg;
	msg.type = UNSOLO_ALL;
	lfqueue_add(engine->graph_q, &msg);
}
