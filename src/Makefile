CC := ccache gcc
CFLAGS_REL := -std=c99 -Wall -pthread -Os -s 
CFLAGS_DBG := -std=c99 -Wall -pthread -O0 -g 

CFLAGS := $(CFLAGS_REL) $(shell sdl2-config --cflags)

LDLIBS := -ljack -lm $(shell sdl2-config --libs)

OBJS := bitset.o lfqueue.o engine.o nodes.o song.o jack_main.o autil_c.o

.PHONY: clean test

main: $(OBJS)
	$(CC) $(CFLAGS) $^ $(LDLIBS) -o main

bitset_test: bitset.o bitset_test.o
	$(CC) $(CFLAGS) $^ -o $@

lfqueue_test: lfqueue.o lfqueue_test.o
	$(CC) $(CFLAGS) $^ -o $@

test: bitset_test lfqueue_test
	./bitset_test
	./lfqueue_test

clean:
	-rm main $(OBJS) bitset_test lfqueue_test bitset_test.o lfqueue_test.o
