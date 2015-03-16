# How much memory will a wavetable take up?

num_waves = 4       # Sin, Saw, Sqr, Tri
num_tables = 128    # one for each note
sample_rate = 48000
sample_sz = 4       # 32bit floats

def table_len(f_s, f_min):
    return int(f_s/f_min + .5)

def midi_freq(note):
    return 440.0*(2.0**((69.0-note)/12.0))

def print_mem_use(b):
    kb = b/1024.0
    mb = kb/1024.0
    print "\t",b, "bytes"
    print "\t",kb, "kilobytes"
    print "\t",mb, "megabytes"

print "One table per note:"
mem = 0
for n in xrange(num_tables):
    mem += sample_sz*table_len(sample_rate, midi_freq(n))
print_mem_use(mem)

print "One table per 8ve:"
octaves = 128/12
mem = 0
for o in xrange(octaves):
    mem += sample_sz*table_len(sample_rate, midi_freq(o*12))
print_mem_use(mem)

