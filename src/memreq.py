nodes = 1024
buffer_size = 1024
channels = 2
song_len = 128
pattern_cell = 5
pattern_len = 128
pattern_channels = 64
pattern_size = pattern_cell*pattern_len*pattern_channels
node_block = 272

audio_buffers = 2*nodes*buffer_size*channels
node_store = node_block*nodes
song = pattern_size*song_len

total = audio_buffers + node_store + song

print "nodes:", nodes
print "buffer_size:", buffer_size
print "channels:", channels
print "song_len:",song_len
print "pattern_cell:",pattern_cell
print "pattern_len:",pattern_len
print "pattern_size:",pattern_size
print "node_block:",node_block
print "------------------------------"
print "audio_buffers:",audio_buffers
print "node_store:",node_store
print "song:",song
print "------------------------------"
print "total:",total
