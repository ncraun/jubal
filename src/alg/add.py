# Addititive Synthesis Efficient Algorithm Comparison
# unfortunately python overhead is taking all the time
# Switch to C

import math
import time

# Direct Sin() approach
f_s = 48000

def direct(amp, buff):
    K = len(amp)
    N = len(buff)
    p = 1.0/f_s

    for n in xrange(N):
        t = p*n
        buff[n] = 0.0
        for k in xrange(K):
            buff[n] += amp[k]*math.sin(k*t)


# For each sample, compute each harmonic
# m1 - Just using the multiplication identity
def m1(amp, buff):
    K = len(amp)
    N = len(buff)
    p = 1.0/f_s

    for n in xrange(N):
        t = p*n
        cos_c = math.cos(t)
        sin_2 = 0
        sin_1 = math.sin(t)
        for k in xrange(K):
            sin = 2*sin_1*cos_c - sin_2
            buff[n] += amp[k]*sin 
            sin_2 = sin_1
            sin_1 = sin

# m2 - using the additive and multiplicative identity
def m2(amp, buff):
    K = len(amp)
    N = len(buff)
    p = 1.0/f_s
    sin_p = math.sin(p)
    cos_p = math.cos(p)
    sin_t = 0.0
    cos_t = 1.0

    for n in xrange(N):
        t = p*n
        sin_2 = 0
        sin_1 = sin_t
        for k in xrange(K):
            sin = 2*sin_1*cos_t - sin_2
            buff[n] += amp[k]*sin 
            sin_2 = sin_1
            sin_1 = sin

        n_sin_t = sin_t*cos_p + sin_p*cos_t
        n_cos_t = cos_t*cos_p - sin_t*sin_p
        cos_t = n_cos_t
        sin_t = n_sin_t

# m3 - Using additive and multiplicative identity

# Use DP table method 1
def dp_m1():
    pass

def bmark(method):
    print method
    amp = [1.0]*1200
    amp[0] = 0.0
    buff = [0]*1024
    t1 = time.clock()
    method(amp, buff)
    t2 = time.clock()
    print t2-t1, "s"

bmark(direct)
bmark(m1)
bmark(m2)
