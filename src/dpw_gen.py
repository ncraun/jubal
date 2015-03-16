import math

def dpw(table_size, order):
	# Integrate
	saw2 = [0]*table_size
	for i in xrange(table_size):
		phasor = (i+0.0)/(table_size-1.0)
		saw = 2*phasor - 1
		saw2[i] = saw*saw

	# Differentiate
	diff_saw = [0]*table_size
	for i in xrange(1,table_size):
		diff_saw[i] = saw2[i] - saw2[i-1]

	# What should diff_saw 0 be?
	diff_saw[0] = diff_saw[1]

	# Sample rate = samples/sec = table_size/seconds = table_size/T0
	# f0 = 1/T0
	# P = fs/f0 = (table_size/T0)/(1/T0) = table_size/T0 * T0 = table_size
	P = float(table_size)

	scale = math.pow(math.pi, order-1)/(math.factorial(order)*math.pow(2*math.sin(math.pi/P), order-1))
	# scale = math.pow(P, order-1)/(math.pow(2, order-1)*math.factorial(order))

	# Scale
	for i in xrange(table_size):
		diff_saw[i] *= scale
		print diff_saw[i]

# DPW equations
# 1: x
# 2: x^2
# 3: x^3 - x
# 4: x^4 - 2x^2
# 5: x^5 - (10/3)x^3 + (7/3)x
# 6: x^6 - 5x^4 + 7x^2

def dpw_online(table_size, order):
	table_size = table_size + order
	# "On-Line" algorithm no caching
	# needs a 1 sample delay
	# or an order sample delay?

	phase_inc = 1.0/(table_size-1.0)
	phasor = 0.0
	P = float(table_size)
	scale = math.pow(math.pi, order-1)/(math.factorial(order)*math.pow(2*math.sin(math.pi/P), order-1))

	saw = 2.0*phasor - 1.0
	phasor += phase_inc

	curr_dpw = [0]*order
	prev_dpw = [0]*order

	prev_dpw[order-1] = math.pow(saw, 6) - 5.0*math.pow(saw, 4) + 7*math.pow(saw, 2)

	# What should the initial value be for prev_dpw[n]?

	# Initial is delay (one sample per order)

	for i in xrange(0, table_size):
		saw = 2.0*phasor - 1.0

		curr_dpw[order-1] = math.pow(saw, 6) - 5.0*math.pow(saw, 4) + 7*math.pow(saw, 2)
		o = order - 2
		while o >= 0:
			curr_dpw[o] = curr_dpw[o+1] - prev_dpw[o+1]
			o -= 1

		output = curr_dpw[0]*scale

		if i >= order:
			print output
		else:
			print 0.0

		phasor += phase_inc
		
		o = 0
		while o < order:
			prev_dpw[o] = curr_dpw[o]
			o += 1


dpw_online(1024, 6)
