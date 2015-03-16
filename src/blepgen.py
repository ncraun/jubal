import math

# sinc
def sinc(x):
	if x == 0.0:
		return 1.0
	else:
		return math.sin(math.pi*x)/(math.pi*x)

def blep(crossings, oversample):

	table_length = crossings*2*oversample + 1
	table = [0]*table_length

	a = -crossings	
	b = crossings

	# Create sinc
	for i in xrange(table_length):
		r = float(i)/float(table_length-1)	
		table[i] = sinc(a + r*(b-a))

	# Window sinc

	# Integrate
	integral = 0
	max = 0
	for i in xrange(table_length):
		integral += table[i]
		table[i] = integral
		if table[i] > max: max = table[i]

	# Normalize to Bipolar
	# and subtract the unit step
	norm = 1.0/table[table_length-1]
	for i in xrange(table_length):
		table[i] = 2*norm*table[i] - 1
		if i >= table_length/2:
			table[i] -= 1

	# Output
	for i in xrange(table_length):
		print table[i]

blep(1, 2048)
