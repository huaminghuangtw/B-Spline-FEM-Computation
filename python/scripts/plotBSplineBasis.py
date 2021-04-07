import pysplinekernel
import numpy
import matplotlib.pyplot as plt


n = 401  # number of sample points
p = 2    # polynomial degree
c = 1    # continuity
s = 5    # number of (non-zero) knot spans


# Construct knot vector
U = [0.0] * (p + 1) 

for i in range( s - 1 ):
	U = U + [(i + 1.0) / s] * (p - c)

U = U + [1.0] * (p + 1)

# Sample coordinates
x = numpy.linspace( 0.0, 1.0, n )

m = len( U )

# Plot knot spans  
for i in range( s + 1 ):
  t = float(i) / s
  plt.plot( [t, t], [0.0, 1.0], 'gray' )
  
# Loop over all basis functions
for i in range( m - p - 1 ):
  y = numpy.zeros( n )

  # Loop over all evaluation points
  for j, xj in enumerate( x ):
    y[j] = pysplinekernel.evaluateBSplineBasis( xj, i, p, U )
    
  plt.plot( x, y )

plt.show( )
