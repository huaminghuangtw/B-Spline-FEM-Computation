import pysplinekernel
import numpy
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.cm as colormap
import matplotlib.pyplot as plt

numberOfSamplesInR = 31;
numberOfSamplesInS = 31;

plotControlPoints = False
plotProjectionOfControlPointsOntoXYPlane = False
plotEvaluationGridLines = True

knotVectorR = [0.0, 0.0, 0.0, 0.5, 1.0, 1.0, 1.0]
knotVectorS = [0.0, 0.0, 0.0, 0.5, 1.0, 1.0, 1.0]

# Triplets of x, y, z
controlPointGrid = [ [ ( 0.0, 0.0, 1.2 ), (-0.3, 1.0, 0.9 ), (-0.4, 2.0, 0.9 ), (-0.2, 3.1, 1.3 ) ],
                     [ ( 1.0,-0.3, 0.6 ), ( 0.7, 0.7, 0.7 ), ( 0.5, 1.6, 0.8 ), ( 0.7, 2.7, 0.9 ) ],
                     [ ( 2.0,-0.2, 1.3 ), ( 1.6, 0.8, 0.7 ), ( 1.4, 1.7, 0.6 ), ( 1.6, 2.8, 1.4 ) ],
                     [ ( 3.0, 0.1, 0.5 ), ( 2.5, 1.1, 0.5 ), ( 2.3, 2.2, 0.7 ), ( 2.4, 3.2, 0.5 ) ] ];

# Transpose from matrix of coordinates vectors to a vector of coordinate matrices (both being a 3D array)
controlPointGrid = numpy.transpose( numpy.array( controlPointGrid ), ( 2, 0, 1 ) )

numberOfSamples = ( numberOfSamplesInR, numberOfSamplesInS )
knotVectors = ( knotVectorR, knotVectorS )

xyz = pysplinekernel.evaluateSurface( knotVectors, controlPointGrid, numberOfSamples )

# Plot surface (https://matplotlib.org/mpl_toolkits/mplot3d/tutorial.html)
fig = plt.figure( )
ax = fig.gca( projection='3d' )

# Compute gradient of z coordinate to add some nice colors to the plot
gradientX, gradientY = numpy.gradient( xyz[2] )
gradientMagnitude = numpy.sqrt( gradientX**2 + gradientY**2 ) 
colorMap = colormap.jet( gradientMagnitude / numpy.max( gradientMagnitude ) )

stepsize = 1

ax.plot_surface( xyz[0], xyz[1], xyz[2], facecolors=colorMap, rstride=stepsize, cstride=stepsize )

if plotEvaluationGridLines:
  ax.plot_wireframe( xyz[0], xyz[1], xyz[2], color="black", linewidth=0.3 )

if plotControlPoints:
  ax.plot_wireframe( controlPointGrid[0], controlPointGrid[1], controlPointGrid[2], color="black" )

if plotProjectionOfControlPointsOntoXYPlane:
  ax.plot_wireframe( controlPointGrid[0], controlPointGrid[1], 0.0 * controlPointGrid[2], color="gray" )

plt.show( )