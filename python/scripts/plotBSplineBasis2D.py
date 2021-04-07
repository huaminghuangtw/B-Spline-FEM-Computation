import pysplinekernel
import numpy
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.cm as colormap
import matplotlib.pyplot as plt

numberOfSamplesInR = 51;
numberOfSamplesInS = 51;

knotVectorR = [0.0, 0.0, 0.0, 0.5, 1.0, 1.0, 1.0]
knotVectorS = [0.0, 0.0, 0.0, 0.5, 1.0, 1.0, 1.0]

polynomialDegreeR = 2
polynomialDegreeS = 2 

basisFunctionIndexToVisualize = ( 3, 1 )

numberOfControlPointsR = len( knotVectorR ) - polynomialDegreeR - 1
numberOfControlPointsS = len( knotVectorS ) - polynomialDegreeS - 1

controlPointGrid = numpy.zeros( ( numberOfControlPointsR, numberOfControlPointsS ) )
controlPointGrid[ basisFunctionIndexToVisualize[0], basisFunctionIndexToVisualize[1] ] = 1.0

functionValues, = pysplinekernel.evaluateSurface( ( knotVectorR, knotVectorS ), [ controlPointGrid ], ( numberOfSamplesInR, numberOfSamplesInS ) )

# Plot surface (https://matplotlib.org/mpl_toolkits/mplot3d/tutorial.html)
fig = plt.figure( )
ax = fig.gca( projection='3d' )

X, Y = numpy.meshgrid( numpy.linspace( 0.0, 1.0, numberOfSamplesInR ), 
                       numpy.linspace( 0.0, 1.0, numberOfSamplesInS ), indexing='ij' )

colorMap = colormap.jet( functionValues / numpy.max( functionValues ) )

stepsize = 1

ax.plot_surface( X, Y, functionValues, facecolors=colorMap, rstride=stepsize, cstride=stepsize )

plt.show( )
