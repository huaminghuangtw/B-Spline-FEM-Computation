import pysplinekernel
import numpy

import finiteelementshelper as fem

polynomialDegrees = (2, 2)
numberOfElements = (10, 10)
lengths = (1.0, 1.0)
origin = (0.0, 0.0)
continuity = tuple( p - 1 for p in polynomialDegrees )

source = lambda x, y : 0.0
#source = lambda x, y : numpy.exp( -( ( x - 0.3 )**2 + ( y - 0.5 )**2 ) * 500 ) * 500

topBoundaryFunction = lambda x : ( x < 0.3 or x > 0.7 ) / 2.0 + 0.5

mesh = pysplinekernel.BSplineFiniteElementPatch( numberOfElements, polynomialDegrees, continuity, 
                                                 lengths, origin, numpy.polynomial.legendre.leggauss )

print( "Assembling linear system ..." ) 

( ( indices, indptr, data ), F ) = mesh.assembleGlobalSystem( source )

K = fem.createSparseMatrix( indices, indptr, data )
F = numpy.array( F )

print( "Imposing boundary conditions ..." )

topBoundaryDofs = mesh.boundaryDofIds( "top" )

topBoundaryValues = [ topBoundaryFunction( x ) for x in numpy.linspace( 0.0, 1.0, len( topBoundaryDofs ) ) ]

fem.imposeDirichletBoundaryCondition( K, F, topBoundaryDofs, topBoundaryValues )
fem.imposeDirichletBoundaryCondition( K, F, mesh.boundaryDofIds( "bottom" ), 0.0 )

print( "Solving linear system ..." )

solutionDofs = fem.solveEquationSystem( K, F )

print( "Postprocessing ... " )

u = mesh.solutionEvaluator( solutionDofs )

numberOfCellsPerElement = [ 1 if p == 1 else p + 1 for p in polynomialDegrees ]

fem.plotScalarFunction( u, lengths, origin, numberOfElements, numberOfCellsPerElement )
