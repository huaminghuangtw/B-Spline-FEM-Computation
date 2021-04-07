import numpy
import sys
import scipy.sparse
import scipy.sparse.linalg

from mpl_toolkits.mplot3d import Axes3D
import matplotlib.cm as colormap
import matplotlib.pyplot as plt

def createSparseMatrix( indices, indptr, data  ):
    numberOfDofs = len( indptr ) - 1;

    K = scipy.sparse.csr_matrix( ( numberOfDofs, numberOfDofs ) ) 
    K.indices, K.indptr, K.data = indices, indptr, data
    K.has_sorted_indices = True

    memory = K.nnz * ( indices.itemsize + data.itemsize ) + ( numberOfDofs + 1 ) * indptr.itemsize

    print( "    number of dofs: " + str( numberOfDofs ) )
    print( "    number of non-zeros: " + str( K.nnz ) )
    print( "    fill ratio: {:.3f}%".format( float( K.nnz ) / ( K.shape[0] * K.shape[1] ) * 100 ) )
    print( "    matrix requires {:.3f} MB memory.".format( memory / 1e6 ) )

    return K

def solveEquationSystem( K, F ):
    def countIterations( _ ):
        countIterations.count += 1
        
        sys.stdout.write( '\r    number of iterations: ' + str( countIterations.count ) )
        sys.stdout.flush( )
        
    countIterations.count = 0

    # diagonal / Jacobi preconditioning
    invDiag = 1.0 / K.diagonal( ) 
    preconditioner = scipy.sparse.linalg.LinearOperator( ( K.shape[0], K.shape[1] ), lambda x : x * invDiag ) 
    solutionDofs, exitCode = scipy.sparse.linalg.bicgstab( K, F, M=preconditioner, tol=1e-12, callback=countIterations )
    
    print( "\n    || K u - f || / || f || = {:.4e}".format( numpy.linalg.norm( K * solutionDofs - F ) / numpy.linalg.norm( F ) ) )

    return solutionDofs

def imposeDirichletBoundaryCondition( K, F, indices, values ):
    K[numpy.array(indices), numpy.array(indices)] = 1e8;
    F[numpy.array(indices)] = numpy.array(values) * 1e8
    
def plotScalarFunction( u, lengths, origin, numberOfElements, numberOfCellsPerElement ):

    numberOfSamples = ( numberOfElements[0] * numberOfCellsPerElement[0] + 1,
                        numberOfElements[1] * numberOfCellsPerElement[1] + 1 )

    solutionEvaluator = numpy.vectorize( u )

    X, Y = numpy.meshgrid( numpy.linspace( origin[0], origin[0] + lengths[0], numberOfSamples[0] ), 
                           numpy.linspace( origin[1], origin[1] + lengths[1], numberOfSamples[1] ), indexing='ij' )

    Z = solutionEvaluator( X, Y )

    colorMap = colormap.jet( Z / numpy.max( Z ) )

    ax = plt.figure( ).gca( projection='3d' )
    ax.plot_surface( X, Y, Z, facecolors=colorMap, rstride=1, cstride=1 )
    ax.plot_wireframe( X, Y, Z, color="black", rstride=numberOfCellsPerElement[0], cstride=numberOfCellsPerElement[1] )

    plt.show( )
