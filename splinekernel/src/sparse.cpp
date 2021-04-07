#include "sparse.hpp"
#include "utilities.hpp"

#include <algorithm>
#include <numeric>

namespace cie
{
namespace splinekernel
{

/* An entry (i, j) in the sparse matrix coming from the finite element method is non-zero if *
 * the two corresponding shape functions Ni and Nj overlap (i and j being global indices).   *
 * The crucial information we need are the location maps, which tell us what shape functions *
 * are supported on an element (supported = non zero). The algorithm to construct a          *
 * compressed sparse row matrix is therefore roughly as follows:                             *
 * 1. For each global shape function index find the elements where it is supported. For      *
 *    example, the shape function with global index 0 might be nonzero in elements 0, 1, 8   *
 *    and 9, the shape function with global index 1 might be nonzero in elements 0, 1, 2, 8, *
 *    9 and 10, and so on.                                                                   *
 * 2. For each global shape function index concatenate/assemble the location maps of all     *
 *    elements that it is supported on. This will give a list of all shape functions that    *
 *    overlap with the current one. As this list has duplicates we need to make it unique.   *
 *    The resulting indices are the non-zero indices of one row of the sparse matrix         *
 * 3. Accumulate the sizes of the vectors from step 2 into the indptr_ vector. If the sizes  *
 *    are for example [2, 2, 2], then the accumulation yields [2, 4, 6]. On top we add a     *
 *    zero at the first index: [0, 2, 4, 6]                                                  *
 * 4. Linearize the vector of vectors coming from step 2 into a single big indices_ vector.  *
 * 5. Resize the data_ vector to the same size as indices_ and initialize the values with 0  */
CompressedSparseRowMatrix::CompressedSparseRowMatrix( const std::vector<LocationMap>& locationMaps )
{
    size_t numberOfElements = locationMaps.size( );
    size_t size = 0;

    // Compute number of dofs
    for( const auto& locationMap : locationMaps )
    {
        size_t maxElementIndex = *std::max_element( locationMap.begin( ), locationMap.end( ) ) + 1;

        if( maxElementIndex > size )
        {
            size = maxElementIndex;
        }
    }
    
    std::vector<std::vector<size_t>> dofToElementCoupling( size );

    // For each dof find the elements that connect to it
    for( size_t iElement = 0; iElement < numberOfElements; ++iElement )
    {
        for( size_t dofIndex : locationMaps[iElement] )
        {
            dofToElementCoupling[dofIndex].push_back( iElement );
        }
    }

    std::vector<std::vector<size_t>> dofToDofCoupling( size );
    std::vector<size_t> tmp;

    // For each dof concatenate the location maps of all connected 
    // elements (from the previous step) and make them unique.
    for( size_t iDof = 0; iDof < size; ++iDof )
    {
        runtime_check( dofToElementCoupling[iDof].size( ) > 0, "Found dof without connected element!" );

        tmp.resize( 0 );

        for( size_t iElement = 0; iElement < dofToElementCoupling[iDof].size( ); ++iElement )
        {
            const auto& locationMap = locationMaps[dofToElementCoupling[iDof][iElement]];

            tmp.insert( tmp.end( ), locationMap.begin( ), locationMap.end( ) );
        }

        std::sort( tmp.begin( ), tmp.end( ) );
        auto end = std::unique( tmp.begin( ), tmp.end( ) );

        dofToDofCoupling[iDof].insert( dofToDofCoupling[iDof].end( ), tmp.begin( ), end );
    }

    // Free memory that we don't need
    dofToElementCoupling.clear( );
    dofToElementCoupling.shrink_to_fit( );

    indptr_.resize( size + 1 );
    indptr_[0] = 0;

    auto vectorSizePredicate = []( const std::vector<size_t>& dofs ) { return static_cast<IndexType>( dofs.size( ) ); };

    // Fill indptr with cumulative sum of the vector sizes of dofToDofCoupling
    std::transform( dofToDofCoupling.begin( ), dofToDofCoupling.end( ), indptr_.begin( ) + 1, vectorSizePredicate );
    std::partial_sum( indptr_.begin( ), indptr_.begin( ) + size + 1, indptr_.begin( ) );

    IndexType nnz = indptr_[size];

    indices_.resize( nnz );

    // Concatenate dof indices
    auto currentRow = indices_.begin( );

    for( const std::vector<size_t>& dofs : dofToDofCoupling )
    {
        std::copy( dofs.begin( ), dofs.end( ), currentRow );

        currentRow += dofs.size( );
    }

    // Free memory that we don't need
    dofToDofCoupling.clear( );
    dofToDofCoupling.shrink_to_fit( );

    // Resize data
    data_.resize( nnz );

    std::fill( data_.begin( ), data_.end( ), 0.0 );
}

std::tuple<CompressedSparseRowMatrix::IndexType*, CompressedSparseRowMatrix::IndexType*, double*> CompressedSparseRowMatrix::dataStructure( )
{
    return { indices_.data( ), indptr_.data( ), data_.data( ) };
}

size_t CompressedSparseRowMatrix::size( ) const
{
    return indptr_.size() - 1;
}

CompressedSparseRowMatrix::IndexType CompressedSparseRowMatrix::nnz( ) const
{
    return indptr_.back();
}

double CompressedSparseRowMatrix::operator( )( size_t i, size_t j ) const
{
    runtime_check( i < this->size( ), "Index out of range." );
    runtime_check( j < this->size( ), "Index out of range." );

    auto begin = indices_.begin( ) + indptr_[i];
    auto end = indices_.begin( ) + indptr_[i + 1];
    
    auto result = std::find( begin, end, j );
    
    if ( result != end )
    {
        return data_[ result - indices_.begin() ];  // non-zero value
        // std::distance( indices_.begin(), result )
        // access the element in data_ with the same index as "result" in indices_
    }
    else
    {
		return 0.0;  // zero value
    }
}

std::vector<double> CompressedSparseRowMatrix::operator*( const std::vector<double>& vector )
{
	runtime_check( vector.size() == this->size(),
		           "Invalid RHS size." );

	size_t numberOfRows = this->size();

	std::vector<double> result(numberOfRows, 0.0);

	for (size_t i = 0; i < numberOfRows; ++i)
	{
		for (size_t j = indptr_[i]; j < indptr_[i + 1]; ++j)
		{
			result[i] += data_[j] * vector[indices_[j]];
		}
	}

	return result;
}

void CompressedSparseRowMatrix::scatter( const linalg::Matrix& elementMatrix, const LocationMap& locationMap )
{
    // Assemble local element matrix to global master matrix!
    runtime_check( elementMatrix.size1() == elementMatrix.size2(), "Non-square matrix." );
    runtime_check( locationMap.size() == elementMatrix.size1(), "Inconsistent size." );

    for (size_t localI = 0; localI < elementMatrix.size1(); ++localI)
    {
        for (size_t localJ = 0; localJ < elementMatrix.size2(); ++localJ)
        {       
			size_t globalI = locationMap[localI];
			size_t globalJ = locationMap[localJ];

            auto begin = indices_.begin() + indptr_[globalI];
            auto end = indices_.begin() + indptr_[globalI + 1];

            /*auto Idx = std::find( begin, end, globalJ );

			runtime_check( Idx != end, "Entry not present in sparsity pattern." );*/

			auto Idx = std::lower_bound( begin, end, globalJ );

			runtime_check( Idx != end && *Idx == globalJ, "Entry not present in sparsity pattern.");

			data_[ Idx - indices_.begin() ] += elementMatrix(localI, localJ);
            // data_[ std::distance( indices_.begin(), Idx ) ] += elementMatrix(i, j);
        }
    }
}

} // splinekernel
} // cie