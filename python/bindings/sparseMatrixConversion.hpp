#pragma once

#include "sparse.hpp"

namespace pybind11
{
namespace detail
{

// Template specialization for pybind11::detail::type_caster, see
// https://pybind11.readthedocs.io/en/stable/advanced/cast/custom.html
template<>
struct type_caster<cie::splinekernel::CompressedSparseRowMatrix>
{
public:

    PYBIND11_TYPE_CASTER( cie::splinekernel::CompressedSparseRowMatrix, _( "cie::splinekernel::CompressedSparseRowMatrix" ) );

    // Convert sparse data structure to three numpy arrays. Note, that this creates
    // a copy of the sparse matrix. It is possible to use dynamic arrays and raw
    // pointers instread of std::vectors and transfer ownership to python without
    // creating a copy. We chose this version at it is easier to understand.
    static pybind11::handle cast( const cie::splinekernel::CompressedSparseRowMatrix& src,
                                  pybind11::return_value_policy policy,
                                  pybind11::handle parent )
    {
        // Don't do this at home
        auto data = const_cast<cie::splinekernel::CompressedSparseRowMatrix&>( src ).dataStructure( );

        auto size = src.size( );
        auto nnz = src.nnz( );
        
        // Construct array with shape ( size1, size2 ) from the given linalg::Matrix
        pybind11::list arrays;

        // Append sparse data structure
        arrays.append( pybind11::array( nnz, std::get<0>( data ) ) );
        arrays.append( pybind11::array( size + 1, std::get<1>( data ) ) );
        arrays.append( pybind11::array( nnz, std::get<2>( data ) ) );

        return arrays.release( );
    }
};

} // detail
} // pybind11
