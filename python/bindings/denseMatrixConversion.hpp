#pragma once

#include "pybind11/numpy.h"
#include "pybind11/functional.h"

namespace pybind11
{
namespace detail
{

// Template specialization for pybind11::detail::type_caster, see
// https://pybind11.readthedocs.io/en/stable/advanced/cast/custom.html
template<>
struct type_caster<cie::linalg::Matrix>
{
public:

    PYBIND11_TYPE_CASTER( cie::linalg::Matrix, _( "cie::linalg::Matrix" ) );

    // Conversion from python numpy array to our C++ linalg::Matrix
    bool load( pybind11::handle src, bool convert )
    {
        if ( convert || pybind11::array_t<double>::check_( src ) )
        {
            auto numpyArray = pybind11::array_t<double, pybind11::array::c_style |
                                                        pybind11::array::forcecast>::ensure( src );

            if ( numpyArray && numpyArray.ndim( ) == 2 )
            {
                size_t size1 = numpyArray.shape( )[0];
                size_t size2 = numpyArray.shape( )[1];

                // value is a member defined by the PYBIND11_TYPE_CASTER macro
                value = cie::linalg::Matrix( size1, size2, 0.0 );

                // Copy matrix
                for( size_t i = 0; i < size1; i++ )
                {
                    for( size_t j = 0; j < size2; j++ )
                    {
                        value( i, j ) = numpyArray.at( i, j );
                    }
                }

                return true; // success
            }
        }
        return false; // failure
    }

    // Conversion from our C++ linalg::Matrix to python numpy array
    static pybind11::handle cast( const cie::linalg::Matrix& src,
                                  pybind11::return_value_policy policy,
                                  pybind11::handle parent )
    {
        // Construct array with shape ( size1, size2 ) from the given linalg::Matrix
        return pybind11::array( std::vector<size_t>{ src.size1( ), src.size2( ) },
                                &const_cast<cie::linalg::Matrix&>( src )( 0, 0 ) ).release( );
    }
};

} // detail
} // pybind11
