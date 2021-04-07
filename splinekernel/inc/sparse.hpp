#pragma once

#include "linalg.hpp"
#include "alias.hpp"

#include <vector>
#include <tuple>

namespace cie
{
namespace splinekernel
{

class CompressedSparseRowMatrix
{
public:
    using IndexType = std::int32_t;

    explicit CompressedSparseRowMatrix( const std::vector<LocationMap>& locationMaps );
    // explicit keyword is used to prevent implicit conversions in one-argument constructor
    
    size_t size( ) const;
    IndexType nnz( ) const;

    double operator()( size_t i, size_t j ) const;
    
    std::vector<double> operator*( const std::vector<double>& vector );
    
    void scatter( const linalg::Matrix& elementMatrix, const LocationMap& locationMap );
    
    std::tuple<IndexType*, IndexType*, double*> dataStructure( );
    
private:
    std::vector<IndexType> indices_, indptr_;
    std::vector<double> data_;
};

} // namespace splinekernel
} // namespace cie