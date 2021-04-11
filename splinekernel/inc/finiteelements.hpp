#pragma once

#include <vector>
#include <tuple>
#include <array>
#include <functional>

#include "linalg.hpp"
#include "alias.hpp"
#include "utilities.hpp"

namespace cie
{
namespace splinekernel
{

/*! Helper class to construct the linear equation system for a finite element method    *
 *  on a B-Spline patch. In this case the mesh of (non-zero) knot spans represents the  *
 *  finite element mesh. For simplicity, the geometry of the computational domain is    *
 *  restricted to an axis-aligned rectangle.                                            */
class BSplineFiniteElementPatch
{
public:
    BSplineFiniteElementPatch( std::array<size_t, 2> numberOfElements,
                               std::array<size_t, 2> polynomialDegrees,
                               std::array<size_t, 2> continuities,
                               std::array<double, 2> lengths,
                               std::array<double, 2> origin,
                               IntegrationPointProvider integrationPointProvider );

    std::vector<double> evaluateActiveBasisAt( std::array<double, 2> globalCoordinates,
                                               std::array<size_t, 2> diffOrders ) const;

    ElementLinearSystem integrateElementSystem( std::array<size_t, 2> elementIndices,
                                                const SpatialFunction& sourceFunction ) const;
    
    GlobalLinearSystem assembleGlobalSystem( const SpatialFunction& sourceFunction ) const;

    std::vector<size_t> boundaryDofIds( const std::string& side ) const;

    SpatialFunction solutionEvaluator( const std::vector<double>& solutionDofs ) const;
    
private:
    std::array<size_t, 2> numberOfElements_, polynomialDegrees_, continuities_;
    std::array<double, 2> lengths_, origin_;

    IntegrationPointProvider integrationPointProvider_;
    
    KnotVectors knotVectors_;
    LocationMaps locationMaps_;
};

namespace detail
{

//! Map from local coordinates [-1, 1] on a given knot span to (x, y)
std::array<double, 2> mapToGlobalCoordinates( std::array<double, 2> localCoordinates,
                                              std::array<size_t, 2> elementIndices,
                                              std::array<double, 2> lengths,
                                              std::array<double, 2> origin,
                                              std::array<size_t, 2> numberOfElements);
    
//! Construct knot vectors in r and s for given data
KnotVectors constructOpenKnotVectors( std::array<size_t, 2> numberOfElements,
                                      std::array<size_t, 2> polynomialDegrees,
                                      std::array<size_t, 2> continuities,
                                      std::array<double, 2> lengths,
                                      std::array<double, 2> origin );

//! Find the indices of the knot span containing the coordinate x, assuming uniform spacing
size_t findKnotSpan( double min, double max, size_t n, double x );

//! Construct the location maps for all elements
LocationMaps constructLocationMaps( std::array<size_t, 2> numberOfElements,
                                    std::array<size_t, 2> polynomialDegrees,
                                    std::array<size_t, 2> continuities );

} // namespace detail
} // namespace splinekernel
} // namespace cie
