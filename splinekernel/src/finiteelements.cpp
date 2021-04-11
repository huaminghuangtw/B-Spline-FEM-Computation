#include "finiteelements.hpp"
#include "basisfunctions.hpp"
#include "sparse.hpp"

#include <cmath>
#include <algorithm>
#include <numeric>

namespace cie
{
namespace splinekernel
{
namespace detail
{

std::array<double, 2> mapToGlobalCoordinates( std::array<double, 2> localCoordinates,
                                              std::array<size_t, 2> elementIndices,
                                              std::array<double, 2> lengths,
                                              std::array<double, 2> origin,
                                              std::array<size_t, 2> numberOfElements )
{
    std::array<double, 2> globalCoordinates;

    for (size_t axis = 0; axis < globalCoordinates.size(); ++axis)
    {
        runtime_check(numberOfElements[axis] > 0, "number of elements cannot be zero!");
        
        double elementLength = lengths[axis] / numberOfElements[axis];

        globalCoordinates[axis] = ( ((localCoordinates[axis] + 1) / 2 + elementIndices[axis])
                                    * elementLength ) + origin[axis];
    }

    return globalCoordinates;
}

KnotVectors constructOpenKnotVectors( std::array<size_t, 2> numberOfElements,
                                      std::array<size_t, 2> polynomialDegrees,
                                      std::array<size_t, 2> continuities,
                                      std::array<double, 2> lengths,
                                      std::array<double, 2> origin )
{
    KnotVectors openKnotVector;

    for (size_t axis = 0; axis < openKnotVector.size(); ++axis)
    {
        // !!!
        runtime_check( continuities[axis] < polynomialDegrees[axis],
                       "Invalid continuity for given polynomial degree!" );
        
        // outer knots (k = p + 1)
        for (size_t iKnot = 0; iKnot <= polynomialDegrees[axis]; ++iKnot)
        {
            openKnotVector[axis].push_back(origin[axis]);

        } // for iKnot

        // inner knots (k = p - c)
        for (size_t iElement = 1; iElement < numberOfElements[axis]; ++iElement)   // No.2 ~ No.(N-1) element
        {
            double elementWidth = lengths[axis] / numberOfElements[axis];
            double internalKnot = iElement * elementWidth + origin[axis];

            size_t numberOfRepetitions = polynomialDegrees[axis] - continuities[axis];

            for (size_t iKnot = 0; iKnot < numberOfRepetitions; ++iKnot)
            {
                openKnotVector[axis].push_back(internalKnot);

            } // for iKnot
        } // for iElement

        // outer knots (k = p + 1)
        for (size_t iKnot = 0; iKnot <= polynomialDegrees[axis]; ++iKnot)
        {
            openKnotVector[axis].push_back(origin[axis] + lengths[axis]);

        } // for iKnot
    } // for axis

    return openKnotVector;
}

size_t findKnotSpan( double min, double max, size_t n, double x )
{
    runtime_check(n > 0, "Number of elements is zero.");

    size_t index;

    if (x == min)
    {
        index = 0;
    }
    else if (x == max)
    {
        index = n - 1;
    }
    else
    {
        double spacing = (max - min) / n;
        index = std::floor( (x - min) / spacing );
    }

    return index;
}

LocationMaps constructLocationMaps( std::array<size_t, 2> numberOfElements,
                                    std::array<size_t, 2> polynomialDegrees,
                                    std::array<size_t, 2> continuities )
{
    runtime_check( numberOfElements[0] * numberOfElements[1] != 0,
                   "Zero number of elements!" );
    runtime_check( polynomialDegrees[0] * polynomialDegrees[1] != 0,
                   "Invalid polynomial degrees!" );
    runtime_check( (continuities[0] < polynomialDegrees[0]) && (continuities[1] < polynomialDegrees[1]),
                   "Invalid continuities!" );

    LocationMaps locationMaps;

    for (size_t iElement = 0; iElement < numberOfElements[0]; ++iElement)
    {
        for (size_t jElement = 0; jElement < numberOfElements[1]; ++jElement)
        {
            LocationMap locationMap;

            for (size_t localI = 0; localI < polynomialDegrees[0] + 1; ++localI)
            {
                for (size_t localJ = 0; localJ < polynomialDegrees[1] + 1; ++localJ)
                {                    
                    size_t strideX = polynomialDegrees[0] - continuities[0];  // i.e., #inner knots = p - c
                    size_t strideY = polynomialDegrees[1] - continuities[1];

                    size_t overlappingY = continuities[1] + 1;

                    size_t numberOfGlobalY = (numberOfElements[1] * strideY) + overlappingY;
                    
                    size_t globalI = localI + ( iElement * strideX );
                    size_t globalJ = localJ + ( jElement * strideY );

                    locationMap.push_back( globalI * numberOfGlobalY + globalJ );
                }
            }

            locationMaps.push_back( locationMap );
        }
    }

    return locationMaps;
}

} // namespace detail

BSplineFiniteElementPatch::BSplineFiniteElementPatch( std::array<size_t, 2> numberOfElements,
                                                      std::array<size_t, 2> polynomialDegrees,
                                                      std::array<size_t, 2> continuities,
                                                      std::array<double, 2> lengths,
                                                      std::array<double, 2> origin,
                                                      IntegrationPointProvider integrationPointProvider ):
    // constructor: member initialization list
    numberOfElements_( numberOfElements ),
    polynomialDegrees_( polynomialDegrees ),
    continuities_( continuities ),
    lengths_( lengths ),
    origin_( origin ),
    integrationPointProvider_( integrationPointProvider )
{ 
    knotVectors_ = detail::constructOpenKnotVectors(numberOfElements, polynomialDegrees, continuities, lengths, origin);
    locationMaps_ = detail::constructLocationMaps(numberOfElements, polynomialDegrees, continuities);
} 

// When we call this function internally we often know the knot span/element indices in which
// globalCoordinates lie already. But this being a public function it would make little sense to
// ask a user to compute the knot span indices which he shouldn't need to care about. 
std::vector<double> BSplineFiniteElementPatch::evaluateActiveBasisAt( std::array<double, 2> globalCoordinates,
                                                                      std::array<size_t, 2> diffOrders ) const
{
    size_t iElement = detail::findKnotSpan(origin_[0], origin_[0] + lengths_[0], numberOfElements_[0], globalCoordinates[0]);
    size_t jElement = detail::findKnotSpan(origin_[1], origin_[1] + lengths_[1], numberOfElements_[1], globalCoordinates[1]);

    size_t index0 = iElement * (polynomialDegrees_[0] - continuities_[0]);   // first basis function with support on iElement (See ppt "Knots&KnotSpans&BasisFunctionSupport")
    size_t index1 = jElement * (polynomialDegrees_[1] - continuities_[1]);   // first basis function with support on jElement
    
    size_t px = polynomialDegrees_[0];
    size_t py = polynomialDegrees_[1];

    std::vector<double> activeBasis;
    activeBasis.reserve( (px + 1) * (py + 1) );   // optional, to increase efficiency
    // std::vector<double> activeBasis( (px + 1) * (py + 1), 0.0 );

    for (size_t i = 0; i < px + 1; ++i)
    {
        for (size_t j = 0; j < py + 1; ++j)
        {
            double Bx = evaluateBSplineDerivative(globalCoordinates[0], index0 + i, px, knotVectors_[0], diffOrders[0]);
            double By = evaluateBSplineDerivative(globalCoordinates[1], index1 + j, py, knotVectors_[1], diffOrders[1]);
            activeBasis.push_back( Bx * By );
            // activeBasis[i * (py + 1) + j] = Bx * By;
        }
    }

    return activeBasis;
}

ElementLinearSystem BSplineFiniteElementPatch::integrateElementSystem( std::array<size_t, 2> elementIndices,
                                                                       const SpatialFunction& sourceFunction ) const
{
    // Note: we are still on the local element-level, not global system-level!
    // Compute number of element dofs!
    size_t numberOfElementDofs = (polynomialDegrees_[0] + 1) * (polynomialDegrees_[1] + 1);

    linalg::Matrix Ke( numberOfElementDofs, numberOfElementDofs, 0.0 );
    std::vector<double> Fe( numberOfElementDofs, 0.0 );

    std::vector<size_t> numberOfIntegrationPoints = { polynomialDegrees_[0] + 1, polynomialDegrees_[1] + 1 };
    
    IntegrationPoints integrationPointsX = integrationPointProvider_(numberOfIntegrationPoints[0]);
    IntegrationPoints integrationPointsY = integrationPointProvider_(numberOfIntegrationPoints[1]);

    // Loop over integration points
    for (size_t iPoint = 0; iPoint < numberOfIntegrationPoints[0]; ++iPoint)
    {
        for (size_t jPoint = 0; jPoint < numberOfIntegrationPoints[1]; ++jPoint)
        {
            std::array<double, 2> localCoordinates = { integrationPointsX[0][iPoint],
                                                       integrationPointsY[0][jPoint] };

            std::array<double, 2> weights = { integrationPointsX[1][iPoint],
                                              integrationPointsY[1][jPoint] };
            double weight = weights[0] * weights[1];

            std::array<double, 2> globalCoordinates = detail::mapToGlobalCoordinates( localCoordinates,
                                                                                      elementIndices,
                                                                                      lengths_,
                                                                                      origin_,
                                                                                      numberOfElements_ );

            auto N = evaluateActiveBasisAt( globalCoordinates, { 0,0 } );
            auto dNdx = evaluateActiveBasisAt( globalCoordinates, { 1,0 } );
            auto dNdy = evaluateActiveBasisAt( globalCoordinates, { 0,1 } );

            double detJ = (1.0 / 4) * (lengths_[0] / numberOfElements_[0]) * (lengths_[1] / numberOfElements_[1]);
            // or double detJ = ((float)1 / 4) * (lengths_[0] / numberOfElements_[0]) * (lengths_[1] / numberOfElements_[1]);
            // or double detJ = (lengths_[0] / numberOfElements_[0]) * (lengths_[1] / numberOfElements_[1]) / 4.0;
            // Do not write (1 / 4) * (lengths_[0] / numberOfElements_[0]) * (lengths_[1] / numberOfElements_[1]);
            // ---> The first term will be truncated to zero!

            double f = sourceFunction( globalCoordinates[0], globalCoordinates[1] );

            // Loop over test function and trial function, i.e., Loop over element dofs
            for (size_t i = 0; i < numberOfElementDofs; ++i)
            {
                for (size_t j = 0; j < numberOfElementDofs; ++j)
                {
                    Ke(i, j) += (dNdx[i] * dNdx[j] + dNdy[i] * dNdy[j]) * weight * detJ;
                }
                Fe[i] += N[i] * f * weight * detJ;
            }
        }
    }

    // ElementLinearSystem elementSystem = std::make_pair(Ke, Fe);
    // return elementSystem; 
    return { Ke, Fe };
}

GlobalLinearSystem BSplineFiniteElementPatch::assembleGlobalSystem(const SpatialFunction& sourceFunction) const
{
    CompressedSparseRowMatrix globalMatrix( locationMaps_ );

    std::vector<double> globalVector( globalMatrix.size(), 0.0 );

    for (size_t iElement = 0; iElement < numberOfElements_[0]; ++iElement)
    {
        for (size_t jElement = 0; jElement < numberOfElements_[1]; ++jElement)
        {
            auto elementSystem = integrateElementSystem( { iElement, jElement }, sourceFunction );

            const LocationMap& locationMap = locationMaps_[iElement * numberOfElements_[1] + jElement];

            globalMatrix.scatter( elementSystem.first, locationMap );
            // globalMatrix.scatter( std::get<0>(elementSystem), locationMap );

            for (size_t iDof = 0; iDof < locationMap.size(); ++iDof)
            {
                globalVector[locationMap[iDof]] += elementSystem.second[iDof];
                // globalVector[locationMap[i]] += std::get<1>(elementSystem)[i];
            }
        }
    }

    return { globalMatrix, globalVector };
}

std::vector<size_t> BSplineFiniteElementPatch::boundaryDofIds(const std::string& side) const
{
    std::vector<size_t> boundaryDofIds;
    
    // size_t numberOfDofsX = numberOfElements_[0] * (polynomialDegrees_[0] - continuities_[0]) + (continuities_[0] + 1);
    size_t numberOfDofsX = knotVectors_[0].size() - polynomialDegrees_[0] - 1;
    // size_t numberOfDofsY = numberOfElements_[1] * (polynomialDegrees_[1] - continuities_[1]) + (continuities_[1] + 1);
    size_t numberOfDofsY = knotVectors_[1].size() - polynomialDegrees_[1] - 1;

    if ( side == "top" )
    {
        for (size_t globalI = 0; globalI < numberOfDofsX; ++globalI)
        {
            boundaryDofIds.push_back( globalI * numberOfDofsY + (numberOfDofsY - 1) );
        }
    }
    else if ( side == "right" )
    {
        for (size_t globalJ = 0; globalJ < numberOfDofsY; ++globalJ)
        {
            boundaryDofIds.push_back( (numberOfDofsX - 1) * numberOfDofsY + globalJ );
        }
    }
    else if ( side == "bottom" )
    {
        for (size_t globalI = 0; globalI < numberOfDofsX; ++globalI)
        {
            boundaryDofIds.push_back( globalI * numberOfDofsY );
        }
    }
    else if ( side == "left" )
    {
        for (size_t globalJ = 0; globalJ < numberOfDofsY; ++globalJ)
        {
            boundaryDofIds.push_back( globalJ );
        }
    }
    else
    {
        throw std::runtime_error("Invalid side string.");
    }


    return boundaryDofIds;
}

SpatialFunction BSplineFiniteElementPatch::solutionEvaluator(const std::vector<double>& solutionDofs) const
{
    return [=]( double x, double y ) -> double
    {
        std::vector<double> N = evaluateActiveBasisAt( { x, y }, { 0, 0 } );

        size_t elementIdxI = detail::findKnotSpan(origin_[0], origin_[0] + lengths_[0], numberOfElements_[0], x);
        size_t elementIdxJ = detail::findKnotSpan(origin_[1], origin_[1] + lengths_[1], numberOfElements_[1], y);
        
        const LocationMap& locationMap = locationMaps_[elementIdxI * numberOfElements_[1] + elementIdxJ];
        
        double value = 0.0;

        for (size_t iDof = 0; iDof < locationMap.size(); ++iDof)
        {
            // basis.size(), i.e., (polynomialDegrees_[0] + 1) * (polynomialDegrees_[1] + 1)
            // is equal to
            // locationMap.size()
            value += N[iDof] * solutionDofs[locationMap[iDof]];
        }

        return value;
    };
}

} // namespace splinekernel
} // namespace cie
