#include "basisfunctions.hpp"
#include "curve.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <string>

namespace cie
{
namespace splinekernel
{

    std::array<std::vector<double>, 2> evaluate2DCurve( const std::vector<double>& tCoordinates,
                                                        const std::vector<double>& xCoordinates,
                                                        const std::vector<double>& yCoordinates,
                                                        const std::vector<double>& knotVector )
    {
        size_t numberOfSamples = tCoordinates.size();   // number of evaluation points/samples
        
        size_t m = knotVector.size();   // number of knots
        size_t n = xCoordinates.size();   // number of control points/vertices
        size_t p = m - n - 1;   // degree of the B-spline curve (m = n + p + 1)

        std::vector<double> curveX(numberOfSamples, 0.0);
        std::vector<double> curveY(numberOfSamples, 0.0);

        // evaluate all n shape functions and multiply with the corresponding n control points
        // at each parametric coordinate tCoordinates
        for (size_t i = 0; i < n; ++i)
        {
            for (size_t j = 0; j < numberOfSamples; ++j)
            {
                double basis = evaluateBSplineBasis(tCoordinates.at(j), i, p, knotVector);
                
                curveX.at(j) += basis * xCoordinates.at(i);
                curveY.at(j) += basis * yCoordinates.at(i);
            }
        }

        /* Alternative:
        for (size_t i = 0; i < numberOfSamples; ++i)
        {
            for (size_t j = 0; j < n; ++j)
            {
                double N = evaluateBSplineBasis(tCoordinates[i], j, p, knotVector);

                curveX[i] += N * xCoordinates[j];
                curveY[i] += N * yCoordinates[j];
            }
        }
        */

        return { curveX, curveY };
    }

} // namespace splinekernel
} // namespace cie
