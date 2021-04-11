#pragma once

#include <vector>
#include "stddef.h"

namespace cie
{
namespace splinekernel
{

double evaluateBSplineBasis( double t,
                             size_t i,
                             size_t p,
                             const std::vector<double>& knotVector );

// evaluates the one-dimensional B-Spline basis functions or their first derivative.
double evaluateBSplineDerivative( double t,
                                  size_t functionIndex,
                                  size_t degree,
                                  const std::vector<double>& knotVector,
                                  size_t diffOrder );

} // splinekernel
} // cie