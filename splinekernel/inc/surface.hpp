#pragma once

#include "linalg.hpp"

namespace cie
{
namespace splinekernel
{

using VectorOfMatrices = std::vector<linalg::Matrix>;

/*
* Evaluates a 2D B-Spline patch.
* @param knotVectors Two knot vectors in r and s directions
* @param controlPoints A vector of matrices for each control point component. The control points of the
*					   patch form a 2D grid. They may have one single value when a scalar field is rep-
*					   resented, or they may have multiple components when, e.g. x, y and z coordinates
*					   are interpolated. There is one matrix for each one component (e.g. x-value).
* @param numberOfSamples The number of times the patch shall be evaluated in each coordinate direction.
*						 For example passing (50, 50) evaluates the entire patch on a 50 x 50 grid.
* @return A vector of matrices, similar to the controlPoints but now with the dimensions specified
*		  by numberOfSamples. In the above example we would return a 50 x 50 matrix of x-values,
*		  a 50 x 50 matrix of y-values, etc.
*/

VectorOfMatrices evaluateSurface( const std::array<std::vector<double>, 2>& knotVectors,
								  const VectorOfMatrices& controlPoints,
								  std::array<size_t, 2> numberOfSamplePoints );

} // namespace splinekernel
} // namespace cie