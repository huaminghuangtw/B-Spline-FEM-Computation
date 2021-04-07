#pragma once

#include <vector>
#include <array>

namespace cie
{
namespace splinekernel
{

/* The comment below is specifically formatted to work with doxygen (http://www.doxygen.nl/), *
 * which is a tool to automatically generate a nice documentation for your source code. We    *
 * won't use it in this course, but it's nice to know about it.                               */

/*! Evaluate B-Spline curve by summing up basis functions times control points.
 *  @param tCoordinates The parametric coordinates at which the curve shall be evaluated
 *  @param xCoordinates The x coordinates of the control points
 *  @param yCoordinates The y coordinates of the control points
 *  @return a vector of x and a vector of y coordinates with one value for each parametric
 *          coordinate tCoordinates
 */
std::array<std::vector<double>, 2> evaluate2DCurve( const std::vector<double>& tCoordinates,
                                                    const std::vector<double>& xCoordinates,
                                                    const std::vector<double>& yCoordinates,
                                                    const std::vector<double>& knotVector );

} // namespace splinekernel
} // namespace cie
