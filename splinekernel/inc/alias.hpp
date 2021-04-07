#pragma once

#include <vector>
#include <tuple>
#include <array>
#include <functional>

#include "linalg.hpp"

namespace cie
{
namespace splinekernel
{

class CompressedSparseRowMatrix;

using GlobalLinearSystem = std::pair<CompressedSparseRowMatrix, std::vector<double>>;
using ElementLinearSystem = std::pair<linalg::Matrix, std::vector<double>>;

using SpatialFunction = std::function<double( double, double )>;

using IntegrationPoints = std::array<std::vector<double>, 2>;
using IntegrationPointProvider = std::function<IntegrationPoints( size_t )>;

using KnotVectors = std::array<std::vector<double>, 2>;

using LocationMap = std::vector<size_t>;
using LocationMaps = std::vector<std::vector<size_t>>;

} // namespace splinekernel
} // namespace cie
