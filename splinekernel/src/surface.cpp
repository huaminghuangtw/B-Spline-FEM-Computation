#include "surface.hpp"
#include "basisfunctions.hpp"

namespace cie
{
namespace splinekernel
{
	
	VectorOfMatrices evaluateSurface( const std::array<std::vector<double>, 2>& knotVectors,
									  const VectorOfMatrices& controlPoints,
									  std::array<size_t, 2> numberOfSamplePoints )
	{
		size_t numberOfControlPointsR = controlPoints[0].size1();
		size_t numberOfControlPointsS = controlPoints[0].size2();

		size_t numberOfFields = controlPoints.size();

		size_t pr = knotVectors[0].size() - numberOfControlPointsR - 1;   // p = m - n - 1
		size_t ps = knotVectors[1].size() - numberOfControlPointsS - 1;

		VectorOfMatrices result(numberOfFields);

		for (size_t iField = 0; iField < numberOfFields; ++iField)
		{
			result[iField] = linalg::Matrix(numberOfSamplePoints[0], numberOfSamplePoints[1], 0.0);

			for (size_t iSampleCoordinate = 0; iSampleCoordinate < numberOfSamplePoints[0]; ++iSampleCoordinate)
			{
				for (size_t jSampleCoordinate = 0; jSampleCoordinate < numberOfSamplePoints[1]; ++jSampleCoordinate)
				{
					double r = iSampleCoordinate / (numberOfSamplePoints[0] - 1.0);  // Normalization
					double s = jSampleCoordinate / (numberOfSamplePoints[1] - 1.0);

					for (size_t iBasisFunction = 0; iBasisFunction < numberOfControlPointsR; ++iBasisFunction)
					{
						for (size_t jBasisFunction = 0; jBasisFunction < numberOfControlPointsS; ++jBasisFunction)
						{
							double Ni = evaluateBSplineBasis(r, iBasisFunction, pr, knotVectors[0]);
							double Nj = evaluateBSplineBasis(s, jBasisFunction, ps, knotVectors[1]);

							double value = Ni * Nj * controlPoints[iField](iBasisFunction, jBasisFunction);
							
							result[iField](iSampleCoordinate, jSampleCoordinate) += value;
						}   // jBasisFunction
					}   // iBasisFunction
				}   // jSampleCoordinate
			}   // iSampleCoordinate
		}   // iField

		return result;
	}

} // namespace splinekernel
} // namespace cie
