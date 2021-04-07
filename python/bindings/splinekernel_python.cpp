#include "pybind11/pybind11.h"
#include "pybind11/stl.h"

#include "basisfunctions.hpp"
#include "curve.hpp"
#include "surface.hpp"
#include "finiteelements.hpp"

#include "denseMatrixConversion.hpp"
#include "sparseMatrixConversion.hpp"

PYBIND11_MODULE( pysplinekernel, m ) 
{
    m.doc( ) = "spline computation kernel"; // optional module doc string

    m.def( "evaluateBSplineBasis", &cie::splinekernel::evaluateBSplineBasis, "Evaluates single b-spline basis function." );
    m.def( "evaluate2DCurve", &cie::splinekernel::evaluate2DCurve, "Evaluate B-Spline curve by summing up basis functions times control points." );
    m.def( "evaluateSurface", &cie::splinekernel::evaluateSurface, "Evaluate B-Spline surface." );

	// Export the class BSplineFiniteElementPatch and corresponding public member functions in pybind11 using the class_ class template.
	pybind11::class_<cie::splinekernel::BSplineFiniteElementPatch> patch( m, "BSplineFiniteElementPatch" );

	patch.def( pybind11::init< std::array<size_t, 2>,
							   std::array<size_t, 2>,
		                       std::array<size_t, 2>,
						       std::array<double, 2>,
						       std::array<double, 2>,
						       cie::splinekernel::IntegrationPointProvider >( ) ); // constructor

	// member functions
	patch.def( "evaluateActiveBasisAt", &cie::splinekernel::BSplineFiniteElementPatch::evaluateActiveBasisAt ); 
	patch.def( "integrateElementSystem", &cie::splinekernel::BSplineFiniteElementPatch::integrateElementSystem );
	patch.def( "assembleGlobalSystem", &cie::splinekernel::BSplineFiniteElementPatch::assembleGlobalSystem );
	patch.def( "boundaryDofIds", &cie::splinekernel::BSplineFiniteElementPatch::boundaryDofIds );
	patch.def( "solutionEvaluator", &cie::splinekernel::BSplineFiniteElementPatch::solutionEvaluator );
}
