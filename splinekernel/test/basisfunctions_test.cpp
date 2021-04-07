#include "catch.hpp"
#include "basisfunctions.hpp" 
#include <vector>

namespace cie
{
namespace splinekernel
{

TEST_CASE("Linear interpolation")
{
  std::vector<double> knotVector{ 0.0, 0.0, 0.5, 1.0, 1.0 };

  const size_t p = 1;

  // First basis function
  REQUIRE_NOTHROW(evaluateBSplineBasis(0.0, 0, p, knotVector));

  CHECK(evaluateBSplineBasis(0.00, 0, p, knotVector) == Approx(1.0));
  CHECK(evaluateBSplineBasis(0.25, 0, p, knotVector) == Approx(0.5));
  CHECK(evaluateBSplineBasis(0.50, 0, p, knotVector) == Approx(0.0));
  CHECK(evaluateBSplineBasis(0.75, 0, p, knotVector) == Approx(0.0));
  CHECK(evaluateBSplineBasis(1.00, 0, p, knotVector) == Approx(0.0));

  //  Second basis function
  REQUIRE_NOTHROW(evaluateBSplineBasis(0.0, 1, p, knotVector));

  CHECK(evaluateBSplineBasis(0.00, 1, p, knotVector) == Approx(0.0));
  CHECK(evaluateBSplineBasis(0.25, 1, p, knotVector) == Approx(0.5));
  CHECK(evaluateBSplineBasis(0.50, 1, p, knotVector) == Approx(1.0));
  CHECK(evaluateBSplineBasis(0.75, 1, p, knotVector) == Approx(0.5));
  CHECK(evaluateBSplineBasis(1.00, 1, p, knotVector) == Approx(0.0));

  //  Third basis function
  REQUIRE_NOTHROW(evaluateBSplineBasis(0.0, 2, p, knotVector));

  CHECK(evaluateBSplineBasis(0.00, 2, p, knotVector) == Approx(0.0));
  CHECK(evaluateBSplineBasis(0.25, 2, p, knotVector) == Approx(0.0));
  CHECK(evaluateBSplineBasis(0.50, 2, p, knotVector) == Approx(0.0));
  CHECK(evaluateBSplineBasis(0.75, 2, p, knotVector) == Approx(0.5));
  CHECK(evaluateBSplineBasis(1.00, 2, p, knotVector) == Approx(1.0));
}

TEST_CASE("Quadratic C1 interpolation")
{
  std::vector<double> knotVector{ 0.0, 0.0, 0.0, 0.5, 1.0, 1.0, 1.0 };

  const size_t p = 2;

  // First basis functiontype ‘struct type_caster’ violates the C++ One Definition Rule
  REQUIRE_NOTHROW(evaluateBSplineBasis(0.0, 0, p, knotVector));

  CHECK(evaluateBSplineBasis(0.00, 0, p, knotVector) == Approx(1.0));
  CHECK(evaluateBSplineBasis(0.25, 0, p, knotVector) == Approx(0.25));
  CHECK(evaluateBSplineBasis(0.50, 0, p, knotVector) == Approx(0.0));
  CHECK(evaluateBSplineBasis(0.75, 0, p, knotVector) == Approx(0.0));
  CHECK(evaluateBSplineBasis(1.00, 0, p, knotVector) == Approx(0.0));

  //  Second basis function
  REQUIRE_NOTHROW(evaluateBSplineBasis(0.0, 1, p, knotVector));

  CHECK(evaluateBSplineBasis(0.00, 1, p, knotVector) == Approx(0.0));
  CHECK(evaluateBSplineBasis(0.25, 1, p, knotVector) == Approx(0.625));
  CHECK(evaluateBSplineBasis(0.50, 1, p, knotVector) == Approx(0.5));
  CHECK(evaluateBSplineBasis(0.75, 1, p, knotVector) == Approx(0.125));
  CHECK(evaluateBSplineBasis(1.00, 1, p, knotVector) == Approx(0.0));

  //  Third basis function
  REQUIRE_NOTHROW(evaluateBSplineBasis(0.0, 2, p, knotVector));

  CHECK(evaluateBSplineBasis(0.00, 2, p, knotVector) == Approx(0.0));
  CHECK(evaluateBSplineBasis(0.25, 2, p, knotVector) == Approx(0.125));
  CHECK(evaluateBSplineBasis(0.50, 2, p, knotVector) == Approx(0.5));
  CHECK(evaluateBSplineBasis(0.75, 2, p, knotVector) == Approx(0.625));
  CHECK(evaluateBSplineBasis(1.00, 2, p, knotVector) == Approx(0.0));

  //  Fourth basis function
  REQUIRE_NOTHROW(evaluateBSplineBasis(0.0, 3, p, knotVector));

  CHECK(evaluateBSplineBasis(0.00, 3, p, knotVector) == Approx(0.0));
  CHECK(evaluateBSplineBasis(0.25, 3, p, knotVector) == Approx(0.0));
  CHECK(evaluateBSplineBasis(0.50, 3, p, knotVector) == Approx(0.0));
  CHECK(evaluateBSplineBasis(0.75, 3, p, knotVector) == Approx(0.25));
  CHECK(evaluateBSplineBasis(1.00, 3, p, knotVector) == Approx(1.0));
}

} // splinekernel
} // cie
