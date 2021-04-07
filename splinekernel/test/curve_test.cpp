#include "catch.hpp"
#include "curve.hpp"

#include <array>
#include <vector>

namespace cie
{
namespace splinekernel
{

TEST_CASE("Linear interpolation curve")
{
    std::vector<double> knotVector{ 0.0, 0.0, 0.5, 1.0, 1.0 };
    std::vector<double> x{ 2.0, 3.0, 0.5 };
    std::vector<double> y{ 1.0, 3.0, 3.0 };
    std::vector<double> t{ 0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0 };
    std::array<std::vector<double>, 2> C;

    REQUIRE_NOTHROW( C = evaluate2DCurve( t, x, y, knotVector ) );

    REQUIRE( C[0].size( ) == t.size( ) );
    REQUIRE( C[1].size( ) == t.size( ) );

    // x-coordinates of curve
    CHECK( C[0][0]  == Approx( 2.0 ) );
    CHECK( C[0][1]  == Approx( 2.2 ) );
    CHECK( C[0][2]  == Approx( 2.4 ) );
    CHECK( C[0][3]  == Approx( 2.6 ) );
    CHECK( C[0][4]  == Approx( 2.8 ) );
    CHECK( C[0][5]  == Approx( 3.0 ) );
    CHECK( C[0][6]  == Approx( 2.5 ) );
    CHECK( C[0][7]  == Approx( 2.0 ) );
    CHECK( C[0][8]  == Approx( 1.5 ) );
    CHECK( C[0][9]  == Approx( 1.0 ) );
    CHECK( C[0][10] == Approx( 0.5 ) );

    // y-coordinates of curve
    CHECK( C[1][0]  == Approx( 1.0 ) );
    CHECK( C[1][1]  == Approx( 1.4 ) );
    CHECK( C[1][2]  == Approx( 1.8 ) );
    CHECK( C[1][3]  == Approx( 2.2 ) );
    CHECK( C[1][4]  == Approx( 2.6 ) );
    CHECK( C[1][5]  == Approx( 3.0 ) );
    CHECK( C[1][6]  == Approx( 3.0 ) );
    CHECK( C[1][7]  == Approx( 3.0 ) );
    CHECK( C[1][8]  == Approx( 3.0 ) );
    CHECK( C[1][9]  == Approx( 3.0 ) );
    CHECK( C[1][10] == Approx( 3.0 ) );
}

} // namespace splinekernel
} // namespace cie
