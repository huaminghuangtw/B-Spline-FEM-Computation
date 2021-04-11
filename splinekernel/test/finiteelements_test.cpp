#include "catch.hpp"
#include "finiteelements.hpp"
#include "sparse.hpp"

#include <vector>

namespace cie
{
namespace splinekernel
{

TEST_CASE("BSplineFiniteElementPatch_mapToGlobalCoordinates_test")
{
    std::array<double, 2> lengths{ 3.0, 4.5 };
    std::array<double, 2> origin{ -1.5, 2.5 };
    std::array<size_t, 2> numberOfElements{ 2, 3 };

    std::array<double, 2> rs1, rs2, rs3, rs4, rs5, rs6;

    // Valid cases
    REQUIRE_NOTHROW(rs1 = detail::mapToGlobalCoordinates({ 0.8,  0.2 }, { 0, 0 }, lengths, origin, numberOfElements));
    REQUIRE_NOTHROW(rs2 = detail::mapToGlobalCoordinates({ -0.3,  0.4 }, { 0, 1 }, lengths, origin, numberOfElements));
    REQUIRE_NOTHROW(rs3 = detail::mapToGlobalCoordinates({ -1.0,  0.6 }, { 0, 2 }, lengths, origin, numberOfElements));
    REQUIRE_NOTHROW(rs4 = detail::mapToGlobalCoordinates({ -0.9, -0.1 }, { 1, 0 }, lengths, origin, numberOfElements));
    REQUIRE_NOTHROW(rs5 = detail::mapToGlobalCoordinates({ 0.7, -0.5 }, { 1, 1 }, lengths, origin, numberOfElements));
    REQUIRE_NOTHROW(rs6 = detail::mapToGlobalCoordinates({ 0.3,  0.9 }, { 1, 2 }, lengths, origin, numberOfElements));

    CHECK(rs1[0] == Approx(-0.15));
    CHECK(rs1[1] == Approx(3.4));

    CHECK(rs2[0] == Approx(-0.975));
    CHECK(rs2[1] == Approx(5.05));

    CHECK(rs3[0] == Approx(-1.5));
    CHECK(rs3[1] == Approx(6.7));

    CHECK(rs4[0] == Approx(0.075));
    CHECK(rs4[1] == Approx(3.175));

    CHECK(rs5[0] == Approx(1.275));
    CHECK(rs5[1] == Approx(4.375));

    CHECK(rs6[0] == Approx(0.975));
    CHECK(rs6[1] == Approx(6.925));

    // Degenerate cases with at least one number of elements being zero
    CHECK_THROWS(detail::mapToGlobalCoordinates({ 0.5, 0.5 }, { 0, 0 }, lengths, origin, { 0, 0 }));
    CHECK_THROWS(detail::mapToGlobalCoordinates({ 0.5, 0.5 }, { 0, 0 }, lengths, origin, { 0, 1 }));
    CHECK_THROWS(detail::mapToGlobalCoordinates({ 0.5, 0.5 }, { 0, 0 }, lengths, origin, { 1, 0 }));

} // BSplineFiniteElementPatch_mapToGlobalCoordinates_test

TEST_CASE("BSplineFiniteElementPatch_constructOpenKnotVector_test")
{
    // Valid case
    auto result = detail::constructOpenKnotVectors({ 5, 7 }, { 2, 3 }, { 1, 2 }, { 3, 5 }, { -2, 4 });

    KnotVectors expected =
    { {
        { -2.0, -2.0, -2.0, -1.4, -0.8, -0.2, 0.4, 1.0, 1.0, 1.0 },
        {  4.0, 4.0, 4.0, 4.0, 4.71429, 5.42857, 6.14286, 6.85714, 7.57143, 8.28571, 9.0, 9.0, 9.0, 9.0 }
    } };

    for (size_t axis = 0; axis < expected.size(); ++axis)
    {
        REQUIRE(result[axis].size() == expected[axis].size());

        for (size_t index = 0; index < expected[axis].size(); ++index)
        {
            CHECK(result[axis][index] == Approx(expected[axis][index]));
        }
    }

    // Degenerate case: Continuity is higher than the polynomial degree
    std::array<size_t, 2> polynomialDegrees = { 2, 2 };

    std::array<size_t, 2> continuities1 = { 2, 1 };
    std::array<size_t, 2> continuities2 = { 1, 2 };
    std::array<size_t, 2> continuities3 = { 2, 2 };

    CHECK_THROWS(detail::constructOpenKnotVectors({ 2, 3 }, polynomialDegrees, continuities1, { 1, 1 }, { 1, 1 }));
    CHECK_THROWS(detail::constructOpenKnotVectors({ 2, 3 }, polynomialDegrees, continuities2, { 1, 1 }, { 1, 1 }));
    CHECK_THROWS(detail::constructOpenKnotVectors({ 2, 3 }, polynomialDegrees, continuities3, { 1, 1 }, { 1, 1 }));
}

TEST_CASE("findKnotSpan_test")
{
    // Spans: [ -2.0, -1.143, -0.286,  0.571,  1.429, 2.286, 3.143, 4.0 ]
    CHECK(detail::findKnotSpan(-2.0, 4.0, 7, -2.0) == 0);
    CHECK(detail::findKnotSpan(-2.0, 4.0, 7, 0.570) == 2);
    CHECK(detail::findKnotSpan(-2.0, 4.0, 7, 0.572) == 3);
    CHECK(detail::findKnotSpan(-2.0, 4.0, 7, 4.0) == 6);

    CHECK_THROWS(detail::findKnotSpan(-2.0, 4.0, 0, -2.0));
}

TEST_CASE("BSplineFiniteElementPatch_constructor_test")
{
    auto dummy_integrator = [](size_t)
    {
        return IntegrationPoints{ };
    };

    CHECK_NOTHROW(BSplineFiniteElementPatch({ 2, 3 }, { 3, 2 }, { 2, 1 }, { 1.0, 1.0 }, { 0.0, 0.0 }, dummy_integrator));
    CHECK_NOTHROW(BSplineFiniteElementPatch({ 2, 3 }, { 1, 1 }, { 0, 0 }, { 1.0, 1.0 }, { 1.0, 1.0 }, dummy_integrator));
    CHECK_NOTHROW(BSplineFiniteElementPatch({ 1, 1 }, { 2, 3 }, { 1, 1 }, { 1.0, 1.0 }, { 1.0, 1.0 }, dummy_integrator));
}

TEST_CASE("BSplineFiniteElementPatch_evaluateActiveBasisAt_test")
{
    // Expected tensor product values for three evaluation points and no derivatives: Ni(r) * Nj(s)
    std::vector<std::vector<double>> expectedBasis_00 =
    {
        { 1.20000000e-02, 1.15703704e-01, 9.21481481e-02, 2.37037037e-03, 3.90000000e-02, 3.76037037e-01,
            2.99481481e-01, 7.70370370e-03, 3.00000000e-03, 2.89259259e-02, 2.30370370e-02, 5.92592593e-04 },
        { 8.53333333e-02, 6.30666667e-01, 2.82666667e-01, 1.33333333e-03, 0.00000000e+00, 0.00000000e+00,
            0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00 },
        { 0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 5.55555556e-02, 0.00000000e+00, 0.00000000e+00,
            0.00000000e+00, 7.22222222e-01, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 2.22222222e-01 }
    };

    // Expected tensor product values for three evaluation points derived with respect to r: dNi/dr(r) * Nj(s)
    std::vector<std::vector<double>> expectedBasis_10 =
    {
        { -6.00000000e-02, -5.78518519e-01, -4.60740741e-01, -1.18518519e-02,  3.00000000e-02,  2.89259259e-01,
            2.30370370e-01,  5.92592593e-03,  3.00000000e-02,  2.89259259e-01,  2.30370370e-01,  5.92592593e-03 },
        { -2.84444444e-01, -2.10222222e+00, -9.42222222e-01, -4.44444444e-03,  2.84444444e-01,  2.10222222e+00,
            9.42222222e-01,  4.44444444e-03,  0.00000000e+00,  0.00000000e+00,  0.00000000e+00,  0.00000000e+00 },
        { -0.00000000e+00, -0.00000000e+00, -0.00000000e+00, -5.55555556e-01, -0.00000000e+00, -0.00000000e+00,
            -0.00000000e+00, -5.55555556e-01,  0.00000000e+00,  0.00000000e+00,  0.00000000e+00,  1.11111111e+00 }
    };

    // Expected tensor product values for three evaluation points derived with respect to s: Ni(r) * dNj/ds(s)
    std::vector<std::vector<double>> expectedBasis_01 =
    {
        { -8.40000000e-02, -1.46222222e-01,  2.05333333e-01,  2.48888889e-02, -2.73000000e-01, -4.75222222e-01,
            6.67333333e-01,  8.08888889e-02, -2.10000000e-02, -3.65555556e-02,  5.13333333e-02,  6.22222222e-03 },
        { -4.48000000e-01, -4.76000000e-01,  8.96000000e-01,  2.80000000e-02, -0.00000000e+00, -0.00000000e+00,
            0.00000000e+00,  0.00000000e+00, -0.00000000e+00, -0.00000000e+00,  0.00000000e+00,  0.00000000e+00 },
        {  0.00000000e+00,  0.00000000e+00, -2.33333333e-01,  2.33333333e-01,  0.00000000e+00,  0.00000000e+00,
            -3.03333333e+00,  3.03333333e+00,  0.00000000e+00,  0.00000000e+00, -9.33333333e-01,  9.33333333e-01 }
    };

    std::array<size_t, 2> numberOfElements{ 5, 7 };
    std::array<size_t, 2> polynomialDegrees{ 2, 3 };
    std::array<size_t, 2> continuities{ 1, 2 };

    std::array<double, 2> lengths{ 3.0, 5.0 };
    std::array<double, 2> origin{ -2.0, 4.0 };

    auto dummy_integrator = [](size_t) { return IntegrationPoints{ }; };

    auto mesh = BSplineFiniteElementPatch(numberOfElements, polynomialDegrees, continuities,
        lengths, origin, dummy_integrator);

    std::vector<std::vector<double>> computedBasis_00(3);
    std::vector<std::vector<double>> computedBasis_10(3);
    std::vector<std::vector<double>> computedBasis_01(3);

    // Evaluation point 1
    REQUIRE_NOTHROW(computedBasis_00[0] = mesh.evaluateActiveBasisAt({ 0.0,  5.0 }, { 0, 0 }));
    REQUIRE_NOTHROW(computedBasis_01[0] = mesh.evaluateActiveBasisAt({ 0.0,  5.0 }, { 0, 1 }));
    REQUIRE_NOTHROW(computedBasis_10[0] = mesh.evaluateActiveBasisAt({ 0.0,  5.0 }, { 1, 0 }));

    // Evaluation point 2
    REQUIRE_NOTHROW(computedBasis_00[1] = mesh.evaluateActiveBasisAt({ -2.0,  7.0 }, { 0, 0 }));
    REQUIRE_NOTHROW(computedBasis_01[1] = mesh.evaluateActiveBasisAt({ -2.0,  7.0 }, { 0, 1 }));
    REQUIRE_NOTHROW(computedBasis_10[1] = mesh.evaluateActiveBasisAt({ -2.0,  7.0 }, { 1, 0 }));

    // Evaluation point 3
    REQUIRE_NOTHROW(computedBasis_00[2] = mesh.evaluateActiveBasisAt({ -1.0, 9.0 }, { 0, 0 }));
    REQUIRE_NOTHROW(computedBasis_01[2] = mesh.evaluateActiveBasisAt({ -1.0, 9.0 }, { 0, 1 }));
    REQUIRE_NOTHROW(computedBasis_10[2] = mesh.evaluateActiveBasisAt({ -1.0, 9.0 }, { 1, 0 }));

    size_t expectedSize = 12; // (p0 + 1) * (p1 + 1)

    for (size_t iPoint = 0; iPoint < 3; ++iPoint)
    {
        REQUIRE(computedBasis_00[iPoint].size() == expectedSize);
        REQUIRE(computedBasis_01[iPoint].size() == expectedSize);
        REQUIRE(computedBasis_10[iPoint].size() == expectedSize);

        for (size_t iFunction = 0; iFunction < expectedSize; ++iFunction)
        {
            CHECK(computedBasis_00[iPoint][iFunction] == Approx(expectedBasis_00[iPoint][iFunction]).margin(1e-12));
            CHECK(computedBasis_01[iPoint][iFunction] == Approx(expectedBasis_01[iPoint][iFunction]).margin(1e-12));
            CHECK(computedBasis_10[iPoint][iFunction] == Approx(expectedBasis_10[iPoint][iFunction]).margin(1e-12));

        } // iFunction
    } // iPoint

} // BSplineFiniteElementPatch_evaluateActiveBasisAt_test

TEST_CASE("BSplineFiniteElementPatch_integrateElementSystem_test")
{
    auto gaussIntegrator = [](size_t size) -> IntegrationPoints
    {
        if (size == 3)
        {
            return { { { -0.77459667, 0.0,        0.77459667 },
                        {  0.55555556, 0.88888889, 0.55555556 } } };
        }
        if (size == 4)
        {
            return { { { -0.86113631, -0.33998104,  0.33998104,  0.86113631 },
                        {  0.34785485,  0.65214515,  0.65214515,  0.34785485 } } };
        }
        else
        {
            throw std::runtime_error("Unexpected integration order.");
        }
    };

    auto testSourceFunction = [](double x, double y)
    {
        return x * y;
    };

    std::vector<std::vector<double>> expectedElementMatrix =
    {
        {  8.26808398e-03,  1.61542329e-02,  7.86167790e-04, -4.06897207e-04,  8.46595810e-03, -2.47711646e-03,
            -1.71930840e-02, -1.19655140e-03, -9.84041956e-04, -6.67711643e-03, -4.59308391e-03, -1.46551399e-04 },
        {  1.61542329e-02,  9.37492510e-02,  6.30478840e-02,  6.59744262e-04, -2.47711646e-03, -2.02746255e-02,
            -5.29906089e-02, -1.10632055e-02, -6.67711643e-03, -4.02246253e-02, -3.68906085e-02, -3.01320547e-03 },
        {  7.86167790e-04,  6.30478840e-02,  1.05437416e-01,  1.26068406e-02, -1.71930840e-02, -5.29906089e-02,
            -2.09853745e-02,  2.29913046e-04, -4.59308391e-03, -3.68906085e-02, -4.47853743e-02, -4.67008694e-03 },
        { -4.06897207e-04,  6.59744262e-04,  1.26068406e-02,  3.67470399e-03, -1.19655140e-03, -1.10632055e-02,
            2.29913046e-04,  3.76264804e-03, -1.46551399e-04, -3.01320547e-03, -4.67008694e-03, -4.37351981e-04 },
        {  8.46595810e-03, -2.47711646e-03, -1.71930840e-02, -1.19655140e-03,  4.60680842e-02,  3.29542330e-02,
            -4.96138324e-02, -4.60689722e-03,  8.46595810e-03, -2.47711646e-03, -1.71930840e-02, -1.19655140e-03 },
        { -2.47711646e-03, -2.02746255e-02, -5.29906089e-02, -1.10632055e-02,  3.29542330e-02,  1.73549252e-01,
            -1.35211593e-03, -3.15402558e-02, -2.47711646e-03, -2.02746255e-02, -5.29906089e-02, -1.10632055e-02 },
        { -1.71930840e-02, -5.29906089e-02, -2.09853745e-02,  2.29913046e-04, -4.96138324e-02, -1.35211593e-03,
            2.00637416e-01,  3.22068407e-02, -1.71930840e-02, -5.29906089e-02, -2.09853745e-02,  2.29913046e-04 },
        { -1.19655140e-03, -1.10632055e-02,  2.29913046e-04,  3.76264804e-03, -4.60689722e-03, -3.15402558e-02,
            3.22068407e-02,  2.04747041e-02, -1.19655140e-03, -1.10632055e-02,  2.29913046e-04,  3.76264804e-03 },
        { -9.84041956e-04, -6.67711643e-03, -4.59308391e-03, -1.46551399e-04,  8.46595810e-03, -2.47711646e-03,
            -1.71930840e-02, -1.19655140e-03,  8.26808398e-03,  1.61542329e-02,  7.86167790e-04, -4.06897207e-04 },
        { -6.67711643e-03, -4.02246253e-02, -3.68906085e-02, -3.01320547e-03, -2.47711646e-03, -2.02746255e-02,
            -5.29906089e-02, -1.10632055e-02,  1.61542329e-02,  9.37492510e-02,  6.30478840e-02,  6.59744262e-04 },
        { -4.59308391e-03, -3.68906085e-02, -4.47853743e-02, -4.67008694e-03, -1.71930840e-02, -5.29906089e-02,
            -2.09853745e-02,  2.29913046e-04,  7.86167790e-04,  6.30478840e-02,  1.05437416e-01,  1.26068406e-02 },
        { -1.46551399e-04, -3.01320547e-03, -4.67008694e-03, -4.37351981e-04, -1.19655140e-03, -1.10632055e-02,
            2.29913046e-04,  3.76264804e-03, -4.06897207e-04,  6.59744262e-04,  1.26068406e-02,  3.67470399e-03 }
    };

    std::vector<double> expectedElementVector =
    {
        { -2.71045920e-02, -1.95578233e-01, -2.10459185e-01, -1.96641158e-02, -9.54081637e-02, -6.88435377e-01,
          -7.40816330e-01, -6.92176874e-02, -2.05994899e-02, -1.48639457e-01, -1.59948980e-01, -1.49447280e-02 }
    };

    auto mesh = BSplineFiniteElementPatch({ 5, 7 }, { 2, 3 }, { 1, 2 }, { 3.0, 5.0 }, { -2.0, 4.0 }, gaussIntegrator);

    auto computedElementSystem = mesh.integrateElementSystem({ 1, 1 }, testSourceFunction);

    // Check the element matrix
    const auto& computedElementMatrix = std::get<0>(computedElementSystem);

    REQUIRE(computedElementMatrix.size1() == expectedElementMatrix.size());
    REQUIRE(computedElementMatrix.size2() == expectedElementMatrix.size());

    for (size_t i = 0; i < expectedElementMatrix.size(); ++i)
    {
        for (size_t j = 0; j < expectedElementMatrix.size(); ++j)
        {
            CHECK(computedElementMatrix(i, j) == Approx(expectedElementMatrix[i][j]));
        }
    }

    // Check the element source vector
    const auto& computedElementRhs = std::get<1>(computedElementSystem);

    REQUIRE(computedElementRhs.size() == expectedElementVector.size());

    for (size_t i = 0; i < expectedElementVector.size(); ++i)
    {
        CHECK(computedElementRhs[i] == Approx(expectedElementVector[i]));
    }
}

TEST_CASE("BSplineFiniteElementPatch_constructLocationMaps_test1")
{
    // n = (3, 2), p = (1, 2), c = (0, 0)
    // 
    // shape function indices:
    // 4   9  14  19
    // 3   8  13  18
    // 2   7  12  17
    // 1   6  11  16
    // 0   5  10  15

    LocationMaps expectedLocationMaps
    {
        {  0,  1,  2,  5,  6,  7 },
        {  2,  3,  4,  7,  8,  9 },
        {  5,  6,  7, 10, 11, 12 },
        {  7,  8,  9, 12, 13, 14 },
        { 10, 11, 12, 15, 16, 17 },
        { 12, 13, 14, 17, 18, 19 }
    };

    auto computedLocationMaps = detail::constructLocationMaps({ 3, 2 }, { 1, 2 }, { 0, 0 });

    REQUIRE(computedLocationMaps.size() == expectedLocationMaps.size());

    for (size_t iElement = 0; iElement < expectedLocationMaps.size(); ++iElement)
    {
        CHECK(computedLocationMaps[iElement] == expectedLocationMaps[iElement]);
    }
}

TEST_CASE("BSplineFiniteElementPatch_constructLocationMaps_test2")
{
    // Zero number of elements
    CHECK_THROWS(detail::constructLocationMaps({ 0, 0 }, { 1, 1 }, { 0, 0 }));
    CHECK_THROWS(detail::constructLocationMaps({ 0, 1 }, { 1, 1 }, { 0, 0 }));
    CHECK_THROWS(detail::constructLocationMaps({ 1, 0 }, { 1, 1 }, { 0, 0 }));

    // Zero polynomial degree
    CHECK_THROWS(detail::constructLocationMaps({ 2, 2 }, { 0, 0 }, { 0, 0 }));
    CHECK_THROWS(detail::constructLocationMaps({ 2, 2 }, { 0, 1 }, { 0, 0 }));
    CHECK_THROWS(detail::constructLocationMaps({ 2, 2 }, { 1, 0 }, { 0, 0 }));

    // Invalid continuity
    CHECK_THROWS(detail::constructLocationMaps({ 2, 2 }, { 3, 3 }, { 3, 3 }));
    CHECK_THROWS(detail::constructLocationMaps({ 2, 2 }, { 3, 3 }, { 3, 0 }));
    CHECK_THROWS(detail::constructLocationMaps({ 2, 2 }, { 3, 3 }, { 0, 3 }));
}

TEST_CASE("BSplineFiniteElementPatch_constructLocationMaps_test3")
{
    // n = (1, 1), p = (2, 3), c = (1, 0)
    auto computedLocationMaps = detail::constructLocationMaps({ 1, 1 }, { 2, 3 }, { 1, 0 });

    REQUIRE(computedLocationMaps.size() == 1);
    REQUIRE(computedLocationMaps[0].size() == 12);

    for (size_t i = 0; i < 12; ++i)
    {
        CHECK(computedLocationMaps[0][i] == i);
    }
}

TEST_CASE("BSplineFiniteElementPatch_constructLocationMaps_test4")
{
    // n = (4, 5), p = (4, 3), c = (2, 1)
    auto computedLocationMaps = detail::constructLocationMaps({ 4, 5 }, { 4, 3 }, { 2, 1 });

    std::vector<std::vector<size_t>> indices(11, std::vector<size_t>(12));

    for (size_t i = 0; i < 11; ++i)
    {
        for (size_t j = 0; j < 12; ++j)
        {
            indices[i][j] = i * 12 + j;
        }
    }

    std::vector<size_t> startX{ 0, 2, 4, 6 };
    std::vector<size_t> startY{ 0, 2, 4, 6, 8 };

    REQUIRE(computedLocationMaps.size() == 20);

    for (size_t iElement = 0; iElement < 4; ++iElement)
    {
        for (size_t jElement = 0; jElement < 5; ++jElement)
        {
            size_t elementIndex = iElement * 5 + jElement;

            REQUIRE(computedLocationMaps[elementIndex].size() == 20);

            for (size_t iDof = 0; iDof < 5; ++iDof)
            {
                for (size_t jDof = 0; jDof < 4; ++jDof)
                {
                    size_t dofIndex = iDof * 4 + jDof;

                    size_t i = startX[iElement] + iDof;
                    size_t j = startY[jElement] + jDof;

                    CHECK(computedLocationMaps[elementIndex][dofIndex] == indices[i][j]);

                } // jDof
            } // iDof

        } // jElement
    } // iElement

} // BSplineFiniteElementPatch_constructLocationMaps_test4

TEST_CASE( "BSplineFiniteElementPatch_boundaryDofIds_test1" )
{
    auto dummy = []( size_t ) -> std::array<std::vector<double>, 2>
    {
        throw std::runtime_error( "Not implemented." );
    };

    // from BSplineFiniteElementPatch_constructLocationMaps_test1 
    BSplineFiniteElementPatch mesh( { 3, 2 }, { 1, 2 }, { 0, 0 }, { 3.0, 2.0 }, { 0.2, 0.4 }, dummy );

    std::vector<size_t> expectedLeftIndices { 0, 1, 2, 3, 4 };
    std::vector<size_t> expectedTopIndices { 4, 9, 14, 19 };
    std::vector<size_t> expectedRightIndices { 15, 16, 17, 18, 19 };
    std::vector<size_t> expectedBottomIndices { 0, 5, 10, 15 };

    CHECK( mesh.boundaryDofIds( "left" ) == expectedLeftIndices );
    CHECK( mesh.boundaryDofIds( "top" ) == expectedTopIndices );
    CHECK( mesh.boundaryDofIds( "right" ) == expectedRightIndices );
    CHECK( mesh.boundaryDofIds( "bottom" ) == expectedBottomIndices );

    CHECK_THROWS( mesh.boundaryDofIds( "" ) );
    CHECK_THROWS( mesh.boundaryDofIds( "wrong" ) );
}

TEST_CASE( "BSplineFiniteElementPatch_boundaryDofIds_test2" )
{
    auto dummy = []( size_t )->std::array<std::vector<double>, 2>
    {
        throw std::runtime_error( "Not implemented." );
    };

    // from BSplineFiniteElementPatch_constructLocationMaps_test4 
    BSplineFiniteElementPatch mesh( { 4, 5 }, { 4, 3 }, { 2, 1 }, { 2.0, 1.4 }, { 2.1, 1.7 }, dummy );

    std::vector<size_t> expectedLeftIndices { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11 };
    std::vector<size_t> expectedTopIndices { 11, 23, 35, 47, 59, 71, 83, 95, 107, 119, 131 };
    std::vector<size_t> expectedRightIndices { 120, 121, 122, 123, 124, 125, 126, 127, 128, 129, 130, 131 };
    std::vector<size_t> expectedBottomIndices { 0, 12, 24, 36, 48, 60, 72, 84, 96, 108, 120 };

    CHECK( mesh.boundaryDofIds( "left" ) == expectedLeftIndices );
    CHECK( mesh.boundaryDofIds( "top" ) == expectedTopIndices );
    CHECK( mesh.boundaryDofIds( "right" ) == expectedRightIndices );
    CHECK( mesh.boundaryDofIds( "bottom" ) == expectedBottomIndices );
}

TEST_CASE("BSplineFiniteElementPatch_assembleGlobalSystem_test")
{
    IntegrationPointProvider provider = [](size_t order)->std::array<std::vector<double>, 2>
    {
        if (order == 3)
        {
            return { { { -0.77459667, 0.0       , 0.77459667 },
                      {  0.55555556, 0.88888889, 0.55555556 } } };
        }
        else if (order == 5)
        {
            return { { {-0.90617985, -0.53846931, 0.0       , 0.53846931, 0.90617985 },
                      { 0.23692689,  0.47862867, 0.56888889, 0.47862867, 0.23692689 }} };
        }
        else
        {
            throw std::runtime_error{ "Integration order not implemented." };
        }
    };

    SpatialFunction source = [](double x, double y)
    {
        return x * x + y * y + 1.0;
    };

    BSplineFiniteElementPatch mesh({ 3, 2 }, { 2, 4 }, { 1, 2 }, { 2.5, 3.5 }, { -1.5, 0.5 }, provider);

    auto system = mesh.assembleGlobalSystem(source);

    auto& computedLhs = system.first;
    auto& computedRhs = system.second;

    REQUIRE(computedLhs.size() == 35);
    REQUIRE(computedRhs.size() == 35);

    REQUIRE(computedLhs.nnz() == 779);

    std::vector<CompressedSparseRowMatrix::IndexType> expectedIndices
    {
        0, 1, 2, 3, 4, 7, 8, 9, 10, 11, 14, 15, 16, 17, 18, 0, 1, 2, 3, 4, 7, 8, 9, 10, 11, 14, 15, 16, 17,
        18, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 0, 1, 2, 3, 4, 5, 6,
        7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13,
        14, 15, 16, 17, 18, 19, 20, 2, 3, 4, 5, 6, 9, 10, 11, 12, 13, 16, 17, 18, 19, 20, 2, 3, 4, 5, 6, 9,
        10, 11, 12, 13, 16, 17, 18, 19, 20, 0, 1, 2, 3, 4, 7, 8, 9, 10, 11, 14, 15, 16, 17, 18, 21, 22, 23,
        24, 25, 0, 1, 2, 3, 4, 7, 8, 9, 10, 11, 14, 15, 16, 17, 18, 21, 22, 23, 24, 25, 0, 1, 2, 3, 4, 5, 6,
        7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 0, 1, 2, 3, 4, 5,
        6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 0, 1, 2, 3, 4,
        5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 2, 3, 4, 5,
        6, 9, 10, 11, 12, 13, 16, 17, 18, 19, 20, 23, 24, 25, 26, 27, 2, 3, 4, 5, 6, 9, 10, 11, 12, 13, 16,
        17, 18, 19, 20, 23, 24, 25, 26, 27, 0, 1, 2, 3, 4, 7, 8, 9, 10, 11, 14, 15, 16, 17, 18, 21, 22, 23,
        24, 25, 28, 29, 30, 31, 32, 0, 1, 2, 3, 4, 7, 8, 9, 10, 11, 14, 15, 16, 17, 18, 21, 22, 23, 24, 25,
        28, 29, 30, 31, 32, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21,
        22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13,
        14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 0, 1, 2, 3, 4, 5,
        6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31,
        32, 33, 34, 2, 3, 4, 5, 6, 9, 10, 11, 12, 13, 16, 17, 18, 19, 20, 23, 24, 25, 26, 27, 30, 31, 32, 33,
        34, 2, 3, 4, 5, 6, 9, 10, 11, 12, 13, 16, 17, 18, 19, 20, 23, 24, 25, 26, 27, 30, 31, 32, 33, 34, 7,
        8, 9, 10, 11, 14, 15, 16, 17, 18, 21, 22, 23, 24, 25, 28, 29, 30, 31, 32, 7, 8, 9, 10, 11, 14, 15, 16,
        17, 18, 21, 22, 23, 24, 25, 28, 29, 30, 31, 32, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20,
        21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18,
        19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16,
        17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 9, 10, 11, 12, 13, 16, 17, 18,
        19, 20, 23, 24, 25, 26, 27, 30, 31, 32, 33, 34, 9, 10, 11, 12, 13, 16, 17, 18, 19, 20, 23, 24, 25, 26,
        27, 30, 31, 32, 33, 34, 14, 15, 16, 17, 18, 21, 22, 23, 24, 25, 28, 29, 30, 31, 32, 14, 15, 16, 17,
        18, 21, 22, 23, 24, 25, 28, 29, 30, 31, 32, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27,
        28, 29, 30, 31, 32, 33, 34, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31,
        32, 33, 34, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 16,
        17, 18, 19, 20, 23, 24, 25, 26, 27, 30, 31, 32, 33, 34, 16, 17, 18, 19, 20, 23, 24, 25, 26, 27, 30,
        31, 32, 33, 34
    };

    std::vector<CompressedSparseRowMatrix::IndexType> expectedIndptr
    {
        0, 15, 30, 51, 72, 93, 108, 123, 143, 163, 191, 219, 247, 267, 287, 312, 337, 372, 407, 442, 467, 492,
        512, 532, 560, 588, 616, 636, 656, 671, 686, 707, 728, 749, 764, 779
    };

    std::vector<double> expectedData
    {
         5.287982e-1,  4.671201e-2, -5.464853e-3, -8.435374e-3, -1.609977e-3, -1.063492e-1, -1.801587e-1, -1.083730e-1,
        -2.269841e-2, -2.420635e-3, -5.963719e-2, -4.795918e-2, -2.675170e-2, -5.147392e-3, -5.045351e-4,  4.671201e-2,
         3.083900e-1,  1.880499e-1,  1.945578e-2, -2.607710e-3, -1.801587e-1, -5.714285e-2, -1.228968e-1, -5.087302e-2,
        -8.928571e-3, -4.795918e-2, -3.356009e-2, -4.247732e-2, -1.393424e-2, -2.069161e-3, -5.464853e-3,  1.880499e-1,
         5.468481e-1,  3.163265e-1,  7.845805e-2, -2.607710e-3, -1.609977e-3, -1.083730e-1, -1.228968e-1, -2.795238e-1,
        -2.154762e-1, -1.023810e-1, -8.928571e-3, -2.420635e-3, -2.675170e-2, -4.247732e-2, -1.040590e-1, -7.363946e-2,
        -3.049887e-2, -2.069161e-3, -5.045351e-4, -8.435374e-3,  1.945578e-2,  3.163265e-1,  4.653061e-1,  3.163265e-1,
         1.945578e-2, -8.435374e-3, -2.269841e-2, -5.087302e-2, -2.154762e-1, -2.619048e-1, -2.154762e-1, -5.087302e-2,
        -2.269841e-2, -5.147392e-3, -1.393424e-2, -7.363946e-2, -9.455782e-2, -7.363946e-2, -1.393424e-2, -5.147392e-3,
        -1.609977e-3, -2.607710e-3,  7.845805e-2,  3.163265e-1,  5.468481e-1,  1.880499e-1, -5.464853e-3, -2.420635e-3,
        -8.928571e-3, -1.023810e-1, -2.154762e-1, -2.795238e-1, -1.228968e-1, -1.083730e-1, -5.045351e-4, -2.069161e-3,
        -3.049887e-2, -7.363946e-2, -1.040590e-1, -4.247732e-2, -2.675170e-2, -2.607710e-3,  1.945578e-2,  1.880499e-1,
         3.083900e-1,  4.671201e-2, -8.928571e-3, -5.087302e-2, -1.228968e-1, -5.714285e-2, -1.801587e-1, -2.069161e-3,
        -1.393424e-2, -4.247732e-2, -3.356009e-2, -4.795918e-2, -1.609977e-3, -8.435374e-3, -5.464853e-3,  4.671201e-2,
         5.287982e-1, -2.420635e-3, -2.269841e-2, -1.083730e-1, -1.801587e-1, -1.063492e-1, -5.045351e-4, -5.147392e-3,
        -2.675170e-2, -4.795918e-2, -5.963719e-2, -1.063492e-1, -1.801587e-1, -1.083730e-1, -2.269841e-2, -2.420635e-3,
         6.739229e-1, -2.585035e-2, -6.170068e-2, -2.294785e-2, -3.424036e-3,  1.878685e-1, -1.328231e-1, -9.772959e-2,
        -2.434240e-2, -2.973356e-3, -2.981860e-2, -2.397959e-2, -1.337585e-2, -2.573696e-3, -2.522676e-4, -1.801587e-1,
        -5.714285e-2, -1.228968e-1, -5.087302e-2, -8.928571e-3, -2.585035e-2,  3.954649e-1,  1.971202e-1,  1.315193e-3,
        -8.049887e-3, -1.328231e-1,  1.138322e-1, -7.633220e-3, -3.417800e-2, -9.197846e-3, -2.397959e-2, -1.678005e-2,
        -2.123866e-2, -6.967120e-3, -1.034581e-3, -1.083730e-1, -1.228968e-1, -2.795238e-1, -2.154762e-1, -1.023810e-1,
        -8.928571e-3, -2.420635e-3, -6.170068e-2,  1.971202e-1,  6.121542e-1,  3.272109e-1,  5.668934e-2, -8.049887e-3,
        -3.424036e-3, -9.772959e-2, -7.633220e-3,  4.592971e-2, -2.049320e-2, -4.790250e-2, -9.197846e-3, -2.973356e-3,
        -1.337585e-2, -2.123866e-2, -5.202948e-2, -3.681973e-2, -1.524943e-2, -1.034581e-3, -2.522676e-4, -2.269841e-2,
        -5.087302e-2, -2.154762e-1, -2.619048e-1, -2.154762e-1, -5.087302e-2, -2.269841e-2, -2.294785e-2,  1.315193e-3,
         3.272109e-1,  5.088435e-1,  3.272109e-1,  1.315193e-3, -2.294785e-2, -2.434240e-2, -3.417800e-2, -2.049320e-2,
         1.802721e-2, -2.049320e-2, -3.417800e-2, -2.434240e-2, -2.573696e-3, -6.967120e-3, -3.681973e-2, -4.727891e-2,
        -3.681973e-2, -6.967120e-3, -2.573696e-3, -2.420635e-3, -8.928571e-3, -1.023810e-1, -2.154762e-1, -2.795238e-1,
        -1.228968e-1, -1.083730e-1, -3.424036e-3, -8.049887e-3,  5.668934e-2,  3.272109e-1,  6.121542e-1,  1.971202e-1,
        -6.170068e-2, -2.973356e-3, -9.197846e-3, -4.790250e-2, -2.049320e-2,  4.592971e-2, -7.633220e-3, -9.772959e-2,
        -2.522676e-4, -1.034581e-3, -1.524943e-2, -3.681973e-2, -5.202948e-2, -2.123866e-2, -1.337585e-2, -8.928571e-3,
        -5.087302e-2, -1.228968e-1, -5.714285e-2, -1.801587e-1, -8.049887e-3,  1.315193e-3,  1.971202e-1,  3.954649e-1,
        -2.585035e-2, -9.197846e-3, -3.417800e-2, -7.633220e-3,  1.138322e-1, -1.328231e-1, -1.034581e-3, -6.967120e-3,
        -2.123866e-2, -1.678005e-2, -2.397959e-2, -2.420635e-3, -2.269841e-2, -1.083730e-1, -1.801587e-1, -1.063492e-1,
        -3.424036e-3, -2.294785e-2, -6.170068e-2, -2.585035e-2,  6.739229e-1, -2.973356e-3, -2.434240e-2, -9.772959e-2,
        -1.328231e-1,  1.878685e-1, -2.522676e-4, -2.573696e-3, -1.337585e-2, -2.397959e-2, -2.981860e-2, -5.963719e-2,
        -4.795918e-2, -2.675170e-2, -5.147392e-3, -5.045351e-4,  1.878685e-1, -1.328231e-1, -9.772959e-2, -2.434240e-2,
        -2.973356e-3,  8.319728e-1, -1.826531e-1, -1.728061e-1, -4.986395e-2, -6.649660e-3,  1.878685e-1, -1.328231e-1,
        -9.772959e-2, -2.434240e-2, -2.973356e-3, -5.963719e-2, -4.795918e-2, -2.675170e-2, -5.147392e-3, -5.045351e-4,
        -4.795918e-2, -3.356009e-2, -4.247732e-2, -1.393424e-2, -2.069161e-3, -1.328231e-1,  1.138322e-1, -7.633220e-3,
        -3.417800e-2, -9.197846e-3, -1.826531e-1,  4.925170e-1,  1.682483e-1, -3.982993e-2, -1.828231e-2, -1.328231e-1,
         1.138322e-1, -7.633220e-3, -3.417800e-2, -9.197846e-3, -4.795918e-2, -3.356009e-2, -4.247732e-2, -1.393424e-2,
        -2.069161e-3, -2.675170e-2, -4.247732e-2, -1.040590e-1, -7.363946e-2, -3.049887e-2, -2.069161e-3, -5.045351e-4,
        -9.772959e-2, -7.633220e-3,  4.592971e-2, -2.049320e-2, -4.790250e-2, -9.197846e-3, -2.973356e-3, -1.728061e-1,
         1.682483e-1,  6.060544e-1,  2.698980e-1, -6.462586e-3, -1.828231e-2, -6.649660e-3, -9.772959e-2, -7.633220e-3,
         4.592971e-2, -2.049320e-2, -4.790250e-2, -9.197846e-3, -2.973356e-3, -2.675170e-2, -4.247732e-2, -1.040590e-1,
        -7.363946e-2, -3.049887e-2, -2.069161e-3, -5.045351e-4, -5.147392e-3, -1.393424e-2, -7.363946e-2, -9.455782e-2,
        -7.363946e-2, -1.393424e-2, -5.147392e-3, -2.434240e-2, -3.417800e-2, -2.049320e-2,  1.802721e-2, -2.049320e-2,
        -3.417800e-2, -2.434240e-2, -4.986395e-2, -3.982993e-2,  2.698980e-1,  4.795918e-1,  2.698980e-1, -3.982993e-2,
        -4.986395e-2, -2.434240e-2, -3.417800e-2, -2.049320e-2,  1.802721e-2, -2.049320e-2, -3.417800e-2, -2.434240e-2,
        -5.147392e-3, -1.393424e-2, -7.363946e-2, -9.455782e-2, -7.363946e-2, -1.393424e-2, -5.147392e-3, -5.045351e-4,
        -2.069161e-3, -3.049887e-2, -7.363946e-2, -1.040590e-1, -4.247732e-2, -2.675170e-2, -2.973356e-3, -9.197846e-3,
        -4.790250e-2, -2.049320e-2,  4.592971e-2, -7.633220e-3, -9.772959e-2, -6.649660e-3, -1.828231e-2, -6.462586e-3,
         2.698980e-1,  6.060544e-1,  1.682483e-1, -1.728061e-1, -2.973356e-3, -9.197846e-3, -4.790250e-2, -2.049320e-2,
         4.592971e-2, -7.633220e-3, -9.772959e-2, -5.045351e-4, -2.069161e-3, -3.049887e-2, -7.363946e-2, -1.040590e-1,
        -4.247732e-2, -2.675170e-2, -2.069161e-3, -1.393424e-2, -4.247732e-2, -3.356009e-2, -4.795918e-2, -9.197846e-3,
        -3.417800e-2, -7.633220e-3,  1.138322e-1, -1.328231e-1, -1.828231e-2, -3.982993e-2,  1.682483e-1,  4.925170e-1,
        -1.826531e-1, -9.197846e-3, -3.417800e-2, -7.633220e-3,  1.138322e-1, -1.328231e-1, -2.069161e-3, -1.393424e-2,
        -4.247732e-2, -3.356009e-2, -4.795918e-2, -5.045351e-4, -5.147392e-3, -2.675170e-2, -4.795918e-2, -5.963719e-2,
        -2.973356e-3, -2.434240e-2, -9.772959e-2, -1.328231e-1,  1.878685e-1, -6.649660e-3, -4.986395e-2, -1.728061e-1,
        -1.826531e-1,  8.319728e-1, -2.973356e-3, -2.434240e-2, -9.772959e-2, -1.328231e-1,  1.878685e-1, -5.045351e-4,
        -5.147392e-3, -2.675170e-2, -4.795918e-2, -5.963719e-2, -2.981860e-2, -2.397959e-2, -1.337585e-2, -2.573696e-3,
        -2.522676e-4,  1.878685e-1, -1.328231e-1, -9.772959e-2, -2.434240e-2, -2.973356e-3,  6.739229e-1, -2.585035e-2,
        -6.170068e-2, -2.294785e-2, -3.424036e-3, -1.063492e-1, -1.801587e-1, -1.083730e-1, -2.269841e-2, -2.420635e-3,
        -2.397959e-2, -1.678005e-2, -2.123866e-2, -6.967120e-3, -1.034581e-3, -1.328231e-1,  1.138322e-1, -7.633220e-3,
        -3.417800e-2, -9.197846e-3, -2.585035e-2,  3.954649e-1,  1.971202e-1,  1.315193e-3, -8.049887e-3, -1.801587e-1,
        -5.714285e-2, -1.228968e-1, -5.087302e-2, -8.928571e-3, -1.337585e-2, -2.123866e-2, -5.202948e-2, -3.681973e-2,
        -1.524943e-2, -1.034581e-3, -2.522676e-4, -9.772959e-2, -7.633220e-3,  4.592971e-2, -2.049320e-2, -4.790250e-2,
        -9.197846e-3, -2.973356e-3, -6.170068e-2,  1.971202e-1,  6.121542e-1,  3.272109e-1,  5.668934e-2, -8.049887e-3,
        -3.424036e-3, -1.083730e-1, -1.228968e-1, -2.795238e-1, -2.154762e-1, -1.023810e-1, -8.928571e-3, -2.420635e-3,
        -2.573696e-3, -6.967120e-3, -3.681973e-2, -4.727891e-2, -3.681973e-2, -6.967120e-3, -2.573696e-3, -2.434240e-2,
        -3.417800e-2, -2.049320e-2,  1.802721e-2, -2.049320e-2, -3.417800e-2, -2.434240e-2, -2.294785e-2,  1.315193e-3,
         3.272109e-1,  5.088435e-1,  3.272109e-1,  1.315193e-3, -2.294785e-2, -2.269841e-2, -5.087302e-2, -2.154762e-1,
        -2.619048e-1, -2.154762e-1, -5.087302e-2, -2.269841e-2, -2.522676e-4, -1.034581e-3, -1.524943e-2, -3.681973e-2,
        -5.202948e-2, -2.123866e-2, -1.337585e-2, -2.973356e-3, -9.197846e-3, -4.790250e-2, -2.049320e-2,  4.592971e-2,
        -7.633220e-3, -9.772959e-2, -3.424036e-3, -8.049887e-3,  5.668934e-2,  3.272109e-1,  6.121542e-1,  1.971202e-1,
        -6.170068e-2, -2.420635e-3, -8.928571e-3, -1.023810e-1, -2.154762e-1, -2.795238e-1, -1.228968e-1, -1.083730e-1,
        -1.034581e-3, -6.967120e-3, -2.123866e-2, -1.678005e-2, -2.397959e-2, -9.197846e-3, -3.417800e-2, -7.633220e-3,
         1.138322e-1, -1.328231e-1, -8.049887e-3,  1.315193e-3,  1.971202e-1,  3.954649e-1, -2.585035e-2, -8.928571e-3,
        -5.087302e-2, -1.228968e-1, -5.714285e-2, -1.801587e-1, -2.522676e-4, -2.573696e-3, -1.337585e-2, -2.397959e-2,
        -2.981860e-2, -2.973356e-3, -2.434240e-2, -9.772959e-2, -1.328231e-1,  1.878685e-1, -3.424036e-3, -2.294785e-2,
        -6.170068e-2, -2.585035e-2,  6.739229e-1, -2.420635e-3, -2.269841e-2, -1.083730e-1, -1.801587e-1, -1.063492e-1,
        -5.963719e-2, -4.795918e-2, -2.675170e-2, -5.147392e-3, -5.045351e-4, -1.063492e-1, -1.801587e-1, -1.083730e-1,
        -2.269841e-2, -2.420635e-3,  5.287982e-1,  4.671201e-2, -5.464853e-3, -8.435374e-3, -1.609977e-3, -4.795918e-2,
        -3.356009e-2, -4.247732e-2, -1.393424e-2, -2.069161e-3, -1.801587e-1, -5.714285e-2, -1.228968e-1, -5.087302e-2,
        -8.928571e-3,  4.671201e-2,  3.083900e-1,  1.880499e-1,  1.945578e-2, -2.607710e-3, -2.675170e-2, -4.247732e-2,
        -1.040590e-1, -7.363946e-2, -3.049887e-2, -2.069161e-3, -5.045351e-4, -1.083730e-1, -1.228968e-1, -2.795238e-1,
        -2.154762e-1, -1.023810e-1, -8.928571e-3, -2.420635e-3, -5.464853e-3,  1.880499e-1,  5.468481e-1,  3.163265e-1,
         7.845805e-2, -2.607710e-3, -1.609977e-3, -5.147392e-3, -1.393424e-2, -7.363946e-2, -9.455782e-2, -7.363946e-2,
        -1.393424e-2, -5.147392e-3, -2.269841e-2, -5.087302e-2, -2.154762e-1, -2.619048e-1, -2.154762e-1, -5.087302e-2,
        -2.269841e-2, -8.435374e-3,  1.945578e-2,  3.163265e-1,  4.653061e-1,  3.163265e-1,  1.945578e-2, -8.435374e-3,
        -5.045351e-4, -2.069161e-3, -3.049887e-2, -7.363946e-2, -1.040590e-1, -4.247732e-2, -2.675170e-2, -2.420635e-3,
        -8.928571e-3, -1.023810e-1, -2.154762e-1, -2.795238e-1, -1.228968e-1, -1.083730e-1, -1.609977e-3, -2.607710e-3,
         7.845805e-2,  3.163265e-1,  5.468481e-1,  1.880499e-1, -5.464853e-3, -2.069161e-3, -1.393424e-2, -4.247732e-2,
        -3.356009e-2, -4.795918e-2, -8.928571e-3, -5.087302e-2, -1.228968e-1, -5.714285e-2, -1.801587e-1, -2.607710e-3,
         1.945578e-2,  1.880499e-1,  3.083900e-1,  4.671201e-2, -5.045351e-4, -5.147392e-3, -2.675170e-2, -4.795918e-2,
        -5.963719e-2, -2.420635e-3, -2.269841e-2, -1.083730e-1, -1.801587e-1, -1.063492e-1, -1.609977e-3, -8.435374e-3,
        -5.464853e-3,  4.671201e-2,  5.287982e-1
    };

    for (size_t i = 0; i < 779; ++i)
    {
        CHECK(std::get<0>(computedLhs.dataStructure())[i] == expectedIndices[i]);
        CHECK(std::get<2>(computedLhs.dataStructure())[i] == Approx(expectedData[i]));
    }

    for (size_t i = 0; i < 35 + 1; ++i)
    {
        CHECK(std::get<1>(computedLhs.dataStructure())[i] == expectedIndptr[i]);
    }

    std::vector<double> expectedRhs
    {
         3.288002e-1, 3.855131e-1, 1.111304e+0, 1.565008e+0, 2.132137e+0, 1.406346e+0, 1.604842e+0, 4.955633e-1,
         6.089892e-1, 1.898534e+0, 2.805941e+0, 3.940201e+0, 2.650656e+0, 3.047647e+0, 5.610532e-1, 7.311921e-1,
         2.483218e+0, 3.844329e+0, 5.545718e+0, 3.793692e+0, 4.389178e+0, 3.740355e-1, 4.874614e-1, 1.655478e+0,
         2.562886e+0, 3.697145e+0, 2.529128e+0, 2.926119e+0, 2.275270e-1, 2.842400e-1, 9.087577e-1, 1.362461e+0,
         1.929591e+0, 1.305073e+0, 1.503569e+0
    };

    for (size_t i = 0; i < 35; ++i)
    {
        CHECK(computedRhs[i] == Approx(expectedRhs[i]));
    }

} // BSplineFiniteElementPatch_assembleGlobalSystem_test

TEST_CASE("BSplineFiniteElementPatch_solutionEvaluator_test")
{
    IntegrationPointProvider provider = [](size_t order)->std::array<std::vector<double>, 2>
    {
        throw std::runtime_error{ "Not implemented." };
    };

    BSplineFiniteElementPatch mesh({ 3, 2 }, { 2, 4 }, { 1, 2 }, { 2.5, 3.5 }, { -1.5, 0.5 }, provider);

    std::vector<double> solution
    {
        0.16394224, 0.31523237, 0.73912341, 0.23852614, 0.52263410,
        0.60238135, 0.19591390, 0.51422093, 0.74025031, 0.50942581,
        0.45894009, 0.50858743, 0.06583951, 0.03543730, 0.43150111,
        0.78337688, 0.25064280, 0.22298908, 0.68154774, 0.20123930,
        0.41964415, 0.01377733, 0.19983029, 0.20277283, 0.61445788,
        0.40233213, 0.15650152, 0.31876451, 0.14404782, 0.43722373,
        0.46045957, 0.28226385, 0.35518153, 0.08915129, 0.42997432
    };

    auto u = mesh.solutionEvaluator(solution);

    std::vector<double> x{ -1.4, -0.6, -0.7, 0.15, 0.2, 0.9 };
    std::vector<double> y{ 0.6, 2.2, 2.3, 3.7 };

    std::vector<std::vector<double>> expectedEvaluations
    {
        { 0.282649, 0.450238, 0.440854, 0.334868 },
        { 0.520361, 0.397263, 0.411312, 0.258750 },
        { 0.528379, 0.414361, 0.425578, 0.230033 },
        { 0.283088, 0.390137, 0.408334, 0.319641 },
        { 0.257205, 0.396974, 0.414373, 0.313951 },
        { 0.174583, 0.370885, 0.368873, 0.282877 }
    };

    for (size_t i = 0; i < x.size(); ++i)
    {
        for (size_t j = 0; j < y.size(); ++j)
        {
            CHECK(u(x[i], y[j]) == Approx(expectedEvaluations[i][j]));
        }
    }
}

} // namespace splinekernel
} // namespace cie
