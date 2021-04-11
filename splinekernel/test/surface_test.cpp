#include "catch.hpp"
#include "surface.hpp"

#include <array>
#include <vector>

namespace cie
{
namespace splinekernel
{

TEST_CASE( "Linear interpolation surface" )
{
    std::vector<double> knotVectorR{ 0.0, 0.0, 0.5, 1.0, 1.0 };
    std::vector<double> knotVectorS{ 0.0, 0.0, 0.5, 1.0, 1.0 };
    
    std::array<std::vector<double>, 2> knotVectors{ knotVectorR, knotVectorS };

    size_t numberOfSamplesR = 11;
    size_t numberOfSamplesS = 9;

    linalg::Matrix xGrid( { -1.0, -1.0, -1.0,
                             0.0,  0.0,  0.0,
                             1.0,  1.0,  1.0 }, 3 );
                             
    linalg::Matrix yGrid( { -1.0,  0.0,  1.0,
                            -1.0,  0.0,  1.0,
                            -1.0,  0.0,  1.0 }, 3 );
                            
    linalg::Matrix zGrid( {  1.0,  1.0,  1.0,
                             1.0,  2.0,  1.0,
                             1.0,  1.0,  1.0 }, 3 );
                             
    VectorOfMatrices controlGrid{ xGrid, yGrid, zGrid };

    // Evaluate
    VectorOfMatrices C;
    REQUIRE_NOTHROW( C = evaluateSurface( knotVectors, controlGrid, 
                     { numberOfSamplesR, numberOfSamplesS } ) );

    // Check sizes
    REQUIRE(C.size() == 3);
    REQUIRE(C[0].size1() == numberOfSamplesR);
    REQUIRE(C[1].size1() == numberOfSamplesR);
    REQUIRE(C[2].size1() == numberOfSamplesR);
    REQUIRE(C[0].size2() == numberOfSamplesS);
    REQUIRE(C[1].size2() == numberOfSamplesS);
    REQUIRE(C[2].size2() == numberOfSamplesS);

    for( size_t r = 0; r < numberOfSamplesR; ++r ) 
    {
        for( size_t s = 0; s < numberOfSamplesS; ++s ) 
        {
            double X = 2.0 * r / ( numberOfSamplesR - 1.0 ) - 1.0;
            double Y = 2.0 * s / ( numberOfSamplesS - 1.0 ) - 1.0;
            double Z = 1 + ( std::abs(X) - 1 ) * ( std::abs(Y) - 1 );

            CHECK( C[0](r, s) == Approx(X) );
            CHECK( C[1](r, s) == Approx(Y) );
            CHECK( C[2](r, s) == Approx(Z) );

            //std::cout << Z << ",\t";        // Uncomment to print correct values
        }
        
        //std::cout << "\n";                // Uncomment to print correct values
    }

} // Linear interpolation surface

TEST_CASE( "Cubic-linear interpolation surface" ) 
{
    std::vector<double> knotVectorR{ 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0 };
    std::vector<double> knotVectorS{ 0.0, 0.0, 0.5, 1.0, 1.0 };
    
    std::array<std::vector<double>, 2> knotVectors{ knotVectorR, knotVectorS };

    size_t numberOfSamplesR = 7;
    size_t numberOfSamplesS = 5;

    linalg::Matrix xGrid( { -3.0, -3.0, -3.0,
                            -1.0, -1.0, -1.0,
                             1.0,  1.0,  1.0,
                             3.0,  3.0,  3.0 }, 4 );
                             
    linalg::Matrix yGrid( { -1.0,  0.0,  1.0,
                            -1.0,  0.0,  1.0,
                            -1.0,  0.0,  1.0,
                            -1.0,  0.0,  1.0 }, 4 );
                            
    linalg::Matrix zGrid( {  1.0,  1.0,  1.0,
                             1.0, 49.0,  1.0,
                             1.0, 49.0,  1.0,
                             1.0,  1.0,  1.0 }, 4 );
                             
    VectorOfMatrices controlGrid{ xGrid, yGrid, zGrid };

    // Evaluate
    VectorOfMatrices C;
    
    REQUIRE_NOTHROW( C = evaluateSurface( knotVectors, controlGrid, 
                     { numberOfSamplesR, numberOfSamplesS } ) );

    // Check sizes
    REQUIRE( C.size( ) == 3 );
    REQUIRE( C[0].size1( ) == numberOfSamplesR );
    REQUIRE( C[1].size1( ) == numberOfSamplesR );
    REQUIRE( C[2].size1( ) == numberOfSamplesR );
    REQUIRE( C[0].size2( ) == numberOfSamplesS );
    REQUIRE( C[1].size2( ) == numberOfSamplesS );
    REQUIRE( C[2].size2( ) == numberOfSamplesS );
    
    double increment = 0.0;

    for( size_t r = 0; r < numberOfSamplesR; ++r ) 
    {
        for( size_t s = 0; s < numberOfSamplesS; ++s ) 
        {
            double X = 6.0 * r / ( numberOfSamplesR - 1.0 ) - 3.0;
            double Y = 2.0 * s / ( numberOfSamplesS - 1.0 ) - 1.0;
            double Z = 1.0 + ( 2.0 - std::abs( s - 2.0 ) ) * increment;

            CHECK( C[0](r, s) == Approx(X) );
            CHECK( C[1](r, s) == Approx(Y) );
            CHECK( C[2](r, s) == Approx(Z) );

            //std::cout << Z << ",\t";                // Uncomment to print correct values
        }
        
        increment += 10 - 4.0 * r;
        
        //std::cout << "\n";                        // Uncomment to print correct values
    }

} // TEST_CASE("Cubic-linear interpolation surface")

} // namespace splinekernel
} // namespace cie