#include "catch.hpp"
#include "sparse.hpp"

#include <vector>
#include <numeric>

namespace cie
{
namespace splinekernel
{

// Elements:       A     B     C
// Nodes:       0-----1-----2-----3
//                    |
//                    |  Element: D
//                    |
//                    4
// Matrix:
//     x, x, 0, 0, 0
//     x, x, x, 0, x 
//     0, x, x, x, 0
//     0, 0, x, x, 0
//     0, x, 0, 0, x
CompressedSparseRowMatrix simpleTestMatrix( )
{
    std::vector<LocationMap> locationMaps
    {
        { 0, 1 },     // Element A
        { 1, 2 },     // Element B
        { 2, 3 },     // Element C
        { 1, 4 }      // Element D
    };

    return CompressedSparseRowMatrix{ locationMaps };
}

TEST_CASE( "SparseMatrix_simpleSizeAndNnz_test" )
{
    auto matrix = simpleTestMatrix( );

    REQUIRE( matrix.size( ) == 5 );
    REQUIRE( matrix.nnz( ) == 13 );

    auto data = matrix.dataStructure( );

    std::vector<CompressedSparseRowMatrix::IndexType> expectedIndptr 
    { 
        0, 2, 6, 9, 11, 13 
    };

    std::vector<CompressedSparseRowMatrix::IndexType> expectedIndices 
    { 
        0, 1, 0, 1, 2, 4, 1, 2, 3, 2, 3, 1, 4 
    };
    
    for( size_t i = 0; i <= 5; ++i )
    {
        CHECK( std::get<1>( data )[i] == expectedIndptr[i] );
    }

    for( size_t i = 0; i < 13; ++i )
    {
        CHECK( std::get<0>( data )[i] == expectedIndices[i] );
        CHECK( std::get<2>( data )[i] == 0.0 );
    }
}

TEST_CASE("SparseMatrix_simpleElementAccess_test")
{
    auto matrix = simpleTestMatrix();

    double* data = std::get<2>(matrix.dataStructure());

    std::iota(data, data + 13, 0.5);

    std::vector<double> expectedValues
    {
        0.5,  1.5,  0.0,  0.0,  0.0,
        2.5,  3.5,  4.5,  0.0,  5.5,
        0.0,  6.5,  7.5,  8.5,  0.0,
        0.0,  0.0,  9.5, 10.5,  0.0,
        0.0, 11.5,  0.0,  0.0, 12.5
    };

    for (size_t i = 0; i < 5; ++i)
    {
        for (size_t j = 0; j < 5; ++j)
        {
            CHECK(matrix(i, j) == Approx(expectedValues[i * 5 + j]));
        }
    }

    REQUIRE_THROWS(matrix(5, 0));
    REQUIRE_THROWS(matrix(0, 5));
    REQUIRE_THROWS(matrix(5, 5));
}

TEST_CASE( "SparseMatrix_simpleMultiplication_test" )
{
    auto matrix = simpleTestMatrix( );

    double* data = std::get<2>( matrix.dataStructure( ) );

    std::iota( data, data + 13, 0.5 );

    std::vector<double> rhs { -1.4, 4.2, 0.7, -3.1, 1.5 };
    
    std::vector<double> expectedResult { 5.6 , 22.6, 6.2 , -25.9, 67.05 };
    std::vector<double> computedResult;

    REQUIRE_NOTHROW( computedResult = matrix * rhs );

    REQUIRE( computedResult.size( ) == 5 );

    for( size_t i = 0; i < 5; ++i )
    {
        CHECK( computedResult[i] == Approx( expectedResult[i] ) );
    }

    REQUIRE_THROWS( matrix * std::vector<double>{ 0.0, 0.0, 0.0, 0.0 } );
}

TEST_CASE( "SparseMatrix_simpleScatter_test" )
{
    auto matrix = simpleTestMatrix( );

    std::vector<LocationMap> locationMaps
    {
        { 0, 1 },     // Element A
        { 1, 2 },     // Element B
        { 2, 3 },     // Element C
        { 1, 4 }      // Element D
    };

    linalg::Matrix elementMatrixA( { { 1.0, -1.0 }, { -1.0, 1.0 } } );
    linalg::Matrix elementMatrixB( { { 2.0, -2.0 }, { -2.0, 2.0 } } );
    linalg::Matrix elementMatrixC( { { 3.0, -3.0 }, { -3.0, 3.0 } } );
    linalg::Matrix elementMatrixD( { { 4.0, -4.0 }, { -4.0, 4.0 } } );

    REQUIRE_NOTHROW( matrix.scatter( elementMatrixA, { 0, 1 } ) );
    REQUIRE_NOTHROW( matrix.scatter( elementMatrixB, { 1, 2 } ) );
    REQUIRE_NOTHROW( matrix.scatter( elementMatrixC, { 2, 3 } ) );
    REQUIRE_NOTHROW( matrix.scatter( elementMatrixD, { 1, 4 } ) );

    auto data = matrix.dataStructure( );

    // Shouldn't change
    std::vector<CompressedSparseRowMatrix::IndexType> expectedIndptr
    {
        0, 2, 6, 9, 11, 13
    };

    // Shouldn't change
    std::vector<CompressedSparseRowMatrix::IndexType> expectedIndices
    {
        0, 1, 0, 1, 2, 4, 1, 2, 3, 2, 3, 1, 4
    };

    std::vector<double> expectedData
    {
        1.0, -1.0, -1.0, 7.0, -2.0, -4.0, -2.0, 5.0, -3.0, -3.0, 3.0, -4.0, 4.0
    };

    for( size_t i = 0; i <= 5; ++i )
    {
        CHECK( std::get<1>( data )[i] == expectedIndptr[i] );
    }

    for( size_t i = 0; i < 13; ++i )
    {
        CHECK( std::get<0>( data )[i] == expectedIndices[i] );
        CHECK( std::get<2>( data )[i] == Approx( expectedData[i] ) );
    }

    REQUIRE_THROWS( matrix.scatter( linalg::Matrix( 2, 1 ), { 1, 2 } ) ); // Non-square matrix
    REQUIRE_THROWS( matrix.scatter( linalg::Matrix( 1, 2 ), { 1, 2 } ) ); // Non-square matrix
    REQUIRE_THROWS( matrix.scatter( elementMatrixB, { 1, 3, 4 } ) );      // Inconsistent size
    REQUIRE_THROWS( matrix.scatter( elementMatrixB, { 1, 3 } ) );         // Entry not in sparsity pattern
}


// Numbers represent nodes
// Letters represent elements
//
//  6---------7---------8
//  |         |         |
//  |    C    |    D    |
//  |         |         |
//  3---------4---------5
//  |         |         |
//  |    A    |    B    |
//  |         |         |
//  0---------1---------2
//
TEST_CASE( "SparseMatrix_combined_test1" )
{
    // Construct the location map as per the above illustration
    std::vector<LocationMap> locationMaps
    {
        { 0, 1, 4, 3 }, // Element A
        { 1, 2, 5, 4 }, // Element B
        { 3, 4, 7, 6 }, // Element C
        { 4, 5, 8, 7 }, // Element D
    };

    // Initialize the sparse matrix
    CompressedSparseRowMatrix sparseMatrix( locationMaps );

    // Initialize the element matrices of each element with arbitrary values
    linalg::Matrix elementMatrixA
    ( {
        {  8.0,  6.0, -8.0, -6.0 },
        {  6.0,  4.5, -6.0, -4.5 },
        { -8.0, -6.0,  8.0,  6.0 },
        { -6.0, -4.5,  6.0,  4.5 }
      } );

    linalg::Matrix elementMatrixB
    ( {
        {  1.5,  0.0, -1.5,  0.0 },
        {  0.0,  0.0,  0.0,  0.5 },
        { -1.5,  0.0,  1.5,  0.0 },
        {  0.0,  0.0,  0.0,  0.0 }
      } );

    linalg::Matrix elementMatrixC
    ( {
        { -2.5,  3.0,  2.5, -3.0 },
        {  0.0, -7.5, -0.0,  7.5 },
        {  2.5, -3.0, -2.5,  3.0 },
        { -0.0,  7.5,  0.0, -7.5 }
      } );

    linalg::Matrix elementMatrixD
    ( {
        {  1.0, -4.0, -1.0,  4.0 },
        { -5.5,  2.0,  5.5, -2.0 },
        { -1.0,  4.0,  1.0, -4.0 },
        {  5.5, -2.0, -5.5,  2.0 }
      } );

    sparseMatrix.scatter( elementMatrixA, locationMaps[0] );
    sparseMatrix.scatter( elementMatrixB, locationMaps[1] );
    sparseMatrix.scatter( elementMatrixC, locationMaps[2] );
    sparseMatrix.scatter( elementMatrixD, locationMaps[3] );

    REQUIRE( sparseMatrix.nnz( ) == 49 );
    REQUIRE( sparseMatrix.size( ) == 9 );

    auto data = sparseMatrix.dataStructure( );

    std::vector<CompressedSparseRowMatrix::IndexType> expectedIndices 
    { 
        0, 1, 3, 4, 0, 1, 2, 3, 4, 5, 1, 2, 4, 5, 0, 1, 3, 
        4, 6, 7, 0, 1, 2, 3, 4, 5, 6, 7, 8, 1, 2, 4, 5, 7, 
        8, 3, 4, 6, 7, 3, 4, 5, 6, 7, 8, 4, 5, 7, 8 
    };

    std::vector<CompressedSparseRowMatrix::IndexType> expectedIndptr
    {
        0, 4, 10, 14, 20, 29, 35, 39, 45, 49
    };

    std::vector<double> expectedData
    { 
        8.0,  6.0, -6.0, -8.0,  6.0,  6.0,  0.0, -4.5, -6.0, -1.5,  0.0,   
        0.0,  0.5,  0.0, -6.0, -4.5,  2.0,  9.0, -3.0,  2.5, -8.0, -6.0,  
        0.0,  6.0,  1.5, -4.0,  7.5,  4.0, -1.0, -1.5,  0.0, -5.5,  3.5, 
       -2.0,  5.5,  0.0,  7.5, -7.5,  0.0,  2.5,  2.5, -2.0, 3.0, -0.5, 
       -5.5, -1.0,  4.0, -4.0,  1.0 
    };
    
    for( size_t i = 0; i <= sparseMatrix.size( ); ++i )
    {
        CHECK( std::get<1>( data )[i] == expectedIndptr[i] );
    }

    for( size_t i = 0; i < sparseMatrix.nnz( ); ++i )
    {
        CHECK( std::get<0>( data )[i] == expectedIndices[i] );
        CHECK( std::get<2>( data )[i] == Approx( expectedData[i] ).margin( 1e-10 ) );
    }

    // Expected dense global matrix
    linalg::Matrix expectedGlobalMatrix
    ( {
        {  8.0,  6.0,  0.0, -6.0, -8.0,  0.0,  0.0,  0.0,  0.0 },
        {  6.0,  6.0,  0.0, -4.5, -6.0, -1.5,  0.0,  0.0,  0.0 },
        {  0.0,  0.0,  0.0,  0.0,  0.5,  0.0,  0.0,  0.0,  0.0 },
        { -6.0, -4.5,  0.0,  2.0,  9.0,  0.0, -3.0,  2.5,  0.0 },
        { -8.0, -6.0,  0.0,  6.0,  1.5, -4.0,  7.5,  4.0, -1.0 },
        {  0.0, -1.5,  0.0,  0.0, -5.5,  3.5,  0.0, -2.0,  5.5 },
        {  0.0,  0.0,  0.0,  0.0,  7.5,  0.0, -7.5,  0.0,  0.0 },
        {  0.0,  0.0,  0.0,  2.5,  2.5, -2.0,  3.0, -0.5, -5.5 },
        {  0.0,  0.0,  0.0,  0.0, -1.0,  4.0,  0.0, -4.0,  1.0 },
    } );

    // Verify the data is inserted into the correct positions
    for( size_t iRow = 0; iRow < sparseMatrix.size( ); ++iRow )
    {
        for( size_t iColumn = 0; iColumn < sparseMatrix.size( ); ++iColumn )
        {
            CHECK( sparseMatrix( iRow, iColumn ) == Approx( expectedGlobalMatrix( iRow, iColumn ) ).margin( 1e-12 ) );
        }
    }

    std::vector<double> rhs { 4.3, 2.2, 3.7, -6.5, 4.2, -0.9, 5.2, -1.3, 2.4 };
    
    std::vector<double> expectedResult { 53.0, 44.4, 2.1, -29.75, -45.3, -13.75, -7.5 , -0.9, -0.2 };
    std::vector<double> computedResult;
    
    REQUIRE_NOTHROW( computedResult = sparseMatrix * rhs );
 
    REQUIRE( computedResult.size( ) == 9 );

    for( size_t i = 0; i < 9; ++i )
    {
        CHECK( computedResult[i] == Approx( expectedResult[i] ) );
    }

} // SparseMatrix_combined_test1


//  ------------ Test quad connectivity -------------
//
//     33----34-----35----36-----37-----------39
//     |            |            |            |
//     |            26  28   30  32           |
//     23    24     |            |            |
//     |            25  27   29  31           |
//     |            |            |            |
//     17----18-----19--20---21--22-----------38
//     |            |            |
//     7   9   11   13    15     |
//     |            |            16
//     6   8   10   12    14     |
//     |            |            |
//     0---1----2---3------4-----5
//
//
// test location maps:
//   [[ 0, 1, 2, 3, 6, 8, 10, 12, 7, 9, 11, 13, 17, 18, 19 ]
//    [ 3, 12, 13, 19, 4, 14, 15, 21, 5, 16, 22, 20 ]
//    [ 17, 23, 33, 34, 35, 26, 25, 19, 18, 24 ]
//    [ 38, 22, 39, 37, 31, 32 ]
//    [ 27, 29, 28, 30, 37, 36, 35, 26, 32, 31, 25, 22, 19, 21, 20 ]]
//
TEST_CASE( "SparseMatrix_combined_test2" )
{
    std::vector<LocationMap> locationMaps
    {
        { 0, 1, 2, 3, 6, 8, 10, 12, 7, 9, 11, 13, 17, 18, 19 },
        { 3, 12, 13, 19, 4, 14, 15, 21, 5, 16, 22, 20 },
        { 17, 23, 33, 34, 35, 26, 25, 19, 18, 24 },
        { 38, 22, 39, 37, 31, 32 },
        { 27, 29, 28, 30, 37, 36, 35, 26, 32, 31, 25, 22, 19, 21, 20 }
    };

    CompressedSparseRowMatrix sparseMatrix( locationMaps );

    REQUIRE( sparseMatrix.size( ) == 40 );

    auto sparseDataStructure = sparseMatrix.dataStructure( );

    auto* indices = std::get<0>( sparseDataStructure );
    auto* indptr = std::get<1>( sparseDataStructure );
    double* data = std::get<2>( sparseDataStructure );

    for( auto& locationMap : locationMaps )
    {
        std::sort( locationMap.begin( ), locationMap.end( ) );
    }

    std::vector<std::vector<size_t>> dofElementMapping =
    {
        { 0 }   , { 0 }   , { 0 }      , { 0, 1 }, { 1 },          // dofs  0 -  4
        { 1 }   , { 0 }   , { 0 }      , { 0 }   , { 0 },          // dofs  5 -  9
        { 0 }   , { 0 }   , { 0, 1 }   , { 0, 1 }, { 1 },          // dofs 10 - 14
        { 1 }   , { 1 }   , { 0, 2 }   , { 0, 2 }, { 0, 1, 2, 4 }, // dofs 15 - 19
        { 1, 4 }, { 1, 4 }, { 1, 3, 4 }, { 2 }   , { 2 },          // dofs 20 - 24
        { 2, 4 }, { 2, 4 }, { 4 }      , { 4 }   , { 4 },          // dofs 25 - 29
        { 4 }   , { 3, 4 }, { 3, 4 }   , { 2 }   , { 2 },          // dofs 30 - 34
        { 2, 4 }, { 4 }   , { 3, 4 }   , { 3 }   , { 3 }           // dofs 35 - 39
    };

    std::vector<size_t> connectingDofSizes;
    std::vector<size_t> expectedIndices;

    for( const std::vector<size_t>& indices : dofElementMapping )
    {
        std::vector<size_t> concatenated;

        for( const auto& elementIndex : indices )
        {
            const auto& locationMap = locationMaps[elementIndex];

            concatenated.insert( concatenated.end( ), 
                                 locationMap.begin( ), 
                                 locationMap.end( ) );
        }

        std::sort( concatenated.begin( ), concatenated.end( ) );

        concatenated.erase( std::unique( concatenated.begin( ), 
                                         concatenated.end( ) ), 
                            concatenated.end( ) );

        connectingDofSizes.push_back( concatenated.size( ) );

        expectedIndices.insert( expectedIndices.end( ), 
                                concatenated.begin( ), 
                                concatenated.end( ) );
    };

    size_t nnz = expectedIndices.size( );

    REQUIRE( indptr[40] == nnz );

    for( size_t iEntry = 0; iEntry < nnz; ++iEntry )
    {
        CHECK( indices[iEntry] == expectedIndices[iEntry] );
        CHECK( data[iEntry] == 0.0 );
    }

    size_t rowPtr = 0;

    for( size_t iDof = 0; iDof < 40; ++iDof )
    {
        CHECK( indptr[iDof] == rowPtr );

        rowPtr += connectingDofSizes[iDof];
    }

    CHECK( indptr[40] == rowPtr );

    std::vector<double> data0( 15 * 15 );
    std::vector<double> data1( 12 * 12 );
    std::vector<double> data2( 10 * 10 );
    std::vector<double> data3(  6 *  6 );
    std::vector<double> data4( 15 * 15 );

    std::iota( data0.begin( ), data0.end( ),  4.7 );
    std::iota( data1.begin( ), data1.end( ), -0.8 );
    std::iota( data2.begin( ), data2.end( ), 13.3 );
    std::iota( data3.begin( ), data3.end( ),  3.2 );
    std::iota( data4.begin( ), data4.end( ), -7.4 );

    linalg::Matrix elementMatrix0( data0, 15 );
    linalg::Matrix elementMatrix1( data1, 12 );
    linalg::Matrix elementMatrix2( data2, 10 );
    linalg::Matrix elementMatrix3( data3,  6 );
    linalg::Matrix elementMatrix4( data4, 15 );

    REQUIRE_NOTHROW( sparseMatrix.scatter( elementMatrix0, locationMaps[0] ) );
    REQUIRE_NOTHROW( sparseMatrix.scatter( elementMatrix1, locationMaps[1] ) );
    REQUIRE_NOTHROW( sparseMatrix.scatter( elementMatrix2, locationMaps[2] ) );
    REQUIRE_NOTHROW( sparseMatrix.scatter( elementMatrix3, locationMaps[3] ) );
    REQUIRE_NOTHROW( sparseMatrix.scatter( elementMatrix4, locationMaps[4] ) );

    REQUIRE( sparseMatrix.nnz( ) == nnz );
    REQUIRE( sparseMatrix.size( ) == 40 );

    CHECK( sparseMatrix( 3, 5 ) == Approx( 1.2 ) );
    CHECK( sparseMatrix( 27, 35 ) == Approx( 94.6 ) );
    CHECK( sparseMatrix( 10, 13 ) == Approx( 135.7 ) );
    CHECK( sparseMatrix( 31, 24 ) == Approx( 0.0 ).margin( 1e-12 ) );
    CHECK( sparseMatrix( 5, 17 ) == Approx( 0.0 ).margin( 1e-12 ) );

    std::vector<double> rhs( 40 );

    std::iota( rhs.begin( ), rhs.end( ), 5.2 );

    std::vector<double> expectedResult
    {
        2881.8,   6091.8,  9301.8, 13831.9, 4036.88, 6753.68, 15721.8, 18931.8,
        22141.8, 25351.8, 28561.8, 31771.8, 44452.3, 50379.1, 14904.1, 17620.9,
        20337.7, 47025.6, 53295.6, 82781.7, 33372.9, 43529.7, 55054.1, 14803.8,
        17863.8, 50845.4, 61345.4, 44801.6, 52241.6, 59681.6, 67121.6, 77310.4,
        86131.6, 27043.8, 30103.8,  122605, 96881.6, 109833,  6892.44, 8273.64
    };

    std::vector<double> computedResult;
    
    REQUIRE_NOTHROW( computedResult = sparseMatrix * rhs );

    for( size_t i = 0; i < 40; ++i )
    {
        CHECK( computedResult[i] == Approx( expectedResult[i] ) );
    }

} // SparseMatrix_combined_test2

} // namespace cie
} // namespace splinekernel
