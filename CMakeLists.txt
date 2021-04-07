cmake_minimum_required( VERSION 3.4.3 )

project( splinecomputation CXX )

# Directories to include header files from
include_directories( external/catch )
include_directories( external/linalg/inc )

set( CMAKE_CXX_STANDARD 14 )

# If compiler is g++: Enable further warnings and treat all warnings as errors. fPIC stands for position independent code.
if( CMAKE_COMPILER_IS_GNUCXX )
    set( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pedantic -Wall -Werror -fPIC" )
endif( CMAKE_COMPILER_IS_GNUCXX )

# -------------- Set up external linalg project ------------------
add_subdirectory( external/linalg )

# ----------------- Set up install paths -------------------------
# Changes the default install path to build/install, assuming build is the project directory
if ( CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT )
    set ( CMAKE_INSTALL_PREFIX "${CMAKE_CURRENT_BINARY_DIR}/install" CACHE PATH "default install path" FORCE )
endif( )

# This will add the install path to the rpath of the installed binaries to make dynamic linking work
SET( CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}" )

# ---------------- Set up splinekernel project --------------------

add_subdirectory( splinekernel )

# ----------------- Set up python bindings ------------------------
set(PYBIND11_CPP_STANDARD -std=c++14)

# Add pybind11 project (this is just the folder available at https://github.com/pybind/pybind11)
add_subdirectory( external/pybind11 )

# Add our python project (containing wrapping code and python scripts)
add_subdirectory( python )
