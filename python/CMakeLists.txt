cmake_minimum_required( VERSION 3.4.3 )

project( python CXX )

# ------------ Set up library for C++ bindings -------------------

set( PYTHON_BINDINGS_FOLDER bindings )

# find all C++ files in bindings folder
file( GLOB BINDING_SOURCES ${PYTHON_BINDINGS_FOLDER}/*.cpp )
file( GLOB BINDING_HEADERS ${PYTHON_BINDINGS_FOLDER}/*.h* )

# Add shared library of python bindings
pybind11_add_module( pysplinekernel ${BINDING_SOURCES} ${BINDING_HEADERS} )

target_link_libraries( pysplinekernel PRIVATE splinekernel linalg )
#target_include_directories( pysplinekernel ../splinekernel/inc )

# specify the relative path the shared library object shall be installed to
install( TARGETS pysplinekernel LIBRARY DESTINATION ${CMAKE_INSTALL_PREFIX} )
  
# -------------- Set up copying python scipts ---------------------
set( PYTHON_SCRIPT_FOLDER scripts )

# find all python files in script folder
file( GLOB PYTHON_SCRIPTS ${PYTHON_SCRIPT_FOLDER}/*.py )

# add all files to install target
foreach( filename ${PYTHON_SCRIPTS} )
  install(FILES ${filename} DESTINATION ${CMAKE_INSTALL_PREFIX} )
endforeach()
