cmake_minimum_required(VERSION 3.15)

if(NOT DEFINED SKBUILD_PROJECT_VERSION)
  set(SKBUILD_PROJECT_VERSION "0.0.0")
endif()

project(
  fesslix
  VERSION ${SKBUILD_PROJECT_VERSION}
  LANGUAGES CXX)

set( FLX_VERSIONb ${PROJECT_VERSION} )

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

set(CMAKE_TOOLCHAIN_FILE "$ENV{CMAKE_TOOLCHAIN_FILE}")
set(VCPKG_ROOT "$ENV{VCPKG_ROOT}")
set(PKG_CONFIG_PATH "$ENV{PKG_CONFIG_PATH}")


###################################################################################
# Configuration Options
###################################################################################

option ( FLX_DEBUG "with debugging stuff .. but maybe more stable" OFF )
option ( FLX_DEBUG_COUT "use std::cout for debugging in your code" OFF )

if("${FLX_DEBUG}" STREQUAL "ON")
set ( FLX_DEBUGb 1 )
else("${FLX_DEBUG}" STREQUAL "ON")
set ( FLX_DEBUGb 0 )
endif("${FLX_DEBUG}" STREQUAL "ON")
if("${FLX_DEBUG_COUT}" STREQUAL "ON")
set ( FLX_DEBUG_COUTb 1 )
else("${FLX_DEBUG_COUT}" STREQUAL "ON")
set ( FLX_DEBUG_COUTb 0 )
endif("${FLX_DEBUG_COUT}" STREQUAL "ON")

option ( FLX_KAHAN_MCI "use Kahan summation for Monte Carlo Integration" ON )
  if("${FLX_KAHAN_MCI}" STREQUAL "ON")
    set ( FLX_KAHAN_MCIb 1 )
  else("${FLX_KAHAN_MCI}" STREQUAL "ON")
    set ( FLX_KAHAN_MCIb 0 )
  endif("${FLX_KAHAN_MCI}" STREQUAL "ON")
option ( FLX_KAHAN_KDE "use Kahan summation for kernel density estimation" ON )
  if("${FLX_KAHAN_KDE}" STREQUAL "ON")
    set ( FLX_KAHAN_KDEb 1 )
  else("${FLX_KAHAN_KDE}" STREQUAL "ON")
    set ( FLX_KAHAN_KDEb 0 )
  endif("${FLX_KAHAN_KDE}" STREQUAL "ON")
option ( FLX_KAHAN_DOT "use Kahan summation for dot-products" OFF )
  if("${FLX_KAHAN_DOT}" STREQUAL "ON")
    set ( FLX_KAHAN_DOTb 1 )
  else("${FLX_KAHAN_DOT}" STREQUAL "ON")
    set ( FLX_KAHAN_DOTb 0 )
  endif("${FLX_KAHAN_DOT}" STREQUAL "ON")
option ( FLX_KAHAN_2NORM "use Kahan summation for calculating the L2-norm" OFF )
  if("${FLX_KAHAN_2NORM}" STREQUAL "ON")
    set ( FLX_KAHAN_2NORMb 1 )
  else("${FLX_KAHAN_2NORM}" STREQUAL "ON")
    set ( FLX_KAHAN_2NORMb 0 )
  endif("${FLX_KAHAN_2NORM}" STREQUAL "ON")
option ( FLX_KAHAN_MTX_FULL "use Kahan summation for internal operations in general full matrices" OFF )
  if("${FLX_KAHAN_MTX_FULL}" STREQUAL "ON")
    set ( FLX_KAHAN_MTX_FULLb 1 )
  else("${FLX_KAHAN_MTX_FULL}" STREQUAL "ON")
    set ( FLX_KAHAN_MTX_FULLb 0 )
  endif("${FLX_KAHAN_MTX_FULL}" STREQUAL "ON")
option ( FLX_KAHAN_MTX_LTRI "use Kahan summation for internal operations in lower triangular matrices" OFF )
  if("${FLX_KAHAN_MTX_LTRI}" STREQUAL "ON")
    set ( FLX_KAHAN_MTX_LTRIb 1 )
  else("${FLX_KAHAN_MTX_LTRI}" STREQUAL "ON")
    set ( FLX_KAHAN_MTX_LTRIb 0 )
  endif("${FLX_KAHAN_MTX_LTRI}" STREQUAL "ON")
option ( FLX_KAHAN_MTX_SLTRI "use Kahan summation for internal operations in sparse lower triangular matrices" OFF )
  if("${FLX_KAHAN_MTX_SLTRI}" STREQUAL "ON")
    set ( FLX_KAHAN_MTX_SLTRIb 1 )
  else("${FLX_KAHAN_MTX_SLTRI}" STREQUAL "ON")
    set ( FLX_KAHAN_MTX_SLTRIb 0 )
  endif("${FLX_KAHAN_MTX_SLTRI}" STREQUAL "ON")
option ( FLX_KAHAN_MTX_SSILU "use Kahan summation for internal operations in sparse matrices of incomplete LU decompositions" OFF )
  if("${FLX_KAHAN_MTX_SSILU}" STREQUAL "ON")
    set ( FLX_KAHAN_MTX_SSILUb 1 )
  else("${FLX_KAHAN_MTX_SSILU}" STREQUAL "ON")
    set ( FLX_KAHAN_MTX_SSILUb 0 )
  endif("${FLX_KAHAN_MTX_SSILU}" STREQUAL "ON")
option ( FLX_KAHAN_MTX_SYM "use Kahan summation for internal operations in full symmetric matrices" OFF )
  if("${FLX_KAHAN_MTX_SYM}" STREQUAL "ON")
    set ( FLX_KAHAN_MTX_SYMb 1 )
  else("${FLX_KAHAN_MTX_SYM}" STREQUAL "ON")
    set ( FLX_KAHAN_MTX_SYMb 0 )
  endif("${FLX_KAHAN_MTX_SYM}" STREQUAL "ON")
option ( FLX_KAHAN_MTX_SSYM "use Kahan summation for internal operations in sparse symmetric matrices" OFF )
  if("${FLX_KAHAN_MTX_SSYM}" STREQUAL "ON")
    set ( FLX_KAHAN_MTX_SSYMb 1 )
  else("${FLX_KAHAN_MTX_SSYM}" STREQUAL "ON")
    set ( FLX_KAHAN_MTX_SSYMb 0 )
  endif("${FLX_KAHAN_MTX_SSYM}" STREQUAL "ON")
option ( FLX_KAHAN_MTX_SSFEM "use Kahan summation for internal operations in SSFEM matrices" OFF )
  if("${FLX_KAHAN_MTX_SSFEM}" STREQUAL "ON")
    set ( FLX_KAHAN_MTX_SSFEMb 1 )
  else("${FLX_KAHAN_MTX_SSFEM}" STREQUAL "ON")
    set ( FLX_KAHAN_MTX_SSFEMb 0 )
  endif("${FLX_KAHAN_MTX_SSFEM}" STREQUAL "ON")
option ( FLX_KAHAN_PFEM "use Kahan summation for pFEM-elements" OFF )
  if("${FLX_KAHAN_PFEM}" STREQUAL "ON")
    set ( FLX_KAHAN_PFEMb 1 )
  else("${FLX_KAHAN_PFEM}" STREQUAL "ON")
    set ( FLX_KAHAN_PFEMb 0 )
  endif("${FLX_KAHAN_PFEM}" STREQUAL "ON")

SET( FLX_KAHAN_CG_STAGE 1 CACHE INTEGER "specify how to use Kahan summation in the conjugate gradient method")

configure_file (
  "${PROJECT_SOURCE_DIR}/src/config.h.in"
  "${PROJECT_BINARY_DIR}/src/config.h"
)
include_directories("${PROJECT_BINARY_DIR}/src")

MARK_AS_ADVANCED(
  FLX_DEBUG_COUT
  FLX_KAHAN_MCI
  FLX_KAHAN_KDE
  FLX_KAHAN_DOT
  FLX_KAHAN_2NORM
  FLX_KAHAN_MTX_FULL
  FLX_KAHAN_MTX_LTRI
  FLX_KAHAN_MTX_SLTRI
  FLX_KAHAN_MTX_SSILU
  FLX_KAHAN_MTX_SYM
  FLX_KAHAN_MTX_SSYM
  FLX_KAHAN_MTX_SSFEM
  FLX_KAHAN_PFEM
  FLX_KAHAN_CG_STAGE
)


###################################################################################
# find dependencies
###################################################################################

find_package(Python REQUIRED COMPONENTS Interpreter Development.Module)
find_package(pybind11 CONFIG REQUIRED)


if(DEFINED ENV{GSL_INCLUDE_DIR})
    set(GSL_INCLUDE_DIR "$ENV{GSL_INCLUDE_DIR}")
    set(GSL_LIBRARY "$ENV{GSL_LIBRARY}")
    set(GSL_CBLAS_LIBRARY "$ENV{GSL_CBLAS_LIBRARY}")
endif()
FIND_PACKAGE(GSL REQUIRED)

if(POLICY CMP0167)
  cmake_policy(SET CMP0167 NEW)
endif()
if(DEFINED ENV{Boost_INCLUDE_DIR})
  set(Boost_INCLUDE_DIR $ENV{Boost_INCLUDE_DIR})
endif()
if(DEFINED ENV{Boost_DIR})
  set(Boost_DIR $ENV{Boost_DIR})
endif()
#find_package(Boost 1.48.0 REQUIRED NO_MODULE)
if(Boost_INCLUDE_DIR)
  include_directories(${Boost_INCLUDE_DIR})
endif()
if(Boost_INCLUDE_DIRS)
  include_directories(${Boost_INCLUDE_DIRS})
endif()

message("202502240926 Boost_INCLUDE_DIR=${Boost_INCLUDE_DIR}")
message("202502240930 Boost_INCLUDE_DIRS=${Boost_INCLUDE_DIRS}")

###################################################################################
# compile
###################################################################################

# Add a library using FindPython's tooling (pybind11 also provides a helper like this)
python_add_library(fesslix MODULE
  src/flxexception.cpp
  src/pdouble.cpp
  src/flxglobal.cpp
  src/flxmath_rnd.cpp
  src/flxmath_interp.cpp
  src/flxmath.cpp
  src/flxVec.cpp
  src/flxio.cpp
  src/python_interface.cpp
  WITH_SOABI)
target_link_libraries(fesslix PRIVATE pybind11::headers PRIVATE GSL::gsl)

target_compile_definitions(fesslix PRIVATE VERSION_INFO=${PROJECT_VERSION})

# The install directory is the output (wheel) directory
install(TARGETS fesslix DESTINATION fesslix)





#find_package(pybind11 REQUIRED)

#find_package (Python3 COMPONENTS Interpreter Development REQUIRED)   # Development
#find_package(Python COMPONENTS Interpreter Development REQUIRED)


#pybind11_add_module(fesslix src/fesslix.cpp)

