# Fesslix - Stochastic Analysis
# Copyright (C) 2010-2025 Wolfgang Betz
#
# Fesslix is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# Fesslix is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with Fesslix.  If not, see <http://www.gnu.org/licenses/>.

cmake_minimum_required(VERSION 3.15)
set(CMAKE_VERBOSE_MAKEFILE ON)

if(NOT DEFINED SKBUILD_PROJECT_VERSION)
  set(SKBUILD_PROJECT_VERSION "0.0.0")
endif()

project(
  fesslix
  VERSION ${SKBUILD_PROJECT_VERSION}
  LANGUAGES CXX)

set( FLX_VERSIONb ${PROJECT_VERSION} )

set(CMAKE_CXX_STANDARD 20)  ## does not compile with manylinux2014
#set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

set(CMAKE_TOOLCHAIN_FILE "$ENV{CMAKE_TOOLCHAIN_FILE}")
set(VCPKG_ROOT "$ENV{VCPKG_ROOT}")
if(DEFINED ENV{TRIPLET})
  set(PKG_CONFIG_PATH "$ENV{VCPKG_ROOT}./installed/$ENV{TRIPLET}/lib/pkgconfig")
else()
  set(PKG_CONFIG_PATH "$ENV{PKG_CONFIG_PATH}")
endif()
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/Modules/")

# add the automatically determined parts of the RPATH
# which point to directories outside the build tree to the install RPATH
SET(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

###################################################################################
# Configuration Options
###################################################################################

if (NOT CMAKE_BUILD_TYPE)
  set (CMAKE_BUILD_TYPE Release CACHE STRING
      "Choose the type of build, options are: None Debug Release."
      FORCE)
endif (NOT CMAKE_BUILD_TYPE)

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

#==========================================================
# pybind11
#==========================================================
find_package(Python REQUIRED COMPONENTS Interpreter Development.Module Development.Embed) #Development.Module
find_package(pybind11 CONFIG REQUIRED)

#==========================================================
# GSL
#==========================================================
if(DEFINED ENV{GSL_INCLUDE_DIR})
    set(GSL_INCLUDE_DIR "$ENV{GSL_INCLUDE_DIR}")
    set(GSL_LIBRARY "$ENV{GSL_LIBRARY}")
    set(GSL_CBLAS_LIBRARY "$ENV{GSL_CBLAS_LIBRARY}")
    set(GSL_STATIC_LIB "$ENV{GSL_LIBRARY}")
    set(GSLCBLAS_STATIC_LIB "$ENV{GSL_CBLAS_LIBRARY}")
else()
    if(DEFINED ENV{TRIPLET})
      set(GSL_INCLUDE_DIR "$ENV{VCPKG_ROOT}./installed/$ENV{TRIPLET}/include")
      if(MSVC)
        set(GSL_LIBRARY "$ENV{VCPKG_ROOT}./installed/$ENV{TRIPLET}/lib/gsl.lib")
        set(GSL_CBLAS_LIBRARY "$ENV{VCPKG_ROOT}./installed/$ENV{TRIPLET}/lib/gslcblas.lib")
      else()
        set(GSL_LIBRARY "$ENV{VCPKG_ROOT}./installed/$ENV{TRIPLET}/lib/libgsl.a")
        set(GSL_CBLAS_LIBRARY "$ENV{VCPKG_ROOT}./installed/$ENV{TRIPLET}/lib/libgslcblas.a")
      endif()
      set(GSL_STATIC_LIB "${GSL_LIBRARY}")
      set(GSLCBLAS_STATIC_LIB "${GSL_CBLAS_LIBRARY}")
    else()
      set(CMAKE_FIND_LIBRARY_SUFFIXES ".a")
      find_library(GSL_STATIC_LIB gsl)
      find_library(GSLCBLAS_STATIC_LIB gslcblas)
      set(CMAKE_FIND_LIBRARY_SUFFIXES ".dylib" ".so" ".a")
    endif()
endif()
if(GSL_STATIC_LIB AND GSLCBLAS_STATIC_LIB)
    message(STATUS "Found static GSL: ${GSL_STATIC_LIB} and ${GSLCBLAS_STATIC_LIB}")
    set(GSL_FOUND_STATIC TRUE)
    if(DEFINED GSL_INCLUDE_DIR)
        include_directories(${GSL_INCLUDE_DIR})
    endif()
else()
    message(STATUS "WARNING: NOT found static GSL")
    FIND_PACKAGE(GSL REQUIRED)
endif()
option ( FLX_USE_GSL "compile with GSL support" ON )

if("${FLX_USE_GSL}" STREQUAL "ON")
  set ( FLX_USE_GSLb 1 )
else("${FLX_USE_GSL}" STREQUAL "ON")
  set ( FLX_USE_GSLb 0 )
endif("${FLX_USE_GSL}" STREQUAL "ON")

MARK_AS_ADVANCED(
  FLX_USE_GSL
)

#==========================================================
# Boost
#==========================================================
if(POLICY CMP0167)
  cmake_policy(SET CMP0167 NEW)
endif()
if(DEFINED ENV{TRIPLET})
    set(Boost_INCLUDE_DIR "$ENV{VCPKG_ROOT}./installed/$ENV{TRIPLET}/include")
    set(Boost_DIR "$ENV{VCPKG_ROOT}./installed/$ENV{TRIPLET}/share/boost")
else()
  if(DEFINED ENV{Boost_INCLUDE_DIR})
    set(Boost_INCLUDE_DIR $ENV{Boost_INCLUDE_DIR})
  endif()
  if(DEFINED ENV{Boost_DIR})
    set(Boost_DIR $ENV{Boost_DIR})
  endif()
endif()
if(NOT DEFINED ENV{Boost_skip_search})
  find_package(Boost 1.48.0 REQUIRED NO_MODULE)
endif()
if(Boost_INCLUDE_DIR)
  include_directories(${Boost_INCLUDE_DIR})
endif()
if(Boost_INCLUDE_DIRS)
  include_directories(${Boost_INCLUDE_DIRS})
endif()

if (Boost_FILESYSTEM_FOUND)
  option ( FLX_USE_BOOST_FS "compile with BOOST::filesystem support" ON )
else()
  option ( FLX_USE_BOOST_FS "compile with BOOST::filesystem support" OFF )
endif()
if (Boost_REGEX_FOUND)
  option ( FLX_USE_BOOST_RX "compile with BOOST::regex support" ON )
else()
  option ( FLX_USE_BOOST_RX "compile with BOOST::regex support" OFF )
endif()

if("${FLX_USE_BOOST_FS}" STREQUAL "ON")
set ( FLX_USE_BOOST_FSb 1 )
else("${FLX_USE_BOOST_FS}" STREQUAL "ON")
set ( FLX_USE_BOOST_FSb 0 )
endif("${FLX_USE_BOOST_FS}" STREQUAL "ON")
if("${FLX_USE_BOOST_RX}" STREQUAL "ON")
set ( FLX_USE_BOOST_RXb 1 )
else("${FLX_USE_BOOST_RX}" STREQUAL "ON")
set ( FLX_USE_BOOST_RXb 0 )
endif("${FLX_USE_BOOST_RX}" STREQUAL "ON")

MARK_AS_ADVANCED(
  FLX_USE_BOOST_FS
  FLX_USE_BOOST_RX
)

#==========================================================
# ARPACK library
#==========================================================

if(UNIX)
  FIND_PACKAGE(ARPACK)
  if (ARPACK_FOUND)
    option ( FLX_USE_ARPACK "compile with arpack support" ON )
  else()
    option ( FLX_USE_ARPACK "compile with arpack support" OFF )
  endif()
else(UNIX)
    set ( FLX_USE_ARPACK "OFF" )
endif(UNIX)

if("${FLX_USE_ARPACK}" STREQUAL "ON")
  set ( FLX_USE_ARPACKb 1 )
else("${FLX_USE_ARPACK}" STREQUAL "ON")
  set ( FLX_USE_ARPACKb 0 )
endif("${FLX_USE_ARPACK}" STREQUAL "ON")


MARK_AS_ADVANCED(
  FLX_USE_ARPACK
)

###################################################################################
# compile
###################################################################################

#==========================================================
# prepare ...
#==========================================================

configure_file (
  "${PROJECT_SOURCE_DIR}/src/config.h.in"
  "${PROJECT_BINARY_DIR}/src/config.h"
)
include_directories(${CMAKE_CURRENT_BINARY_DIR})
include_directories(${PROJECT_BINARY_DIR})
include_directories("${PROJECT_BINARY_DIR}/src")
include_directories("${PROJECT_SOURCE_DIR}/src/flxlib")
include_directories("${PROJECT_SOURCE_DIR}/src/core")

# Determine platform-specific RPATH
if(APPLE)
    set(RPATH_STRING "@loader_path")  # macOS
elseif(UNIX)
    set(RPATH_STRING "$ORIGIN")       # Linux
endif()

if(MSVC)
  include (GenerateExportHeader)
endif(MSVC)

#==========================================================
# flxlib » base library
#==========================================================

# Create a shared library for common functionality
python_add_library(flxlib SHARED
  src/flxlib/flxexception.cpp
  src/flxlib/pdouble.cpp
  src/flxlib/flxMemoryManager.cpp
  src/flxlib/flxglobal.cpp
  src/flxlib/flxmath_rnd.cpp
  src/flxlib/flxmath_interp.cpp
  src/flxlib/flxmath.cpp
  src/flxlib/flxVec.cpp
  src/flxlib/flxio.cpp
  src/flxlib/flxfunction_data.cpp
  src/flxlib/flxfunction_ext.cpp
  src/flxlib/flxfunction_fun.cpp
  src/flxlib/flxfunction_ope.cpp
  src/flxlib/flxfunction.cpp
  src/flxlib/flxrandom_base.cpp
  src/flxlib/flxMtx.cpp
  src/flxlib/flxrbrv.cpp
  src/flxlib/flxrbrv_rvs.cpp
  src/flxlib/flxrbrv_rvs_read.cpp
  src/flxlib/flxmtxfun_data.cpp
  src/flxlib/flxmtxfun_fun.cpp
  src/flxlib/flxmtxfun.cpp
  src/flxlib/flxMtx_Eigen.cpp
  src/flxlib/flxMtx_Eigen_ARP.cpp
  src/flxlib/flxMtx_Eigen_GSL.cpp
  src/flxlib/flxdata.cpp
  src/flxlib/flxstringfun.cpp
  src/flxlib/flxstringfun_fun.cpp
  src/flxlib/flxstring.cpp
  src/flxlib/flxobjects.cpp
  src/flxlib/flxobjmtx.cpp
  src/flxlib/flxobjcommon.cpp
  src/flxlib/flxrandom.cpp
  src/flxlib/flxobjrbrv.cpp
  src/flxlib/flxostools_files.cpp
  src/flxlib/flxBayDA.cpp
  src/flxlib/flxBayUp.cpp
  src/flxlib/flxBayUp_obj.cpp
  src/flxlib/flxobjrandom.cpp
  src/flxlib/flxsensi.cpp
  src/flxlib/flxStatBox.cpp
)
if(MSVC)
  generate_export_header(flxlib)
  target_link_options(flxlib PRIVATE /VERBOSE:LIB)
endif(MSVC)

if(GSL_FOUND_STATIC)
    message(STATUS "GSL: static linking: ${GSL_STATIC_LIB} and ${GSLCBLAS_STATIC_LIB}")
    target_link_libraries(flxlib PRIVATE ${GSL_STATIC_LIB} ${GSLCBLAS_STATIC_LIB})
else()
    target_link_libraries(flxlib PRIVATE GSL::gsl)
endif()
target_link_libraries(flxlib PUBLIC pybind11::headers )
if(DEFINED RPATH_STRING)
  set_target_properties(flxlib PROPERTIES INSTALL_RPATH ${RPATH_STRING} )
endif()

target_compile_definitions(flxlib PRIVATE VERSION_INFO=${PROJECT_VERSION})
target_compile_definitions(flxlib PRIVATE GSL_STATIC)

#==========================================================
# core » main interface of Fesslix
#==========================================================

# Add a library using FindPython's tooling (pybind11 also provides a helper like this)
python_add_library(core MODULE
  src/core/fesslix.cpp
  src/core/python_interface_core.cpp
  WITH_SOABI)
target_link_libraries(core PRIVATE flxlib PRIVATE pybind11::headers)  #  PRIVATE GSL::gsl
if(DEFINED RPATH_STRING)
  set_target_properties(core PROPERTIES INSTALL_RPATH ${RPATH_STRING} )
endif()

target_compile_definitions(core PRIVATE VERSION_INFO=${PROJECT_VERSION})
target_compile_options(core PRIVATE -fvisibility=hidden)

#==========================================================
# install
#==========================================================

# The install directory is the output (wheel) directory
install(TARGETS core DESTINATION fesslix)
install(TARGETS flxlib
    RUNTIME DESTINATION fesslix  # For Windows (.dll)
    LIBRARY DESTINATION fesslix  # For Linux (.so)
    ARCHIVE DESTINATION fesslix  # For static libs (.lib)
)
install(FILES config/gausspoints.dat DESTINATION fesslix)


