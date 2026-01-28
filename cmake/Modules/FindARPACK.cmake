# Fesslix - Stochastic Analysis
# Copyright (C) 2010-2026 Wolfgang Betz
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

# ===========================================================================*\
# This file was originally taken from
# https://graphics.rwth-aachen.de:9000/OpenFlipper-Free/OpenFlipper-Free/blob/31f2b44ab6390f946b27a1501b5eb976c83c97af/cmake/FindARPACK.cmake
# in June 2016
# ===========================================================================*/
# the copyright-notice of the original file is:
# ===========================================================================*\
#                                                                             *
#                               OpenFlipper                                   *
#            Copyright (c) 2001-2015, RWTH-Aachen University                 *
#            Department of Computer Graphics and Multimedia                  *
#                           All rights reserved.                             *
#                             www.openflipper.org                            *
#                                                                            *
# ---------------------------------------------------------------------------*
#  This file is part of OpenFlipper.                                         *
# ---------------------------------------------------------------------------*
#                                                                            *
#  Redistribution and use in source and binary forms, with or without        *
#  modification, are permitted provided that the following conditions        *
#  are met:                                                                  *
#                                                                            *
#  1. Redistributions of source code must retain the above copyright notice, *
#     this list of conditions and the following disclaimer.                  *
#                                                                            *
#  2. Redistributions in binary form must reproduce the above copyright      *
#     notice, this list of conditions and the following disclaimer in the    *
#     documentation and/or other materials provided with the distribution.   *
#                                                                            *
#  3. Neither the name of the copyright holder nor the names of its          *
#     contributors may be used to endorse or promote products derived from   *
#     this software without specific prior written permission.               *
#                                                                            *
#  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS       *
#  "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED *
#  TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A           *
#  PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER *
#  OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,  *
#  EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,       *
#  PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR        *
#  PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF    *
#  LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING      *
#  NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS        *
#  SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.              *
#                                                                            *
# ===========================================================================*/


# - Try to find ARPACK
# Once done this will define
#
#  ARPACK_FOUND        - system has ARPACK
#  ARPACK_INCLUDE_DIR  - the ARPACK include directories
#  ARPACK_LIBRARY      - Link this to use ARPACK
#  ARPACKpp_LIBRARY    - Link this to use ARPACK++
#  ARPACK_LIBRARIES    - Link this to use ARPACK together with ARPACK++


find_path (ARPACK_INCLUDE_DIR NAMES arpack++/argeig.h
  HINTS ENV ARPACK_INCLUDE_DIR
  PATHS /usr/include/arpack++ $ENV{HOME}/opt/arpack++/include $ENV{HOME}/projects/arpack++/include "C:\\libs\\arpack++\\include"
  DOC "ARPACK Include Directory")

IF ( WIN32 )
find_library( ARPACK_LIBRARY arpack.lib
              PATHS "C:\\libs\\arpack++\\lib"
            )
ELSE( WIN32 )
  find_library( ARPACK_LIBRARY arpack
                PATHS /usr/lib /usr/lib64 $ENV{ARPACK_LIBDIR}
              )
  find_library( ARPACKpp_LIBRARY arpack++
                PATHS /usr/lib /usr/lib64 ${ARPACK_INCLUDE_DIR}/../src/.libs
              )

  list( APPEND ARPACK_LIBRARIES  ${ARPACK_LIBRARY} )
  IF(ARPACKpp_LIBRARY)
    list( APPEND ARPACK_LIBRARIES  ${ARPACKpp_LIBRARY} )
  ENDIF(ARPACKpp_LIBRARY)
  
ENDIF( WIN32 )

IF (ARPACK_INCLUDE_DIR AND ARPACK_LIBRARY)
  SET(ARPACK_FOUND TRUE)
ELSE ()
  SET(ARPACK_FOUND FALSE)
ENDIF ()

mark_as_advanced(
  ARPACK_DIR
  ARPACK_INCLUDE_DIR 
  ARPACK_LIBRARY 
  ARPACKpp_LIBRARY
  ARPACK_LIBRARIES
)
