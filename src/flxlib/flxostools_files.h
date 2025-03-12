/* Fesslix - Stochastic Analysis
 * Copyright (C) 2010-2025 Wolfgang Betz
 *
 * Fesslix is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * Fesslix is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Fesslix.  If not, see <http://www.gnu.org/licenses/>. 
 */

#pragma once

#include "flxmath.h"

#include <string>

/**
* @brief copies files recursively from source_dir to dest_dir
* @note the dest_dir must NOT exist!!!
*/
FLXLIB_EXPORT void copyDir( const std::string& source_dir, const std::string& dest_dir );
FLXLIB_EXPORT void copyFile( const std::string& source, const std::string& dest, const bool overwrite=false );
FLXLIB_EXPORT void moveFile( const std::string& source, const std::string& dest );
FLXLIB_EXPORT const bool existsDir( const std::string& dir );
FLXLIB_EXPORT const bool createDir( const std::string& dir );
FLXLIB_EXPORT void getFiles( const std::string& sourcePath, const std::string& sourcePattern, std::vector< std::string >& match_files );

/**
* @brief returns the number of files removed
*/
FLXLIB_EXPORT const tuint delDir( const std::string& dir );
FLXLIB_EXPORT const bool fileIsReadable(const std::string& fileName);


