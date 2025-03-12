/* Fesslix - Stochastic Analysis
 * Copyright (C) 2010-2018 Wolfgang Betz
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

#include "flxMtx.h"



// -------------------------- Arnoldi - ARpack ----------------------------------------------------------------------

/**
* @brief Calculates Eigenvalues of Amtx
*/
void MtxEigenValue_ARP(FlxMtx_baseS& Amtx, const int M, flxVec& EigenValues, std::vector<flxVec>& Eigenvectors);

/**
* @brief Calculates Eigenvalues of a generalized matrix eigenvalue problem Ax = lBx
*/
void MtxEigenValue_ARP(FlxMtx_baseS& Amtx, FlxMtx_baseS& Bmtx, const int M, flxVec& EigenValues, std::vector<flxVec>& Eigenvectors);





