/* Fesslix - Stochastic Analysis
 * Copyright (C) 2010-2026 Wolfgang Betz
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

#include "flxMtx_Eigen_ARP.h"
#include "flxMtx_Eigen_GSL.h"



// -------------------------- General Interface ----------------------------------------------------------------------

// Mode:
//   1 -> ARP
//   2 -> GSL
//   3 -> GLS stabilized


FLXLIB_EXPORT void flxeigen_logInfo(std::ostream& lout);



/**
* @brief Calculates Eigenvalues of Amtx (symmetric)
*/
FLXLIB_EXPORT void MtxEigenValue(FlxMtx_baseS& Amtx, const tuint M, flxVec& EigenValues, std::vector<flxVec>& Eigenvectors, const int Mode);

/**
* @brief Calculates Eigenvalues of Amtx (nonsymmetric)
*/
FLXLIB_EXPORT void MtxEigenValue(FlxMtx& Amtx, const tuint M, flxVec& EigenValues, std::vector<flxVec>& Eigenvectors, const int Mode);


/**
* @brief Calculates Eigenvalues of a generalized matrix eigenvalue problem Ax = lBx
*/
FLXLIB_EXPORT void MtxEigenValue(FlxMtx_baseS& Amtx, FlxMtx_baseS& Bmtx, const tuint M, flxVec& EigenValues, std::vector<flxVec>& Eigenvectors, const int Mode);


FLXLIB_EXPORT void EV_orientation(const tuint M, std::vector< flxVec >& Eigenvectors);



