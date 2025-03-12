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

#ifndef flx_StatBox_H
#define flx_StatBox_H

#include "flxMtx.h"


class FLXLIB_EXPORT FlxStatBox {
  public:
    const tuint N;        // the total number of samples that can be stored
    const tuint M;        // the length (dimension) of each sample
  private:
    tuint Nc;                // the number of samples stored
    /**
    * @brief tp defines the storage vector
    * storage structure
    * -----------------
    *   1,...,N
    *   ... for M times
    */
    tdouble* tp;
  public:
    FlxStatBox(tuint N, tuint M);
    ~FlxStatBox();
    
    /**
    * @returns the number of entries stored
    */
    const tuint get_Nc() const { return Nc; }
    /**
    * @brief empties the list of samples
    */
    void clear() { Nc = 0; }
    /**
    * @brief adds the sample sv to the list
    */
    void add(const flxVec& sv);
    /**
    * @brief computes expectation and standard deviation
    * @note the dimensions of the two vectors must be equal to Nc
    */
    void compute_ExpSd(flxVec& eV, flxVec& sV);
    /**
    * @brief computes the mean
    */
    void compute_mean(flxVec& eV);
    
    /**
    * @brief plots the confidence intervals of the stored data
    * @param ostrm the stream to send the output to
    * @param clevels a vector containing confidence levels (NOTE: the vector will be altered - sorted)
    * @param fixW length of an output-double (see Double2String)
    * 
    * structure of the file:
    *         * expectation
    *          * standard deviation
    *          * all confidence levels stored in clevels
    */
    void plot_confidence(std::ostream& ostrm, flxVec& clevels, const int fixW, flxVec* initCol);
    /**
    * @brief computes the eigen-directions of a fraction of the samples contained
    * @param p fraction of the samples to use
    * @param sde input: the current y-vector to sample around
    *                   output: the standard deviations in the eigen-directions
    * @param covmtx output: the covariance of the samples considered
    * @param hvN vector that is used as a temporary storage (size N)
    * @param hvM vector that is used as a temporary storage (size M)
    */
    void get_smpl_eigen(const tdouble p, const tuint pNmax, flxVec& sde, FlxMtxSym& covmtx, flxVec& hvN1, flxVec& hvN2, flxVec& hvM, iVec& ivN, std::vector<flxVec>& Eigenvectors, FlxMtx& Tinv, const tdouble p_single, const tuint pNmax_single);
    /**
    * @brief computes the eigen-directions of all samples contained in the box
    */
    void get_smpl_eigen(FlxMtxSym& covmtx, flxVec& evs, std::vector<flxVec>& Eigenvectors, FlxMtx& Tinv);
    
};











#endif // flx_StatBox_H
