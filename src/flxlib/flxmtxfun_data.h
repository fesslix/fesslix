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

#include "flxstring.h"


class FLXLIB_EXPORT FlxSMtx {
  private:
    tuint nrows;
    tuint ncols;
    flxVec dptr;
    FlxMtx_base* mtx;
    /**
    * @brief check if rows and columns of matrices match.
    */
    void check_1(const FlxSMtx& rhs) const;
    /**
    * @brief check if types of both matrices match
    */
    void check_2(const FlxSMtx& rhs) const;
  public:
    FlxSMtx(FlxSMtx& rhs);
    FlxSMtx(const tuint nrows, const tuint ncols, const tdouble val);
    /**
    * @brief takes over memory managment of dp
    */
    FlxSMtx(tdouble* dp, const tuint& nrows, const tuint& ncols);
    FlxSMtx(const flxVec& rhs);
    FlxSMtx(FlxMtx_base* mtx) : nrows(0), ncols(0), dptr(0), mtx(mtx) {};
    ~FlxSMtx ();

    const tuint get_nrows() const { return nrows; }
    const tuint get_ncols() const { return ncols; }
    const size_t get_Ncoeff() const { return nrows*ncols; }
    
    /**
    * @brief get coefficient (i,j)  -  has to be within range of matrix (not checked!!!)
    */
    const tdouble operator()(const tuint& i, const tuint& j) const;/**
    * @brief get coefficient (i) from internal array -  has to be within range of matrix (not checked!!!)
    */
    const tdouble operator()(const tuint& i) const;
    /**
    * @brief change a coefficient (i,j)  -  has to be within range of matrix (not checked!!!)
    */
    void insert(const tuint& i, const tuint& j, const tdouble& d);
    FlxSMtx& operator=(const FlxSMtx& rhs);
    FlxSMtx& operator=(const tdouble& rhs);
    FlxSMtx& operator+=(const FlxSMtx& rhs);
    FlxSMtx& operator*=(const tdouble& rhs);
    const tdouble max() const;
    const tdouble min() const;
    const size_t maxID() const;
    const size_t minID() const;
//     virtual const flxVec& get_flxVec() const;
    tdouble* get_internalPtr(const bool throwErr);
    
    /**
    * @brief *this = m1*m2
    */
    void mult(const FlxSMtx& m1, const FlxSMtx& m2);
};

FLXLIB_EXPORT std::ostream& operator<<(std::ostream& os, const FlxSMtx& M);
FLXLIB_EXPORT void SMtxBase_write_fullPrec(std::ostream& os, const FlxSMtx& M);


//=================== Constant-Box ==========================

/**
* @brief A class for storing constant values (the result of the evaluated function is stored).
*/
class FLXLIB_EXPORT FlxConstMtxBox {
  private:
    #if FLX_DEBUG
      /**
      * @brief number of instances of this calss - only one class may exist
      */
      static int Cinst;
    #endif
    std::map<std::string, FlxSMtx*> box;
  public:
    FlxConstMtxBox();
    ~FlxConstMtxBox ();
    /**
    * @brief insert: the reference of a variable is changed by 'insert' (difference to 'FlxConstantBox')
    *                i.e. in case the entry exists, it is overwritten
    */
    void insert ( const std::string& name, FlxSMtx* value);
    /**
    * @brief get: NOTE: you have to COPY the returned reference -> otherwise you will overwrite the stored matrix
    */
    FlxSMtx* get ( const std::string& name, const bool err_if_unknown = false );
    /**
    * @param forceSize true:object must exist and size must match -> otherwise an error is thrown
    */
    FlxSMtx* get ( const std::string& name, const tuint Nrows, const tuint Ncols, const bool forceSize );
    /**
    * @returns reference to vector 'name' with dimension N
    * if the vector does not exist it is defined
    * @param N if N=0: the vector MUST exist! returns N of this vector
    *                  if N=0 and the object does not exist, an error is thrown
    *          if N>0 && forceSize=false: if such an object does not exist, it is created
    * @param forceSize true:vector MUST exist and size must match -> otherwise an error is thrown
    */
    tdouble* get_Vec ( const std::string& name, tuint& N, const bool forceSize=false );
    tdouble* get_Vec ( const tuint N, const std::string& name, const bool forceSize=false);
    const iVec get_Vec_cast2tuint(const std::string& name, const bool No0);
    /**
    * @returns reference to internal storage of 'name' with dimension Nr x Nc
    * if the object does not exist it is defined
    * @param N if Nr=Nc=0: the object MUST exist (otherwise an error is thrown)
    *                      -> returns Nr and Nc of this object
    * @param forceSize true:object must exist and size must match -> otherwise an error is thrown
    */
    tdouble* get_Mtx ( const std::string& name, tuint& Nr, tuint& Nc, const bool forceSize=false );
    tdouble* get_Mtx ( const tuint Nr, const tuint Nc, const std::string& name, const bool forceSize=false);
    void declareC ( const std::string& name );
    /**
    * @brief removes matrix 'name' from the list
    */
    void freeC ( const std::string& name );
};



