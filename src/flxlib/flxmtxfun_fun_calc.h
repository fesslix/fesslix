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

#include "flxmtxfun.h"


//=================== FlxFunction - Part ==========================


/**
* @brief arithmetic class: maximum coefficient in a matrix
*/
class FLXLIB_EXPORT FunMaxMin : public FunBaseFun_MtxConst {
  protected:
    const bool is_max;
  public:
    FunMaxMin (const bool is_max, std::vector< FunBase* > *ParaList_Fun, std::list< FlxMtxConstFun* > *ParaList_Mtx) : FunBaseFun_MtxConst(ParaList_Fun,ParaList_Mtx), is_max(is_max) {};
    const tdouble calc();
    const std::string write_v() { return is_max?"max":"min"; }
    const std::string write();
};

class FunReadFunMaxMin : public FunReadFunBase_MtxConst {
  private:
    const bool is_max;
  public:
    FunReadFunMaxMin(const bool is_max) : is_max(is_max) {}
    FunBase* read ( bool errSerious );
};

// -------------------------------------------------------------------------------

/**
* @brief arithmetic class: index of maximum coefficient in a matrix
*/
class FLXLIB_EXPORT FunMaxMinID : public FunBaseFun_MtxConst {
  protected:
    const bool is_max;
  public:
    FunMaxMinID (const bool is_max, std::vector< FunBase* > *ParaList_Fun, std::list< FlxMtxConstFun* > *ParaList_Mtx) : FunBaseFun_MtxConst(ParaList_Fun,ParaList_Mtx), is_max(is_max) {};
    const tdouble calc();
    const std::string write_v() { return is_max?"max_id":"min_id"; }
    const std::string write();
};

class FunReadFunMaxMinID : public FunReadFunBase_MtxConst {
  private:
    const bool is_max;
  public:
    FunReadFunMaxMinID(const bool is_max) : is_max(is_max) {}
    FunBase* read ( bool errSerious );
};

// -------------------------------------------------------------------------------

/**
* @brief arithmetic class: coefficient of a matrix
*/
class FLXLIB_EXPORT FunMtxCoeff : public FunBaseFun_MtxConst {
  private:
    FunBase* rid;
    FunBase* cid;
  public:
    FunMtxCoeff (std::list< FlxMtxConstFun* > *ParaList_Mtx, FunBase* rid, FunBase* cid) : FunBaseFun_MtxConst(ParaList_Mtx), rid(rid), cid(cid) {};
    ~FunMtxCoeff() { delete rid; delete cid;} 
    const tdouble calc();
    const bool search_circref(FlxFunction* fcr);
    const std::string write();
    const std::string write_v() { return "mtxcoeff"; }
    const bool dependOn_Const(const tdouble* const thenumber);
};

class FunReadFunMtxCoeff : public FunReadFunBase_MtxConst {
  public:
    FunBase* read ( bool errSerious );
};

// -------------------------------------------------------------------------------

/**
* @brief arithmetic class: number of rows of a matrix
*/
class FLXLIB_EXPORT FunMtxRows : public FunBaseFun_MtxConst {
  public:
    FunMtxRows (std::list< FlxMtxConstFun* > *ParaList_Mtx) : FunBaseFun_MtxConst(ParaList_Mtx) {};
    const tdouble calc();
    const std::string write_v() { return "mtxrows"; }
    virtual const bool dependOn_Const(const tdouble* const thenumber);
};

class FunReadFunMtxRows : public FunReadFunBase_MtxConst {
  public:
    FunBase* read ( bool errSerious ) { return new FunMtxRows( read_para(1,errSerious) ); }
};

// -------------------------------------------------------------------------------

/**
* @brief arithmetic class: number of columns of a matrix
*/
class FLXLIB_EXPORT FunMtxCols : public FunBaseFun_MtxConst {
  public:
    FunMtxCols (std::list< FlxMtxConstFun* > *ParaList_Mtx) : FunBaseFun_MtxConst(ParaList_Mtx) {};
    const tdouble calc();
    const std::string write_v() { return "mtxcols"; }
    virtual const bool dependOn_Const(const tdouble* const thenumber);
};

class FunReadFunMtxCols : public FunReadFunBase_MtxConst {
  public:
    FunBase* read ( bool errSerious ) { return new FunMtxCols( read_para(1,errSerious) ); }
};

// -------------------------------------------------------------------------------

/**
* @brief arithmetic class: sum over all coefficients of a matrix
*/
class FLXLIB_EXPORT FunMtxSum : public FunBaseFun_MtxConst {
  public:
    FunMtxSum (std::list< FlxMtxConstFun* > *ParaList_Mtx) : FunBaseFun_MtxConst(ParaList_Mtx) {};
    const tdouble calc();
    const std::string write_v() { return "mtxsum"; }
    virtual const bool dependOn_Const(const tdouble* const thenumber);
};

class FunReadFunMtxSum : public FunReadFunBase_MtxConst {
  public:
    FunBase* read ( bool errSerious ) { return new FunMtxSum( read_para(1,errSerious) ); }
};

// -------------------------------------------------------------------------------

/**
* @brief arithmetic class: multiply all coefficients of a matrix
*/
class FLXLIB_EXPORT FunMtxProd : public FunBaseFun_MtxConst {
  public:
    FunMtxProd (std::list< FlxMtxConstFun* > *ParaList_Mtx) : FunBaseFun_MtxConst(ParaList_Mtx) {};
    const tdouble calc();
    const std::string write_v() { return "mtxprod"; }
    virtual const bool dependOn_Const(const tdouble* const thenumber);
};

class FunReadFunMtxProd : public FunReadFunBase_MtxConst {
  public:
    FunBase* read ( bool errSerious ) { return new FunMtxProd( read_para(1,errSerious) ); }
};

// -------------------------------------------------------------------------------

/**
* @brief arithmetic class: multiply all coefficients of a matrix
*/
class FLXLIB_EXPORT FunMtxMean : public FunBaseFun_MtxConst {
  public:
    FunMtxMean (std::list< FlxMtxConstFun* > *ParaList_Mtx) : FunBaseFun_MtxConst(ParaList_Mtx) {};
    const tdouble calc();
    const std::string write_v() { return "mtxmean"; }
    virtual const bool dependOn_Const(const tdouble* const thenumber);
};

class FunReadFunMtxMean : public FunReadFunBase_MtxConst {
  public:
    FunBase* read ( bool errSerious ) { return new FunMtxMean( read_para(1,errSerious) ); }
};

// -------------------------------------------------------------------------------

/**
* @brief arithmetic class: multiply all coefficients of a matrix
*/
class FLXLIB_EXPORT FunMtxSd : public FunBaseFun_MtxConst {
  public:
    FunMtxSd (std::list< FlxMtxConstFun* > *ParaList_Mtx) : FunBaseFun_MtxConst(ParaList_Mtx) {};
    const tdouble calc();
    const std::string write_v() { return "mtxsd"; }
    virtual const bool dependOn_Const(const tdouble* const thenumber);
};

class FunReadFunMtxSd : public FunReadFunBase_MtxConst {
  public:
    FunBase* read ( bool errSerious ) { return new FunMtxSd( read_para(1,errSerious) ); }
};

// -------------------------------------------------------------------------------

/**
* @brief arithmetic class: 
*/
class FLXLIB_EXPORT FunMtxVec_Norm2 : public FunBaseFun_MtxConst {
  public:
    FunMtxVec_Norm2 (std::list< FlxMtxConstFun* > *ParaList_Mtx) : FunBaseFun_MtxConst(ParaList_Mtx) {};
    const tdouble calc();
    const std::string write_v() { return "vec_norm2"; }
    virtual const bool dependOn_Const(const tdouble* const thenumber);
};

class FunReadFunMtxVec_Norm2 : public FunReadFunBase_MtxConst {
  public:
    FunBase* read ( bool errSerious ) { return new FunMtxVec_Norm2( read_para(1,errSerious) ); }
};

// -------------------------------------------------------------------------------



