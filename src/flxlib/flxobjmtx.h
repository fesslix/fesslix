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

#include "flxobjects.h"


class FLXLIB_EXPORT FlxCreateObjReaders_Mtx : public FlxCreateObjReaders {
  public:
    void createObjReaders (FlxObjectReadBox* objReadBox );
    void createFunReaders (FlxData* dataBox );
};



// -------------------------------------------------------------------------------------------------------------

/**
* @brief object class: extract a coefficient from a matrix
*
* mtxcoeff constName(i,j) = FlxFunction;
*/
class FLXLIB_EXPORT FlxObjMtxCoeff : public FlxObjBase {
  private:
    FlxMtxConstFun* cname;
    FlxFunction *i_fun;
    FlxFunction *j_fun;
    FlxFunction *c_fun;
    void task();
  public:
    FlxObjMtxCoeff ( const bool dolog, FlxMtxConstFun* cname, FlxFunction *i_fun, FlxFunction *j_fun, FlxFunction *c_fun ) : FlxObjBase(dolog), cname(cname), i_fun(i_fun), j_fun(j_fun), c_fun(c_fun) {}
    ~FlxObjMtxCoeff() { delete cname; delete i_fun; delete j_fun; delete c_fun; }
};

class FlxObjReadMtxCoeff : public FlxObjReadOutputBase {
  public:
    FlxObjBase* read ();
};

// -------------------------------------------------------------------------------------------------------------

/**
* @brief object class: output a matrix-constant
*
* mtxcalc FlxMtxFun
*/
class FLXLIB_EXPORT FlxObjMtxCalc : public FlxObjOutputBase {
  private:
    FlxMtxConstFun *fun;
    const bool only_coefs;
    void task();
  public:
    FlxObjMtxCalc ( bool dolog, FlxMtxConstFun *fun, std::string ostreamV, const bool only_coefs ) : FlxObjOutputBase(dolog,ostreamV), fun(fun), only_coefs(only_coefs) {}
    ~FlxObjMtxCalc() { delete fun; }
};

class FlxObjReadMtxCalc : public FlxObjReadOutputBase {
  public:
    FlxObjReadMtxCalc();
    FlxObjBase* read ();
};

// -------------------------------------------------------------------------------------------------------------

/**
* @brief object class: outputs a matrix in Octave(Matlab)-format
*
* transform_mtx2octave FlxMtxConstFun
*/
class FLXLIB_EXPORT FlxObjTransformMtx2Octave: public FlxObjOutputBase {
  private:
    FlxMtxConstFun *fun;
    void task();
  public:
    FlxObjTransformMtx2Octave ( bool dolog, FlxMtxConstFun *fun, std::string ostreamV ) : FlxObjOutputBase(dolog,ostreamV), fun(fun) {}
    ~FlxObjTransformMtx2Octave() { delete fun; }
};

class FlxObjReadTransformMtx2Octave : public FlxObjReadOutputBase {
  public:
    FlxObjBase* read ();
};

// -------------------------------------------------------------------------------------------------------------

/**
* @brief object class: deletes a defined matrix constant
*
* mtxconst_free FlxMtxFun
*/
class FLXLIB_EXPORT FlxObjMtxConst_free: public FlxObjBase {
  private:
    FlxMtxConstFun *funStr;
    void task();
  public:
    FlxObjMtxConst_free ( bool dolog, FlxMtxConstFun *funStr ) : FlxObjBase(dolog), funStr(funStr) {}
    ~FlxObjMtxConst_free() { delete funStr; }
};

class FlxObjReadMtxConst_free : public FlxObjReadBase {
  public:
    FlxObjBase* read ();
};

// -------------------------------------------------------------------------------------------------------------

class FLXLIB_EXPORT FlxObjMtxConstSeq : public FlxObjBase {
  protected:
    FlxMtxConstFun* mcn;
    tdouble* cv;
    FlxFunction* startF;
    FlxFunction* funCond;
    FlxFunction* funConst;
    void task();
  public:
    FlxObjMtxConstSeq ( bool dolog, FlxMtxConstFun* mcn, tdouble* cv, FlxFunction* startF, FlxFunction* funCond, FlxFunction* funConst ) : FlxObjBase(dolog), mcn(mcn), cv(cv), startF(startF), funCond(funCond), funConst(funConst) {}
    virtual ~FlxObjMtxConstSeq();
};
// is part of new

// -------------------------------------------------------------------------------------------------------------

/**
* @brief object class: defines a new mtxconst based on another matrix or a given size
*/
class FLXLIB_EXPORT FlxObjMtxConstNew : public FlxObjBase {
  protected:
    FlxMtxConstFun* mcn;
    FlxMtxConstFun* mtx_right;
    FlxFunction* rows;
    FlxFunction* cols;
    FlxFunction* val;
    void task();
  public:
    FlxObjMtxConstNew ( bool dolog, FlxMtxConstFun* mcn, FlxMtxConstFun* mtx_right, FlxFunction* rows, FlxFunction* cols, FlxFunction* val ) : FlxObjBase(dolog), mcn(mcn), mtx_right(mtx_right), rows(rows), cols(cols), val(val) {}
    virtual ~FlxObjMtxConstNew();
};

/**
* @brief object class: defines a new mtxconst based on a given matrix
*/
class FLXLIB_EXPORT FlxObjMtxConstNewU : public FlxObjBase {
  protected:
    FlxMtxConstFun* mcn;
    std::vector<FlxFunction*> vecV;
    const tuint nrows;
    const tuint ncols;
    void task();
  public:
    FlxObjMtxConstNewU ( bool dolog, FlxMtxConstFun* mcn, std::vector<FlxFunction*>& vecV, const tuint nrows, const tuint ncols ) : FlxObjBase(dolog), mcn(mcn), vecV(vecV), nrows(nrows), ncols(ncols) {}
    virtual ~FlxObjMtxConstNewU();
};

class FlxObjReadMtxConstNew : public FlxObjReadBase {
  public:
    FlxObjBase* read ();
    
    static void read_mtx(std::vector<FlxFunction*>& vecV, tuint& nrows, tuint& ncols);
    static void read_mtx_Matlab(std::vector<FlxFunction*>& vecV, tuint& nrows, tuint& ncols);
    static void read_seq(tdoublePtr& cv, FlxFunctionPtr& startF, FlxFunctionPtr& funCondL, FlxFunctionPtr& funConst);
};

// -------------------------------------------------------------------------------------------------------------

class FLXLIB_EXPORT FlxObjMtxConstOP : public FlxObjBase {
  protected:
    FlxMtxConstFun* mcn;
    const char c;
    FlxFunction* f;
    FlxMtxConstFun* s;
    tdouble* dp;
    void task();
  public:
    FlxObjMtxConstOP ( bool dolog, FlxMtxConstFun* mcn, const char c, FlxFunction* f, FlxMtxConstFun* s, tdouble* dp ) : FlxObjBase(dolog), mcn(mcn), c(c), f(f), s(s), dp(dp) {}
    virtual ~FlxObjMtxConstOP();
};

class FlxObjReadMtxConstOP : public FlxObjReadBase {
  public:
    FlxObjBase* read ();
};

// -------------------------------------------------------------------------------------------------------------

class FLXLIB_EXPORT FlxObjMtxConstSub : public FlxObjBase {
  public:
    enum Meth { cols, rows, seq };
  protected:
    FlxMtxConstFun* mcn_target;                // where to store
    FlxMtxConstFun* mcn_from;                // from where to take the 'sub'-matrix
    Meth methID;                        // store rows,columns,...
    std::vector<FlxFunction*> hv;        // help-vector
    const bool extract;
    
    void task();
  public:
    FlxObjMtxConstSub ( bool dolog, FlxMtxConstFun* mcn_target, FlxMtxConstFun* mcn_from, Meth methID, const std::vector<FlxFunction*> &hvV, const bool extract ) : FlxObjBase(dolog), mcn_target(mcn_target), mcn_from(mcn_from), methID(methID), hv(hvV), extract(extract) {}
    virtual ~FlxObjMtxConstSub();
};

class FlxObjReadMtxConstSub : public FlxObjReadBase {
  private:
     void read_subInfo(FlxObjMtxConstSub::Meth& meth, std::vector<FlxFunction*> &hv);
  public:
    FlxObjBase* read ();
};

// -------------------------------------------------------------------------------------------------------------

class FLXLIB_EXPORT FlxObjMtxConstMult : public FlxObjBase {
  protected:
    FlxMtxConstFun* mcn;
    FlxMtxConstFun* mn1;
    FlxMtxConstFun* mn2;
    void task();
  public:
    FlxObjMtxConstMult ( bool dolog, FlxMtxConstFun* mcn, FlxMtxConstFun* mn1, FlxMtxConstFun* mn2 ) : FlxObjBase(dolog), mcn(mcn), mn1(mn1), mn2(mn2) {}
    virtual ~FlxObjMtxConstMult();
};

class FlxObjReadMtxConstMult : public FlxObjReadBase {
  public:
    FlxObjBase* read ();
};

// -------------------------------------------------------------------------------------------------------------

class FLXLIB_EXPORT FlxObjMtxConstFromFile : public FlxObjBase {
  protected:
    FlxMtxConstFun* mcn_target;                // where to store
    FlxFunction* cols;
    FlxString* strV;

    void task();
  public:
    FlxObjMtxConstFromFile ( bool dolog, FlxMtxConstFun* mcn_target, FlxFunction* cols, FlxString* strV ) : FlxObjBase(dolog), mcn_target(mcn_target), cols(cols), strV(strV) {}
    virtual ~FlxObjMtxConstFromFile();
};

class FlxObjReadMtxConstFromFile : public FlxObjReadBase {
  public:
    FlxObjBase* read ();
};

// -------------------------------------------------------------------------------------------------------------

class FLXLIB_EXPORT FlxObjMtxConstTranspose : public FlxObjBase {
  protected:
    FlxMtxConstFun* mcn_target;                // where to store

    void task();
  public:
    FlxObjMtxConstTranspose ( bool dolog, FlxMtxConstFun* mcn_target ) : FlxObjBase(dolog), mcn_target(mcn_target) {}
    virtual ~FlxObjMtxConstTranspose() { delete mcn_target; }
};

class FlxObjReadMtxConstTranspose : public FlxObjReadBase {
  public:
    FlxObjBase* read ();
};


// -------------------------- FUNCTIONS ----------------------------------------------------------------------






