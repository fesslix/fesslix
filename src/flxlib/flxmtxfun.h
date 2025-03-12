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

#include "flxmtxfun_data.h"
#include <set>


//=================== General Definitions ===========================

class FlxMtxFunBox;

class FLXLIB_EXPORT FlxMtxBoxBase {
  protected:
    static FlxConstMtxBox* ConstMtxBox;
  public:
    virtual ~FlxMtxBoxBase() {};
    static void set_ConstMtxBox(FlxConstMtxBox* ConstMtxBoxV) { ConstMtxBox = ConstMtxBoxV; }
};

class FlxObjBase;
class FLXLIB_EXPORT FlxMtxConstFun : public FlxReaderBase2, public FlxMtxBoxBase {
  private:
    FlxString* mtxName;
    FlxObjBase* block;
    std::string mtxName_str;
    /**
    * @brief instances of this FlxString (used to copy this FlxString)
    */
    tuint* instances;
    tuint intrnl_id;        // = 0: not internal
    static std::set<tuint> internl_lst;
    static tuint intrnl_rqst_id();
    static void intrnl_free_id(const tuint idV);
    static const std::string intrnl_get_id_str(const tuint idV);
    
    void free_mem();
  public:
    FlxMtxConstFun(const bool is_block_allowed);
    FlxMtxConstFun(const char* mtxName_strV, FlxObjBase* block=NULL);
    FlxMtxConstFun(FlxMtxConstFun &rhs);
    ~FlxMtxConstFun();
    void assign ( FlxMtxConstFun* funA );
    
    const std::string& eval();
    const std::string write();
    const bool search_circref(FlxFunction* fcr);
};


//=================== Dealing with FlxFunctions and matrices as parameter ==========================

class FLXLIB_EXPORT FunBaseFun_MtxConst : public FunBaseFun_multPara {
  protected:
    static FlxConstMtxBox* mtxConsts;
    std::list< FlxMtxConstFun* >* ParaList_MtxName;
  public:
    FunBaseFun_MtxConst (std::vector<FunBase*> *ParaListV, std::list< FlxMtxConstFun* >* ParaList_MtxName) : FunBaseFun_multPara(ParaListV), ParaList_MtxName(ParaList_MtxName) {};
    FunBaseFun_MtxConst (std::list< FlxMtxConstFun* >* ParaList_MtxName) : FunBaseFun_multPara(new std::vector<FunBase*>), ParaList_MtxName(ParaList_MtxName) {};
    virtual ~FunBaseFun_MtxConst();
    virtual const bool search_circref(FlxFunction* fcr);
    const std::string write();
    virtual const bool dependOn_Const(const tdouble* const thenumber);
    const bool optimize(FunBasePtr& optf, const Fun_OptimizeInfo &foi) { return false; };
    
    static void set_dataMtxConsts ( FlxConstMtxBox* mtxConstsV ) { mtxConsts = mtxConstsV; }
};

class FLXLIB_EXPORT FunReadFunBase_MtxConst : public FunReadFunBase {
  protected:
    std::list< FlxMtxConstFun* >* read_para( tuint Pnmb, bool errSerious);
  public:
    FunReadFunBase_MtxConst(const bool from_library = false) : FunReadFunBase(from_library) {}
    virtual ~FunReadFunBase_MtxConst() {};
};


