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

#include "flxmtxfun_fun_calc.h"
#include "flxmtxfun_fun.h"

using namespace std;

void flxmtxfun_fun_insert(FlxFunctionBox& FunBox)
{  
  FunBox.insert("max", new FunReadFunMaxMin(true) );
  FunBox.insert("min", new FunReadFunMaxMin(false) );
  FunBox.insert("maxid", new FunReadFunMaxMinID(true) );
  FunBox.insert("minid", new FunReadFunMaxMinID(false) );
  FunBox.insert("mtxcoeff", new FunReadFunMtxCoeff() );
  FunBox.insert("mtxrows", new FunReadFunMtxRows() );
  FunBox.insert("mtxcols", new FunReadFunMtxCols() );
  FunBox.insert("mtxsum", new FunReadFunMtxSum() );
  FunBox.insert("mtxprod", new FunReadFunMtxProd() );
  FunBox.insert("mtxmean", new FunReadFunMtxMean() );
  FunBox.insert("mtxsd", new FunReadFunMtxSd() );
  FunBox.insert("vec_norm2", new FunReadFunMtxVec_Norm2() );
}

const tdouble FunMaxMin::calc()
{
  tdouble cv = ZERO;
  bool b1 = true;
  // search for maximum in function list
    for(std::vector<FunBase*>::const_iterator i=ParaList->begin();i!=ParaList->end();++i) {
      const tdouble t = (*i)->calc();
      if (b1) {
        cv = t;
        b1 = false;
      } else {
        if (is_max) {
          if (t > cv) cv = t;
        } else {
          if (t < cv) cv = t;
        }
      }
    }
  // search for maximum in matrix list
    for (list< FlxMtxConstFun* >::const_iterator i = ParaList_MtxName->begin(); i!=ParaList_MtxName->end(); ++i) {
      const std::string mtxname = (*i)->eval();
      // obtain the local max-value
        FlxSMtx* mtx = mtxConsts->get(mtxname,false);
        tdouble t;
        if (mtx==NULL) {
          std::ostringstream ssV;
          ssV << "A matrix with the name '" << mtxname << "' does not exist.";
          throw FlxException_NeglectInInteractive("FunMaxMin::calc", ssV.str() );
        } else {
          if (is_max) {
            t = mtx->max();
          } else {
            t = mtx->min();
          }
        }
      if (b1) {
        cv = t;
        b1 = false;
      } else {
        if (is_max) {
          if (t > cv) cv = t;
        } else {
          if (t < cv) cv = t;
        }
      }
    } 
  return cv;
}

const std::string FunMaxMin::write()
{
  std::string str1 = write_v() + "(";
  bool b1 = false;
  for(std::vector<FunBase*>::const_iterator i=ParaList->begin();i!=ParaList->end();++i) {
    if (b1) str1 += ",";
    else b1 = true;
    str1 += (*i)->write();
  } 
  for (list< FlxMtxConstFun* >::const_iterator i = ParaList_MtxName->begin(); i!=ParaList_MtxName->end(); ++i) {
    if (b1) str1 += ",";
    else b1 = true;
    str1 += '{';
    str1 += (*i)->write();
    str1 += '}';
  } 
  str1 += ")";
  return str1;
}

const bool FunMtxRows::dependOn_Const(const tdouble*const thenumber)
{
  return false;
}

const tdouble FunMtxRows::calc()
{
  const std::string mtxname = (*ParaList_MtxName->begin())->eval();
  FlxSMtx* mtx = mtxConsts->get(mtxname,true);
  return mtx->get_nrows();
}

const bool FunMtxCols::dependOn_Const(const tdouble*const thenumber)
{
  return false;
}

const tdouble FunMtxCols::calc()
{
  const std::string mtxname = (*ParaList_MtxName->begin())->eval();
  FlxSMtx* mtx = mtxConsts->get(mtxname,true);
  return mtx->get_ncols();
}

const bool FunMtxSum::dependOn_Const(const tdouble*const thenumber)
{
  return false;
}

const tdouble FunMtxSum::calc()
{
  const std::string mtxname = (*ParaList_MtxName->begin())->eval();
  FlxSMtx* mtx = mtxConsts->get(mtxname,true);
  const tdouble* dptr = mtx->get_internalPtr(true);
  const size_t N = mtx->get_Ncoeff();
  qdouble rv(N,false);
  for (size_t i=0;i<N;++i) {
    rv += dptr[i];
  }
  return rv.cast2double();
}

const bool FunMtxProd::dependOn_Const(const tdouble*const thenumber)
{
  return false;
}

const tdouble FunMtxProd::calc()
{
  const std::string mtxname = (*ParaList_MtxName->begin())->eval();
  FlxSMtx* mtx = mtxConsts->get(mtxname,true);
  const tdouble* dptr = mtx->get_internalPtr(true);
  const size_t N = mtx->get_Ncoeff();
  tdouble rv = ONE;
  for (size_t i=0;i<N;++i) {
    rv *= dptr[i];
  }
  return rv;
}

const bool FunMtxVec_Norm2::dependOn_Const(const tdouble*const thenumber)
{
  return false;
}

const tdouble FunMtxVec_Norm2::calc()
{
  const std::string mtxname = (*ParaList_MtxName->begin())->eval();
  tuint N = 0;
  tdouble* vp = mtxConsts->get_Vec(mtxname,N);
  flxVec v(vp,N);
  return v.get_Norm2();
}

FunBase* FunReadFunMtxCoeff::read ( bool errSerious ) 
{
  std::list< FlxMtxConstFun* >* plst = new std::list< FlxMtxConstFun* >();
  FunBase* rid = NULL;
  FunBase* cid = NULL;
  try {
    plst->push_back(new FlxMtxConstFun(false));
    reader->getChar(',');
    rid = FunctionList->read(errSerious);
    if (reader->whatIsNextChar()==',') {
      reader->getChar(',');
      cid = FunctionList->read(errSerious);
    } else {
      cid = new FunNumber(ZERO);
    }
  } catch (FlxException &e) {
    FLXMSG("FunReadFunMtxCoeff::read",1);
    for (list< FlxMtxConstFun* >::const_iterator i = plst->begin(); i!=plst->end(); ++i) {
      delete (*i);
    }
    delete plst;
    if (rid) delete rid;
    if (cid) delete cid;
    throw;
  }
  return new FunMtxCoeff(plst, rid, cid );
}

const tdouble FunMtxCoeff::calc()
{
  const std::string mtxname = (*ParaList_MtxName->begin())->eval();
  FlxSMtx* mtx = mtxConsts->get(mtxname,true);
  const tuint i = tuint_from(rid->calc(),"row index",true,true,rid);
  const tuint j = tuint_from(cid->calc(),"column index",true,true,cid);
  if (i>=mtx->get_nrows() || j>=mtx->get_ncols()) {
    std::ostringstream ssV;
    ssV << "Index (" << i << "," << j << ") not within the valid bounds of [";
    if (mtx->get_nrows()==1) ssV << "0; ";
    else ssV << "0,...," << mtx->get_nrows()-1 << " ; ";
    if (mtx->get_ncols()==1) ssV << "0].";
    else ssV << "0,...," << mtx->get_ncols()-1 << "].";
    throw FlxException("FunMtxCoeff::calc", ssV.str() );
  }
  return mtx->operator()(i,j);
}

const bool FunMtxCoeff::search_circref(FlxFunction* fcr) {
  if (rid->search_circref(fcr)) return true;
  if (cid->search_circref(fcr)) return true;
  return FunBaseFun_MtxConst::search_circref(fcr);
}

const std::string FunMtxCoeff::write() {
  return write_v() + "(" + (*ParaList_MtxName->begin())->write() + "," + rid->write() + "," + cid->write() + ")";
}

const bool FunMtxCoeff::dependOn_Const(const tdouble* const thenumber) {
  if (rid->dependOn_Const(thenumber)) return true;
  if (cid->dependOn_Const(thenumber)) return true;
  return FunBaseFun_MtxConst::dependOn_Const(thenumber);
}

FunBase* FunReadFunMaxMin::read ( bool errSerious ) 
{
  std::list< FunBase* > ParaList_Fun;
  std::list< FlxMtxConstFun* > *ParaList_Mtx = new std::list< FlxMtxConstFun* >;
  try {
    while (true) {
      if (reader->whatIsNextChar()=='{') {        // read MtxConstName
        reader->getChar('{');
        ParaList_Mtx->push_back(new FlxMtxConstFun(true));
        reader->getChar('}');
      } else {        // read Function
        ParaList_Fun.push_back( FunctionList->read(errSerious) );
      }
      if (reader->whatIsNextChar()==',') {
        reader->getChar(',');
      } else {
        break;
      }
    }
    std::vector< FunBase* >* ParaList = new std::vector<FunBase*>();
    ParaList->reserve(ParaList_Fun.size());
    for (list< FunBase* >::const_iterator i = ParaList_Fun.begin(); i!=ParaList_Fun.end(); ++i) {
      ParaList->push_back(*i);
    }
    return new FunMaxMin(is_max,ParaList,ParaList_Mtx);
  } catch (FlxException& e) {
    for (list< FunBase* >::const_iterator i = ParaList_Fun.begin(); i!=ParaList_Fun.end(); ++i) {
      delete (*i);
    }
    for (list< FlxMtxConstFun* >::const_iterator i = ParaList_Mtx->begin(); i!=ParaList_Mtx->end(); ++i) {
      delete (*i);
    }
    throw;
  }
}

const tdouble FunMaxMinID::calc()
{
  if (ParaList_MtxName->size()!=1) throw FlxException_Crude("FunMaxMinID::calc");
  list< FlxMtxConstFun* >::const_iterator i = ParaList_MtxName->begin();
  const std::string mtxname = (*i)->eval();
  FlxSMtx* mtx = mtxConsts->get(mtxname,false);
  if (mtx==NULL) {
    std::ostringstream ssV;
    ssV << "A matrix with the name '" << mtxname << "' does not exist.";
    throw FlxException_NeglectInInteractive("FunMaxMin::calc", ssV.str() );
  } else {
    if (is_max) {
      return mtx->maxID();
    } else {
      return mtx->minID();
    }
  }
}

const std::string FunMaxMinID::write()
{
  std::string str1 = write_v() + "(";
  bool b1 = false;
  for (list< FlxMtxConstFun* >::const_iterator i = ParaList_MtxName->begin(); i!=ParaList_MtxName->end(); ++i) {
    if (b1) str1 += ",";
    else b1 = true;
    str1 += (*i)->write();
  } 
  str1 += ")";
  return str1;
}

FunBase* FunReadFunMaxMinID::read ( bool errSerious ) 
{
  std::vector< FunBase* >* ParaList = new std::vector<FunBase*>();
  std::list< FlxMtxConstFun* > *ParaList_Mtx = new std::list< FlxMtxConstFun* >;
  ParaList_Mtx->push_back(new FlxMtxConstFun(true));
  return new FunMaxMinID(is_max,ParaList,ParaList_Mtx);
}

const bool FunMtxMean::dependOn_Const(const tdouble*const thenumber)
{
  return false;
}

const tdouble FunMtxMean::calc()
{
  const std::string mtxname = (*ParaList_MtxName->begin())->eval();
  FlxSMtx* mtx = mtxConsts->get(mtxname,true);
  const tdouble* dptr = mtx->get_internalPtr(true);
  flxVec tv(dptr,mtx->get_Ncoeff(),false);
  return tv.get_Mean();
}

const bool FunMtxSd::dependOn_Const(const tdouble*const thenumber)
{
  return false;
}

const tdouble FunMtxSd::calc()
{
  const std::string mtxname = (*ParaList_MtxName->begin())->eval();
  FlxSMtx* mtx = mtxConsts->get(mtxname,true);
  const tdouble* dptr = mtx->get_internalPtr(true);
  flxVec tv(dptr,mtx->get_Ncoeff(),false);
  const tdouble vmean = tv.get_Mean();
  return tv.get_sd(vmean);
}



