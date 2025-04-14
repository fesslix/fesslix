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

#include "flxmtxfun.h"
#include "flxobjmtx.h"

using namespace std;

FlxConstMtxBox* FlxMtxBoxBase::ConstMtxBox= NULL;
FlxConstMtxBox* FunBaseFun_MtxConst::mtxConsts = NULL;
std::set<tuint> FlxMtxConstFun::internl_lst;


FlxMtxConstFun::FlxMtxConstFun(const char* mtxName_strV, FlxObjBase* block)
: mtxName(NULL), block(block), mtxName_str(mtxName_strV), instances(new tuint(0)), intrnl_id(0)
{

}

FlxMtxConstFun::FlxMtxConstFun(const bool is_block_allowed)
: mtxName(NULL), block(NULL), instances(new tuint(0)), intrnl_id(0)
{
  try {
    if (is_block_allowed) {        // check for special definitions
      const char rnc = reader->whatIsNextChar();
      if (rnc=='{') {
        std::vector<FlxFunction*> vecV;
        tuint nrows, ncols;
        FlxObjReadMtxConstNew::read_mtx(vecV,nrows,ncols);
        intrnl_id = intrnl_rqst_id();
        mtxName_str = intrnl_get_id_str(intrnl_id);
        block = new FlxObjMtxConstNewU(false,new FlxMtxConstFun(mtxName_str.c_str()),vecV,nrows,ncols);
        return;
      } else if (rnc=='!') {
        reader->getChar();
        const std::string expr = reader->getWord(true);
        if (expr=="unimtx") {
          reader->getChar('(');
          FlxFunction* rows = new FlxFunction(funReader,false);
          FlxFunction* cols = NULL;
          FlxFunction* val = NULL;
          try {
            if (reader->whatIsNextChar()==',') {
              reader->getChar(',',false);
              cols = new FlxFunction(funReader,false);
              if (reader->whatIsNextChar()==',') {
                reader->getChar(',',false);
                val = new FlxFunction(funReader,false);
              }
            }
            reader->getChar(')');
          } catch (FlxException& e) {
            delete rows;
            if (cols) delete cols;
            if (val) delete val;
            throw;
          }
          intrnl_id = intrnl_rqst_id();
          mtxName_str = intrnl_get_id_str(intrnl_id);
          block = new FlxObjMtxConstNew(false,new FlxMtxConstFun(mtxName_str.c_str()),NULL,rows,cols,val);
          return;
        } else if (expr=="seq") {
          tdouble* cv; FlxFunction* startF; FlxFunction* funCondL; FlxFunction* funConst;
          FlxObjReadMtxConstNew::read_seq(cv, startF, funCondL, funConst);
          intrnl_id = intrnl_rqst_id();
          mtxName_str = intrnl_get_id_str(intrnl_id);
          block = new FlxObjMtxConstSeq(false,new FlxMtxConstFun(mtxName_str.c_str()),cv,startF,funCondL,funConst);
          return;
        } else {
          throw FlxException("FlxMtxConstFun::FlxMtxConstFun", "Expression'" + expr + "' is not allowed at this point.");
        }
      }
    }
    mtxName = new FlxString(false,false);
    if (mtxName->is_static()) {
      mtxName_str = mtxName->eval_word(true);
      delete mtxName; mtxName = NULL;
    }
    if (is_block_allowed) {
      if (reader->whatIsNextChar()=='!') {
        reader->getChar('!');
        block = FlxObjReadCodeBlock::read_block(true);
      }
    }
  } catch (FlxException& e) {
    free_mem();
    throw;
  }
}

FlxMtxConstFun::FlxMtxConstFun(FlxMtxConstFun& rhs)
{
  instances = rhs.instances;
  mtxName = rhs.mtxName;
  mtxName_str = rhs.mtxName_str;
  block = rhs.block;
  intrnl_id = rhs.intrnl_id;
  (*instances)++;
}

FlxMtxConstFun::~FlxMtxConstFun()
{
  free_mem();
}
const string& FlxMtxConstFun::eval()
{
  if (block) block->exec();
  if (mtxName) {
    mtxName_str = mtxName->eval_word(true);
  }
  return mtxName_str;
}

const string FlxMtxConstFun::write()
{
  if (intrnl_id==0) {
    std::string str = mtxName?(mtxName->write()):mtxName_str;
    if (block) {
      str += "!{...}";
    }
    return str;
  } else {
    throw FlxException_NotImplemented("FlxMtxConstFun::write");
  }
}

void FlxMtxConstFun::free_mem()
{
  if (instances == NULL) return;
  if ( (*instances) == 0)  {
    if (mtxName) delete mtxName;
    if ( instances ) delete instances;
    if ( block ) delete block;
    intrnl_free_id(intrnl_id);
  } else {
    (*instances)--;
  }
}

const bool FlxMtxConstFun::search_circref(FlxFunction* fcr)
{
  return false;
}

void FlxMtxConstFun::assign(FlxMtxConstFun* funA)
{
  if (this == funA) return;
  if (mtxName && mtxName == funA->mtxName) return;
  if (block && block == funA->block) return;
  if ((*instances)==0) {
    free_mem(); 
  } else {
    (*instances)--;
  }
  mtxName = funA->mtxName;
    funA->mtxName = NULL;
  mtxName_str = funA->mtxName_str;
  instances = funA->instances;
    funA->instances = NULL;
  block = funA->block;
    funA->block = NULL;
  intrnl_id = funA->intrnl_id;
    funA->intrnl_id = 0;
  delete funA;
}

tuint FlxMtxConstFun::intrnl_rqst_id()
{
  tuint count = 1;
  while (true) {
    if (internl_lst.find(count)==internl_lst.end()) {
      internl_lst.insert(count);
      return count;
    }
    ++count;
  }
}

void FlxMtxConstFun::intrnl_free_id(const tuint idV)
{
  if (idV==0) return;
  internl_lst.erase(idV);
  const std::string instr = intrnl_get_id_str(idV);
  if (ConstMtxBox->get(instr,false)) {
    ConstMtxBox->freeC(instr);
  }
}

const string FlxMtxConstFun::intrnl_get_id_str(const tuint idV)
{
  ostringstream ssV;
  ssV << "internal_tmpmtx" << idV;
  return ssV.str();
}


std::list< FlxMtxConstFun* >* FunReadFunBase_MtxConst::read_para(tuint Pnmb, bool errSerious)
{
  std::list< FlxMtxConstFun* >* plst = new std::list< FlxMtxConstFun* >();
  try {
    do {
      plst->push_back(new FlxMtxConstFun(true));
      if ( reader->whatIsNextChar() == ',' ) {
        reader->getChar(',',errSerious);
      } else {
        break;
      }
    } while ( true );
    if (Pnmb > 0 && plst->size() != Pnmb) {
      ostringstream ssV;
      ssV << "Expected '" << Pnmb << "' parameters and not '" << plst->size() << "' parameters.";
      FlxError(errSerious,"FunReadFunBase_MtxConst::read_para_1", ssV.str(), reader->getCurrentPos());
    }
  } catch (FlxException &e) {
    FLXMSG("FunReadFunBase_MtxConst::read_para_2",1);
    for (list< FlxMtxConstFun* >::const_iterator i = plst->begin(); i!=plst->end(); ++i) {
      delete (*i);
    }
    delete plst;
    throw;
  }
  return plst;
}

FunBaseFun_MtxConst::~FunBaseFun_MtxConst()
{
  for (list< FlxMtxConstFun* >::const_iterator i = ParaList_MtxName->begin(); i!=ParaList_MtxName->end(); ++i) {
    delete (*i);
  }
  delete ParaList_MtxName;
}

const std::string FunBaseFun_MtxConst::write()
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

const bool FunBaseFun_MtxConst::search_circref(FlxFunction* fcr)
{
  if ( FunBaseFun_multPara::search_circref(fcr) ) return true;
  for (list< FlxMtxConstFun* >::const_iterator i = ParaList_MtxName->begin(); i!=ParaList_MtxName->end(); ++i) {
    if ( (*i)->search_circref(fcr) ) return true;
  }
  return false;
}

const bool FunBaseFun_MtxConst::dependOn_Const(const tdouble*const thenumber)
{
  throw FlxException_NotImplemented("FunBaseFun_MtxConst::dependOn_Const");
}

FlxMtxFun_const::FlxMtxFun_const(const flxVec& rhs)
: FlxMtxFun_base(rhs.get_N())
{
  res_vec = rhs;
}

void FlxMtxFun_const::eval()
{

}

FlxMtxFun_Py::FlxMtxFun_Py(const tuint N, py::function pyfunc)
: FlxMtxFun_base(N), pyfunc(pyfunc)
{

}

void FlxMtxFun_Py::eval()
{
  py::object result;
  try {
      result = pyfunc();
  } catch (const py::error_already_set &e) {
      std::ostringstream ssV;
      ssV << "Error in evaluating Python expression: " << e.what();
      throw FlxException("FlxMtxFun_Py::eval_01", ssV.str() );
  }
  flxVec rhs = parse_py_obj_as_flxVec(result,"Result of Python function");
  res_vec.assign_save(rhs);
}

FlxMtxFun_base* parse_FlxMtxFun(const tuint N, py::object pyobj, std::string descr)
{
  return get_ReadManager()->parse_FlxMtxFun(N, pyobj, descr);
}

FlxMtxFun_base* parse_py_para_as_flxMtxFun(const tuint N, const std::string& para_name, py::dict config)
{
  if (config.contains(para_name.c_str()) == false) {
    std::ostringstream ssV;
    ssV << "Key '" << para_name << "' not found in Python <dict>.";
    throw FlxException_NeglectInInteractive("parse_py_para_as_flxMtxFun_01", ssV.str());
  }
  return parse_FlxMtxFun(N,config[para_name.c_str()], "key '"+para_name+"' in Python <dict>");

}
