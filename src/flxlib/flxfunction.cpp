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

#define FLXLIB_CPP

#include "flxfunction.h"
#include "flxfunction_fun.h"
#include "flxfunction_ope.h"
#include "flxfunction_ope_calc.h"

using namespace std;

FunReadBase* FunReadBase::startLink = NULL;
FlxConstantBox* FlxBoxBase::ConstantBox = NULL;
FlxVarBox* FlxBoxBase::VarBox = NULL;
FlxFunctionBox* FlxBoxBase::funBox = NULL;
FunReadSTART* FunReadFunBase_0::FunctionList = NULL;
const tdouble* FunPara::ParaList = NULL;
tuint FunPara::ParaListSize = 0;
#if FLX_DEBUG
  int FlxVarBox::Cinst = 0;
  int FlxFunctionBox::Cinst = 0;
  int FlxFunctionReader::Cinst = 0;
#endif
tuint FunReadPara::numbofpara = 0;
tuint FlxFunDeg::default_deg = 0;

void FunBase::child_optimize(FunBasePtr& optf, const Fun_OptimizeInfo& foi)
{
  FunBase* child_t = optf;
  optf = NULL;
  try {
    while (child_t->optimize(optf,foi)) {
      if (optf) {
        delete child_t;
        child_t = optf;
        optf = NULL;
      }
    }
  } catch (FlxException& e) {
    FLXMSG("FunBase::child_optimize",1);
    if (optf) delete optf;
    optf = child_t;
    throw;
  }
  optf = child_t;
}

void FunBase::calc_me(FunBasePtr& numref)
{
  numref = new FunNumber(calc());
}

const bool FunBase::is_number(FunBase* ftc)
{
  return (dynamic_cast<FunNumber*>(ftc) != NULL);
}


FlxVarBox::FlxVarBox()
{
  #if FLX_DEBUG
    ++Cinst;
    if (Cinst > 1) {
      std::ostringstream ssV;
      ssV << "More than one instance of 'FlxVarBox' created ...";
      throw FlxException("FlxVarBox::FlxVarBox", ssV.str() );
    }
  #endif
  FlxBoxBase::set_varBox(this);
}

FlxVarBox::~FlxVarBox() {
  for (map<string, FlxFunction*>::iterator pos = box.begin(); pos != box.end(); ++pos) delete pos->second;
}

void FlxVarBox::insert(const std::string& name, FlxFunction* value) {
  #if FLX_DEBUG
    if (value==NULL) {
      throw FlxException("FlxVarBox::insert", "ERROR" );
    }
  #endif
  pair<std::string, FlxFunction*>Element(name, value);
  if ( ! box.insert(Element).second ) {
    map<std::string, FlxFunction*>::iterator pos;
    pos = box.find(name);
    pos->second->assign(value);
  }
}

const string FlxVarBox::get(FlxFunction* fun)
{
  for (map<string, FlxFunction*>::iterator pos = box.begin(); pos != box.end(); ++pos) {
    if (pos->second == fun) {
      return pos->first;
    }
  }
  std::ostringstream ssV;
  ssV << "Function not part of the list.";
  throw FlxException("FlxVarBox::get", ssV.str() );
}

FlxFunction* FlxVarBox::get(const string& name) {
  map<std::string, FlxFunction*>::iterator pos;
  pos = box.find(name);
  if ( pos != box.end() ) {
    return pos->second;
  } else {
    return NULL;
  }
}

void FlxVarBox::declareV(const string& name) {
  if ( get(name) == NULL ) {
    insert( name, new FlxFunction ( new FunNumber(0) ) );
  }
}

void FlxBoxBase::set_funBox(FlxFunctionBox* funBoxV)
{
  funBox = funBoxV;
}

void FlxBoxBase::set_constBox(FlxConstantBox* ConstantBoxV)
{
  #if FLX_DEBUG
    if (ConstantBox) {
      throw FlxException_Crude("FlxBoxBase::set_constBox");
    }
  #endif
  ConstantBox = ConstantBoxV;
}

void FlxBoxBase::set_varBox(FlxVarBox* VarBoxV)
{
  VarBox = VarBoxV;
}

FlxFunction::~FlxFunction() {
  if (instances == NULL ) return;
  if ( (*instances) == 0)  {
    if ( fun ) delete fun;
    if ( instances ) delete instances;
    if (read_pos) delete read_pos;
  } else {
    (*instances)--;
  }
}

void FlxFunction::check_FlxFunction(const FlxFunction* funtocheck)
{
  const FlxFunction_Combine_Add* fa = dynamic_cast<const FlxFunction_Combine_Add*>(funtocheck);
  if (fa!=NULL) {
    throw FlxException_NotImplemented("FlxFunction::check_FlxFunction");
  }
}

const bool FlxFunction::optimize(FunBasePtr& optf, const Fun_OptimizeInfo& foi)
{
  child_optimize(fun,foi);
  return false;
}

FlxFunction::FlxFunction(FlxFunctionReader* funReader, bool errSerious)
: fun(NULL),instances(new tuint(0)), read_pos(NULL)
{ 
  try {
    read_pos = new FlxReaderPos();
    funReader->get_reader()->getCurrentPos(*read_pos);
    fun = funReader->read(errSerious);
  } catch (...) {
    delete instances;
    if (read_pos) delete read_pos;
    throw;
  }
}

FlxFunction::FlxFunction(const FlxFunction& FlxF)
{
  check_FlxFunction(&FlxF);
  if (FlxF.fun->evalw()) {
    instances = new tuint(0);
    read_pos = NULL;
    fun = new FunNumber(FlxF.fun->calc());
  } else {
    instances = FlxF.instances;
    fun = FlxF.fun;
    read_pos = FlxF.read_pos;
    (*instances)++;
  }
}

FlxFunction& FlxFunction::operator=(const FlxFunction& FlxF)
{
  check_FlxFunction(&FlxF);
  if ( this != &FlxF ) {
    instances = FlxF.instances;
    fun = FlxF.fun;
    read_pos = FlxF.read_pos;
    (*instances)++;
  }
  return *this;
}

void FlxFunction::assign(FlxFunction* funA) {
  check_FlxFunction(funA);
  if (this == funA) return;
  if (this->fun == funA->fun) return;
  if ((*instances)==0) {
    delete fun; 
    delete instances;
    delete read_pos;
  } else {
    (*instances)--;
  }
  fun = funA->fun; 
  funA->fun = NULL;
  instances = funA->instances;
  funA->instances = NULL;
  read_pos = funA->read_pos;
  funA->read_pos = NULL;
  delete funA;
}

const bool FlxFunction::is_zero()
{
  FunNumber* S = dynamic_cast<FunNumber*>(fun);
  if (S==NULL) return false;
  if (fabs(S->get_thenumber()) <= GlobalVar.TOL()) {
    return true;
  } else {
    return false;
  }
}

const tdouble FlxFunction::cast2positive(const bool errSerious)
{
  const tdouble cv = fun->calc();
  if (cv<=ZERO) {
    std::ostringstream ssV; 
    ssV << "Number must not be negative or zero (" << cv << ")." << std::endl;
    ssV << "\tThe function is: " << fun->write();
    if (read_pos) {
      ssV << std::endl << '\t' << ReadStream::write_ReaderPos(*read_pos); 
    }
    FlxError(errSerious,"FlxFunction::cast2positive", "Expected positive double!", ssV.str() );
  }
  return cv;
}

const tdouble FlxFunction::cast2positive_or0(const bool errSerious)
{
  const tdouble cv = fun->calc();
  if (cv<ZERO) {
    std::ostringstream ssV; 
    ssV << "Number must not be negative (" << cv << ")."; 
    if (read_pos) {
      ssV << std::endl << '\t' << ReadStream::write_ReaderPos(*read_pos); 
    }
    FlxError(errSerious,"FlxFunction::cast2positive_or0", "Expected positive double or zero!", ssV.str() );
  }
  return cv;
}

const tnlong FlxFunction::cast2tnlong(const bool errSerious)
{
  const tdouble cv = fun->calc();
  const tdouble rcv  = round_flx(cv);
  if ( rcv <= ZERO) { 
    std::ostringstream ssV; 
    ssV << "Number must not be negative or zero (" << cv << "->" << rcv << ")."; 
    if (read_pos) {
      ssV << std::endl << '\t' << ReadStream::write_ReaderPos(*read_pos); 
    }
    FlxError(errSerious,"FlxFunction::cast2tnlong", "Expected unsigned integer!", ssV.str() );     
  }
  return static_cast<tnlong>(rcv);
}

const tulong FlxFunction::cast2tulong(const bool errSerious)
{
  const tdouble cv = fun->calc();
  const tdouble rcv  = round_flx(cv);
  if ( rcv <= ZERO) { 
    std::ostringstream ssV; 
    ssV << "Number must not be negative or zero (" << cv << "->" << rcv << ").";
    if (read_pos) {
      ssV << std::endl << '\t' << ReadStream::write_ReaderPos(*read_pos); 
    }
    FlxError(errSerious,"FlxFunction::cast2tulong", "Expected unsigned integer!", ssV.str() );     
  }
  return static_cast<tulong>(rcv);
}

const tulong FlxFunction::cast2tulongW0(const bool errSerious)
{
  const tdouble cv = fun->calc();
  const tdouble rcv  = round_flx(cv);
  if ( rcv < ZERO) { 
    std::ostringstream ssV; 
    ssV << "Number must not be negative (" << cv << "->" << rcv << ").";
    if (read_pos) {
      ssV << std::endl << '\t' << ReadStream::write_ReaderPos(*read_pos); 
    }
    FlxError(errSerious,"FlxFunction::cast2tulongW0", "Expected unsigned integer!", ssV.str() );     
  }
  return static_cast<tulong>(rcv);
}

const tuint FlxFunction::cast2tuint(const bool errSerious)
{
  const tdouble cv = fun->calc();
  const tdouble rcv  = round_flx(cv);
  if ( rcv <= ZERO) { 
    std::ostringstream ssV; 
    ssV << "Number must not be negative or zero (" << cv << "->" << rcv << ")."; 
    if (read_pos) {
      ssV << std::endl << '\t' << ReadStream::write_ReaderPos(*read_pos); 
    }
    FlxError(errSerious,"FlxFunction::cast2tuint", "Expected unsigned integer!", ssV.str() );
  }
  return static_cast<tuint>(rcv);
}

const tuint FlxFunction::cast2tuintW0(const bool errSerious)
{
  const tdouble cv = fun->calc();
  const tdouble rcv  = round_flx(cv);
  if ( rcv < ZERO) { 
    std::ostringstream ssV; 
    ssV << "Number must not be negative (" << cv << "->" << rcv << ")."; 
    if (read_pos) {
      ssV << std::endl << '\t' << ReadStream::write_ReaderPos(*read_pos); 
    }
    FlxError(errSerious,"FlxFunction::cast2tuintW0", "Expected unsigned integer!", ssV.str() );     
  }
  return static_cast<tuint>(rcv);
}

const tnlong FlxFunction::cast2tnlongW0(const bool errSerious)
{
  const tdouble cv = fun->calc();
  const tdouble rcv  = round_flx(cv);
  if ( rcv < ZERO) { 
    std::ostringstream ssV; 
    ssV << "Number must not be negative (" << cv << "->" << rcv << ").";
    if (read_pos) {
      ssV << std::endl << '\t' << ReadStream::write_ReaderPos(*read_pos); 
    }
    FlxError(errSerious,"FlxFunction::cast2tnlongW0", "Expected unsigned integer!", ssV.str() );     
  }
  return static_cast<tnlong>(rcv);
}

FlxFunction_Combine_Add::FlxFunction_Combine_Add(FlxFunction* f1V, FlxFunction* f2V)
:f1(f1V), f2(f2V)
{
  instances = new tuint(0);
  fun = new FunAdd(f1->fun,f2->fun);
}

FlxFunction_Combine_Add::~FlxFunction_Combine_Add()
{
  FunAdd* fa = dynamic_cast<FunAdd*>(fun);
  fa->set_childs_zero();
  delete fa;
  fun = NULL;
  
  delete f1; f1=NULL;
  delete f2; f2=NULL;
}

void FlxFunction_Combine_Add::assign(FlxFunction* funA)
{
  throw FlxException_NotImplemented("FlxFunction_Combine_Add::assign");
}

FunBase* FunReadFunDummy::read(bool errSerious)
{
  // an error must be throw at this point, because a user-defined function must not be defined AND called in the same read-environment.
  throw FlxException("FunReadFunDummy::read","This reader is not supposed to be called.");
}

FlxFunctionBox::FlxFunctionBox() {
  #if FLX_DEBUG
    ++Cinst;
    if (Cinst > 1) {
      std::ostringstream ssV;
      ssV << "More than one instance of 'FlxFunctionBox' created ...";
      throw FlxException("FlxFunctionBox::FlxFunctionBox", ssV.str() );
    }
  #endif
  FlxBoxBase::set_funBox(this);
  // initialize random seed
    srand ( tuint(time(NULL)) );
    rv_initialize(true);
  // all objects of predefined functions have to be listed there
    flxfunction_fun_insert(*this);
}

FlxFunctionBox::~FlxFunctionBox() {
  for (map<std::string, FunReadFunBase*>::iterator pos = box.begin(); pos != box.end(); ++pos) {
    FunReadFunBase* tf = pos->second;
    if (tf->is_from_library()==false) delete tf;
  }
}

void FlxFunctionBox::insert(const std::string& name, FunReadFunBase* fR) {
  pair<string, FunReadFunBase*>Element(name, fR);
  if ( ! box.insert(Element).second ) {
    map<string, FunReadFunBase*>::iterator pos;
    pos = box.find(name);
    delete pos->second;
    pos->second = fR;
  }
}

void FlxFunctionBox::declareF(const std::string& name) {
  if ( get(name) == NULL ) {
    insert( name, new FunReadFunDummy() );
  }
}

FunReadFunBase* FlxFunctionBox::get(const string& name) {
  map<string, FunReadFunBase*>::iterator pos;
  pos = box.find(name);
  if ( pos != box.end() ) {
    return pos->second;
  } else {
    return NULL;
  }
}

void FlxFunctionBox::remove(const string& name)
{
  box.erase(name);
}


const bool FunBaseOperat1::optimize(FunBasePtr& optf, const Fun_OptimizeInfo& foi)
{
  child_optimize(child_1, foi);
  if ( is_number(child_1) ) {
    calc_me(optf);
    return true;
  } else {
    return false;
  }
}

const bool FunBaseOperat2::optimize(FunBasePtr& optf, const Fun_OptimizeInfo& foi)
{
  child_optimize(child_1, foi);
  child_optimize(child_2, foi);
  if ( is_number(child_1) && is_number(child_2) ) {
    calc_me(optf);
    return true;
  } else {
    return false;
  }
}

const bool FunBaseOperat3::optimize(FunBasePtr& optf, const Fun_OptimizeInfo& foi)
{
  child_optimize(child_1, foi);
  child_optimize(child_2, foi);
  child_optimize(child_3, foi);
  if ( is_number(child_1) && is_number(child_2) && is_number(child_3) ) {
    calc_me(optf);
    return true;
  } else {
    return false;
  }
}

FunBaseFun_onePara::FunBaseFun_onePara(std::vector<FunBase*>* ParaListV)
:child_1(ParaListV->operator[](0))
{
  for ( size_t i = 1; i < ParaListV->size(); ++i ) 
    delete ParaListV->operator[](i); 
  delete ParaListV;
}

FunBaseFun_onePara::~FunBaseFun_onePara()
{
  delete child_1;
}

const std::string FunBaseFun_onePara::write()
{
  std::string str1 = write_v() + "(";
  str1 += child_1->write();
  str1 += ")";
  return str1;
}

const bool FunBaseFun_onePara::optimize(FunBasePtr& optf, const Fun_OptimizeInfo& foi)
{
  child_optimize(child_1, foi);
  if ( is_number(child_1) ) {
    calc_me(optf);
    return true;
  } else {
    return false;
  }
}

FunBaseFun_multPara::~FunBaseFun_multPara() {
  for ( size_t i = 0; i < ParaList->size(); ++i ) 
    delete ParaList->operator[](i); 
  delete ParaList;
}

const bool FunBaseFun_multPara::search_circref(FlxFunction* fcr) {
  bool b1 = false;
  for ( size_t i = 0; i < ParaList->size(); ++i ) {
    b1 = b1 || ParaList->operator[](i)->search_circref(fcr); 
  }
  return b1;
}

const bool FunBaseFun_multPara::evalw() {
  bool b1 = false;
  for ( size_t i = 0; i < ParaList->size(); ++i ) {
    b1 = b1 || ParaList->operator[](i)->evalw(); 
  }
  return b1;
}

const std::string FunBaseFun_multPara::write()
{
  std::string str1 = write_v() + '(';
  for (size_t i = 0; i < ParaList->size(); ++i) {
    if (i > 0) str1 += ',';
    str1 += ParaList->operator[](i)->write();
  }
  str1 += ')';
  return str1;
}

const bool FunBaseFun_multPara::dependOn_Const(const tdouble*const thenumber)
{
  for ( size_t i = 0; i < ParaList->size(); ++i ) {
    if ( ParaList->operator[](i)->dependOn_Const(thenumber) ) return true;
  }
  return false;
}

const bool FunBaseFun_multPara::optimize(FunBasePtr& optf, const Fun_OptimizeInfo& foi)
{
  bool b1=true;
  for ( size_t i = 0; i < ParaList->size(); ++i ) {
    child_optimize( ParaList->operator[](i), foi );
    b1 = b1 && is_number(ParaList->operator[](i));
  }
  if ( b1 ) {
    calc_me(optf);
    return true;
  } else {
    return false;
  }
}

const tdouble FunBaseFun_Python::calc()
{
  if (ParaList->size()==0) {
    py::object result;
    try {
       result = pyfunc();
    } catch (const py::error_already_set &e) {
        std::ostringstream ssV;
        ssV << "Error in evaluating Python expression: " << e.what();
        throw FlxException("FunBaseFun_Python::calc_01", ssV.str() );
    }
    if (py::isinstance<py::float_>(result)) {
        return result.cast<tdouble>();
    } else {
        throw FlxException("FunBaseFun_Python::calc_02", "Result of Python function has wrong type.");
    }
  } else {
    throw FlxException_NotImplemented("FunBaseFun_Python::calc_03");
  }
}

const std::string FunBaseFun_Python::write()
{
  std::string str1 = write_v() + '(';
  str1 += pyFunName;
  for (size_t i = 0; i < ParaList->size(); ++i) {
    str1 += ',';
    str1 += ParaList->operator[](i)->write();
  }
  str1 += ')';
  return str1;
}

const bool FunBaseFun_Python::optimize(FunBasePtr& optf, const Fun_OptimizeInfo& foi)
{
  for ( size_t i = 0; i < ParaList->size(); ++i ) {
    child_optimize( ParaList->operator[](i), foi );
  }
  return false;
}


FunUser::FunUser(std::vector< FunBase* >* ParaListV, FlxFunction* funV, const string& fname, const tuint numbofpara)
 : FunBaseFun_multPara(ParaListV), fun(funV), fname(fname), numbofpara(numbofpara), tPL(ParaList->size()), tPLp(&tPL[0])
{
  
}

const bool FunUser::dependOn_Const(const tdouble*const thenumber)
{
  if (fun->dependOn_Const(thenumber)) return true;
  return FunBaseFun_multPara::dependOn_Const(thenumber);
}

const bool FunUser::search_circref(FlxFunction* fcr) {
  return FunBaseFun_multPara::search_circref(fcr) || fun->search_circref(fcr);
}

const string FunUser::write()
{
  if (numbofpara>0) {
    return FunBaseFun_multPara::write();
  } else {
    const std::string str1 = write_v() + "()";
    return str1;
  }
}

const tdouble FunUser::calc()
{
  for (size_t i = 0; i < numbofpara; ++i) {
    tPLp[i] = ParaListP[i]->calc();
  }
  const tdouble* const tT = FunPara::ParaList;
  const tuint tTS = FunPara::ParaListSize;
  FunPara::ParaList = tPLp;
  FunPara::ParaListSize = numbofpara;
  const tdouble res = fun->calc();
  FunPara::ParaList = tT;
  FunPara::ParaListSize = tTS;
  return res;
}

const tdouble FunPara::calc()
{
  if (ParaList == NULL) {
    ostringstream ssV_2;
    ssV_2 << "ParaList is not defined.";
    throw FlxException("FunPara::calc_1", ssV_2.str(), "This error should not have been occurred.");  
  }
  if ( ParaListSize >= theindex && theindex >= 1 ) {
    return ParaList[theindex - 1];
  } else {
    ostringstream ssV_2;
    ssV_2 << "A function uses a parameter which is out of range. (index=" << theindex << ")";
    throw FlxException("FunPara::calc_2", ssV_2.str(), "This error is based on faulty function definition.");
  }
}

const std::string FunPara::write()
{
  std::ostringstream ssV;
  ssV << "$" << theindex;
  return ssV.str();;
}

const std::string FunDummyFun::write_v()
{
  ostringstream ssV_2;
  ssV_2 << "ERROR";
  throw FlxException("FunDummyFun::write_v", ssV_2.str());  
}

const bool FunDummyFun::optimize(FunBasePtr& optf, const Fun_OptimizeInfo& foi)
{
  calc_me(optf);
  return true;
}

const std::string FunDummy::write()
{
  ostringstream ssV_2;
  ssV_2 << "ERROR";
  throw FlxException("FunDummy::write", ssV_2.str());  
}

const bool FunDummy::optimize(FunBasePtr& optf, const Fun_OptimizeInfo& foi)
{
  calc_me(optf);
  return true;
}

FunReadBase::FunReadBase (const int Priority, bool isEndNode) : Next(NULL), priority(Priority) {
  if ( (!( isEndNode && priority == -1)) && priority < 0 ) {
    ostringstream ssV_2;
    ssV_2 << "Priority (" << priority << ") not allowed - value has to be greater than '0'.";
    throw FlxException("FunReadBase::FunReadBase_1", ssV_2.str(), "This error is based on faulty source code.");   
  }
}

FunReadBase* FunReadBase::insert ( FunReadBase *TheFunRead ) {
  if ( priority == -1 || TheFunRead->priority < priority ) {
    if ( TheFunRead->priority > 0 ) {
      TheFunRead->Next = this;
      return TheFunRead;
    } else {
      ostringstream ssV_2;
      ssV_2 << "Priority (" << priority << ") not allowed - value has to be greater than '0'.";
      throw FlxException("FunReadBase::insert_1", ssV_2.str(), "This error is based on faulty source code.");
    }
  } else {
    Next = Next->insert(TheFunRead);
    return this;
  }
}

vector< FunBase*>* FunReadFunBase::read_parameters(const int NumbOfPara, const bool errSerious) {
  vector< FunBase*>* paraL = new vector< FunBase*>(0);
  bool bnf = false;
  try {
    if ( reader->whatIsNextChar() != ')' ) {
      do {
        if ( bnf ) { 
          reader->getChar(',',errSerious); 
        } else {
          bnf = true;
        }
        paraL->push_back( FunctionList->read(errSerious) );
      } while ( reader->whatIsNextChar() == ',' );
    }

    if ( NumbOfPara >= 0 && paraL->size() != tuint(NumbOfPara) ) {
      ostringstream ssV_2;
      ssV_2 << "Function expects " << NumbOfPara << " parameters, and not " << paraL->size() << ".";
      FlxError(errSerious,"FunReadFunBase::read_parameters_1", ssV_2.str(), reader->getCurrentPos());
    }
    
    if (NumbOfPara==0) {
      paraL->push_back(new FunDummy());
    }
    
  } catch (FlxException& e) {
    FLXMSG("FunReadFunBase::read_parameters_2",1);
    for (size_t i = 0; i<paraL->size();++i) {
      delete paraL->operator[](i);
    }
    delete paraL;
    throw;
  }

  return paraL;
}

tdoublePtr FunReadFunBase::read_const_var(const bool errSerious,const bool define)
{
  const string strV = reader->getWord(true,errSerious);
  tdouble* d1 = ConstantBox->get(strV,define);
  if ( d1 == NULL ) {
    ostringstream ssV_2;
    ssV_2 << "Const-variable '" << strV << "' does not exist.";
    FlxError(errSerious,"FunReadFunBase::read_const_var", ssV_2.str(), reader->getCurrentPos());
  }
  return d1;
}

FunReadFunUser::~FunReadFunUser() {
  delete fun; 
}



FlxFunctionReader::FlxFunctionReader (  ) {
  #if FLX_DEBUG
    ++Cinst;
    if (Cinst > 1) {
      std::ostringstream ssV;
      ssV << "More than one instance of 'FlxFunctionReader' created ...";
      throw FlxException("FlxFunctionReader::FlxFunctionReader", ssV.str() );
    }
  #endif
  
  FunctionList = new FunReadSTART();

  // thereafter all other FunReadXXXXXX Objects have to be listed (except 'FunReadSTART' and 'FunReadEND')
  FunctionList->insert( new FunReadNumber() );                // 1000
  FunctionList->insert( new FunReadPara() );                // 1020
  
  flxfunction_ope_insert(FunctionList);
}

FunBase* FlxFunctionReader::read (bool errSerious) {
  FunBase* rtn = NULL;
  try {
    rtn = FunctionList->read( errSerious );
    FunBasePtr optf = NULL;
    while (rtn->optimize(optf,Fun_OptimizeInfo())) {
      if (optf) {
        delete rtn;
        rtn = optf;
        optf = NULL;
      }
    }  
  } catch ( FlxException &e) { 
    FLXMSG("FlxFunctionReader::read_1",1);
    if (rtn) delete rtn;
    throw;
  }
  return rtn;
}

ReadStream* FlxFunctionReader::get_reader()
{
  return FunctionList->get_reader();
}


FunBase* FunReadEND::read (bool errSerious) {
  ostringstream ssV_2;
  char c = reader->getChar(false);
  ssV_2 << "Character (" << c << ")[" << int(c) << "] not expected at this point.";
  FlxError(errSerious,"FunReadEND::read_1", ssV_2.str(), reader->getCurrentPos()); return NULL; // dummy-return
}

FunReadSTART::FunReadSTART() : FunReadBase(0) {
  startLink = this;
  FunReadFunBase_0::set_data(this);
  this->Next = new FunReadEND();
}

FunBase* FunReadSTART::read (bool errSerious) {
  return Next->read( errSerious );
}

ReadStream* FunReadSTART::get_reader()
{
  return reader;
}

FunBase* FunReadNumber::read (bool errSerious) {
  if ( reader->nextCanBeNumber() ) {
    bool is_minus = false;
    if (reader->whatIsNextChar()=='-') {
      reader->getChar();
      is_minus = true;
      if ( reader->nextCanBeNumber() == false ) {
        return new FunMult(new FunNumber(-ONE), Next->read(errSerious));
      }
    }
    tdouble d1 = reader->get_Double(errSerious);
    if (is_minus) d1 = -d1;
    return new FunNumber(d1);
  }
  else return Next->read( errSerious );
}

FunBase* FunReadPara::read(bool errSerious) {
  if ( reader->whatIsNextChar() == '$' ) {
    reader->getChar();
    // parameter of functions
      if ( reader->nextCanBeNumber() ) {
        tuint i1 = reader->get_UInt<tuint>(errSerious);
        if ( numbofpara == 0 ) {
          ostringstream ssV;
          ssV << "Parameters of type '$" << i1 << "' are only valid within functions.";
          FlxError(errSerious,"FunReadPara::read_1", ssV.str(), reader->getCurrentPos());
        }
        if ( i1 < 1 || numbofpara < i1 ) {
          ostringstream ssV;
          ssV << "Index '" << i1 << "' out of range [0; " << numbofpara << "].";
          FlxError(errSerious,"FunReadPara::read_2", ssV.str(), reader->getCurrentPos());
        }
        return new FunPara(i1);
      }
    // not allowed
      ostringstream ssV_2;
      ssV_2 << "Character not expeted ('" << reader->whatIsNextChar() << "').";
      FlxError(errSerious,"FunReadPara::read_3", ssV_2.str(), reader->getCurrentPos()); return NULL; // dummy-return
  } else return Next->read( errSerious );
}

void FunReadPara::set_NumbOfPara(const tuint NumbOfPara) {
  numbofpara = NumbOfPara; 
}




const tnlong tnlong_from(const tdouble d, const string Descr, const bool errSerious, const bool zero_is_allowed, FunBase* fun)
{
  const tdouble rd  = round_flx(d);
  if ( (zero_is_allowed && rd<ZERO) || (!zero_is_allowed && rd<= ZERO) ) {
    std::ostringstream ssV;
    ssV << Descr << " must not be a negative number (" << d << "->" << rd << ").";
    if (fun) {
      ssV << " The problem occurred in function: " << fun->write();
    }
    FlxError(errSerious,"tnlong_from", ssV.str() );
  }
  return tnlong(rd);
}


const tulong tulong_from(const tdouble d, const string Descr, const bool errSerious, const bool zero_is_allowed, FunBase* fun)
{
  const tdouble rd  = round_flx(d);
  if ( (zero_is_allowed && rd<ZERO) || (!zero_is_allowed && rd<= ZERO) ) {
    std::ostringstream ssV;
    ssV << Descr << " must not be a negative number (" << d << "->" << rd << ").";
    if (fun) {
      ssV << " The problem occurred in function: " << fun->write();
    }
    FlxError(errSerious,"tulong_from", ssV.str() );
  }
  return tulong(rd);
}

const tuint tuint_from(const tdouble d, const std::string Descr, const bool errSerious, const bool zero_is_allowed, FunBase* fun)
{
  const tdouble rd  = round_flx(d);
  if ( (zero_is_allowed && rd<ZERO) || (!zero_is_allowed && rd<= ZERO) ) {
    std::ostringstream ssV;
    ssV << Descr << " must not be a negative number (" << d << "->" << rd << ").";
    if (fun) {
      ssV << " The problem occurred in function: " << fun->write();
    }
    FlxError(errSerious,"tuint_from", ssV.str() );
  }
  return tuint(rd);
}

FlxFunDeg::FlxFunDeg(ReadStream* reader, FlxFunctionReader* funReader, bool errSerious)
:deg(0),fun(NULL)
{
  FlxFunction* tf = NULL;
  try {
    if (reader->whatIsNextChar()=='[') {
      reader->getChar('[',errSerious);
      fun = new FlxFunction(funReader,errSerious);
      reader->getChar(',',errSerious);
      tf = new FlxFunction(funReader,errSerious);
      deg = tf->cast2tuintW0(errSerious);
      delete tf; tf = NULL;
      reader->getChar(']',errSerious);
    } else {
      deg = default_deg;
      fun = new FlxFunction(funReader,errSerious);
    }
  } catch (FlxException& e) {
    FLXMSG("FlxFunDeg::FlxFunDeg",1);
    if (fun) delete fun;
    if (tf) delete tf;
    throw;
  }
}

FlxFunDeg::FlxFunDeg(FlxFunction* funF)
:deg(default_deg),fun(funF)
{
  
}

FlxFunDeg::FlxFunDeg(const FlxFunDeg& ref)
{
  deg = ref.deg;
  fun = new FlxFunction(*ref.fun);
}

const std::string FlxFunDeg::write()
{
  std::ostringstream ssV;
  ssV << "[" << fun->write() << "," << deg << "]";
  return ssV.str();
}

void FlxFunDeg::assign(FlxFunDeg* funA)
{
  if (this == funA) return;
  deg = funA->deg;
  if (this->fun != funA->fun) {
    delete fun;
    fun = funA->fun;
    funA->fun = NULL;
  }
  delete funA;
}

