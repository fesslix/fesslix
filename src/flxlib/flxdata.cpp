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

#define fesslix_fxldata_CPP

#include "flxdata.h"

#include <fstream>

using namespace std;

FlxData* FlxDataBase::data = nullptr;
#if FLX_DEBUG
  int FlxData::Cinst = 0;
  int FlxReadManager::Cinst = 0;
  int FlxTimerBox::Cinst = 0;
  int FlxIgnoreBox::Cinst = 0;
  int FlxSubBox::Cinst = 0;
  int FlxOstreamBox::Cinst = 0;
  int FlxIstreamBox::Cinst = 0;
#endif

FlxReadManager* readManager_ptr = nullptr;

void set_ReadManager(FlxReadManager* readManager_ptr_)
{
  readManager_ptr = readManager_ptr_;
}


FlxReadManager::FlxReadManager()
{
  #if FLX_DEBUG
    ++Cinst;
    if (Cinst > 1) {
      std::ostringstream ssV;
      ssV << "More than one instance of 'FlxReadManager' created ...";
      throw FlxException("FlxReadManager::FlxReadManager", ssV.str() );
    }
  #endif
}

void FlxReadManager::push(ReadStream* readerV)
{
  if (readerV == NULL) {
    ostringstream ssV_2;
    ssV_2 << "Empty Reader.";
    throw FlxException("FlxReadManager::push_1", ssV_2.str() ); 
  }

  s.push(readerV);
  reader = readerV;
}

void FlxReadManager::pop()
{
  if (s.empty()) {
    ostringstream ssV;
    ssV << "No reader on the stack.";
    throw FlxException("FlxReadManager::pop_1", ssV.str() ); 
  }
  s.pop();
  if ( s.empty() ) {
    reader = NULL;
  } else {
    reader = s.top();
  }
}

FlxFunction* FlxReadManager::parse_function(const std::string& funStr)
{
  ReadStream* rs = new ReadStream(std::string(funStr));
  push(rs);
  FlxFunction* value = NULL;
  try {
    value = new FlxFunction(funReader);        
  } catch (FlxException &e) {
    FLXMSG("FlxReadManager::parse_function_1",1);
    pop();
    delete rs;
    throw;
  }
  pop();
  delete rs;
  return value;
}

FlxFunction * FlxReadManager::parse_function(py::object pyobj)
{
  // ==================================================
  // float
  // ==================================================
  if (py::isinstance<py::float_>(pyobj)) {
    const tdouble val = py::cast<tdouble>(pyobj);
    return new FlxFunction(new FunNumber(val));
  }
  // ==================================================
  // function
  // ==================================================
  if (py::isinstance<py::function>(pyobj)) {
      try {
          // Get the function's __code__ object
          py::object code_obj = pyobj.attr("__code__");
          int arg_count = code_obj.attr("co_argcount").cast<int>();

          // Check if it takes zero arguments
          if (arg_count == 0) {
              py::function pyfunc = py::reinterpret_borrow<py::function>(pyobj);
              vector< FunBase*>* paraL = new vector< FunBase*>(0);
              return new FlxFunction(new FunBaseFun_Python("INTERNAL_CALLABLE", pyfunc, paraL));
          } else {
              std::ostringstream ssV;
              ssV << "The parameter is a callable function, but it requires " << arg_count << " parameters.";
              throw FlxException("FlxReadManager::parse_function_20", ssV.str());
          }
      } catch (const py::error_already_set &e) {
          std::ostringstream ssV;
          ssV << "Error retrieving function signature: " << e.what();
          throw FlxException("FlxReadManager::parse_function_29", ssV.str() );
      }
  }
  // ==================================================
  // string
  // ==================================================
  if (py::isinstance<py::str>(pyobj)) {
    const std::string val = py::cast<std::string>(pyobj);
    return parse_function(val);
  }
  // ==================================================
  // otherwise
  // ==================================================
  throw FlxException("FlxReadManager::parse_function_99", "Unhandled data type of Python object");
}

FlxFunction * parse_function(const std::string& funStr)
{
  return readManager_ptr->parse_function(funStr);
}

FlxFunction * parse_function(py::object pyobj)
{
  return readManager_ptr->parse_function(pyobj);
}


FlxData::FlxData()
{ 
  #if FLX_DEBUG
    ++Cinst;
    if (Cinst > 1) {
      std::ostringstream ssV;
      ssV << "More than one instance of 'FlxData' created ...";
      throw FlxException("FlxData::FlxData", ssV.str() );
    }
  #endif
  FlxDataBase::set_data(this);
  FlxString::set_StrFunBox(&StrFunBox);
}

FlxData::~FlxData()
{
  
}

const ostreamp& FlxOstreamBox::get ( const string& name ) {
  map<std::string, ostreamp*>::iterator pos = box.find(name);
  if ( pos != box.end() ) {
    return (*pos->second);
  } else {
    ostringstream ssV_2;
    ssV_2 << "The output-stream '" << name << "' does not exist.";
    throw FlxException("FlxOstreamBox::get_1", ssV_2.str(), "In oder to use an output-stream, you have to define it first."); 
  }
}

const bool FlxOstreamBox::delete_stream(ostreamp& strm)
{
  // if stream is the logStream
    if ( strm == GlobalVar.get_log() ) return false;
  // if stream is a FileStream
      ofstream* thestream = dynamic_cast<ofstream*> (strm);
      if ( thestream ) {
        thestream->close();
        delete thestream;
        return true;
      }
    // if stream is a DistributorStream
      flxStreamAlloc* thestreamD = dynamic_cast<flxStreamAlloc*> (strm);
      if ( thestreamD ) {
        delete thestreamD;
        return true;
      }
    // if stream is a DummyStream
      flxDummyOstream* thestreamDu = dynamic_cast<flxDummyOstream*> (strm);
      if ( thestreamDu ) {
        delete thestreamDu;
        return true;
      }
    // if stream is a StringStream
      std::ostringstream* thestreamSS = dynamic_cast<std::ostringstream*> (strm);
      if ( thestreamSS ) {
        delete thestreamSS;
        return true;
      }
   return false;
}

FlxOstreamBox::FlxOstreamBox()
{
  #if FLX_DEBUG
    ++Cinst;
    if (Cinst > 1) {
      std::ostringstream ssV;
      ssV << "More than one instance of 'FlxOstreamBox' created ...";
      throw FlxException("FlxOstreamBox::FlxOstreamBox", ssV.str() );
    }
  #endif
  GlobalVar.alert.set_err_stream(GlobalVar.stdcerr);
}

FlxOstreamBox::~FlxOstreamBox() {
  for (map<std::string, ostreamp*>::iterator pos = box.begin(); pos != box.end(); ++pos) {
    if (
      pos->first!="cout" &&
      pos->first!="cerr" &&
      pos->first!="log" &&
      *(pos->second)!=GlobalVar.get_true_cout() && 
      *(pos->second)!=GlobalVar.get_true_cerr() && 
      *(pos->second)!=GlobalVar.get_cout() && 
      *(pos->second)!=GlobalVar.get_cerr() &&
      *(pos->second)!=GlobalVar.get_log() &&
      *(pos->second)!=&GlobalVar.slogcout(0)
    ) delete_stream(*(pos->second));
    if ( 
      pos->second!=&GlobalVar.get_cout() &&
      pos->second!=&GlobalVar.get_cerr() &&
      pos->second!=&GlobalVar.get_true_cerr() &&
      pos->second!=&GlobalVar.get_log()
    )
    delete pos->second;
  }
}

const bool FlxOstreamBox::insert( const std::string& name, ostreamp value) {
  map<std::string, ostreamp*>::iterator pos = box.find(name);
  if ( pos != box.end() ) {
    if (name=="cout" || name=="log" || name=="cerr" ) {
      ostringstream ssV;
      ssV << "The output-stream '" << name << "' is not allowed to be redefined.";
      throw FlxException("FlxOstreamBox::insert_1", "Stream already exists", ssV.str() ); 
    }
    ostreamp& ot = *(pos->second);
    if (delete_stream(ot)==false) {
      ostringstream ssV;
      ssV << "The output-stream '" << name << "' was already defined; it could not be closed.";
      throw FlxException("FlxOstreamBox::insert_2", "Stream already exists", ssV.str() ); 
    }
    ot = value;
    return false;
  } else {
    ostreamp* valuep = NULL;
    // check for some special cases
      if (name=="cout") {
        if ( value == GlobalVar.get_cout() ) {
          valuep = &GlobalVar.stdcout;
        }
      } else if (name=="cerr") {
        if ( value == GlobalVar.get_cerr() ) {
          valuep = &GlobalVar.stdcerr;
        }
      } else if (name=="log") {
        if ( value == GlobalVar.get_log() ) {
          valuep = &GlobalVar.slogP;
        }
      }
    // the not-so-special case
      if (valuep==NULL) {
        valuep = new ostreamp(value);
      }
    pair<std::string, ostreamp*> Element(name, valuep);
    box.insert(Element);
    return true;
  }
}

void FlxOstreamBox::close( const std::string& name, const bool err ) 
{
  map<std::string, ostreamp*>::iterator pos = box.find(name);
  if ( pos != box.end() ) {
    ostreamp& ot = *(pos->second);
    if ( (dynamic_cast<flxDummyOstream*> (ot)) != NULL) {
      if (err) {
        ostringstream ssV;
        ssV << "'" << name << "' is already closed.";
        GlobalVar.alert.alert("FlxOstreamBox::close_1",ssV.str());
      }
      return;
    }
    if (name!="cout" && name!="log" && name!="cerr" && delete_stream(ot)) {
      ot = new flxDummyOstream();
    } else {
      ostringstream ssV;
      ssV << "'" << name << "' cannot be closed.";
      GlobalVar.alert.alert("FlxOstreamBox::close_2",ssV.str());
    }
  } else {
    if (err) {
      ostringstream ssV;
      ssV << "The output-stream '" << name << "' does not exist.";
      throw FlxException("FlxOstreamBox::close_3", "Stream does not exist", ssV.str()); 
    }
  }
}

void FlxIstreamBox::insert(const string& name, FlxIstream* value, bool errSerious)
{
  pair<std::string, FlxIstream*> Element(name, value);
  map<std::string, FlxIstream*>::iterator pos = box.find(name);
  if ( pos != box.end() ) {
    pos->second->copyStream(value,errSerious);
    return;
  }
  box.insert(Element);
}

FlxIstream& FlxIstreamBox::get(const string& name)
{
  map<std::string, FlxIstream*>::iterator pos = box.find(name);
  if ( pos != box.end() ) {
    return (*pos->second);
  } else {
    ostringstream ssV_2;
    ssV_2 << "The input-stream '" << name << "' does not exist.";
    throw FlxException("FlxIstreamBox::get_1", ssV_2.str(), "In oder to use an input-stream, you have to define it first."); 
  }
}

FlxIstream_vector* FlxIstreamBox::get_isVector(const string& name)
{
  map<std::string, FlxIstream*>::iterator pos = box.find(name);
  if (pos == box.end()) return NULL;
  FlxIstream* tp = pos->second;
  return dynamic_cast<FlxIstream_vector*> (tp);
}

FlxIstreamBox::FlxIstreamBox()
{
  #if FLX_DEBUG
    ++Cinst;
    if (Cinst > 1) {
      std::ostringstream ssV;
      ssV << "More than one instance of 'FlxIstreamBox' created ...";
      throw FlxException("FlxIstreamBox::FlxIstreamBox", ssV.str() );
    }
  #endif
}

FlxIstreamBox::~FlxIstreamBox()
{
  for (map<std::string, FlxIstream*>::iterator pos = box.begin(); pos != box.end(); ++pos) {
    delete pos->second;
  }
}


FlxSubBox::FlxSubBox()
{
  #if FLX_DEBUG
    ++Cinst;
    if (Cinst > 1) {
      std::ostringstream ssV;
      ssV << "More than one instance of 'FlxSubBox' created ...";
      throw FlxException("FlxSubBox::FlxSubBox", ssV.str() );
    }
  #endif
}

FlxSubBox::~FlxSubBox() {
  for (map<std::string, FlxObjBase*>::iterator pos = box.begin(); pos != box.end(); ++pos) delete pos->second;
}

void FlxSubBox::insert(const std::string& name, FlxObjBase* value) {
  pair<std::string, FlxObjBase*>Element(name, value);
  if ( ! box.insert(Element).second ) {
    delete value;
    ostringstream ssV;
    ssV << "Procedure '" << name << "' is already defined.";
    throw FlxException("FlxSubBox::insert_1", ssV.str() );
  }
}

FlxObjBase* FlxSubBox::get(const string& name) {
  map<std::string, FlxObjBase*>::iterator pos = box.find(name);
  if ( pos != box.end() ) {
    return pos->second;
  } else {
    return NULL;
  }
}


FlxIgnoreBox::FlxIgnoreBox(): iL_recur_active(0) 
{
  #if FLX_DEBUG
    ++Cinst;
    if (Cinst > 1) {
      std::ostringstream ssV;
      ssV << "More than one instance of 'FlxIgnoreBox' created ...";
      throw FlxException("FlxIgnoreBox::FlxIgnoreBox", ssV.str() );
    }
  #endif
}

const bool FlxIgnoreBox::isOnIgnoreList_recur( const string& objName ) {
  vector<string>::const_iterator iter = find ( ignoreList_recur.begin(), ignoreList_recur.end(), objName);
  return (iter != ignoreList_recur.end() )?true:false; 
}



void FlxTimerBox::insert(const string& name, FlxTimer* value) {
  pair<std::string, FlxTimer*>Element(name, value);
  if ( ! box.insert(Element).second ) {;
    delete value;
    ostringstream ssV;
    ssV << "Timer '" << name << "' is already defined.";
    throw FlxException("FlxTimer::insert_1", ssV.str() );
  }
}


FlxTimer* FlxTimerBox::get(const string& name) {
  map<std::string, FlxTimer*>::iterator pos = box.find(name);
  if ( pos != box.end() ) {
    return pos->second;
  } else {
    ostringstream ssV;
    ssV << "Timer '" << name << "' does not exist.";
    throw FlxException("FlxTimer::get_1", ssV.str() );
  }
}


void FlxTimerBox::deleteEl(const string& name)
{
  FlxTimer* tt = get(name);
  delete tt;
  box.erase(name);
}

FlxTimerBox::FlxTimerBox()
{
  #if FLX_DEBUG
    ++Cinst;
    if (Cinst > 1) {
      std::ostringstream ssV;
      ssV << "More than one instance of 'FlxTimerBox' created ...";
      throw FlxException("FlxTimerBox::FlxTimerBox", ssV.str() );
    }
  #endif
}

FlxTimerBox::~FlxTimerBox() {
  for (map<std::string, FlxTimer*>::iterator pos = box.begin(); pos != box.end(); ++pos) delete pos->second;
}

void FlxObjBase::exec() {
  if (NOTdolog) {
    NOTdolog_start();
    try {
      task();
      NOTdolog_stop();
    } catch (...) {
      FLXMSG("FlxTimerBox::~FlxTimerBox",1);
      NOTdolog_stop();
      throw;
    }
  } else {
    task();
  }
  #if FLX_DEBUG
    GlobalVar.slog_flush();
  #endif
  if ( Next ) {
    Next->exec();
  }
}

void FlxObjBase::attach_obj(FlxObjBase* NextV) {
  if ( Next == NULL ) {
    Next = NextV;
  } else {
    Next->attach_obj(NextV);
  }
}


FlxObjBase::~FlxObjBase() {
  if (Next) delete Next;
}


void FlxCodeBlock::exec()
{
  try {
    loop_block_exec_1();
  } catch (FlxReturnBreakContinue_baseE &e) {
    loop_block_exec_2();
    throw;
  }
  loop_block_exec_2();
}

void FlxCodeBlock::loop_block_exec_1()
{
  // remember the value of the variables marked as 'internal'
    const size_t N = cvec.size();
    if (N>0 && dvec.size()==0) dvec.resize(N);
    for (size_t i=0;i<N;++i) dvec[i] = *(cvec[i]);
  // executed all objects inside this code-block
    try {        // catch return-statement
        FlxObjBase::exec();
    } catch ( FlxReturnE &e) { 
      if (!catch_return) throw;
    } catch ( FlxContinueE &e) { 
      if (!catch_continue) throw;
    }
}

void FlxCodeBlock::loop_block_exec_2()
{
  const size_t N = cvec.size();
  // restore the value of all variables marked as 'internal'
    for (size_t i=0;i<N;++i) *(cvec[i]) = dvec[i];
}


void FlxCodeBlock::add_internal_const(tdouble* icp)
{
  bool add = true;
  for (size_t i=0;i<cvec.size();++i) {
    if (cvec[i]==icp) {
      add = false;
      break;
    }
  }
  if (add) {
    cvec.push_back(icp);
  }
}

void flxStrConstBox::insert( const std::string& name, const std::string& value)
{
  std::map<std::string, std::string>::iterator pos = box.find(name);
  if ( pos != box.end() ) {
    pos->second = value;
    return;
  }
  std::pair<std::string,std::string> Element(name,value);
  box.insert(Element);
}

const std::string& flxStrConstBox::get(const std::string& name)
{
  std::map<std::string,std::string>::iterator pos = box.find(name);
  if ( pos != box.end() ) {
    return (pos->second);
  } else {
    std::ostringstream ssV_2;
    ssV_2 << "The string-constant '" << name << "' does not exist.";
    throw FlxException("flxStrConstBox::get_1", ssV_2.str(), "In oder to use an string-constant, you have to define it first."); 
  }
}

std::string& flxStrConstBox::get_ref(const std::string& name)
{
  std::map<std::string,std::string>::iterator pos = box.find(name);
  if ( pos != box.end() ) {
    return (pos->second);
  } else {
    return box[name];
  }
}

