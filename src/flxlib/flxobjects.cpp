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

#include "flxobjects.h"
#include "flxobjcommon.h"

using namespace std;

FlxDefParaBox* FlxDefParaBoxBase::AllDefParaBox = NULL;
FlxEvaluateCmd* FlxEvaluateCmdBase::EvaluateCmd = NULL;
#if FLX_DEBUG
  int FlxDefParaBox::Cinst = 0;
  int FlxEvaluateCmd::Cinst = 0;
#endif

void FlxOptionalParaBase::read(bool errSerious) {
  set( read_value(errSerious) );
  is_set = true;
}

void FlxOptionalParaFun::set_default(void* defV)
{
  FlxFunction* valueV = static_cast<FlxFunction*>(defV);
  if (defv!=NULL) {
    delete defv;
  }
  defv = new FlxFunction(*valueV); 
  GlobalVar.slog(3) << "default: set '" << pName << "' to '" << defv->write() << "'." << std::endl;
}

FlxOptionalParaFun::FlxOptionalParaFun(const tdouble defV, const string pName)
: FlxOptionalParaBase(pName), defv(NULL), value(NULL)
{
  defv = new FlxFunction(new FunNumber(defV));
}

FlxOptionalParaFun::~FlxOptionalParaFun()
{
  if ( value != NULL ) {
    delete value; 
  }
  if ( defv != NULL ) {
    delete defv;
  }
}

void FlxOptionalParaFun::free_value(void* valueP)
{
  FlxFunction* valueV = static_cast<FlxFunction*>(valueP);
  if ( valueP ) { delete valueV; }
}

void FlxOptionalParaFun::set(void* valueP)
{
  FlxFunction* valueV = static_cast<FlxFunction*>(valueP);
  if ( value==NULL) {
    value = valueV;
  } else {
    value->assign(valueV);
  }
}

void* FlxOptionalParaFun::get()
{
  if ( value == NULL || is_set == false ) {
    is_set = false;
    return new FlxFunction(*defv);
  } else {
    is_set = false;
    return new FlxFunction(*value);
  }
}

FlxFunction* FlxOptionalParaFun::get_ref()
{
  if ( value == NULL || is_set == false ) {
    is_set = false;
    return defv;
  } else {
    is_set = false;
    return value;
  }
}


void FlxOptionalParaFlxString::set(void* valueP)
{
  FlxString* valueV = static_cast<FlxString*>(valueP);
  if ( value==NULL) {
    value = valueV;
  } else {
    value->assign(valueV);
  }
}

void* FlxOptionalParaFlxString::get()
{
  if ( value == NULL || is_set == false ) {
    is_set = false;
    return new FlxString(*defv);
  } else {
    is_set = false;
    return new FlxString(*value);
  }
}

FlxString* FlxOptionalParaFlxString::get_ref()
{
  if ( value == NULL || is_set == false ) {
    is_set = false;
    return defv;
  } else {
    is_set = false;
    return value;
  }
}

void FlxOptionalParaFlxString::set_default(void* defV)
{
  FlxString* valueV = static_cast<FlxString*>(defV);
  if (defv) {
    delete defv;
  }
  defv = new FlxString(*valueV); 
  GlobalVar.slog(3) << "default: set '" << pName << "' to '" << defv->write() << "'." << std::endl;
}

FlxOptionalParaFlxString::FlxOptionalParaFlxString(string defV, string pName, const bool is_Word)
: FlxOptionalParaBase(pName), defv(new FlxString(new FlxString_String(defV,is_Word),false)), value(NULL)
{

}

FlxOptionalParaFlxString::~FlxOptionalParaFlxString()
{
  if ( value ) {
    delete value; 
  }
  if ( defv ) {
    delete defv;
  }
}

void FlxOptionalParaFlxString::free_value(void* valueP)
{
  FlxString* valueV = static_cast<FlxString*>(valueP);
  if ( valueP != NULL) { delete valueV; }
}

void FlxOptionalParaMtxFun::set_default(void* defV)
{
  FlxMtxConstFun* valueV = static_cast<FlxMtxConstFun*>(defV);
  if (defv) {
    delete defv;
  }
  defv = new FlxMtxConstFun(*valueV);
  GlobalVar.slog(3) << "default: set '" << pName << "' to '...'." << std::endl;
}

void* FlxOptionalParaMtxFun::get()
{
  if ( value == NULL || is_set == false ) {
    is_set = false;
    return new FlxMtxConstFun(*defv);
  } else {
    is_set = false;
    return new FlxMtxConstFun(*value);
  }
}

void FlxOptionalParaMtxFun::free_value(void* valueP)
{
  FlxMtxConstFun* valueV = static_cast<FlxMtxConstFun*>(valueP);
  if ( valueP ) { delete valueV; }
}

FlxOptionalParaMtxFun::~FlxOptionalParaMtxFun()
{
  if ( value ) {
    delete value; 
  }
  if ( defv ) {
    delete defv;
  }
}

void FlxOptionalParaMtxFun::set(void* valueP)
{
  FlxMtxConstFun* valueV = static_cast<FlxMtxConstFun*>(valueP);
  if ( value==NULL) {
    value = valueV;
  } else {
    value->assign(valueV);
  }
}


FlxOptionalParaStream::~FlxOptionalParaStream() {
  if ( value ) 
    delete value;
}

void* FlxOptionalParaStream::get() {
  if ( value == NULL || is_set == false ) {
    is_set = false;
    return &defv;
  } else {
    is_set = false;
    return value;
  }
}

void* FlxOptionalParaStream::read_value(const bool errSerious) {
  FlxString tmpStr(false,errSerious);
  return ( new string(tmpStr.eval_word(true)) );
}

void FlxOptionalParaStream::set(void* valueP) {
  string* valueV = static_cast<string*>(valueP);
  if ( value == NULL ) {
    value = new string( *valueV );
  } else {
    *value = *valueV;
  }
  delete valueV;
}

void FlxOptionalParaStream::free_value(void* valueP) {
  string* valueV = static_cast<string*>(valueP);
  if ( valueP != NULL) { delete valueV; }
}


void FlxOptionalParaStream::set_default(void* defV) {
  defv = *(static_cast<string*>(defV));
  transform(defv.begin(), defv.end(), defv.begin(), (int(*)(int)) std::tolower);
  GlobalVar.slog(3) << "default: set '" << pName << "' to '" << defv << "'." << std::endl;
}

void* FlxOptionalParaString::read_value(const bool errSerious) {
  FlxString tmpStr(false,errSerious);
  return ( new string(tmpStr.eval_word(false)) );
}

void FlxOptionalParaString::set_default(void* defV) {
  defv = *(static_cast<string*>(defV));
  GlobalVar.slog(3) << "default: set '" << pName << "' to '" << defv << "'." << std::endl;
}

void* FlxOptionalParaText::read_value(bool errSerious)
{
  FlxString tmpStr(false,errSerious);
  return ( new string(tmpStr.eval(false)) );
}

FlxOptionalParaBool::~FlxOptionalParaBool() {
  if ( value != NULL ) {
    delete value;
  }
}

void* FlxOptionalParaBool::read_value(bool errSerious) {
  FlxFunction* t = new FlxFunction(funReader,errSerious);
  const bool bt = (t->cast2int()!=0);
  delete t;
  return ( new bool(bt) );
}

void FlxOptionalParaBool::set(void* valueP) {
  bool* valueV = static_cast<bool*>(valueP);
  if ( value == NULL ) {
    value = new bool( *valueV );
  } else {
    *value = *valueV;
  }
  delete valueV;
}

void* FlxOptionalParaBool::get() {
  if ( value == NULL || is_set == false ) {
    is_set = false;
    return &defv;
  } else {
    is_set = false;
    return value;
  }
}

void FlxOptionalParaBool::free_value(void* value) {
  bool* valueV = static_cast<bool*>(value);
  if ( value != NULL) { delete valueV; }
}





void FlxOptionalParaBox::insert(string name, string defname) {
  transform(name.begin(), name.end(), name.begin(), (int(*)(int)) std::tolower);
  transform(defname.begin(), defname.end(), defname.begin(), (int(*)(int)) std::tolower);
  FlxOptionalParaBase* para = AllDefParaBox->get(defname);
  #if FLX_DEBUG
    if (para==NULL) {
      throw FlxException_Crude("FlxOptionalParaBox::insert");
    }
  #endif
  pair<std::string, FlxOptionalParaBase*>Element(name, para);
  if ( ! box.insert(Element).second ) {
    ostringstream ssV;
    ssV << "Optional parameter '" << name << "(" << defname << ")' could not be inserted into ParaBox.";
    throw FlxException("FlxOptionalParaBox::insert_1", ssV.str() );
  }
}

FlxOptionalParaBase* FlxOptionalParaBox::get(std::string name) {
  transform(name.begin(), name.end(), name.begin(), (int(*)(int)) std::tolower);
  map<std::string, FlxOptionalParaBase*>::iterator pos;
  pos = box.find(name);
  if ( pos != box.end() ) {
    if (pos->second==NULL) {
      throw FlxException_Crude("FlxOptionalParaBox::get_1");
    }
    return pos->second;
  } else {
    std::ostringstream ssV;
    ssV << "An optional parameter with name '" << name << "' does not exist.";
    throw FlxException("FlxOptionalParaBox::get_2", ssV.str() );
  }
}

void FlxOptionalParaBool::set_default(void* defV) {
  defv = *(static_cast<bool*>(defV));
  GlobalVar.slog(3) << "default: set '" << pName << "' to '" << (defv?"true":"false") << "'." << std::endl;
}

FlxDefParaBox::FlxDefParaBox()
{
  #if FLX_DEBUG
    ++Cinst;
    if (Cinst > 1) {
      std::ostringstream ssV;
      ssV << "More than one instance of 'FlxDefParaBox' created ...";
      throw FlxException("FlxDefParaBox::FlxDefParaBox", ssV.str() );
    }
  #endif
  FlxDefParaBoxBase::set_AllDefParaBox(this);
}

FlxDefParaBox::~FlxDefParaBox() {
  for (map<std::string, FlxOptionalParaBase*>::iterator pos = box.begin(); pos != box.end(); ++pos) delete pos->second;
}

void FlxDefParaBox::insert(FlxOptionalParaBase* value) {
  pair<std::string, FlxOptionalParaBase*>Element(value->get_name(), value);
  if ( ! box.insert(Element).second ) {;
    delete value;
  }
}

FlxOptionalParaBase* FlxDefParaBox::get(std::string name) {
  transform(name.begin(), name.end(), name.begin(), (int(*)(int)) std::tolower);
  map<std::string, FlxOptionalParaBase*>::iterator pos = box.find(name);
  if ( pos != box.end() ) {
    return pos->second;
  } else {
    return NULL;
  }
}



void FlxEvaluateCmd::check_ending() {
  char ch = reader->getChar();
  if ( ch != ';' ) {
    if ( int(ch) != -1 ) {         // cheap and nasty - this produces a bug. (';' is not necessary at the end of the file)
      ostringstream ssV;
      ssV << "Expected ';' (and NOT '" << ch << "' [" << int(ch) << "])";
      throw FlxException_NeglectInInteractive("FlxEvaluateCmd::check_ending_1", ssV.str(), reader->getCurrentPos());
    }
  }
}

FlxEvaluateCmd::FlxEvaluateCmd()
{
  #if FLX_DEBUG
    ++Cinst;
    if (Cinst > 1) {
      std::ostringstream ssV;
      ssV << "More than one instance of 'FlxEvaluateCmd' created ...";
      throw FlxException("FlxEvaluateCmd::FlxEvaluateCmd", ssV.str() );
    }
  #endif
  FlxEvaluateCmdBase::set_EvaluateCmd(this);   
}

FlxObjBase* FlxEvaluateCmd::evaluateCmd( ) {
  // Detect multiple ;;;'s
  while (reader->whatIsNextChar() == ';') {
    reader->getChar(true,false);
    if (reader->whatIsNextChar() == '}') {        // make sure we are not inside a block ...
      return new FlxObjDummy();
    }
  }
  if (reader->whatIsNextString(2)=="#!") {
    const std::string ts = reader->getNextLine();
    return new FlxObjEcho(true, new FlxString(new FlxString_String(ts,false),false), "log", true );
  }
  const string str1 = reader->getWord(true,false);
  FlxObjBase* obj1 = NULL;

  // Are there some objects which are not allowed ( for example "fun" in loops? ) 
    if ( data->IgnoreBox.isActive_iL_recur() ) {
      if ( data->IgnoreBox.isOnIgnoreList_recur(str1) ) {
        ostringstream ssV;
        ssV << "'" << str1 << "'not allowed at this point.";
        throw FlxException_NeglectInInteractive("FlxEvaluateCmd::evaluateCmd_1", ssV.str(), reader->getCurrentPos());
      }
    }

  // check if there is a keyword for an object
    FlxObjReadBase *objR = ObjReadBox.get(str1);
    if (objR != NULL) {
      obj1 = NULL;
      try {
        obj1 = objR->read();
        check_ending();
      } catch ( FlxException &e) {
        FLXMSG("FlxEvaluateCmd::evaluateCmd_2",1);
        if (obj1!=NULL) delete obj1;
        throw;
      };
      GlobalVar.prelog_write();
      return obj1;
    } 

  // check if the keyword refers to a procedure
  obj1 = data->SubBox.get(str1);
  if ( obj1 ) {
    FlxObjReadBase *objR = ObjReadBox.get("procedure");
    FlxObjReadSubDo* objRSD = dynamic_cast<FlxObjReadSubDo*>( objR );
    try {
      obj1 = objRSD->read(obj1, str1);
      check_ending();
    } catch ( FlxException &e) {
      FLXMSG("FlxEvaluateCmd::evaluateCmd_3",1);
      if (obj1) delete obj1;
      throw;
    };
    GlobalVar.prelog_write();
    return obj1;
  }

  // if no valid keyword was found, check if it is an assignment of a 'const' variable
    if ( reader->whatIsNextChar() != '=' ) {
      ostringstream ssV;
      ssV << "Expected keyword ('" << str1 << "' is not valid).";
      throw FlxException_NeglectInInteractive("FlxEvaluateCmd::evaluateCmd_4", ssV.str(), reader->getCurrentPos());
    } else {        // assign 'const' variables
      FlxObjReadConst* objRC = dynamic_cast<FlxObjReadConst*> ( ObjReadBox.get("const") );
      obj1 = objRC->read(str1);
      try {
        check_ending();
      } catch ( FlxException &e) {
        FLXMSG("FlxEvaluateCmd::evaluateCmd_4",1);
        delete obj1;
        throw;
      };
      GlobalVar.prelog_write();
      return obj1;
    }

}

void FlxObjOutputBase::write(const tdouble val, std::ostream& sout_)
{
  if (boost_str.empty()) {
    sout_ << GlobalVar.Double2String(val,checkTOL,prec,fixW);
  } else {
    try {
      sout_ << std::vformat(boost_str, std::make_format_args(val));
    } catch (...) {
      std::ostringstream ssV;
      ssV << "The output-format string '" << boost_str << "' caused an error.";
      throw FlxException_NeglectInInteractive("FlxObjOutputBase::write",ssV.str());
    }
  }
}



FlxObjReadBase* FlxObjectReadBox::get(std::string name) {
  map<string, FlxObjReadBase*>::iterator pos;
  std::transform(name.begin(), name.end(), name.begin(), (int(*)(int)) std::tolower);
  pos = box.find(name);
  if ( pos != box.end() ) {
    return pos->second;
  } else {
    return NULL;
  }
}

void FlxObjectReadBox::remove(const string& name)
{
  box.erase(name);
}

void FlxObjectReadBox::insert(std::string name, FlxObjReadBase* value) {
  pair<string, FlxObjReadBase*> Element(name, value);
  if ( ! box.insert(Element).second ) {
    ostringstream ssV;
    ssV << "Error during inserting " << name << " in FlxObjectReadBox.";
    throw FlxException("FlxObjectReadBox::insert_1", ssV.str());
  }
}

FlxObjectReadBox::~FlxObjectReadBox() {
  for (map<string, FlxObjReadBase*>::iterator pos = box.begin(); pos != box.end(); ++pos) {
    FlxObjReadBase* to = pos->second;
    if (to->is_from_library()==false) delete to;
  }
}

FlxObjectReadBox::FlxObjectReadBox() {
  insert("const", new FlxObjReadConst());
  insert("sub", new FlxObjReadSub());
  insert("procedure", new FlxObjReadSubDo());
  
}

FlxObjReadBase::FlxObjReadBase(const bool from_library)
: from_library(from_library)
{
  // Default do-log-setting
    AllDefParaBox->insert(new FlxOptionalParaBool(true,"flxlog::dolog"));
    ParaBox.insert("dolog", "flxlog::dolog" );
}

const bool FlxObjReadBase::get_doLog()
{
  return *(static_cast<bool*>(ParaBox.get("dolog")->get() ));
}

void FlxObjReadBase::read_optionalPara(bool errSerious) {
  if ( reader->whatIsNextChar() != '{' ) { return; }
  string strV;
  reader->getChar('{',errSerious);
  while ( reader->getNextType() == ReadStream::STRING ) {
    strV = reader->getWord(true,errSerious);
    FlxOptionalParaBase* para = ParaBox.get(strV);
    if ( para == NULL ) {
      ostringstream ssV;
      ssV << "Unknown optional Parameter '" << strV << "'.";
      FlxError(errSerious,"FlxObjReadBase::read_optionalPara_1", ssV.str());
    }
    reader->getChar('=',errSerious);
    para->read(errSerious);
    reader->getChar(';',errSerious);
  }
  reader->getChar('}',errSerious);
}

FlxFunction* FlxObjReadBase::get_optPara_FlxFunction(const string& str)
{
  void *v1 = ParaBox.get(str)->get();
  FlxFunction* fun = static_cast<FlxFunction*>( v1 ); 
  if (fun == NULL) {
    throw FlxException_Crude("FlxObjReadBase::get_optPara_FlxFunction");
  }
  return fun;
}

const tdouble FlxObjReadBase::get_optPara_tdouble_from_FlxFunction(const string& str, const bool only_positiveORzero, const bool errSerious)
{
  FlxOptionalParaBase* v1 = ParaBox.get(str);
  FlxOptionalParaFun* ffun = dynamic_cast<FlxOptionalParaFun*>( v1 );
  if (ffun==NULL) throw FlxException_Crude("FlxObjReadBase::get_optPara_tuint_from_FlxFunction");
  if (only_positiveORzero) {
    return ffun->get_ref()->cast2positive_or0(errSerious);
  } else {
    return ffun->get_ref()->calc();
  }
}

const int FlxObjReadBase::get_optPara_int_from_FlxFunction(const string& str)
{
  FlxOptionalParaBase* v1 = ParaBox.get(str);
  FlxOptionalParaFun* ffun = dynamic_cast<FlxOptionalParaFun*>( v1 );
  if (ffun==NULL) throw FlxException_Crude("FlxObjReadBase::get_optPara_tuint_from_FlxFunction");
  return ffun->get_ref()->cast2int();
}

const tuint FlxObjReadBase::get_optPara_tuint_from_FlxFunction(const string& str, const bool zero_is_allowed, const bool errSerious)
{
  FlxOptionalParaBase* v1 = ParaBox.get(str);
  FlxOptionalParaFun* ffun = dynamic_cast<FlxOptionalParaFun*>( v1 );
  if (ffun==NULL) throw FlxException_Crude("FlxObjReadBase::get_optPara_tuint_from_FlxFunction");
  if (zero_is_allowed) {
    return ffun->get_ref()->cast2tuintW0(errSerious);
  } else {
    return ffun->get_ref()->cast2tuint(errSerious);
  }
}

FlxFunDeg* FlxObjReadBase::get_optPara_FlxFunDeg(const string& str)
{
  void *v1 = ParaBox.get(str)->get();
  FlxFunDeg* fun = static_cast<FlxFunDeg*>( v1 ); 
  if (fun == NULL) {
    throw FlxException_Crude("FlxObjReadBase::get_optPara_FlxFunDeg");
  }
  return fun;
}

FlxMtxConstFun* FlxObjReadBase::get_optPara_FlxMtxFun(const string& str)
{
  void *v1 = ParaBox.get(str)->get();
  FlxMtxConstFun* fun = static_cast<FlxMtxConstFun*>( v1 ); 
  if (fun == NULL) {
    throw FlxException_Crude("FlxObjReadBase::get_optPara_FlxMtxFun");
  }
  return fun;
}

FlxString* FlxObjReadBase::get_optPara_FlxString(const string& str)
{
  void *v1 = ParaBox.get(str)->get();
  FlxString* fun = static_cast<FlxString*>( v1 ); 
  if (fun == NULL) {
    throw FlxException_Crude("FlxObjReadBase::get_optPara_FlxString");
  }
  return fun;
}

const string FlxObjReadBase::get_optPara_string_from_FlxString(const string& str, const bool lowercase)
{
  FlxOptionalParaBase *v1 = ParaBox.get(str);
  FlxOptionalParaFlxString* fstr = dynamic_cast<FlxOptionalParaFlxString*>( v1 );
  if (fstr==NULL) throw FlxException_Crude("FlxObjReadBase::get_optPara_string_from_FlxString");
  return fstr->get_ref()->eval(lowercase);
}

const std::string FlxObjReadBase::get_optPara_word_from_FlxString(const string& str, const bool lowercase, const bool emptyAllow, const bool numbegallow)
{
  FlxOptionalParaBase *v1 = ParaBox.get(str);
  FlxOptionalParaFlxString* fstr = dynamic_cast<FlxOptionalParaFlxString*>( v1 );
  if (fstr==NULL) throw FlxException_Crude("FlxObjReadBase::get_optPara_word_from_FlxString");
  return fstr->get_ref()->eval_word(lowercase,emptyAllow,numbegallow);
}

const string& FlxObjReadBase::get_optPara_string(const string& str, const bool lowercase)
{
  void *v1 = ParaBox.get(str)->get();
  std::string* fun = static_cast<std::string*>( v1 ); 
  if (fun == NULL) {
    throw FlxException_Crude("FlxObjReadBase::get_optPara_string");
  }
  if (lowercase) {
    std::transform(fun->begin(), fun->end(), fun->begin(), (int(*)(int)) std::tolower);
  }
  return *fun;
}

const bool FlxObjReadBase::get_optPara_bool(const string& str)
{
  bool* b1 = static_cast<bool*>( ParaBox.get(str)->get() ); 
  if (b1 == NULL) {
    throw FlxException_Crude("FlxObjReadBase::get_optPara_bool");
  }
  return *b1;
}

FlxObjReadLogBase::FlxObjReadLogBase(const bool from_library)
: FlxObjReadBase(from_library)
{
  // Default verbose-setting
    AllDefParaBox->insert(new FlxOptionalParaBool(true,"flxlog::verbose"));
    ParaBox.insert("vlog", "flxlog::verbose" );
}

const bool FlxObjReadLogBase::get_verboseLog()
{
  return *(static_cast<bool*>(ParaBox.get("vlog")->get() ));
}

FlxObjReadOutputBase::FlxObjReadOutputBase(const bool from_library)
: FlxObjReadBase(from_library) {
  // Default stream
    AllDefParaBox->insert(new FlxOptionalParaStream("cout","flxoutput::stream"));
    ParaBox.insert("stream", "flxoutput::stream" );
  // Default verbose-setting
    AllDefParaBox->insert(new FlxOptionalParaBool(true,"flxoutput::verbose"));
    ParaBox.insert("verbose", "flxoutput::verbose" );
  // Default check TOL for numerical values
    AllDefParaBox->insert(new FlxOptionalParaBool(false,"flxoutput::checktol"));
    ParaBox.insert("checktol", "flxoutput::checktol" );
  // prec
    AllDefParaBox->insert(new FlxOptionalParaFun(-ONE,"flxoutput::prec"));
    ParaBox.insert("prec", "flxoutput::prec" );  
  // fixW
    AllDefParaBox->insert(new FlxOptionalParaFun(ZERO,"flxoutput::fixw"));
    ParaBox.insert("fixw", "flxoutput::fixw" );
  // boost_str
    AllDefParaBox->insert(new FlxOptionalParaText("","flxoutput::boost_str"));
    ParaBox.insert("boost_str", "flxoutput::boost_str" );
}

const bool FlxObjReadOutputBase::get_verbose()
{
  return get_optPara_bool("verbose");
}

const bool FlxObjReadOutputBase::get_checkTOL()
{
  return get_optPara_bool("checktol");
}

const int FlxObjReadOutputBase::get_prec()
{
  FlxFunction* precF = get_optPara_FlxFunction("prec");
  try {
    const int prec = precF->cast2int();
    delete precF;
    return prec;
  } catch (...) {
    FLXMSG("FlxObjReadOutputBase::get_prec",1);
    delete precF;
    throw;
    return -1;        // dummy return
  }
}

const int FlxObjReadOutputBase::get_fixW()
{
  FlxFunction* fixWF = get_optPara_FlxFunction("fixw");
  try {
    const int fixW = fixWF->cast2int();
    delete fixWF;
    return fixW;
  } catch (...) {
    FLXMSG("FlxObjReadOutputBase::get_fixW",1);
    delete fixWF;
    throw;
    return -1;        // dummy return
  }
}

const string FlxObjReadOutputBase::get_boost_str()
{
  return get_optPara_string("boost_str");
}

const std::string FlxObjReadOutputBase::get_stream() {
  return *(static_cast<string*>( ParaBox.get("stream")->get() ));
}

void FlxObjConst::task()
{
  const tdouble d = fun->calc();
  if (opc=='=') {
    *cptr = d;
  } else if (opc=='+') {
    *cptr += d;
  } else if (opc=='-') {
    *cptr -= d;
  } else if (opc=='*') {
    *cptr *= d;
  } else if (opc=='/') {
    *cptr /= d;
  } else {
    throw FlxException_Crude("FlxObjConst::task");
  }
}

void FlxObjReadFCVbase::isdefined(const std::string& thename, const char theFCVindex, const bool errSerious) const {
  if ( data->FunBox.get(thename) != NULL && theFCVindex != 'F' ) {
    ostringstream ssV;
    ssV << "A function with the name ('" << thename << "') is already defined.";
    FlxError(errSerious,"FlxObjReadFCVbase::isdefined_1", ssV.str(), reader->getCurrentPos());
  }
  if ( data->ConstantBox.get(thename) != NULL && theFCVindex != 'C' ) {
    ostringstream ssV;
    ssV << "A 'const' variable with the name ('" << thename << "') is already defined.";
    FlxError(errSerious,"FlxObjReadFCVbase::isdefined_2", ssV.str(), reader->getCurrentPos());
  }
  if ( data->VarBox.get(thename) != NULL && theFCVindex != 'V' ) {
    ostringstream ssV;
    ssV << "A 'var' variable with the name ('" << thename << "') is already defined.";
    FlxError(errSerious,"FlxObjReadFCVbase::isdefined_3", ssV.str(), reader->getCurrentPos());
  }
  if ( data->ConstMtxBox.get(thename) != NULL && theFCVindex != 'M' ) {
    ostringstream ssV;
    ssV << "A 'mtxconst' matrix-variable with the name ('" << thename << "') is already defined.";
    FlxError(errSerious,"FlxObjReadFCVbase::isdefined_5", ssV.str(), reader->getCurrentPos());
  }
}

FlxObjReadLoops::FlxObjReadLoops() {
  // optional parameters
  AllDefParaBox->insert(new FlxOptionalParaFun(ZERO,"loops::maxpasses"));
  ParaBox.insert("maxpasses", "loops::maxpasses" );
}

const tuint FlxObjReadLoops::get_maxpasses() {
  void *v1 = ParaBox.get("maxpasses")->get();
  FlxFunction* fun = static_cast<FlxFunction*>( v1 ); 
  const tuint i1 = fun->cast2tuintW0();
  delete fun;
  return i1;
}

FlxObjReadConst::FlxObjReadConst()
{
  AllDefParaBox->insert(new FlxOptionalParaBool(false,"const::only_init"));
  ParaBox.insert("only_init", "const::only_init" );
}

FlxObjBase* FlxObjReadConst::read() {   
  const string cname = reader->getWord(true,false);
  return read(cname,true);
}

FlxObjBase* FlxObjReadConst::read(const string& cname, const bool allow_operators) {
  isdefined(cname, 'C', false);
  char opc = '=';
  if (allow_operators) {
    opc = reader->getChar(false);
    if (opc!='=') {
      if (opc!='+'&&opc!='-'&&opc!='*'&&opc!='/') {
        ostringstream ssV;
        ssV << "Character '" << opc << "' not allowed at this point.";
        throw FlxException("FlxObjReadConst::read", ssV.str(), reader->getCurrentPos() );
      }
      reader->getChar('=',false);
    } else {
      reader->setNext();
    }
  } else {
    reader->getChar('=',false);
  }
  FlxFunction *fun = new FlxFunction(funReader,false);
  try {
    read_optionalPara(false);
    if (get_optPara_bool("only_init")) {
      if (data->ConstantBox.get(cname,false)) {
        delete fun;
        return new FlxObjDummy();
      }
    }
    return new FlxObjConst(get_doLog(), cname, fun,opc );
  } catch (FlxException &e) {
    delete fun;
    throw;
  }
}

void FlxObjSub::task()
{
  try {
    data->SubBox.insert(sub_name, InternListSub);
    InternListSub = NULL;
  } catch (FlxException &e) {
    InternListSub = NULL;
    throw;
  }
  GlobalVar.slog(4) << "sub: defined procedure '" << sub_name << "()'." << std::endl;
}

void FlxObjSubDo::task() {
  if (InternListSub==NULL) {
    InternListSub = data->SubBox.get(sub_name);
    if (InternListSub==NULL) {
      ostringstream ssV;
      ssV << "A procedure with the name ('" << sub_name << "') does not exist.";
      throw FlxException("FlxObjReadSubDo::read_01", ssV.str());
    }
  }
  InternListSub->exec();
}

FlxCodeBlock* FlxObjReadCodeBlock::read_block(const bool UseIgnoreList, const bool errSerious) {
  if (UseIgnoreList) {
    data->IgnoreBox.activate_iL_recur();        // activate IgnoreList
  }
  reader->getChar('{',errSerious);
  FlxCodeBlock* obR = new FlxCodeBlock(true);
  try {
    while ( reader->whatIsNextChar() != '}' ) {
      FlxObjBase* ob = EvaluateCmd->evaluateCmd( );
      obR->attach_obj(ob);
    }
    reader->getChar('}',errSerious);
    // check for list of 'internal' variables
      if (reader->whatIsNextChar()==':') {
        reader->getChar(':',errSerious);
        reader->getChar(':',errSerious);
        reader->getChar('{',errSerious);
        while (reader->whatIsNextChar()!='}') {
          tdouble *icp = data->ConstantBox.get(reader->getWord(true,errSerious),true);
          obR->add_internal_const(icp);
          if (reader->whatIsNextChar()!='}') reader->getChar(',',errSerious);
        }
        reader->getChar('}',errSerious);
      }
  } catch ( FlxException &e ) {
    FLXMSG("FlxObjReadCodeBlock::read_block",1);
    if (obR) delete obR;
    throw;
  }
  if (UseIgnoreList) {
    data->IgnoreBox.deactivate_iL_recur();        // deactivate IgnoreList
  }
  return obR;
}

FlxObjBase* FlxObjReadSub::read() {
  const string sub_name = reader->getWord(true,false);
  reader->getChar('(',false);
  reader->getChar(')',false);
  FlxObjBase* loopblock = FlxObjReadCodeBlock::read_block(true,false);
  read_optionalPara(false);
  return new FlxObjSub( get_doLog(), sub_name, loopblock);
}

FlxObjBase* FlxObjReadSubDo::read()
{
  const string str1 = reader->getWord(true,false);
  FlxObjBase* obj1 = data->SubBox.get(str1);
  return read(obj1,str1);
}

FlxObjBase* FlxObjReadSubDo::read(FlxObjBase* objProcedure, const std::string& sub_name)
{
  reader->getChar('(',false);
  reader->getChar(')',false);
  read_optionalPara(false);
  return new FlxObjSubDo(get_doLog(),objProcedure,sub_name);
}

FunLSF_callback::FunLSF_callback(flx_lsf_callback lsfp, const std::string& lsf_name, RBRV_constructor* RndBox, std::string rbrvsetn)
: lsfp(lsfp), lsf_name(lsf_name), RndBox(RndBox), NOX(RndBox->get_NOX()), rvv(NOX), rbrvsetn(rbrvsetn)
{
  
}

FunLSF_callback::~FunLSF_callback()
{
  delete RndBox;
}

const tdouble FunLSF_callback::calc()
{
  RndBox->get_x_Vec(rvv.get_tmp_vptr());
  return lsfp(rvv.get_tmp_vptr_const(),NOX);
}

const std::string FunLSF_callback::write()
{
  std::string str1 = write_v() + '(' + rbrvsetn + ')';
  return str1;
}

FunReadFunLSF_callback::FunReadFunLSF_callback(flx_lsf_callback lsfp, string lsf_nameV, const bool from_library)
: FunReadFunBase(from_library), lsfp(lsfp)
{
  std::transform(lsf_nameV.begin(), lsf_nameV.end(), lsf_nameV.begin(), (int(*)(int)) std::tolower);
  lsf_name = lsf_nameV;
}

FunBase* FunReadFunLSF_callback::read( bool errSerious )
{
  // check if we need a user-defined rbrv-set
    std::string rbrvsetn = "nataf";
    if (reader->whatIsNextChar()!=')') {
      FlxString st(false,errSerious);
      rbrvsetn = st.eval(true);
    }
  // get relevant variables
    const std::vector<std::string> set_str_vec = parse_strseq_as_vec(rbrvsetn);
    RBRV_constructor* RndBox = new RBRV_constructor(set_str_vec,data->rbrv_box);
  try {
    return new FunLSF_callback(lsfp, lsf_name,RndBox, rbrvsetn ); 
  } catch (FlxException &e) {
    FLXMSG("FunReadFunLSF_callback::read_2",1);
    delete RndBox;
    throw;
  }
}
