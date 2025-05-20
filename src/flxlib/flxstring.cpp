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

#include "flxstring.h"

// #include <boost/algorithm/string.hpp>

using namespace std;

FlxFunctionReader* FlxReaderBase2::funReader = NULL;
FlxStringFunBox* FlxString::StrFunBox = NULL;



const string& trim(string& str)
{
  boost::algorithm::trim(str);
  return str;
}

const std::string& find_and_replace_all(std::string& source_str, const std::string& find_str, const std::string& replace_str)
{
  for(string::size_type i = 0; (i = source_str.find(find_str, i)) != string::npos;)  {
    source_str.replace(i, find_str.length(), replace_str);
    i += replace_str.length();
  }
  return source_str;
}




FlxString::FlxString(FlxString_Base* strB, const bool multiline) :
strLst(new std::list<FlxString_Base*>),instances(new tuint(0)), multiline(multiline)
{
  strLst->push_back(strB);
}

FlxString::FlxString(const bool multiline, const bool errSerious) :
strLst(new std::list<FlxString_Base*>),instances(new tuint(0)), multiline(multiline)
{
  bool b1 = false;
  FlxString_Base* strB = NULL;
  try {
    do {
      if (b1) reader->getChar('&',errSerious);
      else b1 = true;
      strB = NULL;
      // --------------------------------------------------------------------------------
      // Text
      // --------------------------------------------------------------------------------
      if (reader->whatIsNextChar() == '"' ) {
        strB = new FlxString_String(reader->getText(multiline,errSerious),false);
      // --------------------------------------------------------------------------------
      // FlxFunction
      // --------------------------------------------------------------------------------
      } else if ( reader->whatIsNextChar() == '{' ) {
        reader->getChar('{',errSerious);
        FlxFunction* fs = new FlxFunction(funReader,errSerious);
        // type of output?
          FlxString_Fun::otype id = FlxString_Fun::ot_dbl;
          std::string boost_str;
          try {
            reader->getChar('}',errSerious);
            if (reader->whatIsNextChar()=='$') {
              reader->getChar('$',errSerious);
              const std::string ids = reader->getWord(true,errSerious);
              id = FlxString_Fun::parse_ot(ids);
              if (id == FlxString_Fun::ot_boost) {
                reader->getChar(':');
                boost_str = reader->getText();
              }
            }
          } catch (FlxException &e) {
            FLXMSG("FlxString::FlxString_2",1);
            delete fs;
            throw;
          }
        strB = new FlxString_Fun(fs,id,boost_str);
      // --------------------------------------------------------------------------------
      // String-Function
      // --------------------------------------------------------------------------------
      } else if ( reader->whatIsNextChar() == '$' ) {
        reader->getChar('$',errSerious);
        strB = new FlxString_StrFun(StrFunBox->read(reader,errSerious));
      // --------------------------------------------------------------------------------
      // Word
      // --------------------------------------------------------------------------------
      } else {
        strB = new FlxString_String(reader->getWord(false,errSerious,true,true),true);
      }
      strLst->push_back(strB);
    } while (reader->whatIsNextChar() == '&');
  } catch (FlxException &e) {
    FLXMSG("FlxString::FlxString_4",1);
    if (strB) delete strB;
    std::list<FlxString_Base*>::iterator iter;
    for (iter = strLst->begin(); iter != strLst->end(); ++iter) {
      delete (*iter);
    }
    delete strLst;
    delete instances;
    throw;
  }
}

FlxString::FlxString(const FlxString& FlxStrF)
: strLst(FlxStrF.strLst),instances(FlxStrF.instances), multiline(FlxStrF.multiline)
{
  (*instances)++;
}

FlxString& FlxString::operator=(const FlxString& FlxStrF)
{
  if ( this != &FlxStrF ) {
    free_mem();
    instances = FlxStrF.instances;
    strLst = FlxStrF.strLst;
    (*instances)++;
    multiline = FlxStrF.multiline;
  }
  return *this;
}

void FlxString::free_mem()
{
  if (instances == NULL) return;
  if ( (*instances) == 0)  {
    if (strLst!=NULL) {
      std::list<FlxString_Base*>::iterator iter;
      for (iter = strLst->begin(); iter != strLst->end(); ++iter) {
        delete (*iter);
      }
      delete strLst;
    }
    if ( instances != NULL ) delete instances;
  } else {
    (*instances)--;
  }
}

FlxString::~FlxString()
{
  free_mem();
}

void FlxString::eval(ostream& os)
{
  std::list<FlxString_Base*>::iterator iter;
  for (iter = strLst->begin(); iter != strLst->end(); ++iter) {
    (*iter)->eval(os);
  }
}

const std::string FlxString::eval(const bool lowercase)
{
  std::ostringstream ssV;
  eval(ssV);
  std::string str1 = ssV.str();
  if (!multiline) {
    if (str1.find('\n')!=std::string::npos) {
      throw FlxException("FlxString::eval","The evaluated string-expression must not contain a line break.","The FlxString is: " + write());
    }
  }
  if (lowercase) {
    std::transform(str1.begin(), str1.end(), str1.begin(), (int(*)(int)) std::tolower);
  }
  return str1;
}

const std::string FlxString::eval_word(const bool lowercase, const bool emptyAllow, const bool numbegallow)
{
  std::string strV = eval(lowercase);
  if (strV.length() == 0) {
    if (emptyAllow) {
      return "";
    } else {
      std::ostringstream ssV_2;
      ssV_2 << "Word must not be empty.";
      throw FlxException("FlxString::eval_word_1", ssV_2.str()); 
    }
  }
  tuint i_start = 0;
  if (!numbegallow) {
    if ( ReadStream::getType( strV[0] ) != ReadStream::STRING ) {
      std::ostringstream ssV_2;
      ssV_2 << "Evaluated string '" << strV << "' is not of type 'word'.";
      throw FlxException("FlxString::eval_word_2", ssV_2.str()); 
    }
    i_start = 1;
  }
  for ( std::string::size_type i = i_start; i < strV.length(); ++i) {
    ReadStream::InpType iT = ReadStream::getType( strV[i] );
    if ( iT != ReadStream::STRING && iT != ReadStream::NUMERAL ) {
      std::ostringstream ssV_2;
      ssV_2 << "Evaluated string '" << strV << "' is not of type 'word'.";
      throw FlxException("FlxString::eval_word_3", ssV_2.str()); 
    }
  }
  return strV;
}

const std::string FlxString::write()
{
  std::string strV = "";
  bool b1 = false;
  std::list<FlxString_Base*>::iterator iter;
  for (iter = strLst->begin(); iter != strLst->end(); ++iter) {
    if (b1) {
      strV+="&";
    } else {
      b1 = true;
    };
    strV+= (*iter)->write();
  }
  return strV;
}

void FlxString::assign(FlxString* strA)
{
  if (this == strA) return;
  if (this->strLst == strA->strLst) return;
  if ((*instances)==0) {
    std::list<FlxString_Base*>::iterator iter;
    for (iter = strLst->begin(); iter != strLst->end(); ++iter) {
      delete (*iter);
    }
    delete strLst; 
    delete instances;
  } else {
    (*instances)--;
  }
  strLst = strA->strLst; 
  strA->strLst = NULL;
  instances = strA->instances;
  strA->instances = NULL;
  multiline = strA->multiline;
  delete strA;
}

const bool FlxString::is_static() const
{
  if (strLst->size()==1) {
    FlxString_Base* vp = strLst->front();
    FlxString_String* sp = dynamic_cast<FlxString_String*>(vp);
    if (sp) return true;
  }
  return false;
}

const bool FlxString::search_circref(FlxFunction* fcr)
{
  std::list<FlxString_Base*>::iterator iter;
  for (iter = strLst->begin(); iter != strLst->end(); ++iter) {
    if ((*iter)->search_circref(fcr)) {
      return true;
    }
  }
  return false;
}

const std::string FlxString_String::write()
{
  if (isWord) {
    return strV;
  } else {
    return "\"" + strV + "\"";
  }
}

void FlxString_Fun::eval(ostream& os)
{
  #if FLX_DEBUG
    if (id==ot_boost && boost_str.empty()) throw FlxException_Crude("FlxString_Fun::eval_1");
    if (id!=ot_boost && boost_str.empty()==false) throw FlxException_Crude("FlxString_Fun::eval_2");
  #endif
  switch (id) {
    case ot_dbl:
      os << GlobalVar.Double2String(fun->calc());
      break;
    case ot_int:
      os << round_flx(fun->calc());
      break;
    case ot_boost:
    {
      const tdouble res = fun->calc();
      try {
      os << format(boost_str) % res;
      } catch (...) {
        std::ostringstream ssV;
        ssV << "The output-format string '" << boost_str << "' caused an error.";
        throw FlxException_NeglectInInteractive("FlxString_Fun::eval_4",ssV.str());
      }
      break;
    }
    #if FLX_DEBUG
    default:
      throw FlxException_Crude("FlxString_Fun::eval");
    #endif
  }
}

const string FlxString_Fun::write()
{
  switch (id) {
    case ot_dbl:
      return "{" + fun->write() + "}";
    case ot_int:
      return "{" + fun->write() + "}$int";
    case ot_boost:
      return "{" + fun->write() + "}$boost:\"" + boost_str + "\"";
  }
  throw FlxException_Crude("FlxString_Fun::write");
}

const FlxString_Fun::otype FlxString_Fun::parse_ot(const string& str)
{
  otype res;
  if (str=="dbl") res = ot_dbl;
  else if (str=="int") res = ot_int;
  else if (str=="udef") res = ot_boost;
  else {
    std::ostringstream ssV;
    ssV << "Unknown keyword '" << str << "'.";
    throw FlxException("FlxString_Fun::parse_ot", ssV.str()); 
  }
  return res;
}

const string FlxString_Fun::parse_ot(const FlxString_Fun::otype ot)
{
  switch (ot) {
    case ot_dbl:
      return "dbl";
    case ot_int:
      return "int";
    case ot_boost:
      return "udef";
    default:
      throw FlxException_Crude("FlxString_Fun::parse_ot");
  }
}


