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

#include "flxstringfun_fun.h"
#include "flxdata.h"
#include "flxostools_files.h"

#include <ctime>
#include <fstream>
#include <streambuf>

// ------------------------------------------------------------------------------------------------


void flxString_fun_insert(FlxStringFunBox& StrFunBox)
{
  StrFunBox.insert("trim",new FunReadFlxStringFunTrim());
  StrFunBox.insert("replace_all",new FunReadFlxStringFunReplaceAll());
  StrFunBox.insert("pwd",new FunReadFlxStringFunPWD());
  StrFunBox.insert("filename",new FunReadFlxStringFunFileName());
  StrFunBox.insert("datefromtoday",new FunReadFlxStringFunDateFromToday());
  StrFunBox.insert("strconst",new FunReadFlxStringFunStrConst());
  StrFunBox.insert("substr",new FunReadFlxStringFunSubStr());
  StrFunBox.insert("strfromfile",new FunReadFlxStringFunStrFromFile());
  StrFunBox.insert("stringstream",new FunReadFlxStringFunStringStream());
  StrFunBox.insert("file_list",new FunReadFlxStringFunFileList());
  StrFunBox.insert("equwrite",new FunReadFlxStringFunEquWrite());
  StrFunBox.insert("varwrite",new FunReadFlxStringFunVarWrite());
  StrFunBox.insert("randstr",new FunReadFlxStringFunRandStr());
}


void FlxCreateObjReaders_FlxString::createObjReaders(FlxObjectReadBox* objReadBox) {
  objReadBox->insert("strconst", new FlxObjReadStrConst());
}

void FlxCreateObjReaders_FlxString::createFunReaders(FlxData* dataBox)
{
  dataBox->FunBox.insert("numberfromstring", new FunReadStringFun_NumberFromString() );
  dataBox->FunBox.insert("strlen", new FunReadStringFun_StrLen() );
  dataBox->FunBox.insert("strequal", new FunReadStringFun_StrEqual() );
  dataBox->FunBox.insert("strcontains", new FunReadStringFun_StrContains() );
}

// ------------------------------------------------------------------------------------------------



void FlxObjStrConst::task()
{
  const std::string cname_ = cname->eval_word(true);
  const std::string fun_ = fun->eval(false);
  if (append) {
    std::string& astr = data->strConstBox.get_ref(cname_);
    astr += fun_;
  } else {
    data->strConstBox.insert(cname_,fun_);
  }
}

FlxObjBase* FlxObjReadStrConst::read() {   
  FlxString* cname = new FlxString(false,false);
  FlxString* fun = NULL;
  bool append = false;
  try {
    if (reader->whatIsNextChar()=='+') {
      reader->getChar('+') ;
      append = true;
    }
    reader->getChar('=');
    fun = new FlxString(true,false);
    read_optionalPara(false);
    return new FlxObjStrConst(get_doLog(),cname,fun,append);
  } catch (FlxException& e) {
    delete cname;
    if (fun) delete fun;
    throw;
  }
}

void StringFunStrConst::eval(std::ostream& os)
{
  const std::string cname_ = cname->eval_word(true);
  os << data->strConstBox.get(cname_);
}

const std::string StringFunStrConst::write()
{
  return "strconst(" + cname->write() + ")";
}

StringFunBase* FunReadFlxStringFunStrConst::read(bool errSerious)
{
  return new StringFunStrConst(new FlxString(false,errSerious));
}


// ------------------------------------------------------------------------------------------------

void StringFunStrStringStream::getContent(FlxString* strmID, std::string& res)
{
  const std::string sname_ = strmID->eval_word(true);
  std::ostringstream* thestreamSS = dynamic_cast<std::ostringstream*> (data->OstreamBox.get(sname_));
  if ( thestreamSS ) {
    res = thestreamSS->str();
    thestreamSS->str("");
    thestreamSS->clear();
  } else {
    throw FlxException("StringFunStrStringStream::getContent","The stream '"+sname_+"' is not a string-stream.");
  }
}

void StringFunStrStringStream::eval(std::ostream& os)
{
  std::string res;
  getContent(sname,res);
  os << res;
}

const std::string StringFunStrStringStream::write()
{
  return "stringstream(" + sname->write() + ")";
}

StringFunBase* FunReadFlxStringFunStringStream::read(bool errSerious)
{
  return new StringFunStrStringStream(new FlxString(false,errSerious));
}

// ------------------------------------------------------------------------------------------------

void StringFunStrFileList::eval(std::ostream& os)
{
  const std::string dirStr = dir->eval(false);
  std::vector< std::string > match_files;
  getFiles(dirStr,pattern->eval(false),match_files);
  std::sort(match_files.begin(),match_files.end());
  const size_t N = match_files.size();
  for (size_t i=0;i<N;++i) {
    if (i>0) os << sep;
    os << match_files[i];
  }
}

const std::string StringFunStrFileList::write()
{
  return "file_list(" + dir->write() + ",pattern=\"" + pattern->write() + "\",sep=\"" + sep + "\")";
}

StringFunBase* FunReadFlxStringFunFileList::read(bool errSerious)
{
  FlxString* dir = new FlxString(false,false);
  FlxString* pattern = NULL;
  try {
    std::string sep = ";";
    while (reader->whatIsNextChar()==',') {
      reader->getChar(',');
      const std::string pid = reader->getWord(true,false);
      reader->getChar('=');
      if (pid=="pattern") {
        pattern = new FlxString(false,false); 
      } else if (pid=="sep") {
        sep = reader->getText(false,false);
      } else {
        throw FlxException("FunReadFlxStringFunFileList::read","Unknown parameter '"+pid+"'.");
      }
    }
    if (pattern==NULL) pattern = new FlxString(new FlxString_String(".+",false),false);
    return new StringFunStrFileList(dir,pattern,sep);
  } catch (FlxException& e) {
    delete dir;
    if (pattern) delete pattern;
    throw;
  }
}

// ------------------------------------------------------------------------------------------------

void StringFunTrim::eval(std::ostream& os)
{
  std::string str = expr->eval(false);
  os << trim(str);
}

const std::string StringFunTrim::write()
{
  std::string ostr = "trim(" + expr->write() + ')';
  return ostr;
}

StringFunBase* FunReadFlxStringFunTrim::read(bool errSerious)
{
  FlxString* expr = new FlxString(true,false);
  return new StringFunTrim(expr);
}

// ------------------------------------------------------------------------------------------------

void StringFunReplaceAll::eval(std::ostream& os)
{
  std::string str = expr->eval(false);
  os << find_and_replace_all(str,find_str->eval(false),replace_str->eval(false));
}

const std::string StringFunReplaceAll::write()
{
  std::string ostr = "replace_all(" + expr->write() + "," + find_str->write() + "," + replace_str->write() + ')';
  return ostr;
}

StringFunBase* FunReadFlxStringFunReplaceAll::read(bool errSerious)
{
  FlxString* expr = new FlxString(true,false);
  reader->getChar(',',errSerious);
  FlxString* find_str = new FlxString(true,false);
  reader->getChar(',',errSerious);
  FlxString* replace_str = new FlxString(true,false);
  return new StringFunReplaceAll(expr,find_str,replace_str);
}

// ------------------------------------------------------------------------------------------------

void StringFunPWD::eval(std::ostream& os)
{
  os << GlobalVar.pwd;
}

const std::string StringFunPWD::write()
{
  return "pwd()";
}

StringFunBase* FunReadFlxStringFunPWD::read(bool errSerious)
{
  return new StringFunPWD();
}

// ------------------------------------------------------------------------------------------------


void StringFunFileName::eval(std::ostream& os)
{
  os << filename;
}

const std::string StringFunFileName::write()
{
  return "filename()";
}

StringFunBase* FunReadFlxStringFunFileName::read(bool errSerious)
{
  return new StringFunFileName(reader->get_FileName());
}

// ------------------------------------------------------------------------------------------------

StringFunDateFromToday::~StringFunDateFromToday()
{
  delete daydiff_fun; 
  delete format;
}

void StringFunDateFromToday::eval(std::ostream& os)
{
  const int ddiff = daydiff_fun->cast2int();
  // current date/time based on current system
    time_t now = (date_specified?bdate:time(0));
    now +=  ddiff*(60*60*24);
    tm *ltm = localtime(&now);
  // output
    char buffer [80];
    const std::string fstr = format->eval(false);
    std::strftime (buffer,80,fstr.c_str(),ltm);
  os << buffer;
}

const std::string StringFunDateFromToday::write()
{
  std::ostringstream ssV;
  ssV << "datefromtoday(" << daydiff_fun->write() << "," << format->write();
  if (date_specified) {
    char buffer [80];
      tm *ltm = localtime(&bdate);
      std::strftime (buffer,80,"%F",ltm);
      ssV << "," << buffer;
  }
  ssV << ")";
  return ssV.str();
}

StringFunBase* FunReadFlxStringFunDateFromToday::read(bool errSerious)
{
  FlxFunction* ddiff = new FlxFunction(funReader,false);
  FlxString* format = NULL;
  bool date_specified = false;
  time_t bdate = time(0);
  try {
    if (reader->whatIsNextChar()==',') {
      reader->getChar(',');
      format = new FlxString(false,false);
    } else {
      format = new FlxString(new FlxString_String("%d.%m.%Y",false),false);
    }
    if (reader->whatIsNextChar()==',') {
      reader->getChar(',');
      bdate = reader->getDate(false);
      date_specified = true;
    }
  } catch (FlxException& e) {
    delete ddiff;
    if (format) delete format;
    throw;
  }
  return new StringFunDateFromToday(ddiff,format,date_specified,bdate);
}

// ------------------------------------------------------------------------------------------------

StringFunSubStr_search::StringFunSubStr_search()
: fV(NULL), cV(' ')
{
  const char c = reader->getChar();
  switch (c) {
    case 'p':
      actio = pos;
      break;
    case 'l':
      actio = length;
      break;
    case 'c':
      actio = fchar;
      break;
    case 's':
      actio = fstring;
      break;
    default:
    {
      std::ostringstream ssV;
      ssV << "Character '" << c << "' not expected.";
      throw FlxException("StringFunSubStr_search::StringFunSubStr_search",ssV.str());
    }
  }
  reader->getChar(':');
  switch (actio) {
    case pos:
    case length:
      fV = new FlxFunction(funReader,false);
      break;
    case fchar:
      cV = reader->getChar();
      break;
    case fstring:
      sV = reader->getText(true);
      break;
    default:
      throw FlxException_Crude("StringFunSubStr_search::StringFunSubStr_search_1");
  }
}

const std::string StringFunSubStr_search::write()
{
  switch (actio) {
    case pos:
      return "p:" + fV->write();
    case length:
      return "l:" + fV->write();
    case fchar:
      return "c:" + cV;
    case fstring:
      return "s:\"" + sV + '\"';
    default:
      throw FlxException_Crude("StringFunSubStr_search::StringFunSubStr_search_1");
  }
}

const size_t StringFunSubStr_search::get_pos(const std::string& expr, const size_t start_pos)
{
  switch (actio) {
    case pos:
    {
      const size_t prop_pos = fV->cast2tuintW0();
      if (prop_pos<start_pos) {
        std::ostringstream ssV;
        ssV << "The position (" << prop_pos << ") must not be smaller than the starting position (" << start_pos << ").";
        throw FlxException("StringFunSubStr_search::get_pos_1",ssV.str());
      }
      return prop_pos;
    }
    case length:
      return start_pos + fV->cast2tuint();
    case fchar:
    {
      const size_t prop_pos = expr.find(cV,start_pos);
      if (prop_pos == std::string::npos) {
        std::ostringstream ssV;
        ssV << "The character '" << cV << "' was not found in the string (" << expr.substr(start_pos) << ").";
        throw FlxException("StringFunSubStr_search::get_pos_2",ssV.str());
      }
      return prop_pos;
    }
    case fstring:
    {
      const size_t prop_pos = expr.find(sV,start_pos);
      if (prop_pos == std::string::npos) {
        throw FlxException("StringFunSubStr_search::get_pos_3","The expression '" + sV + "' was not found in the string (" + expr.substr(start_pos) + ").");
      }
      return prop_pos;
    }
    default:
      throw FlxException_Crude("StringFunSubStr_search::get_pos_4");
  }
}

void StringFunSubStr::eval(std::ostream& os)
{
  const std::string str_expr = expr->eval(false);
  const size_t pos_start = search_start->get_pos(str_expr,0);
  if (pos_start>=str_expr.length()) {
    std::ostringstream ssV;
    ssV << "The starting position (" << pos_start << ") must be smaller than the length of the expression (" << str_expr.length() << ").";
    throw FlxException("StringFunSubStr::eval_1",ssV.str());
  }
  if (search_end) {
    const size_t pos_end = search_end->get_pos(str_expr,pos_start);
    if (pos_end>str_expr.length()) {
      std::ostringstream ssV;
      ssV << "The ending position (" << pos_start << ") must be smaller or equal than the length of the expression (" << str_expr.length() << ").";
      throw FlxException("StringFunSubStr::eval_2",ssV.str());
    }
    os << str_expr.substr(pos_start,pos_end-pos_start);
  } else {
    os << str_expr.substr(pos_start);
  }
}

const std::string StringFunSubStr::write()
{
  std::string ostr = "substr(" + expr->write() + ',' + search_start->write();
  if (search_end) {
    ostr += ',' + search_end->write();
  }
  ostr += ')';
  return ostr;
}

StringFunSubStr::~StringFunSubStr()
{
  delete expr;
  delete search_start;
  if (search_end) delete search_end;
}

StringFunBase* FunReadFlxStringFunSubStr::read(bool errSerious)
{
  FlxString* expr = new FlxString(true,false);
  StringFunSubStr_search* search_start = NULL;
  StringFunSubStr_search* search_end = NULL;
  try {
    reader->getChar(',');
    search_start = new StringFunSubStr_search();
    if (reader->whatIsNextChar()==',') {
      reader->getChar(',');
      search_end = new StringFunSubStr_search();
    }
    return new StringFunSubStr(expr,search_start,search_end);
  } catch (FlxException& e) {
    FLXMSG("FunReadFlxStringFunSubStr_1",1);
    delete expr;
    if (search_start) delete search_start;
    if (search_end) delete search_end;
    throw;
  }
}

// ------------------------------------------------------------------------------------------------

void StringFunStrFromFile::eval(std::ostream& os)
{
  const std::string fname_ = fname->eval(false);
  // open the file
    std::ifstream t(fname_.c_str());
    if (t.is_open()==false) {
      throw FlxException_NeglectInInteractive("StringFunStrFromFile::eval", "The file '" + fname_ + "' could not be opened.");
    }
  // read the file
    std::string str;
    t.seekg(0, std::ios::end);   
    str.reserve(t.tellg());
    t.seekg(0, std::ios::beg);
    str.assign((std::istreambuf_iterator<char>(t)),
            std::istreambuf_iterator<char>());
  os << str;
}

const std::string StringFunStrFromFile::write()
{
  return "strfromfile(" + fname->write() + ")";
}

StringFunBase* FunReadFlxStringFunStrFromFile::read(bool errSerious)
{
  return new StringFunStrFromFile(new FlxString(false,errSerious));
}

// ------------------------------------------------------------------------------------------------

const tdouble FunStringFun_NumberFromString::calc()
{
  // determine string-expression to read from
    const std::string cname = strConstName->eval_word(true);
    std::string expr = data->strConstBox.get(cname);
    ReadStream myreader(expr);
  // find possible start for number in string
    while (true) {
      if (myreader.getNextType()==ReadStream::ENDOFFILE) break;
      if (myreader.nextCanBeNumber()) break;
      myreader.getChar();
    }
  // ensure that end has not yet been reached
    if (myreader.getNextType()==ReadStream::ENDOFFILE) {
      throw FlxException_NeglectInInteractive("FunStringFun_NumberFromString::calc_01","A number could not be retrieved");
    }
  // actually retrieve the number from the string
    const tdouble res = myreader.get_Double(false);
  // retrieve remaining string after number
    std::string expr2;
    while (myreader.getNextType()!=ReadStream::ENDOFFILE) {
      expr2 += myreader.getNextLine(false);
      expr2 += '\n';
    }
  #if FLX_DEBUG
    const std::string str3 = expr2.substr(expr2.size()-(ReadStream_String_End.size()+1),ReadStream_String_End.size());
    if (str3!=ReadStream_String_End) {
      throw FlxException_Crude("FunStringFun_NumberFromString::calc_02");
    }
  #endif
  expr = expr2.substr(0,expr2.size()-(ReadStream_String_End.size()+1));
  data->strConstBox.insert(cname,expr);
  return res;
}

const bool FunStringFun_NumberFromString::dependOn_Const(const tdouble*const thenumber)
{
  throw FlxException_NotImplemented("FunStringFun_NumberFromString::dependOn_Const");
}

const std::string FunStringFun_NumberFromString::write()
{
  std::ostringstream ssV;
  ssV << "numberfromstring(" << strConstName->write() << ")";
  return ssV.str();
}

FunBase* FunReadStringFun_NumberFromString::read(bool errSerious)
{
  FlxString* strConstName = new FlxString(false,errSerious);
  return new FunStringFun_NumberFromString(strConstName);
}

const tdouble FunStringFun_StrLen::calc()
{
  const std::string expr = strExpr->eval(false);
  return expr.length();
}

const bool FunStringFun_StrLen::dependOn_Const(const tdouble*const thenumber)
{
  throw FlxException_NotImplemented("FunStringFun_StrLen::dependOn_Const");
}

const std::string FunStringFun_StrLen::write()
{
  std::ostringstream ssV;
  ssV << "strlen(" << strExpr->write() << ")";
  return ssV.str();
}


FunBase* FunReadStringFun_StrLen::read(bool errSerious)
{
  FlxString* strExpr = new FlxString(true,errSerious);
  return new FunStringFun_StrLen(strExpr);
}

const tdouble FunStringFun_StrEqual::calc()
{
  const std::string expr1 = strExpr1->eval(false);
  const std::string expr2 = strExpr2->eval(false);
  return (expr1==expr2)?ONE:ZERO;
}

const bool FunStringFun_StrEqual::dependOn_Const(const tdouble*const thenumber)
{
  throw FlxException_NotImplemented("FunStringFun_StrEqual::dependOn_Const");
}

const std::string FunStringFun_StrEqual::write()
{
  std::ostringstream ssV;
  ssV << "strequal(" << strExpr1->write() << "," << strExpr2->write() << ")";
  return ssV.str();
}

FunBase* FunReadStringFun_StrEqual::read(bool errSerious)
{
  FlxString* strExpr1 = new FlxString(true,errSerious);
  FlxString* strExpr2 = NULL;
  try {
    reader->getChar(',',errSerious);
    strExpr2 = new FlxString(true,errSerious);
    return new FunStringFun_StrEqual(strExpr1, strExpr2);
  } catch (FlxException& e) {
    delete strExpr1;
    if (strExpr2) delete strExpr2;
    throw;
  }
}

void StringFunEquWrite::eval(std::ostream& os)
{
  os << fun->write();
}

const std::string StringFunEquWrite::write()
{
  return "equwrite(" + fun->write() + ")";
}

StringFunBase* FunReadFlxStringFunEquWrite::read(bool errSerious)
{
  FlxFunction* fun = new FlxFunction(funReader,false);
  return new StringFunEquWrite(fun);
}

void StringFunVarWrite::eval(std::ostream& os)
{
  FlxFunction* fun = data->VarBox.get(var_name);
  if (fun) {
    os << fun->write();
  } else {
    throw FlxException("StringFunVarWrite::eval", "var-variable '" + var_name + "' does not exist.");
  }
}

void StringFunRandStr::eval(std::ostream& os)
{ // code adapted from https://stackoverflow.com/questions/440133/how-do-i-create-a-random-alpha-numeric-string-in-c
  // length of string
    const tuint N = fun->cast2tuint();
  // allowed characters
    auto randchar = []() -> char
      {
          const char charset[] =
          "0123456789"
          "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
          "abcdefghijklmnopqrstuvwxyz";
          const size_t L = (sizeof(charset)-1);
          const tuint rid = data->RndCreator.gen_smp_index(L);
          return charset[ rid ];
      };
  // sample string
    std::string str(N,0);
    std::generate_n( str.begin(), N, randchar );
  os << str;
}

const std::string StringFunRandStr::write()
{
  return "randstr(" + fun->write() + ")";
}

StringFunBase* FunReadFlxStringFunRandStr::read(bool errSerious)
{
  FlxFunction* fun = new FlxFunction(funReader,false);
  return new StringFunRandStr(fun);
}

const tdouble FunStringFun_StrContains::calc()
{
  const std::string expr = strExpr->eval(false);
  const std::string search = strSearch->eval(false);
  size_t pos = 0;
  if (position) {
    pos = position->cast2tulongW0(false);
  }
  const size_t prop_pos = expr.find(search,pos);
  if (prop_pos == std::string::npos) {
    return -ONE;
  } else {
    return prop_pos;
  }
}

const std::string FunStringFun_StrContains::write()
{
  std::ostringstream ssV;
  ssV << "strcontains(" << strExpr->write() << "," << strSearch->write();
  if (position) {
    ssV << "," << position->write(); 
  }
  ssV << ")";
  return ssV.str();
}

const bool FunStringFun_StrContains::dependOn_Const(const tdouble*const thenumber)
{
  throw FlxException_NotImplemented("FunStringFun_StrContains::dependOn_Const");
}

FunBase* FunReadStringFun_StrContains::read(bool errSerious)
{
  FlxString* strExpr = new FlxString(true,errSerious);
  FlxString* strSearch = NULL;
  FlxFunction* position = NULL;
  try {
  reader->getChar(',',errSerious);
  FlxString* strSearch = new FlxString(true,errSerious);
  if (reader->whatIsNextChar()==',') {
    reader->getChar(',');
    position = new FlxFunction(funReader,false);
  }
  return new FunStringFun_StrContains(strExpr, strSearch,position);  
  } catch (FlxException &e) {
    delete strExpr;
    if (strSearch) delete strSearch;
    if (position) delete position;
    throw;
  }
}

const std::string StringFunVarWrite::write()
{
  return "varwrite(" + var_name + ")";
}

StringFunBase* FunReadFlxStringFunVarWrite::read(bool errSerious)
{
  return new StringFunVarWrite(reader->getWord(true,false));
}




