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

#include "flxstringfun.h"
#include "flxdata.h"
#include "flxobjects.h"

void flxString_fun_insert(FlxStringFunBox& StrFunBox);

class FLXLIB_EXPORT FlxCreateObjReaders_FlxString : public FlxCreateObjReaders {
  public:
    void createObjReaders (FlxObjectReadBox* objReadBox );
    void createFunReaders (FlxData* dataBox );
};


// ------------------------------------------------------------------------------------------------


/**
* @brief object class: defines String-Constants
*
* strconst varName = FlxString;
*/
class FlxObjStrConst : public FlxObjBase {
  private:
    FlxString* cname;
    FlxString *fun;
    const bool append;
    void task();
  public:
    FlxObjStrConst ( const bool dolog, FlxString* cname, FlxString *fun, const bool append ) : FlxObjBase(dolog), cname(cname), fun(fun), append(append) {}
    ~FlxObjStrConst() { delete cname; delete fun; }
};

/**
* @brief object read class: FlxObjStrConst
*/
class FlxObjReadStrConst: public FlxObjReadBase {
  public:
    FlxObjBase* read ();
};

class StringFunStrConst : public StringFunBase, public FlxDataBase {
  private:
    FlxString* cname;
  public:
    StringFunStrConst(FlxString* cname) : cname(cname) {}
    ~StringFunStrConst() { delete cname; }
    virtual void eval(std::ostream& os);
    virtual const std::string write();
};

class FunReadFlxStringFunStrConst : public FunReadFlxStringFunBase {
  public:
    virtual StringFunBase* read ( bool errSerious );
};

// ------------------------------------------------------------------------------------------------

class FLXLIB_EXPORT StringFunStrStringStream: public StringFunBase, public FlxDataBase {
  private:
    FlxString* sname;
  public:
    StringFunStrStringStream(FlxString* sname) : sname(sname) {}
    ~StringFunStrStringStream() { delete sname; }
    virtual void eval(std::ostream& os);
    virtual const std::string write();
    
    static void getContent(FlxString* strmID, std::string& res);
};

class FunReadFlxStringFunStringStream : public FunReadFlxStringFunBase {
  public:
    virtual StringFunBase* read ( bool errSerious );
};

// ------------------------------------------------------------------------------------------------

class FLXLIB_EXPORT StringFunStrFileList: public StringFunBase, public FlxDataBase {
  private:
    FlxString* dir;
    FlxString* pattern;
    const std::string sep;
  public:
    StringFunStrFileList(FlxString* dir, FlxString* pattern, const std::string& sep) : dir(dir), pattern(pattern), sep(sep) {}
    ~StringFunStrFileList() { delete dir; delete pattern; }
    virtual void eval(std::ostream& os);
    virtual const std::string write();
};

class FunReadFlxStringFunFileList : public FunReadFlxStringFunBase {
  public:
    virtual StringFunBase* read ( bool errSerious );
};
// ------------------------------------------------------------------------------------------------

class StringFunTrim : public StringFunBase {
  private:
    FlxString* expr;
  public:
    StringFunTrim(FlxString* expr) : expr(expr) {}
    ~StringFunTrim() { delete expr; }
    virtual void eval(std::ostream& os);
    virtual const std::string write();
};

class FunReadFlxStringFunTrim : public FunReadFlxStringFunBase {
  public:
    virtual StringFunBase* read ( bool errSerious );
};

// ------------------------------------------------------------------------------------------------

class StringFunReplaceAll : public StringFunBase {
  private:
    FlxString* expr;
    FlxString* find_str;
    FlxString* replace_str;
  public:
    StringFunReplaceAll(FlxString* expr, FlxString* find_str, FlxString* replace_str) : expr(expr), find_str(find_str), replace_str(replace_str) {}
    ~StringFunReplaceAll() { delete expr; delete find_str; delete replace_str; }
    virtual void eval(std::ostream& os);
    virtual const std::string write();
};

class FunReadFlxStringFunReplaceAll : public FunReadFlxStringFunBase {
  public:
    virtual StringFunBase* read ( bool errSerious );
};
// ------------------------------------------------------------------------------------------------

class StringFunPWD : public StringFunBase {
  public:
    virtual void eval(std::ostream& os);
    virtual const std::string write();
};

class FunReadFlxStringFunPWD : public FunReadFlxStringFunBase {
  public:
    virtual StringFunBase* read ( bool errSerious );
};

// ------------------------------------------------------------------------------------------------

class StringFunFileName : public StringFunBase {
  private:
    const std::string filename;
  public:
    StringFunFileName(const std::string filename) : filename(filename) {}
    virtual void eval(std::ostream& os);
    virtual const std::string write();
};

class FunReadFlxStringFunFileName : public FunReadFlxStringFunBase {
  public:
    virtual StringFunBase* read ( bool errSerious );
};

// ------------------------------------------------------------------------------------------------

class StringFunDateFromToday : public StringFunBase {
  protected:
    FlxFunction* daydiff_fun;
    FlxString* format;
    const bool date_specified;
    const time_t bdate;
  public:
    StringFunDateFromToday(FlxFunction* daydiff_fun, FlxString* format, const bool date_specified, const time_t bdate) : daydiff_fun(daydiff_fun), format(format), date_specified(date_specified), bdate(bdate) {}
    ~StringFunDateFromToday();
    virtual void eval(std::ostream& os);
    virtual const std::string write();
};

class FunReadFlxStringFunDateFromToday : public FunReadFlxStringFunBase, public FlxReaderBase2  {
  public:
    virtual StringFunBase* read ( bool errSerious );
};

// ------------------------------------------------------------------------------------------------

class StringFunSubStr_search  : public FlxReaderBase2, public FlxBoxBase {
  public:
    enum flxSearchAction { pos, length, fchar, fstring };
  private:
    flxSearchAction actio;
    FlxFunction* fV;
    char cV;
    std::string sV;
  public:
    StringFunSubStr_search();
    ~StringFunSubStr_search() { if (fV) delete fV; }
    const size_t get_pos(const std::string& expr,const size_t start_pos);
    const std::string write();
};

class StringFunSubStr : public StringFunBase {
  private:
    FlxString* expr;
    StringFunSubStr_search* search_start;
    StringFunSubStr_search* search_end;
  public:
    StringFunSubStr(FlxString* expr, StringFunSubStr_search* search_start, StringFunSubStr_search* search_end) : expr(expr), search_start(search_start), search_end(search_end) {}
    ~StringFunSubStr();
    virtual void eval(std::ostream& os);
    virtual const std::string write();
};

class FunReadFlxStringFunSubStr : public FunReadFlxStringFunBase {
  public:
    virtual StringFunBase* read ( bool errSerious );
};

// ------------------------------------------------------------------------------------------------

class StringFunStrFromFile : public StringFunBase {
  private:
    FlxString* fname;
  public:
    StringFunStrFromFile(FlxString* cname) : fname(cname) {}
    ~StringFunStrFromFile() { delete fname; }
    virtual void eval(std::ostream& os);
    virtual const std::string write();
};

class FunReadFlxStringFunStrFromFile : public FunReadFlxStringFunBase {
  public:
    virtual StringFunBase* read ( bool errSerious );
};

// ------------------------------------------------------------------------------------------------

/**
* @brief arithmetic class: numerical integration of a function
*/
class FunStringFun_NumberFromString : public FunBase, public FlxDataBase {
  private:
    FlxString* strConstName;
  public:
    FunStringFun_NumberFromString(FlxString* strConstName) : strConstName(strConstName)  {}
    ~FunStringFun_NumberFromString() { delete strConstName; }
    const tdouble calc();
    const bool search_circref(FlxFunction* fcr) { return strConstName->search_circref(fcr); }
    const std::string write();
    const tuint precedence() { return 0; }
    const bool dependOn_Const(const tdouble* const thenumber);
    const bool optimize(FunBasePtr& optf, const Fun_OptimizeInfo &foi) { return false; }
};

class FunReadStringFun_NumberFromString : public FunReadFunBase {
  public:
    FunBase* read ( bool errSerious );
};

// ------------------------------------------------------------------------------------------------

/**
* @brief arithmetic class: numerical integration of a function
*/
class FunStringFun_StrLen : public FunBase, public FlxDataBase {
  private:
    FlxString* strExpr;
  public:
    FunStringFun_StrLen(FlxString* strExpr) : strExpr(strExpr)  {}
    ~FunStringFun_StrLen() { delete strExpr; }
    const tdouble calc();
    const bool search_circref(FlxFunction* fcr) { return strExpr->search_circref(fcr); }
    const std::string write();
    const tuint precedence() { return 0; }
    const bool dependOn_Const(const tdouble* const thenumber);
    const bool optimize(FunBasePtr& optf, const Fun_OptimizeInfo &foi) { return false; }
};

class FunReadStringFun_StrLen : public FunReadFunBase {
  public:
    FunBase* read ( bool errSerious );
};

// ------------------------------------------------------------------------------------------------

/**
* @brief arithmetic class: numerical integration of a function
*/
class FunStringFun_StrEqual : public FunBase, public FlxDataBase {
  private:
    FlxString* strExpr1;
    FlxString* strExpr2;
  public:
    FunStringFun_StrEqual(FlxString* strExpr1, FlxString* strExpr2) : strExpr1(strExpr1), strExpr2(strExpr2)  {}
    ~FunStringFun_StrEqual() { delete strExpr1; delete strExpr2; }
    const tdouble calc();
    const bool search_circref(FlxFunction* fcr) { return strExpr1->search_circref(fcr) || strExpr2->search_circref(fcr); }
    const std::string write();
    const tuint precedence() { return 0; }
    const bool dependOn_Const(const tdouble* const thenumber);
    const bool optimize(FunBasePtr& optf, const Fun_OptimizeInfo &foi) { return false; }
};

class FunReadStringFun_StrEqual : public FunReadFunBase {
  public:
    FunBase* read ( bool errSerious );
};

// ------------------------------------------------------------------------------------------------

/**
* @brief arithmetic class: numerical integration of a function
*/
class FunStringFun_StrContains : public FunBase, public FlxDataBase {
  private:
    FlxString* strExpr;    // the string in which to search
    FlxString* strSearch;  // the string for which we search
    FlxFunction* position; // the position at which to start (0 by default)
  public:
    FunStringFun_StrContains(FlxString* strExpr, FlxString* strSearch, FlxFunction* position) : strExpr(strExpr), strSearch(strSearch), position(position)  {}
    ~FunStringFun_StrContains() { delete strExpr; delete strSearch; if (position) delete position; }
    const tdouble calc();
    const bool search_circref(FlxFunction* fcr) { return strExpr->search_circref(fcr) || strSearch->search_circref(fcr) || position->search_circref(fcr); }
    const std::string write();
    const tuint precedence() { return 0; }
    const bool dependOn_Const(const tdouble* const thenumber);
    const bool optimize(FunBasePtr& optf, const Fun_OptimizeInfo &foi) { return false; }
};

class FunReadStringFun_StrContains : public FunReadFunBase, public FlxReaderBase2 {
  public:
    FunBase* read ( bool errSerious );
};

// ------------------------------------------------------------------------------------------------

class StringFunEquWrite : public StringFunBase {
  private:
    FlxFunction* fun;
  public:
    StringFunEquWrite(FlxFunction* fun) : fun(fun) {}
    ~StringFunEquWrite() { delete fun; }
    virtual void eval(std::ostream& os);
    virtual const std::string write();
};

class FunReadFlxStringFunEquWrite : public FunReadFlxStringFunBase, public FlxReaderBase2 {
  public:
    virtual StringFunBase* read ( bool errSerious );
};

// ------------------------------------------------------------------------------------------------

class StringFunVarWrite : public StringFunBase, public FlxDataBase {
  private:
    const std::string var_name;
  public:
    StringFunVarWrite(const std::string& var_name) : var_name(var_name) {}
    virtual void eval(std::ostream& os);
    virtual const std::string write();
};

class FunReadFlxStringFunVarWrite : public FunReadFlxStringFunBase {
  public:
    virtual StringFunBase* read ( bool errSerious );
};

// ------------------------------------------------------------------------------------------------

class StringFunRandStr : public StringFunBase, public FlxDataBase {
  private:
    FlxFunction* fun;
  public:
    StringFunRandStr(FlxFunction* fun) : fun(fun) {}
    ~StringFunRandStr() { delete fun; }
    virtual void eval(std::ostream& os);
    virtual const std::string write();
};

class FunReadFlxStringFunRandStr : public FunReadFlxStringFunBase, public FlxReaderBase2 {
  public:
    virtual StringFunBase* read ( bool errSerious );
};

// ------------------------------------------------------------------------------------------------


