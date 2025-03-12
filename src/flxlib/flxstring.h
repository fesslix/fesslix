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

#include <list>

FLXLIB_EXPORT const std::string& trim(std::string& str);
FLXLIB_EXPORT const std::string& find_and_replace_all(std::string& source_str, const std::string& find_str, const std::string& replace_str);


class FLXLIB_EXPORT FlxReaderBase2 : public FlxReaderBase {
  protected:
    static FlxFunctionReader *funReader;
  public:
    virtual ~FlxReaderBase2() {};
};


class FLXLIB_EXPORT FlxString_Base {
  public:
    virtual ~FlxString_Base() {};
    virtual void eval(std::ostream& os) = 0;
    virtual const bool search_circref(FlxFunction* fcr) { return false; }
    virtual const std::string write() = 0;
};

class FLXLIB_EXPORT FlxString_String : public FlxString_Base {
  private:
    const std::string strV;
    const bool isWord;
  public:
    FlxString_String(const std::string strV, const bool isWord) : strV(strV), isWord(isWord) {}
    void eval(std::ostream& os) { os << strV;}
    const std::string write();
};

class FLXLIB_EXPORT FlxString_Fun : public FlxString_Base {
  public:
    enum otype { ot_dbl, ot_int, ot_boost };
  private:
    FlxFunction* fun;
    const otype id;
    const std::string boost_str;
  public:
    FlxString_Fun(FlxFunction* fun, const otype id, const std::string& boost_str) : fun(fun), id(id), boost_str(boost_str) {}
    ~FlxString_Fun() { delete fun; }
    void eval(std::ostream& os);
    const bool search_circref(FlxFunction* fcr) { return fun->search_circref(fcr); }
    const std::string write();
    
    static const otype parse_ot(const std::string& str);
    static const std::string parse_ot(const otype ot);
};

class FLXLIB_EXPORT FlxString_StrFun : public FlxString_Base {
  private:
    StringFunBase* strFun;
  public:
    FlxString_StrFun(StringFunBase* strFun) : strFun(strFun) {}
    ~FlxString_StrFun() { delete strFun; }
    void eval(std::ostream& os) { strFun->eval(os); }
    const std::string write() { return "$" + strFun->write(); }
};

class FLXLIB_EXPORT FlxString : public FlxReaderBase2 {
  private:
    std::list<FlxString_Base*>* strLst;
    /**
    * @brief instances of this FlxString (used to copy this FlxString)
    */
    tuint* instances;
    bool multiline;
    
    void free_mem();
    
    static FlxStringFunBox* StrFunBox;
  public:
    FlxString(FlxString_Base* strB, const bool multiline);
    FlxString(const bool multiline, const bool errSerious);
    FlxString(const FlxString &FlxStrF);
    FlxString& operator=(const FlxString &FlxStrF);
    ~FlxString();
    const std::string eval(const bool lowercase);
    void eval(std::ostream& os);
    const std::string eval_word(const bool lowercase, const bool emptyAllow=false,const bool numbegallow=false);
    const bool search_circref(FlxFunction* fcr);
    const std::string write();
    /**
    * @brief set expression of strA to this class; and deletes strA;
    */
    void assign ( FlxString *strA );
    /**
    * @brief checks if the string is a function or a constant
    */
    const bool is_static() const;
    
    static void set_StrFunBox(FlxStringFunBox* StrFunBoxV) { StrFunBox = StrFunBoxV; }
};

