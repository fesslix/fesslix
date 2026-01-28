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

#pragma once

#include "flxdata.h"



class FLXLIB_EXPORT FlxOptionalParaBase : public FlxReaderBase2 {
  protected:
    std::string pName;
    bool is_set;
    virtual void set(void* valueP) = 0;
  public:
    FlxOptionalParaBase(std::string pName) : pName(pName), is_set(false) {}
    std::string get_name() const { return pName; };
    virtual ~FlxOptionalParaBase() {};
    virtual void read(bool errSerious = true);
    virtual void* read_value(bool errSerious) = 0;
    virtual void* get() = 0;
    virtual void set_default(void* defV) = 0;
    virtual void free_value(void* valueP) = 0;
};


class FLXLIB_EXPORT FlxOptionalParaFun : public FlxOptionalParaBase {
  private:
    /**
    * @brief The default value; this value will be deleted by the destructor
    */
    FlxFunction* defv;
    /**
    * @brief The value; this value will be deleted by the destructor
    */
    FlxFunction *value;
    /**
    * @brief Set a new function;
    */
    void set(void* valueP);
  public:
    FlxOptionalParaFun(FlxFunction* defV, const std::string pName) : FlxOptionalParaBase(pName), defv(defV), value(NULL) {}
    FlxOptionalParaFun(const tdouble defV, const std::string pName);
    ~FlxOptionalParaFun();
    /**
    * @brief Read the function - and allocate memory;
    */
    void* read_value(bool errSerious) { return ( new FlxFunction(funReader,errSerious)); }
    /**
    * @brief Return the function - the function is copied and has to be deleted by some other function;
    */
    void* get();
    FlxFunction* get_ref();
    /**
    * @brief Set a new default function;
    */
    void set_default(void* defV);
    /** 
    * @brief Delete valueP (used by flxObjDefault);
    */
    void free_value(void* valueP);
};

class FLXLIB_EXPORT FlxOptionalParaFlxString : public FlxOptionalParaBase {
  private:
    /**
    * @brief The default value; this value will be deleted by the destructor
    */
    FlxString* defv;
    /**
    * @brief The value; this value will be deleted by the destructor
    */
    FlxString *value;
    /**
    * @brief Set a new FlxString;
    */
    void set(void* valueP);
  public:
    FlxOptionalParaFlxString(std::string defV, std::string pName, const bool is_Word);
    FlxOptionalParaFlxString(FlxString* defV, std::string pName) : FlxOptionalParaBase(pName), defv(defV), value(NULL) {}
    
    ~FlxOptionalParaFlxString();
    /**
    * @brief Read the FlxString - and allocate memory;
    */
    void* read_value(bool errSerious) { return ( new FlxString(false,errSerious)); }
    /**
    * @brief Return the FlxString - the FlxString is copied and has to be deleted by some other function;
    */
    void* get();
    FlxString* get_ref();
    /**
    * @brief Set a new default FlxString;
    */
    void set_default(void* defV);
    /** 
    * @brief Delete valueP (used by flxObjDefault);
    */
    void free_value(void* valueP);
};

class FLXLIB_EXPORT FlxOptionalParaMtxFun : public FlxOptionalParaBase {
  private:
    /**
    * @brief The default value; this value will be deleted by the destructor
    */
    FlxMtxConstFun* defv;
    /**
    * @brief The value; this value will be deleted by the destructor
    */
    FlxMtxConstFun *value;
    /**
    * @brief Set a new function;
    */
    void set(void* valueP);
  public:
    FlxOptionalParaMtxFun(FlxMtxConstFun* defV, std::string pName) : FlxOptionalParaBase(pName), defv(defV), value(NULL) {}
    ~FlxOptionalParaMtxFun();
    /**
    * @brief Read the function - and allocate memory;
    */
    void* read_value(bool errSerious) { return ( new FlxMtxConstFun(true) ); }
    /**
    * @brief Return the function - the function is copied and has to be deleted by some other function;
    */
    void* get();
    /**
    * @brief Set a new default function;
    */
    void set_default(void* defV);
    /** 
    * @brief Delete valueP (used by flxObjDefault);
    */
    void free_value(void* valueP);
};

class FLXLIB_EXPORT FlxOptionalParaStream : public FlxOptionalParaBase {
  protected:
    std::string defv;
    std::string *value;
    void set(void* valueP);
  public:
    /**
    * @param defV MUST BE LOWERCASE
    */
    FlxOptionalParaStream(const std::string& defV, const std::string& pName) : FlxOptionalParaBase(pName), defv(defV), value(NULL) { }
    virtual ~FlxOptionalParaStream();
    virtual void* read_value(const bool errSerious);
    virtual void* get();
    virtual void set_default(void* defV);
    virtual void free_value(void* valueP);
};

class FLXLIB_EXPORT FlxOptionalParaBool : public FlxOptionalParaBase {
  private:
    bool defv;
    bool *value;
    void set(void* valueP);
  public:
    FlxOptionalParaBool(bool defV, std::string pName) : FlxOptionalParaBase(pName), defv(defV), value(NULL) { }
    ~FlxOptionalParaBool();
    void* read_value(bool errSerious);
    void* get();
    void set_default(void* defV);
    void free_value(void* valueP);
};

class FLXLIB_EXPORT FlxOptionalParaString : public FlxOptionalParaStream {
  public:
    FlxOptionalParaString(std::string defV, std::string pName) : FlxOptionalParaStream(defV,pName) {};
    virtual void* read_value(const bool errSerious);
    virtual void set_default(void* defV);
};

class FLXLIB_EXPORT FlxOptionalParaText : public FlxOptionalParaStream {
  public:
    FlxOptionalParaText(std::string defV, std::string pName) : FlxOptionalParaStream(defV,pName) {};
    virtual void* read_value(bool errSerious);
};

class FLXLIB_EXPORT FlxDefParaBox {
  private:
    #if FLX_DEBUG
      /**
      * @brief number of instances of this calss - only one class may exist
      */
      static int Cinst;
    #endif
    std::map<std::string, FlxOptionalParaBase*> box;
  public:
    FlxDefParaBox();
    void insert (FlxOptionalParaBase* value);
    FlxOptionalParaBase* get(std::string name);
    ~FlxDefParaBox();
};

class FLXLIB_EXPORT FlxDefParaBoxBase {
  protected:
    static FlxDefParaBox* AllDefParaBox;
  public:
    virtual ~FlxDefParaBoxBase() {};
    static void set_AllDefParaBox ( FlxDefParaBox* AllDefParaBoxV ) { AllDefParaBox = AllDefParaBoxV; }
};

class FLXLIB_EXPORT FlxOptionalParaBox : public FlxDefParaBoxBase {
  private:
    std::map<std::string, FlxOptionalParaBase*> box;
  public:
    void insert ( std::string name, std::string defname);
    FlxOptionalParaBase* get ( std::string name );
};



/**
* @brief The base class for all other 'object read classes'
*/
class FLXLIB_EXPORT FlxObjReadBase : public FlxReaderBase2, public FlxDataBase, public FlxDefParaBoxBase {
  private:
    /**
    * @brief true: Fesslix will not delete this object - this has to be done be the corresponding external library
    */
    const bool from_library;
  protected:
    FlxOptionalParaBox ParaBox;
    void read_optionalPara(bool errSerious = true);
    const bool get_doLog();
    FlxFunction* get_optPara_FlxFunction(const std::string& str);
    const tdouble get_optPara_tdouble_from_FlxFunction(const std::string& str, const bool only_positiveORzero, const bool errSerious);
    const int get_optPara_int_from_FlxFunction(const std::string& str);
    const tuint get_optPara_tuint_from_FlxFunction(const std::string& str, const bool zero_is_allowed, const bool errSerious=true);
    FlxFunDeg* get_optPara_FlxFunDeg(const std::string& str);
    FlxMtxConstFun* get_optPara_FlxMtxFun(const std::string& str);
    const bool get_optPara_bool(const std::string& str);
    FlxString* get_optPara_FlxString(const std::string& str);
    const std::string get_optPara_string_from_FlxString(const std::string& str, const bool lowercase);
    const std::string get_optPara_word_from_FlxString(const std::string& str, const bool lowercase, const bool emptyAllow=false,const bool numbegallow=false);
    const std::string& get_optPara_string(const std::string& str, const bool lowercase=false);
  public:
    FlxObjReadBase(const bool from_library = false);
    virtual ~FlxObjReadBase () {};
    /**
    * @brief This function has to be overriden by all 'object read classes' (that is the hart)
    */
    virtual FlxObjBase* read() = 0;
    const bool is_from_library() const { return from_library; }
};

/**
* @brief object read class: base class for objects using output functionality
*/
class FLXLIB_EXPORT FlxObjReadLogBase : public FlxObjReadBase {
  protected:
    const bool get_verboseLog();
  public:
    FlxObjReadLogBase(const bool from_library = false);
};

/**
* @brief object read class: base class for objects using output functionality
*/
class FLXLIB_EXPORT FlxObjReadOutputBase : public FlxObjReadBase {
  protected:
    const std::string get_stream();
    const bool get_verbose();
    const bool get_checkTOL();
    const int get_prec();
    const int get_fixW();
    const std::string get_boost_str();
  public:
    FlxObjReadOutputBase(const bool from_library = false);
};


/**
* @brief A class for storing all FlxReadObjects (keywords in parameter-files)
*/
class FLXLIB_EXPORT FlxObjectReadBox {
  private:
    std::map<std::string, FlxObjReadBase*> box;
  public:
    FlxObjectReadBox();
    ~FlxObjectReadBox ();
    /**
    * @brief Get a FlxObject by name.
    * @param name The name of the FlxObject to return.
    * @return FlxObjBase*
    */
    FlxObjReadBase* get ( std::string name );
    /**
    * @brief removes a reader - without deleting it
    * @note be carefull - this is intented for libraries that link to Fesslix
    */
    void remove( const std::string& name );
    /**
    * @brief Insert a Flx-Read-Object.
    * @param name An unique name
    * @param value The FlxReadObject to store
    */
    void insert ( std::string name, FlxObjReadBase* value);
};


/**
* @brief This class evaluates a single command in the parameter file
*/
class FLXLIB_EXPORT FlxEvaluateCmd : public FlxDataBase, public FlxReaderBase2 {
  private:
    #if FLX_DEBUG
      /**
      * @brief number of instances of this calss - only one class may exist
      */
      static int Cinst;
    #endif
    FlxObjectReadBox ObjReadBox;
    void check_ending();
  public:
    FlxEvaluateCmd();
    FlxObjBase* evaluateCmd( );
    FlxObjectReadBox* get_ObjReadBox() { return &ObjReadBox; }
};

class FLXLIB_EXPORT FlxEvaluateCmdBase {
  protected:
    static FlxEvaluateCmd* EvaluateCmd;
  public:
    ~FlxEvaluateCmdBase() {};
    static void set_EvaluateCmd (FlxEvaluateCmd *EvaluateCmdV) {EvaluateCmd = EvaluateCmdV;};
    static FlxEvaluateCmd* get_EvaluateCmd() { return EvaluateCmd; }
};


class FLXLIB_EXPORT FlxCreateObjReaders {
  public:
    virtual void createObjReaders (FlxObjectReadBox* objReadBox ) = 0;
    virtual void createFunReaders (FlxData* dataBox ) {}
    virtual ~FlxCreateObjReaders() {}
};

/**
* @brief object class: base class for output objects
*/
class FLXLIB_EXPORT FlxObjLogBase : public FlxObjBase {
  protected:
    const bool verboseLog;
  public:
    FlxObjLogBase(bool dolog, bool verboseLog) : FlxObjBase(dolog), verboseLog(verboseLog) {}
    virtual ~FlxObjLogBase() { } 
};

/**
* @brief object class: base class for output objects
*/
class FLXLIB_EXPORT FlxObjOutputBase : public FlxObjBase {
  private:
    const std::string ostreamv;
    const bool verbose;
  protected:
    const bool checkTOL;
    const int prec;
    const int fixW;
    const std::string boost_str;        // how to format output
  public:
    FlxObjOutputBase( const bool dolog, const std::string ostreamV, const bool verboseV=false, const bool checkTOLv=true, const int prec=-1, const int fixW=0, const std::string& boost_str="") : FlxObjBase(dolog), ostreamv(ostreamV), verbose(verboseV),checkTOL(checkTOLv), prec(prec), fixW(fixW), boost_str(boost_str) {}
    virtual ~FlxObjOutputBase() { }
    std::ostream& sout () { return *(data->OstreamBox.get(ostreamv)); }    
    const bool is_verbose() const {return verbose;}
    void write(const tdouble val, std::ostream& sout_);
};

/**
* @brief object class: define a constant
*
* const constName = FlxFunction
*/
class FLXLIB_EXPORT FlxObjConst : public FlxObjBase {
  private:
    const std::string cname;
    FlxFunction *fun;
    tdouble* cptr;
    const char opc;
    void task();
  public:
    FlxObjConst ( const bool dolog, const std::string& cnameV, FlxFunction* funV, const char opc ) : FlxObjBase(dolog), cname(cnameV), fun(funV), cptr(data->ConstantBox.get(cname,true)), opc(opc) {}
    ~FlxObjConst() { delete fun; }
};



/**
* @brief object read class: for FlxObjReadFCVbase
*/
class FLXLIB_EXPORT FlxObjReadFCVbase : public FlxObjReadOutputBase {
  protected:
    /**
    * @brief checks if a function or variable with the name 'thename' already exists - if it does, an exception is thrown
    * @param thename name to check (MUST BE LOWERCASE)
    * @param theFCVindex 'F'=function; 'C'='const'variable 'V'='var'variable
    */
    void isdefined(const std::string& thename, const char theFCVindex, const bool errSerious = true) const;
};

/**
* @brief object read class: for FlxObjConst
*/
class FLXLIB_EXPORT FlxObjReadConst : public FlxObjReadFCVbase {
  public:
    FlxObjReadConst();
    FlxObjBase* read ();
    /**
    * @note cname must be lowercase!
    */
    FlxObjBase* read (const std::string& cname, const bool allow_operators=false);
};

/**
* @brief object class: procedure - sub
*
* sub name { ... };
*/
class FLXLIB_EXPORT FlxObjSub : public FlxObjBase {
  private:
    FlxObjBase* InternListSub;
    const std::string sub_name;
    void task ();
  public:
    FlxObjSub ( const bool dolog, const std::string& sub_nameV, FlxObjBase* InternListSubV ) 
      : FlxObjBase(dolog), InternListSub(InternListSubV), sub_name(sub_nameV) {}
    ~FlxObjSub() { if ( InternListSub != NULL) delete InternListSub; }    
};

/**
* @brief object class: procedure - sub (execution)
*
* sub_name(); 
* procedure sub_name();
*/
class FLXLIB_EXPORT FlxObjSubDo : public FlxObjBase {
  private:
    FlxObjBase* InternListSub;
    const std::string sub_name;
    void task();
  public:
    FlxObjSubDo ( bool dolog, FlxObjBase* InternListSubV, const std::string &sub_name )
      : FlxObjBase(dolog), InternListSub(InternListSubV), sub_name(sub_name) {}
    ~FlxObjSubDo() { }    
};

class FLXLIB_EXPORT FlxObjReadCodeBlock : public FlxObjReadBase, public FlxEvaluateCmdBase {
  protected:
    FlxObjReadCodeBlock(const bool from_library = false) : FlxObjReadBase(from_library) {}
  public:
    static FlxCodeBlock* read_block (const bool UseIgnoreList, const bool errSerious=true);
};

/**
* @brief object read class: for FlxSub
*/
class FLXLIB_EXPORT FlxObjReadSub : public FlxObjReadBase {
  public:
    FlxObjReadSub() { data->IgnoreBox.set_iL_recur("sub"); }
    FlxObjBase* read();
};

/**
* @brief object read class: for FlxSub
*/
class FLXLIB_EXPORT FlxObjReadSubDo : public FlxObjReadBase {
  public:
    FlxObjBase* read();
    FlxObjBase* read(FlxObjBase* objProcedure, const std::string& sub_name);
};

/**
* @brief object read base class: for Loops
*/
class FLXLIB_EXPORT FlxObjReadLoops : public FlxObjReadBase {
  protected:
    const tuint get_maxpasses();
  public:
    FlxObjReadLoops ();
};


/**
* @brief object class: output a string
*
* echo "text"
*/
class FLXLIB_EXPORT FlxObjEcho : public FlxObjOutputBase {
  private:
    FlxString* strV;
    const bool newline;
    void task() { strV->eval(sout()); if (newline) sout() << std::endl; }
  public:
    FlxObjEcho ( bool dolog, FlxString* strV, std::string ostreamV, const bool newline ) : FlxObjOutputBase(dolog,ostreamV), strV(strV), newline(newline) {}
    ~FlxObjEcho() { delete strV; }
};


typedef const tdouble (*flx_lsf_callback)(const tdouble* const, const tuint);
/**
* @brief arithmetic class: a callback limit-state function (can be defined by the user)
*/
class FLXLIB_EXPORT FunLSF_callback : public FunBaseFun, public FlxDataBase {
  private:
    flx_lsf_callback lsfp;
    const std::string lsf_name;
    RBRV_constructor* RndBox;
    const tuint NOX;
    flxVec rvv;
    std::string rbrvsetn;
  public:
    FunLSF_callback (flx_lsf_callback lsfp, const std::string& lsf_name, RBRV_constructor* RndBox, std::string rbrvsetn);
    ~FunLSF_callback();
    const tdouble calc();
    const bool search_circref(FlxFunction* fcr) { return false; }
    const std::string write_v() { return lsf_name; }
    const std::string write();
    const bool dependOn_Const(const tdouble* const thenumber) { return false; }
    const bool optimize(FunBasePtr& optf, const Fun_OptimizeInfo &foi) { return false; }
};

class FLXLIB_EXPORT FunReadFunLSF_callback : public FunReadFunBase, public FlxDataBase {
  private:
    const tdouble (*lsfp)(const tdouble* const, const tuint);
    std::string lsf_name;
  public:
    FunReadFunLSF_callback(flx_lsf_callback lsfp, std::string lsf_nameV, const bool from_library);
    FunBase* read ( bool errSerious );
};



