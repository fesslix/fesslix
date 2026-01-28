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

#ifndef fesslix_flxfunction_H
#define fesslix_flxfunction_H

class FlxFunction;

#include "flxrandom_base.h"


/**
* @brief Returns a value of type tnlong from tdouble - Descr: what must not be a negative number?
*/
FLXLIB_EXPORT const tnlong tnlong_from(const tdouble d, const std::string Descr, const bool errSerious=true, const bool zero_is_allowed=true, FunBase* fun=NULL);
FLXLIB_EXPORT const tulong tulong_from(const tdouble d, const std::string Descr, const bool errSerious=true, const bool zero_is_allowed=true, FunBase* fun=NULL);
FLXLIB_EXPORT const tuint tuint_from(const tdouble d, const std::string Descr, const bool errSerious=true, const bool zero_is_allowed=true, FunBase* fun=NULL);

class FlxFunctionBox;

/**
* @brief A class for storing variables (the function to evaluate is stored).
*/
class FLXLIB_EXPORT FlxVarBox {
  private:
    #if FLX_DEBUG
      /**
      * @brief number of instances of this calss - only one class may exist
      */
      static int Cinst;
    #endif
    std::map<std::string, FlxFunction*> box;
  public:
    FlxVarBox();
    ~FlxVarBox ();
    void insert ( const std::string& name, FlxFunction* fun);
    /**
    * @brief returns the function stored in the variable (NULL if the variable does not exist)
    */
    FlxFunction* get ( const std::string& name );
    /**
    * @brief returns the name of the variable based on the address of the function
    */
    const std::string get ( FlxFunction* fun);
    void declareV ( const std::string& name );
};

class FLXLIB_EXPORT FlxBoxBase : public FlxBoxBaseR {
  protected:
    static FlxConstantBox *ConstantBox;
    static FlxVarBox *VarBox;
    static FlxFunctionBox* funBox;
  public:
    virtual ~FlxBoxBase() {};
    static void set_funBox( FlxFunctionBox* funBoxV );
    static void set_constBox( FlxConstantBox *ConstantBoxV );
    static void set_varBox( FlxVarBox *VarBoxV );
};

class Fun_OptimizeInfo {
  public:
    Fun_OptimizeInfo() : computeConstants(false) {}
    Fun_OptimizeInfo(const bool computeConstants) : computeConstants(computeConstants) {}
    const bool computeConstants;
};

class FunBase;
typedef FunBase* FunBasePtr;
/**
* @brief The base class for all other arithmetic classes.
*/
class FLXLIB_EXPORT FunBase : public FlxBoxBaseR {
  protected:
    /**
    * @brief calculates the numerical value of this function (used within optimization)
    * @param numref returns the pointer to a FunNumber object
    */
    void calc_me(FunBasePtr& numref);
    /**
    * @brief check if a function is a double value (number)
    * @return true if ftc is of type FunNumber
    * @param ftc the function to check
    */
    const bool is_number(FunBase* ftc);
  public:
    virtual const tdouble calc() = 0;
    virtual ~FunBase () {};
    virtual const bool search_circref(FlxFunction* fcr) = 0;
    virtual const std::string write() = 0;
    virtual const tuint precedence() = 0;
    /**
    * @brief Determines if the expression depends on a specific const-variable
    */
    virtual const bool dependOn_Const(const tdouble* const thenumber) = 0;
    /**
    * @brief Determines on which random variables the expression depends
    * @return true, if the function has been optimized
    * @param optf replace optimized function with optf if it is not NULL
    */
    virtual const bool optimize(FunBasePtr& optf, const Fun_OptimizeInfo &foi) = 0;
    
    /**
    * @brief optimizes a sub-FunBase-Function
    * @param optf a reference to the pointer of the function to optimize
    */
    static void child_optimize(FunBasePtr& optf, const Fun_OptimizeInfo &foi);
    virtual const bool evalw() { return false; }
};

class FLXLIB_EXPORT FunDummy : public FunBase {
  public:
    const tdouble calc() { return ZERO; }
    const bool search_circref(FlxFunction* fcr) {return false;};
    const std::string write();
    const tuint precedence() { return 0; }
    const bool dependOn_Const(const tdouble* const thenumber) { return false; }
    const bool optimize(FunBasePtr& optf, const Fun_OptimizeInfo &foi);
};

/**
* @brief The base class for all arithmetic classes with one operand.
*/
class FLXLIB_EXPORT FunBaseOperat1 : public FunBase {
  protected:
    FunBase *child_1;
  public:
    FunBaseOperat1 (FunBase *child1) : child_1(child1) {};
    virtual ~FunBaseOperat1 () { delete child_1; }
    const bool search_circref(FlxFunction* fcr) { return child_1->search_circref(fcr); }
    const bool dependOn_Const(const tdouble* const thenumber) { return child_1->dependOn_Const(thenumber); }
    const bool optimize(FunBasePtr& optf, const Fun_OptimizeInfo &foi);
    virtual const bool evalw() { return child_1->evalw(); }
};

/**
* @brief The base class for all arithmetic classes with two operands.
*/
class FLXLIB_EXPORT FunBaseOperat2 : public FunBase {
  protected:
    FunBase *child_1;
    FunBase *child_2;
  public:
    FunBaseOperat2 (FunBase *child1, FunBase *child2) : child_1(child1), child_2(child2) {};
    virtual ~FunBaseOperat2 () {
      if (child_1!=NULL) delete child_1;
      if (child_2!=NULL) delete child_2;
    }
    const bool search_circref(FlxFunction* fcr) { return child_1->search_circref(fcr) || child_2->search_circref(fcr); }
    const bool dependOn_Const(const tdouble* const thenumber) { return child_1->dependOn_Const(thenumber) || child_2->dependOn_Const(thenumber); }
    const bool optimize(FunBasePtr& optf, const Fun_OptimizeInfo &foi);
    /**
    * @brief do not use this function unless you are sure what you are doing
    */
    void set_childs_zero() { child_1 = NULL; child_2 = NULL; }
    void reset_childs(FunBase *child1V, FunBase *child2V) { child_1 = child1V; child_2 = child2V; }
    virtual const bool evalw() { return child_1->evalw() || child_2->evalw(); }
};

/**
* @brief The base class for all arithmetic classes with tree operands.
*/
class FLXLIB_EXPORT FunBaseOperat3 : public FunBase {
  protected:
    FunBase *child_1;
    FunBase *child_2;
    FunBase *child_3;
  public:
    FunBaseOperat3 (FunBase *child1, FunBase *child2, FunBase *child3) : child_1(child1), child_2(child2), child_3(child3) {};
    virtual ~FunBaseOperat3 () {
      delete child_1;
      delete child_2;
      delete child_3;
    }
    const bool search_circref(FlxFunction* fcr) { return child_1->search_circref(fcr) || child_2->search_circref(fcr) || child_3->search_circref(fcr); }
    const bool dependOn_Const(const tdouble* const thenumber) { return child_1->dependOn_Const(thenumber) || child_2->dependOn_Const(thenumber) || child_3->dependOn_Const(thenumber); }
    const bool optimize(FunBasePtr& optf, const Fun_OptimizeInfo &foi);
    virtual const bool evalw() { return child_1->evalw() || child_2->evalw() || child_3->evalw(); }
};

/**
* @brief The base class for all arithmetic classes of -functions-.
*/
class FLXLIB_EXPORT FunBaseFun : public FunBase {
  public:
    virtual ~FunBaseFun () {}
    virtual const std::string write_v() = 0;
    const tuint precedence() { return 0; }
};

/**
* @brief The base class for all arithmetic classes of -functions- with one parameter
*/
class FLXLIB_EXPORT FunBaseFun_onePara : public FunBaseFun {
  protected:
    FunBase *child_1;
  public:
    FunBaseFun_onePara (FunBase *child_1) : child_1(child_1) {};
    FunBaseFun_onePara (std::vector<FunBase*> *ParaListV);
    virtual ~FunBaseFun_onePara();
    virtual const bool search_circref(FlxFunction* fcr) { return child_1->search_circref(fcr); }
    const std::string write();
    const bool dependOn_Const(const tdouble* const thenumber) { return child_1->dependOn_Const(thenumber); }
    virtual const bool optimize(FunBasePtr& optf, const Fun_OptimizeInfo &foi);
    virtual const bool evalw() { return child_1->evalw(); }
};

/**
* @brief The base class for all arithmetic classes of -functions- with multiple parameters.
*/
class FLXLIB_EXPORT FunBaseFun_multPara : public FunBaseFun {
  protected:
    FunBase** ParaListP;
    std::vector<FunBase*> *ParaList;
  public:
    FunBaseFun_multPara (std::vector<FunBase*> *ParaListV);
    virtual ~FunBaseFun_multPara();
    virtual const bool search_circref(FlxFunction* fcr);
    const std::string write();
    const bool dependOn_Const(const tdouble* const thenumber);
    virtual const bool optimize(FunBasePtr& optf, const Fun_OptimizeInfo &foi);
    virtual const bool evalw();
};

/**
* @brief The base class for all arithmetic classes of -functions- with multiple parameters.
*/
class PYBIND11_EXPORT FunBaseFun_Python : public FunBaseFun_multPara {
  protected:
    const std::string pyFunName;
    py::function pyfunc;
    const tuint ParaN;

  public:
    FunBaseFun_Python (const std::string& pyFunName, py::function pyfunc, const tuint ParaN);
    virtual ~FunBaseFun_Python() {}
    const tdouble calc();
    const std::string write_v() { return "pyfun"; }
    const std::string write();
    virtual const bool optimize(FunBasePtr& optf, const Fun_OptimizeInfo &foi);
    virtual const bool evalw() { return false; }
};

/**
* @brief arithmetic class: stores a number
*/
class FLXLIB_EXPORT FunNumber : public FunBase {
  private:
    tdouble thenumber;
  public:
    FunNumber (tdouble TheNumber) : thenumber(TheNumber) {};
    const tdouble calc() { return thenumber; }
    const bool search_circref(FlxFunction* fcr) { return false; }
    const std::string write() { return GlobalVar.Double2String(thenumber); }
    const tuint precedence() { return 0; }
    const bool dependOn_Const(const tdouble* const thenumber) { return false; }
    const bool optimize(FunBasePtr& optf, const Fun_OptimizeInfo &foi) { return false; }
    
    const tdouble& get_thenumber() { return thenumber; }
    /**
    * @brief do not use this function unless you are sure what you are doing
    */
    void set_thenumber(const tdouble &tn) { thenumber = tn; }
};

/**
* @brief arithmetic class: stores a index of a paramter to use
*/
class FLXLIB_EXPORT FunPara : public FunBase {
  private:
    const tuint theindex;
  public:
    static const tdouble* ParaList;
    static tuint ParaListSize;
    FunPara ( const tuint TheIndex) : theindex(TheIndex) {};
    const tdouble calc();
    const bool search_circref(FlxFunction* fcr) { return false; };
    const std::string write();
    const tuint precedence() { return 0; }
    const bool dependOn_Const(const tdouble* const thenumber) { return false; }
    const bool optimize(FunBasePtr& optf, const Fun_OptimizeInfo &foi) { return false; }
};

/**
* @brief arithmetic class: user defined function
*/
class FLXLIB_EXPORT FunUser : public FunBaseFun_multPara {
  private:
    FlxFunction* fun;
    const std::string& fname;
    const tuint numbofpara;
    tVec tPL;        // constains the values of the parameters
    tdouble* const tPLp;
  public:
    FunUser (std::vector<FunBase*> *ParaListV, FlxFunction* funV, const std::string& fname, const tuint numbofpara);
    const tdouble calc();
    const bool search_circref(FlxFunction* fcr);
    const std::string write_v() { return fname; }
    const std::string write();
    const bool dependOn_Const(const tdouble* const thenumber);
    const bool optimize(FunBasePtr& optf, const Fun_OptimizeInfo &foi) { return false; }
};

class FLXLIB_EXPORT FunDummyFun : public FunBaseFun_multPara {
  public:
    FunDummyFun (std::vector<FunBase*> *ParaListV) : FunBaseFun_multPara(ParaListV) {};
    const tdouble calc() { return ZERO; }
    const std::string write_v();
    const bool optimize(FunBasePtr& optf, const Fun_OptimizeInfo &foi);
};


/**
* @brief Each arithmetic class needs a reader-class - this is the base class for all other arithmetic read classes.
*/
class FLXLIB_EXPORT FunReadBase : public FlxReaderBase, public FlxBoxBase {
  protected:
    FunReadBase *Next;
    const int priority;
    FunReadBase* insert ( FunReadBase *TheFunRead );
    static FunReadBase *startLink;
  public: 
    FunReadBase (const int Priority, bool isEndNode = false);
    virtual ~FunReadBase () { delete Next; };
    virtual FunBase* read ( bool errSerious ) = 0;
};

/**
* @brief if there is no suitable arithmetic read class found, this 'read'-function will be used (end of the list of arithmetic read classes)
*/
class FLXLIB_EXPORT FunReadEND : public FunReadBase {
  protected:
    FunBase* read ( bool errSerious );
  public:
    FunReadEND () : FunReadBase(-1, true) {}
};

/**
* @brief head of the list of arithmetic read classes
*/
class FLXLIB_EXPORT FunReadSTART : public FunReadBase {
  public:
    FunReadSTART (  );
    void insert ( FunReadBase *TheFunRead ) { FunReadBase::insert(TheFunRead); }
    FunBase* read ( bool errSerious );
    ReadStream* get_reader();
};

/**
* @brief arithmetic read class of FunNumber
*/
class FLXLIB_EXPORT FunReadNumber : public FunReadBase {
  private:
    FunBase* read ( bool errSerious );
  public:
    FunReadNumber () : FunReadBase(1000) {}
};

/**
* @brief arithmetic read class of FunPara
*/
class FLXLIB_EXPORT FunReadPara : public FunReadBase {
  private:
    FunBase* read ( bool errSerious );
    static tuint numbofpara;
  public:
    FunReadPara () : FunReadBase(1020) {}
    static void set_NumbOfPara(const tuint NumbOfPara);
};



/**
* @brief This class builds an arithmetic function out of a stream
*/
class FLXLIB_EXPORT FlxFunctionReader {
  private:
    #if FLX_DEBUG
      /**
      * @brief number of instances of this calss - only one class may exist
      */
      static int Cinst;
    #endif
    /**
    * @brief A list containing all arithmetic read classes.
    */
    FunReadSTART* FunctionList;
  public:
    /**
    * @brief Create a FlxFunctionReader class.
    *
    * An instance of each arithmetic read class has to be created !!!
    *
    * @param StrmReader The Stream-Reader class
    */
    FlxFunctionReader ( );                 // FunReadXXXXXX Objects have to be listed there
    virtual ~FlxFunctionReader () {
      delete FunctionList;
    }

    /**
    * @brief builds an arithmetic function (a ReadStream class (see constructor) is used for input)
    */
    FunBase* read (bool errSerious = true);
    ReadStream* get_reader();
};

class FLXLIB_EXPORT FunReadFunBase_0 {
  protected:
    static FunReadSTART* FunctionList;
  public:
    virtual ~FunReadFunBase_0() {};
    static void set_data ( FunReadSTART* FunctionListV) { FunctionList = FunctionListV; }
};

class FLXLIB_EXPORT FunReadFunBase : public FlxReaderBase, public FunReadFunBase_0, public FlxBoxBase {
  private:
    /**
    * @brief true: Fesslix will not delete this object - this has to be done be the corresponding external library
    */
    const bool from_library;
  protected:
    /**
    * @brief reads the passed arguments of the function
    * @param NumbOfPara -1: read as many as passed; >=0: number of passed arguments has to match!
    */
    std::vector<FunBase*>* read_parameters( const int NumbOfPara, const bool errSerious );
    tdoublePtr read_const_var(const bool errSerious, const bool define=false);
  public:
    FunReadFunBase(const bool from_library = false) : from_library(from_library) {}
    virtual ~FunReadFunBase() {};
    virtual FunBase* read ( bool errSerious ) = 0;
    
    /**
    * @brief is executed when the function is added to the function list
    */
    virtual void initialize() {}
    const bool is_from_library() const { return from_library; }
};

class FLXLIB_EXPORT FunReadFunUser : public FunReadFunBase {
  private:
    const std::string fname;
    FlxFunction* fun;
    const tuint numbofpara;
  public:
    FunReadFunUser ( const std::string fname, FlxFunction* funV, int NumbOfPara) : fname(fname), fun(funV), numbofpara(NumbOfPara) {};
    ~FunReadFunUser();
    FunBase* read ( bool errSerious ) { return new FunUser( read_parameters( int(numbofpara), errSerious ), fun, fname, numbofpara ); }
    
    FlxFunction* get_fun_ptr() { return fun; }
    const tuint get_numbofpara() const { return numbofpara; }
};

class FLXLIB_EXPORT FunReadFunDummy : public FunReadFunBase {
  public:
    FunBase* read ( bool errSerious );
};


/**
* @brief A class for storing pre-defined functions
*/
class FLXLIB_EXPORT FlxFunctionBox {
  private:
    #if FLX_DEBUG
      /**
      * @brief number of instances of this calss - only one class may exist
      */
      static int Cinst;
    #endif
    std::map<std::string, FunReadFunBase*> box;
  public:
    FlxFunctionBox ();
    ~FlxFunctionBox ();
    void insert ( const std::string& name, FunReadFunBase* fR);
    FunReadFunBase* get ( const std::string& name);
    /**
    * @brief removes a reader - without deleting it
    * @note be carefull - this is intented for libraries that link to Fesslix
    */
    void remove( const std::string& name );
    void declareF ( const std::string& name );
};



/**
* @brief Using (not developing) this library, this is the class you will have to work with.
*
* FlxFunction *fun1 = new FlxFunction(funReader);
* 
*/
class FLXLIB_EXPORT FlxFunction : public FunBase {
  protected:
    /**
    * @brief the arithmetic function
    */
    FunBase *fun;
    /**
    * @brief instances of this function (used to copy this function)
    */
    tuint* instances;
    /**
    * @brief position of start of function
    */
    FlxReaderPos* read_pos;
    /**
    * @brief of no particular meaning ...
    */
    const tuint precedence() { return 0; }
    
    FlxFunction() : fun(NULL),instances(NULL) {};
    /**
    * @brief check if it is really a true FlxFunction
    */
    void check_FlxFunction(const FlxFunction *funtocheck);
    const bool optimize(FunBasePtr& optf, const Fun_OptimizeInfo &foi);
  public:
    FlxFunction (FlxFunctionReader *funReader, bool errSerious = true );
    FlxFunction (FunBase* funV) : fun(funV),instances(new tuint(0)),read_pos(NULL) {}
    FlxFunction (const FlxFunction& FlxF);
    FlxFunction& operator=( const FlxFunction& FlxF);
    virtual ~FlxFunction ();
    /**
    * @brief calculates the arithmetic function
    */
    const tdouble calc() { return fun->calc(); }
    /**
    * @brief set arithmetic expression of funA to this class; and deletes funA;
    */
    virtual void assign ( FlxFunction *funA );
    /**
    * @brief returns the arithmetic function stored as a string
    */
    const std::string write() { return fun->write(); }
    /**
    * @brief check if a 'var' definition causes a circular reference
    * @return true - if a circular reference was found
    */
    const bool search_circref(FlxFunction* fcr) { if (fcr==NULL) {return false; } else { return fun->search_circref(fcr);} }
    /**
    * @brief Determines if the expression depends on a specific const-variable
    */
    const bool dependOn_Const(const tdouble* const thenumber) { return fun->dependOn_Const(thenumber); }
    
    /**
    * @brief evaluates the result of the function and returns a positiv integer (zero is not allowed)
    */
    const tuint cast2tuint(const bool errSerious=true) ;
    /**
    * @brief evaluates the result of the function and returns a positiv integer (zero is allowed)
    */
    const tuint cast2tuintW0(const bool errSerious=true) ;
    /**
    * @brief evaluates the result of the function and returns a positiv integer (zero is not allowed)
    */
    const tnlong cast2tnlong(const bool errSerious=true);
    /**
    * @brief evaluates the result of the function and returns a positiv integer (zero is allowed)
    */
    const tnlong cast2tnlongW0(const bool errSerious=true);
    /**
    * @brief evaluates the result of the function and returns a positiv integer (zero is not allowed)
    */
    const tulong cast2tulong(const bool errSerious=true);
    /**
    * @brief evaluates the result of the function and returns a positiv integer (zero is allowed)
    */
    const tulong cast2tulongW0(const bool errSerious=true);
    /**
    * @brief cast the result to a positive value (>0.0)
    * @return a value t with t > 0.0;
    */
    const tdouble cast2positive(const bool errSerious=true);
    const tdouble cast2positive_or0(const bool errSerious=true);
    const int cast2int( ) { return static_cast<int>(fun->calc()); }
    /**
    * @brief check if the function is zero
    */
    const bool is_zero();
    /**
    * @brief returns a reference to the inner funBase-Element
    */
    FunBase& get_ptrFun() { return *fun; }
    
    friend class FlxFunction_Combine_Add;
};

typedef FlxFunction* FlxFunctionPtr;
inline FlxFunctionPtr double2flx(const tdouble d) { return new FlxFunction(new FunNumber(d)); }

class FLXLIB_EXPORT FlxFunction_Combine_Add : public FlxFunction {
  private:
    FlxFunction *f1;
    FlxFunction *f2;
    // dummys
      FlxFunction_Combine_Add(const FlxFunction_Combine_Add& FlxF);
      FlxFunction_Combine_Add& operator=( const FlxFunction_Combine_Add& FlxF);
  public:
    /**
    * @brief combines two functions by addition. 
    * NOTE: do not delete f1V and f2V
    * NOTE: do not use both functions afterwards
    */
    FlxFunction_Combine_Add ( FlxFunction* f1V, FlxFunction* f2V);
    ~FlxFunction_Combine_Add();
    virtual void assign ( FlxFunction *funA );
};

/**
* @brief stores a function and its approximate polynomial degree
* (FlxFunction,Degree)
*/
class FLXLIB_EXPORT FlxFunDeg {
  private:
    tuint deg;
    FlxFunction* fun;
    static tuint default_deg;
  public:
    FlxFunDeg(ReadStream *reader, FlxFunctionReader *funReader, bool errSerious=true);
    FlxFunDeg(FlxFunction* funF);
    FlxFunDeg(const FlxFunDeg& ref);
    ~FlxFunDeg() { if (fun) delete fun; }
    const tuint get_deg() { return deg; }
    FlxFunction* get_fun() { return fun; }
    const std::string write();
    void assign(FlxFunDeg* funA);
    const bool is_zero() { return fun->is_zero(); }
    static void set_deg(const tuint degV) { default_deg = degV; }
};

#endif // fesslix_flxfunction_H
