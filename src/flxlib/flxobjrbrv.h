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

#include "flxobjects.h"
#include "flxrbrv_rvs_read.h"


class FLXLIB_EXPORT FlxCreateObjReaders_RBRV : public FlxCreateObjReaders {
  public:
    void createObjReaders (FlxObjectReadBox* objReadBox );
    void createFunReaders (FlxData* dataBox );
};

//------------------------------------------------------------

class FLXLIB_EXPORT RBRV_vfset : public RBRV_set_parents, public FlxDataBase {
  protected:
    const tuint Nentries;        // number of entries in the set
    flxVec x_of_set;                // current x-vector of the set
    FlxMtxConstFun* vecfun;        //the function to evaluate
    #if FLX_DEBUG
      bool valid;
    #endif
  public:
    RBRV_vfset(const bool internal, const std::string& name, const bool noID, const tuint Nentries, FlxMtxConstFun* vecfun, const tuint Nparents, RBRV_set_base** const parents);
    virtual ~RBRV_vfset() { delete vecfun; }
    virtual const tuint get_NRV() const { return 0; }
    virtual const tuint get_NOX() const { return Nentries; }
    virtual const tuint get_NRV_only_this() const { return 0; }
    virtual const tuint get_NOX_only_this() const { return Nentries; }
    #if FLX_DEBUG
    virtual void set_is_valid(const bool is_valid) { valid = is_valid; }
    #else
    virtual void set_is_valid(const bool is_valid) {}
    #endif
    virtual void transform_y2x();
    virtual const bool allow_x2y() const { return false; }
    virtual void set_x(const tdouble* const x_vec);
    virtual void set_x_only_this(const tdouble* const x_vec) { set_x(x_vec); }
    virtual void get_x(tdouble* const x_vec);
    virtual void get_x_only_this(tdouble* const x_vec) { get_x(x_vec); }
    virtual const bool check_xVec(const tdouble* xp) { return true; }
    virtual void get_mean(tdouble* const m_vec);
    virtual void get_mean_only_this(tdouble* const m_vec) { get_mean(m_vec); };
    virtual void get_sd(tdouble* const s_vec);
    virtual void get_sd_only_this(tdouble* const s_vec) { get_sd(s_vec); };
    void print(std::ostream& sout, const std::string prelim, tuint& counter, const bool printID);
};
//------------------------------------------------------------

class FLXLIB_EXPORT RBRV_dirichlet : public RBRV_set_parents, public FlxDataBase {
  protected:
    const tuint Nentries;        // number of entries in the set
    flxVec x_of_set;                // current x-vector of the set
    flxVec alpha_vec;                // current x-vector of the set
    FlxMtxConstFun* vecfun;        //the function to evaluate
    #if FLX_DEBUG
      bool valid;
    #endif

    virtual void get_pars();
  public:
    RBRV_dirichlet(const bool internal, const std::string& name, const bool noID, const tuint Nentries, FlxMtxConstFun* vecfun, const tuint Nparents, RBRV_set_base** const parents, flxVec* avec=NULL, const tuint idim=0);
    virtual ~RBRV_dirichlet();
    virtual const tuint get_NRV() const { return sRV; }
    virtual const tuint get_NOX() const { return Nentries; }
    virtual const tuint get_NRV_only_this() const { return get_NRV(); }
    virtual const tuint get_NOX_only_this() const { return get_NOX(); }
    #if FLX_DEBUG
    virtual void set_is_valid(const bool is_valid) { valid = is_valid; }
    #else
    virtual void set_is_valid(const bool is_valid) {}
    #endif
    virtual void transform_y2x();
    virtual const bool allow_x2y() const { return false; }
    virtual void set_x(const tdouble* const x_vec);
    virtual void set_x_only_this(const tdouble* const x_vec) { set_x(x_vec); }
    virtual void get_x(tdouble* const x_vec);
    virtual void get_x_only_this(tdouble* const x_vec) { get_x(x_vec); }
    virtual const bool check_xVec(const tdouble* xp);
    virtual void get_mean(tdouble* const m_vec);
    virtual void get_mean_only_this(tdouble* const m_vec) { get_mean(m_vec); };
    virtual void get_sd(tdouble* const s_vec);
    virtual void get_sd_only_this(tdouble* const s_vec) { get_sd(s_vec); };
    void print(std::ostream& sout, const std::string prelim, tuint& counter, const bool printID);
};

//------------------------------------------------------------

class FLXLIB_EXPORT RBRV_multinomial : public RBRV_dirichlet {
  protected:
    const tuint Ntrials;

    virtual void get_pars();
  public:
    RBRV_multinomial(const bool internal, const std::string& name, const bool noID, const tuint Nentries, FlxMtxConstFun* vecfun, const tuint Nparents, RBRV_set_base** const parents, const tuint Ntrials, flxVec* avec=NULL);
    virtual ~RBRV_multinomial();
    virtual void transform_y2x();
    virtual const bool check_xVec(const tdouble* xp);
    virtual void get_mean(tdouble* const m_vec);
    virtual void get_sd(tdouble* const s_vec);
    void print(std::ostream& sout, const std::string prelim, tuint& counter, const bool printID);
};

//------------------------------------------------------------



/**
* @brief object class: defines a RBRV_set_create
*/
class FlxObjRBRV_set_new : public FlxObjBase {
  private:
    FlxString* set_name;
    std::vector<FlxString*> set_parents;
    const bool allow_x2y;
    const bool is_Nataf;
    const bool is_Nataf_evalOnce;
    
    void task();
  public:
    FlxObjRBRV_set_new ( const bool dolog, FlxString* set_name, std::vector<FlxString*> set_parents, const bool allow_x2y, const bool is_Nataf, const bool is_Nataf_evalOnce );
    ~FlxObjRBRV_set_new();
};

/**
* @brief object read class: for FlxObjRBRV_set_new
*/
class FlxObjReadRBRV_set_new : public FlxObjReadBase {
  public:
    FlxObjReadRBRV_set_new();
    FlxObjBase* read ();
};


class FlxObjRBRV_set_addRV : public FlxObjBase {
  private:
    FlxString* set_name;
    RBRV_entry_read_base* entry;
    
    void task();
  public:
    FlxObjRBRV_set_addRV ( const bool dolog, FlxString* set_name, RBRV_entry_read_base* entry ) : FlxObjBase(dolog), set_name(set_name), entry(entry) {}
    ~FlxObjRBRV_set_addRV() { delete set_name; delete entry; }
};

/**
* @brief object read class: for FlxObjRBRV_set
*/
class FlxObjReadRBRV_set_addRV : public FlxObjReadBase {
  public:
    FlxObjBase* read ();
};


class FlxObjRBRV_set_addCorr : public FlxObjBase {
  private:
    FlxString* set_name;
    FlxString* name_rv1;
    FlxString* name_rv2;
    FlxFunction* corrVal;
    const bool corr_approx;
    const bool rhogauss;
    
    void task();
  public:
    FlxObjRBRV_set_addCorr ( const bool dolog, FlxString* set_name, FlxString* name_rv1, FlxString* name_rv2, FlxFunction* corrVal, const bool corr_approx, const bool rhogauss ) : FlxObjBase(dolog), set_name(set_name), name_rv1(name_rv1), name_rv2(name_rv2), corrVal(corrVal), corr_approx(corr_approx), rhogauss(rhogauss) {}
    ~FlxObjRBRV_set_addCorr() { delete set_name; delete name_rv1; delete name_rv2; delete corrVal; }
};

/**
* @brief object read class: for FlxObjRBRV_set_addCorr
*/
class FlxObjReadRBRV_set_addCorr : public FlxObjReadBase {
  public:
    FlxObjReadRBRV_set_addCorr();
    FlxObjBase* read ();
};


/**
* @brief object class: defines a RBRV_set
*/
class FlxObjRBRV_set_create : public FlxObjBase {
  private:
    FlxString* set_name;
    
    void task();
  public:
    FlxObjRBRV_set_create ( const bool dolog, FlxString* set_name ) : FlxObjBase(dolog), set_name(set_name) {}
    ~FlxObjRBRV_set_create() { delete set_name; }
};

/**
* @brief object read class: for FlxObjRBRV_set
*/
class FlxObjReadRBRV_set_create : public FlxObjReadBase {
  public:
    FlxObjBase* read ();
};




/**
* @brief object class: defines a RBRV_set
*/
class FlxObjRBRV_set : public FlxObjBase {
  private:
    FlxString* set_name;
    std::vector<FlxString*> set_parents;
    std::vector<RBRV_entry_read_base*> set_entries;
    const bool allow_x2y;
    
    void task();
  public:
    FlxObjRBRV_set ( const bool dolog, FlxString* set_name, std::vector<FlxString*> set_parents, std::vector<RBRV_entry_read_base*> set_entries, const bool allow_x2y );
    ~FlxObjRBRV_set();
};

/**
* @brief object read class: for FlxObjRBRV_set
*/
class FlxObjReadRBRV_set : public FlxObjReadBase {
  public:
    FlxObjReadRBRV_set();
    FlxObjBase* read ();
};

//------------------------------------------------------------


/**
* @brief object class: defines a RBRV_noise
*/
class FlxObjRBRV_noise : public FlxObjBase {
  protected:
    FlxString* set_name;
    std::vector<FlxString*> set_parents;
    FlxFunction* Nfun;
    RBRV_entry_read_base* transf;
    
    void task();
  public:
    FlxObjRBRV_noise ( const bool dolog, FlxString* set_name, std::vector<FlxString*> set_parents, FlxFunction* Nfun, RBRV_entry_read_base* transf );
    ~FlxObjRBRV_noise();
};

/**
* @brief object read class: for FlxObjRBRV_noise
*/
class FlxObjReadRBRV_noise : public FlxObjReadBase {
  public:
    FlxObjBase* read ();
};


/**
* @brief object class: defines a RBRV_proc
*/
class FlxObjRBRV_proc : public FlxObjRBRV_noise {
  private:
    FlxFunction* rho;
    FlxFunction* dxf;
    const tuint M;
    const tuint evtype;
    const bool only_once;
    const bool rhoGauss;
    
    void task();
  public:
    FlxObjRBRV_proc ( const bool dolog, FlxString* set_name, std::vector<FlxString*> set_parents, FlxFunction* Nfun, RBRV_entry_read_base* transf, FlxFunction* rho, FlxFunction* dxf, const tuint M, const tuint evtype, const bool only_once, const bool rhoGauss );
    ~FlxObjRBRV_proc();
};

/**
* @brief object read class: for FlxObjRBRV_proc
*/
class FlxObjReadRBRV_proc : public FlxObjReadBase {
  public:
    FlxObjReadRBRV_proc();
    FlxObjBase* read ();
};


/**
* @brief object class: defines a RBRV_mvn (as a convolution integral)
*/
class FlxObjRBRV_mvn_conv : public FlxObjBase {
  private:
    FlxString* set_name;
    FlxString* set1_name;
    FlxString* set2_name;
    const tuint M;
    const tuint evtype;
    
    void task();
  public:
    FlxObjRBRV_mvn_conv ( const bool dolog, FlxString* set_name, FlxString* set1_name, FlxString* set2_name, const tuint M, const tuint evtype );
    ~FlxObjRBRV_mvn_conv();
};

/**
* @brief object class: defines a RBRV_mvn (as a posterior distribution)
*/
class FlxObjRBRV_mvn_post : public FlxObjBase {
  private:
    FlxString* set_name;
    FlxString* set1_name;
    FlxString* set2_name;
    std::string ov_name;
    const bool only_obsv;
    const tuint M;
    const tuint evtype;
    flxVec* helpV;
    FlxMtxSym* helpM;
    
    void task();
  public:
    FlxObjRBRV_mvn_post ( const bool dolog, FlxString* set_name, FlxString* set1_name, FlxString* set2_name, std::string ov_name, const bool only_obsv, const tuint M, const tuint evtype );
    ~FlxObjRBRV_mvn_post();
};

/**
* @brief object read class: for FlxObjRBRV_mvn
*/
class FlxObjReadRBRV_mvn : public FlxObjReadBase {
  public:
    FlxObjReadRBRV_mvn();
    FlxObjBase* read ();
};


/**
* @brief object class: defines a RBRV_mvn (as a posterior distribution)
*/
class FlxObjRBRV_mvn_cond_obsv : public FlxObjBase {
  private:
    FlxString* set_name;
    FlxString* tv_name;
    
    void task();
  public:
    FlxObjRBRV_mvn_cond_obsv ( const bool dolog, FlxString* set_name, FlxString* tv_name );
    ~FlxObjRBRV_mvn_cond_obsv();
};

/**
* @brief object read class: for FlxObjRBRV_mvn_cond_obsv
*/
class FlxObjReadRBRV_mvn_cond_obsv : public FlxObjReadBase {
  public:
    FlxObjBase* read ();
};


/**
* @brief object class: defines a RBRV_psd
*/
class FlxObjRBRV_psd : public FlxObjBase {
  protected:
    FlxString* set_name;
    std::vector<FlxString*> set_parents;
    FlxFunction* Nfun;
    FlxFunction* psd_fun;
    FlxFunction* lbfun;
    FlxFunction* ubfun;
    
    void task();
  public:
    FlxObjRBRV_psd ( const bool dolog, FlxString* set_name, std::vector<FlxString*> set_parents, FlxFunction* Nfun, FlxFunction* psd_fun, FlxFunction* lbfun, FlxFunction* ubfun );
    ~FlxObjRBRV_psd();
};

/**
* @brief object read class: for FlxObjRBRV_psd
*/
class FlxObjReadRBRV_psd : public FlxObjReadBase {
  public:
    FlxObjBase* read ();
};


/**
* @brief object class: defines a RBRV_sphere
*/
class FlxObjRBRV_sphere : public FlxObjBase {
  protected:
    FlxString* set_name;
    std::vector<FlxString*> set_parents;
    FlxFunction* Nfun;
    FlxFunction* r;
    
    void task();
  public:
    FlxObjRBRV_sphere ( const bool dolog, FlxString* set_name, std::vector<FlxString*> set_parents, FlxFunction* Nfun, FlxFunction* r );
    ~FlxObjRBRV_sphere();
};

/**
* @brief object read class: for FlxObjRBRV_sphere
*/
class FlxObjReadRBRV_sphere : public FlxObjReadBase {
  public:
    FlxObjBase* read ();
};


//------------------------------------------------------------


/**
* @brief object class: defines a RBRV_sphere
*/
class FlxObjRBRV_vfset : public FlxObjBase {
  protected:
    const tuint type;
    FlxString* set_name;
    std::vector<FlxString*> set_parents;
    FlxFunction* Nfun;
    FlxMtxConstFun* vecfun;
    FlxFunction* Ntrials;
    
    void task();
  public:
    FlxObjRBRV_vfset ( const bool dolog, FlxString* set_name, std::vector<FlxString*> set_parents, FlxFunction* Nfun, FlxMtxConstFun* vecfun, FlxFunction* Ntrials, const tuint type=0 );
    ~FlxObjRBRV_vfset();
};

/**
* @brief object read class: for FlxObjRBRV_sphere
*/
class FlxObjReadRBRV_vfset : public FlxObjReadBase {
  protected:
    const tuint type;
      // 0: vfset
      // 1: dirichlet variable
  public:
    FlxObjReadRBRV_vfset(const tuint type) : type(type) {}
    FlxObjBase* read ();
};

//------------------------------------------------------------

/**
* @brief object class: defines a RBRV_print
*/
class FlxObjRBRV_print : public FlxObjOutputBase {
  private:
    FlxString* rbstr;
    
    void task();
  public:
    FlxObjRBRV_print ( const bool dolog, std::string ostreamV, const bool verbose, FlxString* rbstr ) : FlxObjOutputBase(dolog,ostreamV,verbose), rbstr(rbstr) {}
    ~FlxObjRBRV_print() { if (rbstr) delete rbstr; }
};

/**
* @brief object read class: for FlxObjRBRV_print
*/
class FlxObjReadRBRV_print : public FlxObjReadOutputBase {
  public:
    FlxObjBase* read ();
};

/**
* @brief object class: defines a RBRV_info
*/
class FlxObjRBRV_info : public FlxObjOutputBase {
  private:
    RBRV_entry_RV_base* rep;
    
    void task();
  public:
    FlxObjRBRV_info ( const bool dolog, std::string ostreamV, RBRV_entry_RV_base* rep ) : FlxObjOutputBase(dolog,ostreamV), rep(rep) {}
    ~FlxObjRBRV_info() { if (rep) delete rep; }
};

/**
* @brief object read class: for FlxObjRBRV_info
*/
class FlxObjReadRBRV_info : public FlxObjReadOutputBase {
  public:
    FlxObjBase* read ();
};



/**
* @brief object class:
*/
class FlxObjRBRV_vec_get : public FlxObjBase {
  public:
    enum RBRVvecGetType { x, y, mean, sd };
  private:
    FlxMtxConstFun* MtxConstStr;                // -> eval -> = vecName
    FlxString* rbstr;                        // -> eval -> = rbstr || constr
    RBRV_constructor* constr;                // only_this=false: constructor
    RBRV_set_base* rbrvSet;                // only_this=true:  set
    const bool only_this;
    tuint NOX;                                // number of random variables (original sapce)
    tuint NRV;                                // number of random variables (standard normal space)
    std::string vecName;                // name of the vector to assign
    const RBRVvecGetType gType;
    
    void task();
  public:
    FlxObjRBRV_vec_get ( const bool dolog, FlxMtxConstFun* MtxConstStr, FlxString* rbstr, const bool only_this, const RBRVvecGetType gType ) : FlxObjBase(dolog), MtxConstStr(MtxConstStr), rbstr(rbstr), constr(NULL), rbrvSet(NULL), only_this(only_this), NOX(0), NRV(0), vecName(""), gType(gType) {}
    ~FlxObjRBRV_vec_get();
};

/**
* @brief object read class: for FlxObjRBRV_vec_get
*/
class FlxObjReadRBRV_vec_get : public FlxObjReadBase {
  public:
    FlxObjReadRBRV_vec_get();
    FlxObjBase* read ();
};


/**
* @brief object class:
*/
class FlxObjRBRV_vec_set : public FlxObjBase {
  public:
    enum RBRVvecSetType { x, y };
  private:
    FlxMtxConstFun* MtxConstStr;                // -> eval -> = vecName
    FlxString* rbstr;                        // -> eval -> = rbstr || constr
    RBRV_constructor* constr;                // only_this=false: constructor
    RBRV_set_base* rbrvSet;                // only_this=true:  set
    const bool only_this;
    tuint NOX;                                // number of random variables (original sapce)
    tuint NRV;                                // number of random variables (standard normal space)
    std::string vecName;                // name of the vector to assign
    const RBRVvecSetType sType;
    
    void task();
  public:
    FlxObjRBRV_vec_set ( const bool dolog, FlxMtxConstFun* MtxConstStr, FlxString* rbstr, const bool only_this, const RBRVvecSetType sType ) : FlxObjBase(dolog), MtxConstStr(MtxConstStr), rbstr(rbstr), constr(NULL), rbrvSet(NULL), only_this(only_this), NOX(0), NRV(0), vecName(""), sType(sType) {}
    ~FlxObjRBRV_vec_set();
};

/**
* @brief object read class: for FlxObjRBRV_vec_set
*/
class FlxObjReadRBRV_vec_set : public FlxObjReadBase {
  public:
    FlxObjReadRBRV_vec_set();
    FlxObjBase* read ();
};



//------------------------------------------------------------


/**
* @brief arithmetic class: returns the realization of a RBRV
*/
class FunRBRV : public FunBase, public FlxDataBase {
  private:
    RBRV_entry* thenumber;
    const std::string rbrv_name;
    const bool is_log;
  public:
    FunRBRV (const std::string& rbrv_name, const bool is_log) : thenumber(NULL), rbrv_name(rbrv_name), is_log(is_log) {};
    const tdouble calc();
    const bool search_circref(FlxFunction* fcr) { return false; }
    const tuint precedence() { return 0; }
    const std::string write();
    const bool dependOn_Const(const tdouble* const theconst);
    const bool optimize(FunBasePtr& optf, const Fun_OptimizeInfo &foi) { return false; }
};

class FunRBRV_fast : public FunBase, public FlxDataBase {
  private:
    FlxString* rbrv_name;
    const bool is_log;
  public:
    FunRBRV_fast (FlxString* rbrv_name, const bool is_log) : rbrv_name(rbrv_name), is_log(is_log) {};
    ~FunRBRV_fast();
    const tdouble calc();
    const bool search_circref(FlxFunction* fcr) { return rbrv_name->search_circref(fcr); }
    const tuint precedence() { return 0; }
    const std::string write();
    const bool dependOn_Const(const tdouble* const theconst);
    const bool optimize(FunBasePtr& optf, const Fun_OptimizeInfo &foi) { return false; }
};


class FunReadFunRBRV : public FunReadFunBase {
  public:
    FunBase* read ( bool errSerious );
};


/**
* @brief arithmetic class: returns log-transform of joint PDF of realization
*/
class FunRBRV_prob : public FunBase, public FlxDataBase {
  private:
    FlxString* MtxConstStr;                // -> eval -> = vecName
    FlxString* rbstr;                        // -> eval -> = rbstr
    RBRV_set_base* rbrvSet;                // the set
    tuint N;                                // number of random variables
    std::string vecName;                // name of the vector to assign
    
  public:
    FunRBRV_prob (FlxString* MtxConstStr, FlxString* rbstr) : MtxConstStr(MtxConstStr), rbstr(rbstr), rbrvSet(NULL), N(0) {};
    ~FunRBRV_prob();
    const tdouble calc();
    const bool search_circref(FlxFunction* fcr);
    const tuint precedence() { return 0; }
    const std::string write();
    const bool dependOn_Const(const tdouble* const theconst);
    const bool optimize(FunBasePtr& optf, const Fun_OptimizeInfo &foi) { return false; }
};

class FunReadFunRBRV_prob : public FunReadFunBase {
  public:
    FunBase* read ( bool errSerious );
};


class FunRBRV_rp : public FunBaseFun_onePara {
  protected:
    virtual const std::string& get_name() = 0;
  public:
    FunRBRV_rp(std::vector<FunBase*> *ParaListV) : FunBaseFun_onePara(ParaListV) {}
    const std::string write();
    const std::string write_v() { return "rbrv_rp"; }
};

class FunRBRV_rp_psd : public FunRBRV_rp {
  protected:
    RBRV_set_psd* rp;
    virtual const std::string& get_name() { return rp->name; }
  public:
    FunRBRV_rp_psd(RBRV_set_psd* rp, std::vector<FunBase*> *ParaListV) : FunRBRV_rp(ParaListV), rp(rp) {}
    virtual const tdouble calc();
};

/**
* @brief arithmetic read class of FunRBRV_rp
*/
class FunReadFunRBRV_rp : public FunReadFunBase, public FlxReaderBase2, public FlxDataBase {
  public:
    FunBase* read ( bool errSerious );
};


/**
* @brief arithmetic class: pdf
*/
class FunPDF : public FunBaseFun {
  protected:
    FunBase* fun;
    RBRV_entry_RV_base *rep;
    const bool free_rep;
  public:
    FunPDF (FunBase* fun, RBRV_entry_RV_base *rep, const bool free_rep) : fun(fun), rep(rep), free_rep(free_rep) {};
    virtual ~FunPDF() { if (fun) delete fun; if (free_rep) delete rep; }
    virtual const tdouble calc();
    virtual const std::string write_v() { return "pdf";}
    virtual const bool search_circref(FlxFunction* fcr);
    virtual const std::string write();
    virtual const bool dependOn_Const(const tdouble* const thenumber);
    virtual const bool optimize(FunBasePtr& optf, const Fun_OptimizeInfo &foi) { return false; }
    virtual const bool evalw() { return fun->evalw(); }
};

class FunReadFunPDF : public FunReadFunBase {
  protected:
    const int methID;        // 0:pdf, 1:cdf, 2:cdf_inv; 3:entropy; 4:mean; 5:sd; 6:coeffofvar; 7:y2x; 8:pdf_ln; 9:rnd_sample; 10:x2y; 11:hpd; 12:median, 13:mode, 14:sf
  public:
    FunReadFunPDF() : methID(0) {}
    FunReadFunPDF(const int methID) : methID(methID) {}
    virtual FunBase* read ( bool errSerious );
};


/**
* @brief arithmetic class: pdf_log
*/
class FunPDF_log : public FunPDF {
  public:
    FunPDF_log (FunBase* fun, RBRV_entry_RV_base *rep, const bool free_rep) : FunPDF(fun,rep,free_rep) {};
    virtual const tdouble calc();
    virtual const std::string write_v() { return "pdf_ln";}
};

class FunReadFunPDF_log : public FunReadFunPDF {
  public:
    FunReadFunPDF_log() : FunReadFunPDF(8) {}
};


/**
* @brief arithmetic class: cdf
*/
class FunCDF : public FunPDF {
  public:
    FunCDF (FunBase* fun, RBRV_entry_RV_base *rep, const bool free_rep) : FunPDF(fun,rep,free_rep) {};
    virtual const tdouble calc();
    virtual const std::string write_v() { return "cdf";}
};

class FunReadFunCDF : public FunReadFunPDF {
  public:
    FunReadFunCDF() : FunReadFunPDF(1) {}
};

/**
* @brief arithmetic class: cdf
*/
class FunSF : public FunPDF {
  public:
    FunSF (FunBase* fun, RBRV_entry_RV_base *rep, const bool free_rep) : FunPDF(fun,rep,free_rep) {};
    virtual const tdouble calc();
    virtual const std::string write_v() { return "sf";}
};

class FunReadFunSF : public FunReadFunPDF {
  public:
    FunReadFunSF() : FunReadFunPDF(14) {}
};


/**
* @brief arithmetic class: cdf_inv
*/
class FunCDF_inv : public FunPDF {
  public:
    FunCDF_inv (FunBase* fun, RBRV_entry_RV_base *rep, const bool free_rep) : FunPDF(fun,rep,free_rep) {};
    virtual const tdouble calc();
    virtual const std::string write_v() { return "cdf_inv";}
};

class FunReadFunCDF_inv : public FunReadFunPDF {
  public:
    FunReadFunCDF_inv() : FunReadFunPDF(2) {}
};



/**
* @brief arithmetic class: cdf
*/
class FunHPD : public FunPDF {
  public:
    FunHPD (FunBase* fun, RBRV_entry_RV_base *rep, const bool free_rep) : FunPDF(fun,rep,free_rep) {};
    virtual const tdouble calc();
    virtual const std::string write_v() { return "hpd";}
};

class FunReadFunHPD : public FunReadFunPDF {
  public:
    FunReadFunHPD() : FunReadFunPDF(11) {}
};


/**
* @brief arithmetic class: rnd_y2x
*/
class FunRBRV_y2x : public FunPDF {
  public:
    FunRBRV_y2x (FunBase* fun, RBRV_entry_RV_base *rep, const bool free_rep) : FunPDF(fun,rep,free_rep) {};
    virtual const tdouble calc();
    virtual const std::string write_v() { return "rnd_y2x";}
};

class FunReadFunRBRV_y2x : public FunReadFunPDF {
  public:
    FunReadFunRBRV_y2x() : FunReadFunPDF(7) {}
};

/**
* @brief arithmetic class: rnd_x2y
*/
class FunRBRV_x2y : public FunPDF {
  public:
    FunRBRV_x2y (FunBase* fun, RBRV_entry_RV_base *rep, const bool free_rep) : FunPDF(fun,rep,free_rep) {};
    virtual const tdouble calc();
    virtual const std::string write_v() { return "rnd_x2y";}
};

class FunReadFunRBRV_x2y : public FunReadFunPDF {
  public:
    FunReadFunRBRV_x2y() : FunReadFunPDF(10) {}
};


/**
* @brief arithmetic class: entropy
*/
class FunEntropy : public FunPDF {
  public:
    FunEntropy (FunBase* fun, RBRV_entry_RV_base *rep, const bool free_rep) : FunPDF(fun,rep,free_rep) {};
    virtual const tdouble calc();
    virtual const std::string write_v() { return "entropy";}
};

class FunReadFunEntropy : public FunReadFunPDF {
  public:
    FunReadFunEntropy() : FunReadFunPDF(3) {}
};

/**
* @brief arithmetic class: entropy
*/
class FunRndSample : public FunPDF {
  public:
    FunRndSample (FunBase* fun, RBRV_entry_RV_base *rep, const bool free_rep) : FunPDF(fun,rep,free_rep) {};
    virtual const tdouble calc();
    virtual const std::string write_v() { return "rnd_sample";}
};

class FunReadFunRndSample : public FunReadFunPDF {
  public:
    FunReadFunRndSample() : FunReadFunPDF(9) {}
};

/**
* @brief arithmetic class: mean
*/
class FunRBRV_mean : public FunPDF {
  public:
    FunRBRV_mean (FunBase* fun, RBRV_entry_RV_base *rep, const bool free_rep) : FunPDF(fun,rep,free_rep) {};
    virtual const tdouble calc();
    virtual const std::string write_v() { return "mean";}
};

class FunReadFunRBRV_mean : public FunReadFunPDF {
  public:
    FunReadFunRBRV_mean() : FunReadFunPDF(4) {}
};

/**
* @brief arithmetic class: stddev
*/
class FunRBRV_sd : public FunPDF {
  public:
    FunRBRV_sd (FunBase* fun, RBRV_entry_RV_base *rep, const bool free_rep) : FunPDF(fun,rep,free_rep) {};
    virtual const tdouble calc();
    virtual const std::string write_v() { return "stddev";}
};

class FunReadFunRBRV_sd : public FunReadFunPDF {
  public:
    FunReadFunRBRV_sd() : FunReadFunPDF(5) {}
};

/**
* @brief arithmetic class: C.o.V.
*/
class FunRBRV_coeffofvar : public FunPDF {
  public:
    FunRBRV_coeffofvar (FunBase* fun, RBRV_entry_RV_base *rep, const bool free_rep) : FunPDF(fun,rep,free_rep) {};
    virtual const tdouble calc();
    virtual const std::string write_v() { return "coeffofvar";}
};

class FunReadFunRBRV_coeffofvar : public FunReadFunPDF {
  public:
    FunReadFunRBRV_coeffofvar() : FunReadFunPDF(6) {}
};

/**
* @brief arithmetic class: median
*/
class FunRBRV_median : public FunPDF {
  public:
    FunRBRV_median (FunBase* fun, RBRV_entry_RV_base *rep, const bool free_rep) : FunPDF(fun,rep,free_rep) {};
    virtual const tdouble calc();
    virtual const std::string write_v() { return "median";}
};

class FunReadFunRBRV_median : public FunReadFunPDF {
  public:
    FunReadFunRBRV_median() : FunReadFunPDF(12) {}
};

/**
* @brief arithmetic class: mode
*/
class FunRBRV_mode : public FunPDF {
  public:
    FunRBRV_mode (FunBase* fun, RBRV_entry_RV_base *rep, const bool free_rep) : FunPDF(fun,rep,free_rep) {};
    virtual const tdouble calc();
    virtual const std::string write_v() { return "mode";}
};

class FunReadFunRBRV_mode : public FunReadFunPDF {
  public:
    FunReadFunRBRV_mode() : FunReadFunPDF(13) {}
};

/**
* @brief arithmetic class: expectation of an expression
*/
class FunExpectation_1d : public FunBase, public FlxDataBase {
  private:
    FunBase* fun;
    RBRV_entry_RV_base* thenumber;
    FlxString* rbrv_name;
    FunBase* ni;
    FunBase* ns;
    FunBase* rd;
    FunBase* lb;
    FunBase* ub;
  public:
    FunExpectation_1d (FunBase* fun, FlxString* rbrv_name, FunBase* ni, FunBase* ns, FunBase* rd, FunBase* lb, FunBase* ub) : fun(fun), thenumber(NULL), rbrv_name(rbrv_name), ni(ni), ns(ns), rd(rd), lb(lb), ub(ub) {};
    ~FunExpectation_1d();
    const tdouble calc();
    const std::string write() { return "expectation_1d(...)";}
    const bool search_circref(FlxFunction* fcr);
    const tuint precedence() { return 0; }
    const bool dependOn_Const(const tdouble* const theconst);
    const bool optimize(FunBasePtr& optf, const Fun_OptimizeInfo &foi);
};

class FunReadFunExpectation_1d : public FunReadFunBase, public FlxDataBase {
  public:
    FunBase* read ( bool errSerious );
};

/**
* @brief arithmetic class: expectation of an expression
*/
class FunExpectation_mci : public FunBase, public FlxDataBase {
  private:
    FunBase* fun;
    RBRV_constructor* rndBox;
    FlxString* rbrv_sets;
    FunBase* ni;
    FunBase* nsi;
    FunBase* nsr;
    FunBase* rd;
    FunBase* bound;
  public:
    FunExpectation_mci (FunBase* fun, FlxString* rbrv_sets, FunBase* ni, FunBase* nsi, FunBase* nsr, FunBase* rd, FunBase* bound) : fun(fun), rndBox(NULL), rbrv_sets(rbrv_sets), ni(ni), nsi(nsi), nsr(nsr), rd(rd), bound(bound) {};
    ~FunExpectation_mci ();
    const tdouble calc();
    const std::string write() { return "expectation_mci(...)";}
    const bool search_circref(FlxFunction* fcr);
    const tuint precedence() { return 0; }
    const bool dependOn_Const(const tdouble* const theconst);
    const bool optimize(FunBasePtr& optf, const Fun_OptimizeInfo &foi);
};

class FunReadFunExpectation_mci : public FunReadFunBase, public FlxDataBase {
  public:
    FunBase* read ( bool errSerious );
};


// -------------------------------------------------------------------------------
// more functions
// -------------------------------------------------------------------------------


class FLXLIB_EXPORT FunRBRV_calc_R_for_rhoPrime : public FunBase, public FlxDataBase {
  private:
    const tdouble& eps_x;
    const tdouble& eps_y;
    const tdouble& ubound;
    const tdouble& intn;
    RBRV_entry_RV_base* rv1;
    RBRV_entry_RV_base* rv2;
    FlxFunction* rhoF;
    tdouble rho;        // correlation in original space
    tdouble y1,y2,r;        // variables needed for iteration; r: correlation in standard normal space
    tdouble m1, m2, s1, s2;        // mean and standard deviations
    FunBase* fun;
    const bool corr_approx;
    
    bool last_num;
  public:
    /**
    * @brief calculates the correlation of the underlying standard normal random variables
    *   rv1 and rv2 has to be managed by the parent!!!
    */
    FunRBRV_calc_R_for_rhoPrime(RBRV_entry_RV_base* rv1, RBRV_entry_RV_base* rv2, FlxFunction* rhoF, const bool corr_approx);
    virtual ~FunRBRV_calc_R_for_rhoPrime () { delete rhoF; delete fun; }
    
    virtual const tdouble calc();
    const tdouble calc_(const bool throwErr);
    const tuint precedence() { return 0; }
    virtual const std::string write();
    virtual const bool dependOn_Const(const tdouble* const thenumber);
    virtual const bool optimize(FunBasePtr& optf, const Fun_OptimizeInfo &foi) { return false; }
    virtual const bool search_circref(FlxFunction* fcr);
    
    const bool last_was_numerical() const { return last_num; }
};





