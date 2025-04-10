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

#include "flxparse.h"


class RBRV_set;

class PYBIND11_EXPORT RBRV_entry {
  protected:
    #if FLX_DEBUG
      bool valid;                // true, if value has been set
    #endif
    tdouble value;        // the 'realization' of the entry
    RBRV_set* parent;
    
  public:
    const std::string name;
  
    RBRV_entry(const std::string& name);
    virtual ~RBRV_entry() {}
  
    virtual const std::string get_type() const = 0;
    virtual void eval_para();
    /**
    * @brief assigns value
    * @param y_vec a vector of random variables of the parent-set
    */
    virtual void transform_y2x(const tdouble* const y_vec) = 0;
    virtual const tdouble transform_x2y(const tdouble& x_val);
    virtual const tdouble calc_pdf_x(const tdouble& x_val, const bool safeCalc=false);
    virtual const tdouble calc_pdf_x_log(const tdouble& x_val, const bool safeCalc=false);
    virtual const tdouble calc_cdf_x(const tdouble& x_val, const bool safeCalc=false);
    virtual const tdouble calc_sf_x(const tdouble& x_val, const bool safeCalc=false);
    virtual const tdouble calc_entropy();
    /**
    * @brief sets the current x-value of the random variable
    * @note the physical bounds are NOT checked
    */
    void set_x(const tdouble& xV);
    #if FLX_DEBUG
    void set_is_valid(const bool is_valid) { valid = is_valid; }
    #endif
    void set_parent(RBRV_set* parentV);
    #if FLX_DEBUG
      const tdouble& get_value() const;
    #else
      const tdouble& get_value() const { return value; }
    #endif
    virtual const tdouble get_value_log() const { return log(get_value()); }
    const tdouble* get_value_addr() const { return &value; }
    virtual const tdouble get_mean_current_config() = 0;
    virtual const tdouble get_sd_current_config() = 0;
    virtual const tdouble get_median_current_config();
    virtual const tdouble get_mode_current_config();
    virtual const bool check_x(const tdouble xV) = 0;
    virtual const bool search_circref(FlxFunction* fcr) = 0;
    
    void print(std::ostream& sout, const std::string prelim, tuint& counter, const bool printID);      
    virtual py::dict info();
};

class PYBIND11_EXPORT RBRV_entry_fun : public RBRV_entry {
  protected:
    FlxFunction* fun;
    
    tdouble fun_val;
  public:
    RBRV_entry_fun(const std::string& name, FlxFunction* fun);
    virtual ~RBRV_entry_fun() { delete fun; }
        
    const std::string get_type() const { return "fun"; }
    virtual void eval_para();
    virtual void transform_y2x(const tdouble* const y_vec);
    
    virtual const tdouble get_mean_current_config() { return fun_val; }
    virtual const tdouble get_sd_current_config() { return ZERO; }
    virtual const tdouble get_median_current_config() { return fun_val; }
    virtual const tdouble get_mode_current_config() { return fun_val; }
    virtual const bool check_x(const tdouble xV);
    virtual const bool search_circref(FlxFunction* fcr) { return fun->search_circref(fcr); }
};

class PYBIND11_EXPORT RBRV_entry_RV_base : public RBRV_entry {
  protected:
    const tuint iID;        // the ID of the RV in the set
    // correlation with a single RBRV in the set (IF allowed)
      RBRV_entry_RV_base* corr_rv;
      FunBase* corr_valF;                // if !NULL: corr_val needs to be re-evaluated every time!
      tdouble corr_val;                        // correlation in underlying standard normal space
      bool throwErrors;

      void init();
  public:
    RBRV_entry_RV_base(const std::string& name, const tuint iID);
    virtual ~RBRV_entry_RV_base();
    
    virtual void transform_y2x(const tdouble* const y_vec);
    virtual const tdouble transform_y2x(const tdouble y_val) = 0;
    const bool allow_x2y() const { return (corr_rv==NULL); }
    /**
    * @brief returns the lower quantile value of the HPD (highest probability density) interval of the distribution
    */
    virtual const tdouble get_HPD(const tdouble p);
    
    // handle correlation (with a a single RBRV)
      void set_corr(RBRV_entry_RV_base* corr_rv_, FlxFunction* corrVal_, bool corrFixed_, const bool throwErrorsV);
      void eval_corr();
    const tuint get_iID() const { return iID; }
};


//------------------------------------------------------------

class RBRV_set_box;

class FLXLIB_EXPORT RBRV_set_base : public FlxBoxBaseR {
  protected:
    const tuint ID;                // unique identifier (in order of definition -> Nataf: 0)
    const bool internal;        // true: can't be used by other sets  (exceptions!!!)
    const tuint sRV;                // number of random variables in the set
    flxVec y_of_set;                // current y-vector of the set
    
    static tuint static_ID;
    
    const tdouble get_pdf_y_eval_log() const;
  public:
    const std::string name;
    
    RBRV_set_base(const bool internal, const tuint sRV, const std::string& name, const bool noID);
    virtual ~RBRV_set_base() {}
    
    virtual const tuint get_NRV() const = 0;
    virtual const tuint get_NOX() const = 0;
    /**
    * @brief number of rv's only of this set - not including rv's in sub-sets
    */
    virtual const tuint get_NRV_only_this() const = 0;
    virtual const tuint get_NOX_only_this() const = 0;
    const tuint get_ID() const { return ID; }
    const bool is_internal() const { return internal; }
    virtual void set_is_valid(const bool is_valid) = 0;
    virtual const flxVec& propose_y();
    virtual void transform_y2x() = 0;
    virtual void transform_x2y();
    virtual void transform_y2w(const tdouble* const y_vec, tdouble* const w_vec);
    virtual const bool allow_x2y() const = 0;
    
    /**
    * @brief sets y without setting x
    */
    virtual void set_y(const tdouble* const y_vec);
    virtual void set_y_only_this(const tdouble* const y_vec) { set_y(y_vec); }
    virtual void get_y(tdouble* const y_vec);
    virtual void get_y_only_this(tdouble* const y_vec) { get_y(y_vec); };
    virtual const flxVec& get_y() { return y_of_set; }
    /**
    * @brief sets x without setting y first
    */
    virtual void set_x(const tdouble* const x_vec) = 0;
    virtual void set_x_only_this(const tdouble* const x_vec) = 0;
    virtual void get_x(tdouble* const x_vec) = 0;
    virtual void get_x_only_this(tdouble* const x_vec) = 0;
    virtual const bool check_xVec(const tdouble* xp) = 0;
    virtual void get_mean(tdouble* const m_vec) = 0;
    virtual void get_mean_only_this(tdouble* const m_vec) = 0;
    virtual void get_sd(tdouble* const s_vec) = 0;
    virtual void get_sd_only_this(tdouble* const s_vec) = 0;
    /**
    * @returns matrix pointer that must not be deallocated
    */
    virtual const FlxMtxSparsLTri* calc_Jinv_1();
    virtual void calc_Jinv_2(tdouble* dxdw);
    
    /**
    * @brief the sets must add all the sets it depends on ... AND itself
    */
    virtual void find_dependent_sets(std::vector<RBRV_set_base*>& setvec ) = 0;
    /**
    * @brief some sets are managed by other sets ... these need to be removed by the managing set
    * @param setvec the vector that contains all the sets that are (still) in the list
    * @param pos_this the index of 'this' set in the vector
    * @param returns the number of sets that where removed
    */
    virtual const tuint group_dependent_sets(std::vector<RBRV_set_base*>& setvec, const tuint pos_this ) { return 0; }
    
    virtual void print(std::ostream& sout, const std::string prelim, tuint& counter, const bool printID) = 0;
    
    virtual const tdouble get_pdf_x_eval_log();
    /**
    * @brief adds the covariance matrix of this set to cm
    *   the dimension of the matrix cm must match
    */
    virtual void add_covMTX(FlxMtxSym& cm);
    
    virtual const std::string get_rv_name(const tuint index);
    
    static const tuint static_get_ID() { return static_ID; }
};
typedef RBRV_set_base** RBRV_set_baseDPtr;


class FLXLIB_EXPORT RBRV_set_parents : public RBRV_set_base {
  protected:
    const tuint Nparents;        // number of parents of this set
    RBRV_set_base** const parents;        // the parents of this set
  public:
    RBRV_set_parents(const bool internal, const tuint sRV, const std::string& name, const tuint Nparents, RBRV_set_base** const parents, const bool noID);
    virtual ~RBRV_set_parents();
    virtual void find_dependent_sets(std::vector<RBRV_set_base*>& setvec );
    void print_parents(std::ostream& sout);
    
};

class FLXLIB_EXPORT RBRV_set : public RBRV_set_parents {
  protected:
    const tuint Nentries;        // number of entries in the set
    RBRV_entry** const entries;        // the entries
    bool x2y_allowed;        // true, if x2y-transformation is allowed
  public:
    RBRV_set(const bool internal, const tuint sRV, const std::string& name, const bool noID, const tuint Nentries, RBRV_entry** const entries,
             const tuint Nparents, RBRV_set_base** const parents, const bool x2y_allowedV);
    virtual ~RBRV_set();
    virtual const tuint get_NRV() const { return sRV; }
    virtual const tuint get_NOX() const { return Nentries; }
    virtual const tuint get_NRV_only_this() const { return sRV; }
    virtual const tuint get_NOX_only_this() const { return Nentries; }
    virtual void set_is_valid(const bool is_valid);
    virtual void transform_y2x();
    virtual void transform_x2y();
    virtual void transform_y2w(const tdouble* const y_vec, tdouble* const w_vec);
    virtual const bool allow_x2y() const { return x2y_allowed; }
    virtual void set_x(const tdouble* const x_vec);
    virtual void set_x_only_this(const tdouble* const x_vec) { set_x(x_vec); }
    virtual void get_x(tdouble* const x_vec);
    virtual void get_x_only_this(tdouble* const x_vec) { get_x(x_vec); }
    virtual const bool check_xVec(const tdouble* xp);
    virtual void get_mean(tdouble* const m_vec);
    virtual void get_mean_only_this(tdouble* const m_vec) { get_mean(m_vec); };
    virtual void get_sd(tdouble* const s_vec);
    virtual void get_sd_only_this(tdouble* const s_vec) { get_sd(s_vec); };
    virtual const FlxMtxSparsLTri* calc_Jinv_1();
    virtual void calc_Jinv_2(tdouble* dxdw);
    void print(std::ostream& sout, const std::string prelim, tuint& counter, const bool printID);
    virtual const tdouble get_pdf_x_eval_log();
    virtual const std::string get_rv_name(const tuint index);
};

class FLXLIB_EXPORT RBRV_set_Nataf : public RBRV_set_base {
  protected:
    const tuint Nentries;        // number of entries in the set
    flxVec w_of_set;                // current x-vector of the set
    RBRV_entry** const entries;        // the entries
    /**
    * @brief Cholesky-Decomposition of rhoPrime (calculates the inverse multiplication implicitly)
    */
    FlxMtxSparsLTri* L;
    
  public:
    RBRV_set_Nataf(const bool internal, const tuint sRV, const std::string& name, const bool noID, const tuint Nentries, RBRV_entry** const entries, FlxMtxSparsLTri* L);
    virtual ~RBRV_set_Nataf();
    
    virtual const tuint get_NRV() const { return sRV; }
    virtual const tuint get_NOX() const { return get_NRV(); }
    virtual const tuint get_NRV_only_this() const { return sRV; }
    virtual const tuint get_NOX_only_this() const { return get_NRV_only_this(); }
    virtual void set_is_valid(const bool is_valid);
    virtual void transform_y2x();
    virtual void transform_x2y();
    virtual void transform_y2w(const tdouble* const y_vec, tdouble* const w_vec);
    virtual const bool allow_x2y() const { return true; }
    virtual void set_x(const tdouble* const x_vec);
    virtual void set_x_only_this(const tdouble* const x_vec) { set_x(x_vec); }
    virtual void get_x(tdouble* const x_vec);
    virtual void get_x_only_this(tdouble* const x_vec) { get_x(x_vec); }
    virtual const bool check_xVec(const tdouble* xp);
    virtual void get_mean(tdouble* const m_vec);
    virtual void get_mean_only_this(tdouble* const m_vec) { get_mean(m_vec); };
    virtual void get_sd(tdouble* const s_vec);
    virtual void get_sd_only_this(tdouble* const s_vec) { get_sd(s_vec); };
    virtual const FlxMtxSparsLTri* calc_Jinv_1();
    virtual void calc_Jinv_2(tdouble* dxdw);
    void print(std::ostream& sout, const std::string prelim, tuint& counter, const bool printID);
    virtual const tdouble get_pdf_x_eval_log();
    virtual const std::string get_rv_name(const tuint index);
    virtual void find_dependent_sets(std::vector<RBRV_set_base*>& setvec );
};


class FLXLIB_EXPORT RBRV_set_noise : public RBRV_set_parents {
  protected:
    flxVec x_of_set;                // current x-vector of the set
    RBRV_entry* const transf;        // the transformation
    const bool is_stdN;
    #if FLX_DEBUG
      bool valid;
    #endif
      
  public:
    RBRV_set_noise(const bool internal, const tuint sRV, const std::string& name, const bool noID, RBRV_entry* const transf, const tuint Nparents, RBRV_set_base** const parents);
    virtual ~RBRV_set_noise() { delete transf; }

    virtual const tuint get_NRV() const { return sRV; }
    virtual const tuint get_NOX() const { return  get_NRV(); }
    virtual const tuint get_NRV_only_this() const { return sRV; }
    virtual const tuint get_NOX_only_this() const { return get_NRV_only_this(); }
    #if FLX_DEBUG
    virtual void set_is_valid(const bool is_valid) { valid = is_valid; }
    #else
    virtual void set_is_valid(const bool is_valid) {}
    #endif
    virtual void transform_y2x();
    virtual void transform_x2y();
    virtual void transform_y2w(const tdouble* const y_vec, tdouble* const w_vec);
    virtual const bool allow_x2y() const;
    virtual void set_x(const tdouble* const x_vec);
    virtual void set_x_only_this(const tdouble* const x_vec) { set_x(x_vec); }
    virtual void get_x(tdouble* const x_vec);
    virtual void get_x_only_this(tdouble* const x_vec) { get_x(x_vec); }
    virtual const bool check_xVec(const tdouble* xp);
    virtual void get_mean(tdouble* const m_vec);
    virtual void get_mean_only_this(tdouble* const m_vec) { get_mean(m_vec); };
    virtual void get_sd(tdouble* const s_vec);
    virtual void get_sd_only_this(tdouble* const s_vec) { get_sd(s_vec); };
    virtual const FlxMtxSparsLTri* calc_Jinv_1();
    virtual void calc_Jinv_2(tdouble* dxdw);
    virtual void print(std::ostream& sout, const std::string prelim, tuint& counter, const bool printID);
    virtual const tdouble get_pdf_x_eval_log();
    virtual void add_covMTX(FlxMtxSym& cm);
};


class FLXLIB_EXPORT RBRV_set_proc : public RBRV_set_parents {
  protected:
    flxVec x_of_set;                // current x-vector of the set
    RBRV_entry* const transf;        // the transformation
    FlxFunction* rho;                // correlation coefficient function (original space)
    const tdouble dx;                // lag between points
    const tuint N;
    const tuint M;
    const int ev_solver;
    const bool rho_Gauss;
    tdouble eole_err;
    tdouble Jdet;                // log-transformed Jacobian
    const bool only_once;
    #if FLX_DEBUG
      bool valid;
    #endif
    
    flxVec* Eigenvalues;
    std::vector<flxVec> Eigenvectors;        // Eigenvectors (DIM*Ndof)
    FlxMtxLTri* Lt;
    flxVec* xhelp;
      
    void assemble_rhoPrime(FlxMtxSym& rhoPrime);
    void assemble_system();
  public:
    RBRV_set_proc(const bool internal, const tuint Nv, const tuint Mv, const std::string& name, const bool noID, RBRV_entry* const transf, FlxFunction* rho, const tdouble dx, const tuint Nparents, RBRV_set_base** const parents, const int ev_solver, const bool only_once, const bool rho_Gauss);
    virtual ~RBRV_set_proc();

    virtual const tuint get_NRV() const { return sRV; }
    virtual const tuint get_NOX() const { return  N; }
    virtual const tuint get_NRV_only_this() const { return sRV; }
    virtual const tuint get_NOX_only_this() const { return N; }
    #if FLX_DEBUG
    virtual void set_is_valid(const bool is_valid) { valid = is_valid; }
    #else
    virtual void set_is_valid(const bool is_valid) {}
    #endif
    virtual void transform_y2x();
    virtual void transform_x2y();
    virtual const bool allow_x2y() const;
    virtual void set_x(const tdouble* const x_vec);
    virtual void set_x_only_this(const tdouble* const x_vec) { set_x(x_vec); }
    virtual void get_x(tdouble* const x_vec);
    virtual void get_x_only_this(tdouble* const x_vec) { get_x(x_vec); }
    virtual const bool check_xVec(const tdouble* xp);
    virtual void get_mean(tdouble* const m_vec);
    virtual void get_mean_only_this(tdouble* const m_vec) { get_mean(m_vec); };
    virtual void get_sd(tdouble* const s_vec);
    virtual void get_sd_only_this(tdouble* const s_vec) { get_sd(s_vec); };
    virtual void print(std::ostream& sout, const std::string prelim, tuint& counter, const bool printID);
    virtual const tdouble get_pdf_x_eval_log();
    virtual void add_covMTX(FlxMtxSym& cm);
};

class FLXLIB_EXPORT RBRV_set_MVN : public RBRV_set_parents {
  protected:
    flxVec x_of_set;                // current x-vector of the set
    const tuint N;
    const tuint M;
    flxVec* mu;
    FlxMtxSym* CovM;
    const int ev_solver;
    
    tdouble eole_err;
    tdouble Jdet;                // log-transformed Jacobian
    #if FLX_DEBUG
      bool valid;
    #endif
    
    flxVec* Eigenvalues;
    std::vector<flxVec> Eigenvectors;        // Eigenvectors (DIM*Ndof)
    FlxMtxLTri* Lt;
    flxVec* xhelp;
      
    void assemble_system();
    void deallocate();
  public:
    RBRV_set_MVN(const bool internal, const tuint Nv, const tuint Mv, const std::string& name, const bool noID, flxVec* mu, FlxMtxSym* CovM, const int ev_solver);
    virtual ~RBRV_set_MVN() { deallocate(); }

    virtual const tuint get_NRV() const { return sRV; }
    virtual const tuint get_NOX() const { return  N; }
    virtual const tuint get_NRV_only_this() const { return sRV; }
    virtual const tuint get_NOX_only_this() const { return N; }
    #if FLX_DEBUG
    virtual void set_is_valid(const bool is_valid) { valid = is_valid; }
    #else
    virtual void set_is_valid(const bool is_valid) {}
    #endif
    virtual void transform_y2x();
    virtual void transform_x2y();
    virtual const bool allow_x2y() const;
    virtual void set_x(const tdouble* const x_vec);
    virtual void set_x_only_this(const tdouble* const x_vec) { set_x(x_vec); }
    virtual void get_x(tdouble* const x_vec);
    virtual void get_x_only_this(tdouble* const x_vec) { get_x(x_vec); }
    virtual const bool check_xVec(const tdouble* xp) { return true; }
    virtual void get_mean(tdouble* const m_vec);
    virtual void get_mean_only_this(tdouble* const m_vec) { get_mean(m_vec); };
    virtual void get_sd(tdouble* const s_vec);
    virtual void get_sd_only_this(tdouble* const s_vec) { get_sd(s_vec); };
    virtual void print(std::ostream& sout, const std::string prelim, tuint& counter, const bool printID);
    virtual const tdouble get_pdf_x_eval_log();
    virtual void add_covMTX(FlxMtxSym& cm);
    
    FlxMtxSym* get_CovM() { return CovM; }
    flxVec* get_mu() { return mu; }
    void update_EVP();
};

class FLXLIB_EXPORT RBRV_set_MVN_cond : public RBRV_set_parents {
  protected:
    flxVec x_of_set;                // current x-vector of the set (of dimension Nrnd)
    flxVec x_obsv;                // y-transform of observations (of dimension Nobsv)
    flxVec y_obsv;                // y-transform of observations (of dimension Nobsv)
    const tuint Nrnd;
    const tuint Nobsv;
    const tuint Ntotal;                // = Nrnd+Nobsv
    flxVec* mu;                        // of dimension Ntotal        
    FlxMtxSym* CovM;                // of dimension Ntotal

    tdouble Jdet;                // log-transformed Jacobian
    #if FLX_DEBUG
      bool valid;
    #endif

    FlxMtxLTri* Lt;
    flxVec xhelp;                // of dimension Ntotal
    flxVec yhelp;                // of dimension Ntotal
      
    void assemble_system();
    void deallocate();
    void comp_yobsv();                // transforms x_obsv to yobsv
  public:
    RBRV_set_MVN_cond(const bool internal, const tuint NrndV, const tuint NobsvV, const std::string& name, const bool noID, flxVec* mu, FlxMtxSym* CovM, const flxVec& x_obsvV);
    virtual ~RBRV_set_MVN_cond() { deallocate(); }

    virtual const tuint get_NRV() const { return Nrnd; }
    virtual const tuint get_NOX() const { return  Nrnd; }
    virtual const tuint get_NRV_only_this() const { return Nrnd; }
    virtual const tuint get_NOX_only_this() const { return Nrnd; }
    #if FLX_DEBUG
    virtual void set_is_valid(const bool is_valid) { valid = is_valid; }
    #else
    virtual void set_is_valid(const bool is_valid) {}
    #endif
    virtual void transform_y2x();
    virtual void transform_x2y();
    virtual const bool allow_x2y() const;
    virtual void set_x(const tdouble* const x_vec);
    virtual void set_x_only_this(const tdouble* const x_vec) { set_x(x_vec); }
    virtual void get_x(tdouble* const x_vec);
    virtual void get_x_only_this(tdouble* const x_vec) { get_x(x_vec); }
    virtual const bool check_xVec(const tdouble* xp) { return true; }
    virtual void get_mean(tdouble* const m_vec);
    virtual void get_mean_only_this(tdouble* const m_vec) { get_mean(m_vec); };
    virtual void get_sd(tdouble* const s_vec);
    virtual void get_sd_only_this(tdouble* const s_vec) { get_sd(s_vec); };
    virtual void print(std::ostream& sout, const std::string prelim, tuint& counter, const bool printID);
    virtual const tdouble get_pdf_x_eval_log();
    
    const tuint get_Nobsv() const { return Nobsv; }
    void set_x_obsv(const flxVec& x_obsvV);
    FlxMtxSym* get_CovM() { return CovM; }
    flxVec* get_mu() { return mu; }
    void update_EVP();
};

class FLXLIB_EXPORT RBRV_set_psd : public RBRV_set_parents {
  protected:
    const tuint N;                // number of discretization intervals
    FlxFunction* psd_fun;        // power spectral density function
    const tdouble lb;                // lower bound
    const tdouble ub;                // upper bound
    tdouble& wp;
    #if FLX_DEBUG
      bool valid;
    #endif
      
  public:
    RBRV_set_psd(const bool internal, const std::string& name, const tuint Nv, FlxFunction* psd_fun, const tdouble lb, const tdouble ub, const tuint Nparents, RBRV_set_base** const parents, tdouble& wp);
    virtual ~RBRV_set_psd() { delete psd_fun; }

    virtual const tuint get_NRV() const { return sRV; }
    virtual const tuint get_NOX() const { return  0; }
    virtual const tuint get_NRV_only_this() const { return sRV; }
    virtual const tuint get_NOX_only_this() const { return 0; }
    #if FLX_DEBUG
    virtual void set_is_valid(const bool is_valid) { valid = is_valid; }
    #else
    virtual void set_is_valid(const bool is_valid) {}
    #endif
    virtual void transform_y2x() {}
    virtual const bool allow_x2y() const { return false; }
    virtual void set_x(const tdouble* const x_vec) {}
    virtual void set_x_only_this(const tdouble* const x_vec) { set_x(x_vec); }
    virtual void get_x(tdouble* const x_vec) {}
    virtual void get_x_only_this(tdouble* const x_vec) { get_x(x_vec); }
    virtual const bool check_xVec(const tdouble* xp) { return true; }
    virtual void get_mean(tdouble* const m_vec);
    virtual void get_mean_only_this(tdouble* const m_vec) { get_mean(m_vec); };
    virtual void get_sd(tdouble* const s_vec);
    virtual void get_sd_only_this(tdouble* const s_vec) { get_sd(s_vec); };
    virtual void print(std::ostream& sout, const std::string prelim, tuint& counter, const bool printID);
    
    const tdouble eval_realization(const tdouble t);
};


class FLXLIB_EXPORT RBRV_set_sphere : public RBRV_set_parents {
  protected:
    flxVec x_of_set;                // current x-vector of the set
    FlxFunction* r;                // radius of hyper-sphere
    #if FLX_DEBUG
      bool valid;
    #endif
      
  public:
    RBRV_set_sphere(const bool internal, const tuint sRV, const std::string& name, const bool noID, const tuint Nparents, RBRV_set_base** const parents, FlxFunction* r);
    virtual ~RBRV_set_sphere() { delete r; }

    virtual const tuint get_NRV() const { return sRV; }
    virtual const tuint get_NOX() const { return  get_NRV(); }
    virtual const tuint get_NRV_only_this() const { return sRV; }
    virtual const tuint get_NOX_only_this() const { return get_NRV_only_this(); }
    #if FLX_DEBUG
    virtual void set_is_valid(const bool is_valid) { valid = is_valid; }
    #else
    virtual void set_is_valid(const bool is_valid) {}
    #endif
    virtual void transform_y2x();
    virtual void transform_x2y();
    virtual const bool allow_x2y() const { return true; }
    virtual void set_x(const tdouble* const x_vec);
    virtual void set_x_only_this(const tdouble* const x_vec) { set_x(x_vec); }
    virtual void get_x(tdouble* const x_vec);
    virtual void get_x_only_this(tdouble* const x_vec) { get_x(x_vec); }
    virtual const bool check_xVec(const tdouble* xp);
    virtual void get_mean(tdouble* const m_vec);
    virtual void get_mean_only_this(tdouble* const m_vec) { get_mean(m_vec); };
    virtual void get_sd(tdouble* const s_vec);
    virtual void get_sd_only_this(tdouble* const s_vec) { get_sd(s_vec); };
    virtual void print(std::ostream& sout, const std::string prelim, tuint& counter, const bool printID);
    virtual const tdouble get_pdf_x_eval_log();
};


//------------------------------------------------------------


class FLXLIB_EXPORT RBRV_constructor {
  private:
    std::vector<RBRV_set_base*> setvec;        // a list with all the relevant sets (the order does matter)
    tuint NRV;                // the total number of random variables
    tuint NOX;                // the total number of entries (includes functions and random variables)
    tuint Nsets;                // number of entries of setvec(c)
    bool allow_x2y;                // is an x2y transformation allowed?
    
    /**
    * @brief transforms y to original space
    */
    void transform_y2x();
    void transform_x2y();
  public:
    
    RBRV_constructor(const std::vector<RBRV_set_base*>& setvec);
    RBRV_constructor(const std::vector<std::string>& set_str_vec, RBRV_set_box &rbrv_box);
    
    const tuint get_NRV() const { return NRV; }
    const tuint get_NOX() const { return NOX; }
    
    /**
    * @brief propose a new y-realization .. without a y2x transformation
    */
    void propose_y(flxVec& yV);
    void propose_y();
    /**
    * @brief generates a new set of observations (y and y2x-transformation)
    */
    void gen_smp();
    /**
    * @brief sets y=y_new and transfroms y2x
    */
    void set_smp(const flxVec& y_new);
    /**
    * @brief sets y=y_new ... without! a y2x transformation!
    */
    void set_smp_y(const flxVec& y_new);
    /**
    * @brief sets x=x_new ... without! setting y first
    */
    void set_smp_x(const flxVec& x_new);
    /**
    * @brief sets x=x_new ... + and x2y transformation
    */
    void set_smp_x_transform(const flxVec& x_new);
    void set_is_valid(const bool is_valid);
    void get_y_Vec(tdouble* const y_vec) const;
    void get_x_Vec(tdouble* const x_vec) const;
    void get_mean_Vec(tdouble* const m_vec) const;
    void get_sd_Vec(tdouble* const s_vec) const;
    /**
    * @brief checks if x is within physical bounds
    */
    const bool check_xVec(const flxVec& xV) const;
    void calc_Jinv(FlxMtxLTri& dxdy);
    void transform_y2w(const tdouble* const y_vec, tdouble* const w_vec);
    
    void print_info(std::ostream& sout, const std::string prelim="  ");
    /**
    * @returns the name of the ith random variable in the set; counting starts with zero
    */
    const std::string get_rv_name(const tuint index);
    
    /**
    * @brief finds all dependent sets of a list of sets 
    * @param setstr list of sets - separator ','
    * @param setvec returns vector with dependent sets -> the order does matter!
    */
    static void find_dependent_sets(const std::vector<std::string>& set_str_vec, std::vector<RBRV_set_base*>& setvec, RBRV_set_box& RBRVbox);
    /**
    * @brief counts the number of random variables in the set
    */
    static const tuint count_NRV(const std::vector<RBRV_set_base*>& setvec);
    static const tuint count_NRV_long(const std::vector<RBRV_set_base*>& setvec);
    /**
    * @brief counts the number of values in original space (includes functions as well)
    */
    static const tuint count_NOX(const std::vector<RBRV_set_base*>& setvec);
    static const tuint count_NOX_long(const std::vector<RBRV_set_base*>& setvec);
};


//------------------------------------------------------------

class FLXLIB_EXPORT RBRV_set_box {
  private:
    std::map<std::string,RBRV_set_base*> set_box;
    std::map<std::string,RBRV_entry*> entry_box;
    std::vector<RBRV_set_base*> set_vec;
  public:
    ~RBRV_set_box();
    
    void register_set(RBRV_set_base* sta);
    RBRV_set_base* get_set(const std::string& name, const bool throwErr) const;
    const tuint get_set_N() const { return set_vec.size(); }
    
    void register_entry(RBRV_entry* eta);
    RBRV_entry* get_entry(const std::string& name, const bool throwErr) const;
    
    /**
    * @brief print out all defined RBRV sets (and their parents)
    */    
    void print_sets(std::ostream& sout ,const std::string prelim);
    
    const std::vector<RBRV_set_base*>& get_entire_set_vec() { return set_vec; }
};





