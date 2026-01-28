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

#ifndef fesslix_flxgp_kernel_H
#define fesslix_flxgp_kernel_H

#include "flxfunction.h"

#if FLX_USE_NLOPT == 0
  #include <gsl/gsl_multimin.h>
#endif

// ------------------------------------------------------------------------------------------------

class FLXLIB_EXPORT flxGP_mean_base {
  protected:
    const tuint Ndim;
    flxVec pVec;
  public:
    flxGP_mean_base(const tuint Npara, const tuint Ndim) : Ndim(Ndim), pVec(Npara)  {}
    virtual ~flxGP_mean_base() {}

    const tuint get_Npara() const { return pVec.get_N(); }
    virtual void initialize_pVec(const tdouble mean_output) = 0;
    void initialize_from_template(const flxVec& vt);
    virtual void reset_pVec();
    virtual const tdouble eval_mean_f(const tdouble* pV) = 0;
    virtual void assemble_F(FlxMtx& F, const tdouble* dm_ptr) = 0;
    virtual void assemble_f_vec(flxVec& f, const tdouble* dm_ptr) = 0;
    flxVec& get_paraVec() { return pVec; }
    virtual py::dict info() = 0;
};

class FLXLIB_EXPORT flxGP_mean_0 : public flxGP_mean_base {
  public:
    flxGP_mean_0(const tuint Ndim);

    virtual void initialize_pVec(const tdouble mean_output) {}
    virtual const tdouble eval_mean_f(const tdouble* pV) { return ZERO; }
    virtual void assemble_F(FlxMtx& F, const tdouble* dm_ptr);
    virtual void assemble_f_vec(flxVec& f, const tdouble* dm_ptr);
    virtual py::dict info();
};

class FLXLIB_EXPORT flxGP_mean_const : public flxGP_mean_base {
  private:
    const tdouble val;
  public:
    flxGP_mean_const(const tdouble val, const tuint Ndim);

    virtual void initialize_pVec(const tdouble mean_output) {}
    virtual const tdouble eval_mean_f(const tdouble* pV) { return val; }
    virtual void assemble_F(FlxMtx& F, const tdouble* dm_ptr);
    virtual void assemble_f_vec(flxVec& f, const tdouble* dm_ptr);
    virtual py::dict info();
};

class FLXLIB_EXPORT flxGP_mean_ordinary : public flxGP_mean_base {
  protected:
    FlxFunction* mean_f;    // just a pointer, memory is not managed by class
  public:
    flxGP_mean_ordinary(FunReadFunUser* mean_f_, const tuint Ndim);
    flxGP_mean_ordinary(FlxFunction* mean_f, const tuint Ndim);

    virtual void initialize_pVec(const tdouble mean_output);
    virtual const tdouble eval_mean_f(const tdouble* pV);
    virtual void assemble_F(FlxMtx& F, const tdouble* dm_ptr);
    virtual void assemble_f_vec(flxVec& f, const tdouble* dm_ptr);
    virtual py::dict info();
};

class FLXLIB_EXPORT flxGP_mean_universal : public flxGP_mean_base {
  protected:
    const tuint type;
    tdouble normalizef;
  public:
    flxGP_mean_universal(const tuint type_, const tuint Ndim);

    virtual void initialize_pVec(const tdouble mean_output);
    virtual const tdouble eval_mean_f(const tdouble* pV);
    virtual void assemble_F(FlxMtx& F, const tdouble* dm_ptr);
    virtual void assemble_f_vec(flxVec& f, const tdouble* dm_ptr);
    virtual py::dict info();
};


// ------------------------------------------------------------------------------------------------

class FLXLIB_EXPORT flxGP_kernel_base {
  protected:
    const tuint Ndim;
  public:
    flxGP_kernel_base(const tuint Ndim);
    virtual ~flxGP_kernel_base() {}

    virtual const tuint get_Npara() const = 0;
    virtual void initialize_pVec(const tdouble sd_output, const flxVec& sd_vec) = 0;
    virtual void initialize_from_template(const flxVec& vt) = 0;
    virtual void reset_pVec() = 0;
    virtual const tdouble eval_kernel_f(const tdouble* pV) = 0;
    virtual const tdouble eval_kernel_corrl_f(const tdouble* pV);
    virtual const tdouble eval_kernel_sd();
    virtual void set_sd(const tdouble sd);
    virtual flxVec* get_paraVec() { return NULL; }
    virtual py::dict info() = 0;
};

class FLXLIB_EXPORT flxGP_kernel_user : public flxGP_kernel_base {
  protected:
    FlxFunction* kernel_f;    // just a pointer, memory is not managed by class
  public:
    flxGP_kernel_user(FunReadFunUser* kernel_f_, const tuint Ndim);
    flxGP_kernel_user(FlxFunction* kernel_f, const tuint Ndim);

    virtual const tuint get_Npara() const { return 0; }
    virtual void initialize_pVec(const tdouble sd_output, const flxVec& sd_vec) {}
    virtual void initialize_from_template(const flxVec& vt) {}
    virtual void reset_pVec() {}
    virtual const tdouble eval_kernel_f(const tdouble* pV);
    virtual py::dict info();
};

class FLXLIB_EXPORT flxGP_kernel_auto : public flxGP_kernel_base {
  protected:
    flxVec pVec;    // parameter vector
    flxVec nVec;    // scaling vector
    iVec tVec;     // type vector
  public:
    flxGP_kernel_auto(const iVec& tVec);
    flxGP_kernel_auto(std::vector<std::string> &kernel_type_lst);

    virtual const tuint get_Npara() const { return pVec.get_N(); }
    virtual void initialize_pVec(const tdouble sd_output, const flxVec& sd_vec);
    virtual void initialize_from_template(const flxVec& vt);
    virtual void reset_pVec();
    virtual const tdouble eval_kernel_f(const tdouble* pV);
    virtual const tdouble eval_kernel_corrl_f(const tdouble* pV);
    virtual const tdouble eval_kernel_sd();
    virtual void set_sd(const tdouble sd);
    virtual flxVec* get_paraVec() { return &pVec; }
    virtual py::dict info();
};

// ------------------------------------------------------------------------------------------------

class FLXLIB_EXPORT flxGP_data_base {
  public:
    flxGP_data_base() {}
    virtual ~flxGP_data_base() {}

    virtual flxGP_data_base* copy() const = 0;
    virtual const tdouble* get_data_ptr_mtx(tuint &N_obsv, const tuint N_dim) = 0;
    virtual const tdouble* get_data_ptr_vec(const tuint N_obsv) = 0;
};

// ------------------------------------------------------------------------------------------------
#if FLX_USE_NLOPT
  double gp_likeli_f (unsigned n, const double *v, double *grad, void *params);
#else
  double gp_likeli_f (const gsl_vector *v, void *params);
#endif

class FLXLIB_EXPORT flxGPProj_base {
  protected:
    const std::string name;   // name of this project
    const tuint Ndim;
  public:
    flxGPProj_base(const std::string& name, const tuint Ndim);
    virtual ~flxGPProj_base() {}

    virtual const tdouble register_observation(const flxGP_data_base& d_info_, const bool initialize_pVec, const bool optimize_noise_val) = 0;
    virtual void register_noise(const tdouble noise_sd) = 0;
    virtual void unassemble() = 0;
    virtual const tdouble optimize(const tuint iterMax, const bool opt_noise) = 0;

    virtual const tdouble get_log_likeli_obsv() = 0;
    /**
    * @brief evaluate posterior mean and variance at x_vec
    *
    * @note can be called in parallel !!!
    */
    virtual void predict_mean_var(const flxVec& x_vec, const bool predict_noise, tdouble& res_mean, tdouble& res_var) = 0;
    virtual const tdouble eval_trend(const flxVec& x_vec, const bool predict_noise) = 0;
    virtual const tdouble eval_kernel(const flxVec& x_vec, const bool predict_noise) = 0;
    virtual py::dict info() = 0;

    const std::string& get_name() const { return name; }
    const tuint get_Ndim() const { return Ndim; }

};

class FLXLIB_EXPORT flxGPProj : public flxGPProj_base {
  private:
    flxGP_mean_base* gp_mean;    // just a pointer, memory is not managed by class
    flxGP_kernel_base* gp_kernel;  // just a pointer, memory is not managed by class
    const bool use_LSE;

    bool obsv_changed;      // true, if observations changed since last assembly
    tuint N_obsv;           // number of data-points in observation
    const tdouble* dm_ptr;  // pointer to observed model input data; matrix of size N_obsv x Ndim
    const tdouble* do_ptr;  // pointer to observed model output data; vector of size N_obsv
    
    /**
    * @brief Covariance/Correlation matrix
    *
    * for use_LSE=true: Q
    * otherwise:        sigma_Z^2 R + sigma_eps^2 I
    */
    FlxMtxSym* cov_mtx_obsv;    // covariance matrix of observed data-points
    FlxMtxLTri* lt_mtx_obsv;    // cholesky decomposition of cov_mtx_obsv
    // matrices required for LSE
      FlxMtxLTri* ltinv_mtx_obsv;
      FlxMtxSym* covinv_mtx_obsv;    // inverse covariance matrix of observed data-points
      FlxMtxSym* FRinvF_mtx;
      FlxMtxLTri* FRinvF_cdc_mtx;
    tdouble noise_val;
    tdouble lt_mtx_obsv_ldet;   // determinant of matrix lt_mtx_obsv (log-transform)
    tdouble lpr_obsv;           // log-likelihood of the observation
    /**
    * @brief trend-corrected observation vector » y - Ba
    */
    flxVec* do_0mean;
    /**
    * @brief alpha = (sigma_Z^2 R + sigma_eps^2 I)^-1 * (y-Ba) = (sigma_Y^2*Q)^-1 * (y-Ba)
    */
    flxVec* alpha;               // = L^T^-1 * L^-1 * do_0mean
    /**
    * @brief matrix B (only assembled if use_LSE==true)
    */
    FlxMtx* Fmtx;               // shape functions of trend
    
    // parameters for optimizing the noise value
      std::stringstream noise_opt_stream;
      tdouble noise_best_guess_val;
      tdouble noise_best_guess_sdZ;
      tdouble noise_best_guess_logl;
    // parameters for optimizing the process parameters
      std::stringstream para_opt_stream;
      bool keep_LSE_results;
    flxGP_data_base* d_info;    // contains information about the observed model data

    /**
    * @brief actually assemble the Gaussian process ond data and (free) process parameters
    *
    * If use_LSE==true, this method modifies the standard deviation of noise (noise_val) and the standrd deviation associated with the kernel.
    *
    * @returns log-likelihood of observation (conditional on (free) process parameters
    */
    const tdouble assemble_observations_help();
    /**
    * @brief ensures that observations are assembled
    */
    void assemble_observations(const bool initialize_pVec, const bool optimize_noise_val);
    /**
    * @brief evaluates the covariance between x_vec and the observed points
    *
    * the returned prior_var is WITHOUT noise
    */
    void eval_covar_point(flxVec& K_star, const flxVec& x_vec, tdouble& prior_mean, tdouble& prior_var);
    const tdouble optimize_help(const tdouble step_size, const tuint iterMax, const bool opt_noise, const bool output_ini, std::ostream& ostrm);
  public:
    flxGPProj(const std::string& name, const tuint Ndim, flxGP_mean_base* gp_mean, flxGP_kernel_base* gp_kernel, const bool use_LSE);
    virtual ~flxGPProj();
    
    virtual const tdouble register_observation(const flxGP_data_base& d_info_, const bool initialize_pVec, const bool optimize_noise_val);
    virtual void register_noise(const tdouble noise_sd);
    /**
    * @brief unassemble the Gaussian process
    *
    * call this function whenever either the process parameters or the data change
    */
    virtual void unassemble() { obsv_changed = true; }
    virtual const tdouble optimize(const tuint iterMax, const bool opt_noise);
    
    virtual const tdouble get_log_likeli_obsv();
    virtual void predict_mean_var(const flxVec& x_vec, const bool predict_noise, tdouble& res_mean, tdouble& res_var);
    virtual const tdouble eval_trend(const flxVec& x_vec, const bool predict_noise);
    virtual const tdouble eval_kernel(const flxVec& x_vec, const bool predict_noise);
    virtual py::dict info();

    void initialize_from_template_trend(const flxGPProj& ref);
    void initialize_from_template_kernel(const flxGPProj& ref);
    const tuint get_AIC_penalty() const;
    
    #if FLX_USE_NLOPT
      friend double gp_likeli_f (unsigned n, const double *v, double *grad, void *params);
    #else
      friend double gp_likeli_f (const gsl_vector *v, void *params);
    #endif
    friend double gp_likeli_f_nv (const tdouble nv, void *params);
};


class FLXLIB_EXPORT flxGP_avgModel : public flxGPProj_base {
  private:
    const tuint max_poM;              // maximum polynomial order of mean function
    const tuint max_NKernel;          // maximum number of kernels to investigate
    std::vector<flxGPProj*> model_lst;
    flxVec modProbVec;
    const tuint iterMax;
    // parameters for optimizing the process parameters
      std::stringstream para_opt_stream;

    const bool increase_kernel_switch(iVec& tVec) const;
  public:
    flxGP_avgModel(const std::string& name, const tuint Ndim, const tuint iterMax, const bool use_LSE);
    virtual ~flxGP_avgModel();

    virtual const tdouble register_observation(const flxGP_data_base& d_info_, const bool initialize_pVec, const bool optimize_noise_val);
    virtual void register_noise(const tdouble noise_sd);
    virtual void unassemble();
    virtual const tdouble optimize(const tuint iterMax, const bool opt_noise);

    virtual const tdouble get_log_likeli_obsv();
    virtual void predict_mean_var(const flxVec& x_vec, const bool predict_noise, tdouble& res_mean, tdouble& res_var);
    virtual const tdouble eval_trend(const flxVec& x_vec, const bool predict_noise);
    virtual const tdouble eval_kernel(const flxVec& x_vec, const bool predict_noise);
    virtual py::dict info();
};



// ------------------------------------------------------------------------------------------------

class FLXLIB_EXPORT flxGPBox {
  private:
    #if FLX_DEBUG
      /**
      * @brief number of instances of this class - only one class may exist
      */
      static int Cinst;
    #endif
    std::map<std::string, flxGPProj_base*> box;
    
  public:
    flxGPBox();
    ~flxGPBox();
    
    void insert (const std::string name, const tuint Ndim, flxGP_mean_base* gp_mean, flxGP_kernel_base* gp_kernel, const bool use_LSE);
    void insert (const std::string name, flxGPProj_base* gp_model);
    flxGPProj_base& get ( const std::string& name );
    
};












#endif // fesslix_flxgp_kernel_H

