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

#include "flxrbrv.h"

#include <boost/math/interpolators/pchip.hpp>



class PYBIND11_EXPORT RBRV_entry_RV_normal : public RBRV_entry_RV_base {
  protected:
    int pid;        // 0: mean,std.dev; 1: quantile value; 2: C.o.V, quantile value; 3: std.dev, quantile value
    FlxFunction* p1;
    FlxFunction* p2;
    FlxFunction* p3;
    FlxFunction* p4;
    bool eval_once;
    tdouble mu;
    tdouble sigma;
    
  public:
    RBRV_entry_RV_normal(const std::string& name, const tuint iID, const int pid, FlxFunction* p1v, FlxFunction* p2v, FlxFunction* p3v, FlxFunction* p4v, const bool eval_once);
    RBRV_entry_RV_normal(const std::string& name, const tuint iID, py::dict config);
    virtual ~RBRV_entry_RV_normal();
    
    const std::string get_type() const { return "normal"; }
    virtual void eval_para();
    virtual const tdouble transform_y2x(const tdouble y_val);
    virtual const tdouble transform_x2y(const tdouble& x_val);
    virtual const tdouble calc_pdf_x(const tdouble& x_val, const bool safeCalc=false);
    virtual const tdouble calc_pdf_x_log(const tdouble& x_val, const bool safeCalc=false);
    virtual const tdouble calc_cdf_x(const tdouble& x_val, const bool safeCalc=false);
    virtual const tdouble calc_sf_x(const tdouble& x_val, const bool safeCalc=false);
    virtual const tdouble calc_entropy();
    virtual const tdouble get_mean_current_config() { return mu; }
    virtual const tdouble get_sd_current_config() { return sigma; }
    virtual const tdouble get_median_current_config() { return get_mean_current_config(); }
    virtual const tdouble get_mode_current_config() { return get_mean_current_config(); }
    virtual const bool check_x(const tdouble xV) { return true; }
    virtual const bool search_circref(FlxFunction* fcr);
    virtual py::dict info();
    virtual const tdouble get_HPD(const tdouble p);
    
    // compute distribution based on two given quantile values
    static void get_para_from_quantile(tdouble& meanV, tdouble& sdV, const tdouble x1, const tdouble p1, const tdouble x2, const tdouble p2);
    // compute distribution based on quantile value and C.o.V.
    static void get_para_from_quantile2(tdouble& meanV, tdouble& sdV, const tdouble x1, const tdouble p1, const tdouble delta);
    // compute distribution based on quantile value and standard deviation
    static void get_para_from_quantile3(tdouble& meanV, const tdouble x1, const tdouble p1, const tdouble sdV);
};


class PYBIND11_EXPORT RBRV_entry_RV_stdN : public RBRV_entry_RV_base {
  public:
    RBRV_entry_RV_stdN(const std::string& name, const tuint iID) : RBRV_entry_RV_base(name,iID) { this->init(); }
    RBRV_entry_RV_stdN(const std::string& name, const tuint iID, py::dict config) : RBRV_entry_RV_base(name,iID) { this->init(); }
    virtual ~RBRV_entry_RV_stdN() {}
    
    const std::string get_type() const { return "stdn"; }
    virtual const tdouble transform_y2x(const tdouble y_val);
    virtual const tdouble transform_x2y(const tdouble& x_val);
    virtual const tdouble calc_pdf_x(const tdouble& x_val, const bool safeCalc=false);
    virtual const tdouble calc_pdf_x_log(const tdouble& x_val, const bool safeCalc=false);
    virtual const tdouble calc_cdf_x(const tdouble& x_val, const bool safeCalc=false);
    virtual const tdouble calc_sf_x(const tdouble& x_val, const bool safeCalc=false);
    virtual const tdouble calc_entropy();
    virtual const tdouble get_mean_current_config() { return ZERO; }
    virtual const tdouble get_sd_current_config() { return ONE; }
    virtual const tdouble get_median_current_config() { return get_mean_current_config(); }
    virtual const tdouble get_mode_current_config() { return get_mean_current_config(); }
    virtual const bool check_x(const tdouble xV) { return true; }
    virtual const bool search_circref(FlxFunction* fcr) { return false; }
    virtual py::dict info();
    virtual const tdouble get_HPD(const tdouble p);
};

class PYBIND11_EXPORT RBRV_entry_RV_lognormal : public RBRV_entry_RV_base {
  protected:
    int pid;                // ... see get_paras()
    FlxFunction* p1;
    FlxFunction* p2;
    FlxFunction* p3;
    FlxFunction* p4;
    FlxFunction* epsilon;
    bool eval_once;
    tdouble lambda;
    tdouble zeta;
    tdouble eps;
    
  public:
    RBRV_entry_RV_lognormal(const std::string& name, const tuint iID, const int pid, FlxFunction* p1, FlxFunction* p2, FlxFunction* p3, FlxFunction* p4, FlxFunction* epsilon, const bool eval_once);
    RBRV_entry_RV_lognormal(const std::string& name, const tuint iID, py::dict config);
    ~RBRV_entry_RV_lognormal();

    const std::string get_type() const { return "logn"; }
    virtual void eval_para();
    virtual const tdouble transform_y2x(const tdouble y_val);
    virtual const tdouble transform_x2y(const tdouble& x_val);
    virtual const tdouble calc_pdf_x(const tdouble& x_val, const bool safeCalc=false);
    virtual const tdouble calc_pdf_x_log(const tdouble& x_val, const bool safeCalc=false);
    virtual const tdouble calc_cdf_x(const tdouble& x_val, const bool safeCalc=false);
    virtual const tdouble calc_sf_x(const tdouble& x_val, const bool safeCalc=false);
    virtual const tdouble calc_entropy();
    virtual const tdouble get_mean_current_config();
    virtual const tdouble get_sd_current_config();
    virtual const tdouble get_median_current_config();
    virtual const tdouble get_mode_current_config();
    virtual const bool check_x(const tdouble xV);
    virtual const bool search_circref(FlxFunction* fcr);
    virtual py::dict info();
    
    const tdouble get_CoeffOfVar_withoutEpsilon();
};

class PYBIND11_EXPORT RBRV_entry_RV_uniform : public RBRV_entry_RV_base {
  protected:
    FlxFunction* a;
    FlxFunction* b;
    bool eval_once;
    tdouble av;
    tdouble bv;
    
  public:
    RBRV_entry_RV_uniform(const std::string& name, const tuint iID, FlxFunction* a, FlxFunction* b, const bool eval_once);
    RBRV_entry_RV_uniform(const std::string& name, const tuint iID, py::dict config);
    ~RBRV_entry_RV_uniform();

    const std::string get_type() const { return "uniform"; }
    virtual void eval_para();
    virtual const tdouble transform_y2x(const tdouble y_val);
    virtual const tdouble transform_x2y(const tdouble& x_val);
    virtual const tdouble calc_pdf_x(const tdouble& x_val, const bool safeCalc=false);
    virtual const tdouble calc_cdf_x(const tdouble& x_val, const bool safeCalc=false);
    virtual const tdouble calc_icdf_x(const tdouble p_val);
    virtual const tdouble calc_sf_x(const tdouble& x_val, const bool safeCalc=false);
    virtual const tdouble calc_entropy();
    virtual const tdouble get_mean_current_config();
    virtual const tdouble get_sd_current_config();
    virtual const tdouble get_median_current_config() { return get_mean_current_config(); }
    virtual const bool check_x(const tdouble xV);
    virtual const bool search_circref(FlxFunction* fcr);
    virtual const tdouble get_HPD(const tdouble p);
    
    const tdouble Inv_cdf_x(const tdouble p);
};

class PYBIND11_EXPORT RBRV_entry_RV_Gumbel : public RBRV_entry_RV_base {
  protected:
    int methID;                // 0: p1=location, p2=scale; 1: p1=mean and p2=sd; 2: P(p1)=p2, P(p3)=p4
    FlxFunction* p1;
    FlxFunction* p2;
    FlxFunction* p3;
    FlxFunction* p4;
    bool eval_once;
    tdouble u;
    tdouble alpha;
    
  public:
    RBRV_entry_RV_Gumbel(const std::string& name, const tuint iID, const int methID, FlxFunction* p1, FlxFunction* p2, FlxFunction* p3, FlxFunction* p4, const bool eval_once);
    RBRV_entry_RV_Gumbel(const std::string& name, const tuint iID, py::dict config);
    ~RBRV_entry_RV_Gumbel();

    const std::string get_type() const { return "gumbel"; }
    virtual void eval_para();
    virtual const tdouble transform_y2x(const tdouble y_val);
    virtual const tdouble transform_x2y(const tdouble& x_val);
    virtual const tdouble calc_pdf_x(const tdouble& x_val, const bool safeCalc=false);
    virtual const tdouble calc_pdf_x_log(const tdouble& x_val, const bool safeCalc=false);
    virtual const tdouble calc_cdf_x(const tdouble& x_val, const bool safeCalc=false);
    virtual const tdouble calc_entropy();
    virtual const tdouble get_mean_current_config();
    virtual const tdouble get_sd_current_config();
    virtual const tdouble get_median_current_config();
    virtual const tdouble get_mode_current_config();
    virtual const bool check_x(const tdouble xV) { return true; }
    virtual const bool search_circref(FlxFunction* fcr);
    virtual py::dict info();
};

class PYBIND11_EXPORT RBRV_entry_RV_normal_trunc : public RBRV_entry_RV_base {
  protected:
    FlxFunction* m;
    FlxFunction* s;
    FlxFunction* a;
    FlxFunction* b;
    const bool eval_once;
    tdouble mV;
    tdouble sV;
    tdouble aV;
    tdouble bV;
    tdouble alpha, beta, q;
    
  public:
    RBRV_entry_RV_normal_trunc(const std::string& name, const tuint iID, FlxFunction* m, FlxFunction* s, FlxFunction* a, FlxFunction* b, const bool eval_once);
    ~RBRV_entry_RV_normal_trunc();

    const std::string get_type() const { return "normal_trunc"; }
    virtual void eval_para();
    virtual const tdouble transform_y2x(const tdouble y_val);
    virtual const tdouble transform_x2y(const tdouble& x_val);
    virtual const tdouble calc_pdf_x(const tdouble& x_val, const bool safeCalc=false);
    virtual const tdouble calc_pdf_x_log(const tdouble& x_val, const bool safeCalc=false);
    virtual const tdouble calc_cdf_x(const tdouble& x_val, const bool safeCalc=false);
    virtual const tdouble calc_sf_x(const tdouble& x_val, const bool safeCalc=false);
    virtual const tdouble calc_entropy();
    virtual const tdouble get_mean_current_config();
    virtual const tdouble get_sd_current_config();
    virtual const tdouble get_median_current_config();
    virtual const tdouble get_mode_current_config();
    virtual const bool check_x(const tdouble xV);
    virtual const bool search_circref(FlxFunction* fcr);
    virtual py::dict info();
};

class PYBIND11_EXPORT RBRV_entry_RV_beta : public RBRV_entry_RV_base {
  protected:
    bool is_mean;                // true if p1=mean and p2=sd
    FlxFunction* p1;
    FlxFunction* p2;
    FlxFunction* a;
    FlxFunction* b;
    bool eval_once;
    tdouble alpha;
    tdouble beta;
    tdouble av;
    tdouble bv;
  public:
    RBRV_entry_RV_beta(const std::string& name, const tuint iID, const bool is_mean, FlxFunction* p1, FlxFunction* p2, FlxFunction* a, FlxFunction* b, const bool eval_once);
    RBRV_entry_RV_beta(const std::string& name, const tuint iID, py::dict config);
    ~RBRV_entry_RV_beta();

    const std::string get_type() const { return "beta"; }
    virtual void eval_para();
    virtual const tdouble transform_y2x(const tdouble y_val);
    virtual const tdouble transform_x2y(const tdouble& x_val);
    virtual const tdouble calc_pdf_x(const tdouble& x_val, const bool safeCalc=false);
    virtual const tdouble calc_pdf_x_log(const tdouble& x_val, const bool safeCalc=false);
    virtual const tdouble calc_cdf_x(const tdouble& x_val, const bool safeCalc=false);
    virtual const tdouble calc_sf_x(const tdouble& x_val, const bool safeCalc=false);
    virtual const tdouble calc_entropy();
    virtual const tdouble get_mean_current_config();
    virtual const tdouble get_sd_current_config();
    virtual const tdouble get_median_current_config();
    virtual const tdouble get_mode_current_config();
    virtual const bool check_x(const tdouble xV);
    virtual const bool search_circref(FlxFunction* fcr);
    virtual py::dict info();
    
    const tdouble Inv_cdf_x(const tdouble p);
};

class PYBIND11_EXPORT RBRV_entry_RV_exponential : public RBRV_entry_RV_base {
  protected:
    FlxFunction* lambda;
    FlxFunction* epsilon;
    tdouble lambdaV;
    tdouble eps;

  public:
    RBRV_entry_RV_exponential(const std::string& name, const tuint iID, FlxFunction* lambda, FlxFunction* epsilon);
    RBRV_entry_RV_exponential(const std::string& name, const tuint iID, py::dict config);
    ~RBRV_entry_RV_exponential();

    const std::string get_type() const { return "exponential"; }
    virtual void eval_para();
    virtual const tdouble transform_y2x(const tdouble y_val);
    virtual const tdouble transform_x2y(const tdouble& x_val);
    virtual const tdouble calc_pdf_x(const tdouble& x_val, const bool safeCalc=false);
    virtual const tdouble calc_pdf_x_log(const tdouble& x_val, const bool safeCalc=false);
    virtual const tdouble calc_cdf_x(const tdouble& x_val, const bool safeCalc=false);
    virtual const tdouble calc_sf_x(const tdouble& x_val, const bool safeCalc=false);
    virtual const tdouble calc_entropy();
    virtual const tdouble get_mean_current_config();
    virtual const tdouble get_sd_current_config();
    virtual const tdouble get_median_current_config();
    virtual const tdouble get_mode_current_config();
    virtual const bool check_x(const tdouble xV);
    virtual const bool search_circref(FlxFunction* fcr);
    virtual py::dict info();
};

class PYBIND11_EXPORT RBRV_entry_RV_gamma : public RBRV_entry_RV_base {
  protected:
    const bool is_mean;                // true if p1=mean and p2=sd
    FlxFunction* p1;
    FlxFunction* p2;
    FlxFunction* epsilon;
    const bool eval_once;
    tdouble k;
    tdouble lambda;
    tdouble eps;
    
  public:
    RBRV_entry_RV_gamma(const std::string& name, const tuint iID, const bool is_mean, FlxFunction* p1, FlxFunction* p2, FlxFunction* epsilon, const bool eval_once);
    ~RBRV_entry_RV_gamma();

    const std::string get_type() const { return "gamma"; }
    virtual void eval_para();
    virtual const tdouble transform_y2x(const tdouble y_val);
    virtual const tdouble transform_x2y(const tdouble& x_val);
    virtual const tdouble calc_pdf_x(const tdouble& x_val, const bool safeCalc=false);
    virtual const tdouble calc_pdf_x_log(const tdouble& x_val, const bool safeCalc=false);
    virtual const tdouble calc_cdf_x(const tdouble& x_val, const bool safeCalc=false);
    virtual const tdouble calc_sf_x(const tdouble& x_val, const bool safeCalc=false);
    virtual const tdouble calc_entropy();
    virtual const tdouble get_mean_current_config();
    virtual const tdouble get_sd_current_config();
    virtual const tdouble get_median_current_config();
    virtual const tdouble get_mode_current_config();
    virtual const bool check_x(const tdouble xV);
    virtual const bool search_circref(FlxFunction* fcr);
    virtual py::dict info();
};

class PYBIND11_EXPORT RBRV_entry_RV_Poisson : public RBRV_entry_RV_base {
  protected:
    FlxFunction* mean;
    tdouble meanV;
  public:
    RBRV_entry_RV_Poisson(const std::string& name, const tuint iID, FlxFunction* mean);
    ~RBRV_entry_RV_Poisson();

    const std::string get_type() const { return "poisson"; }
    virtual void eval_para();
    virtual const tdouble transform_y2x(const tdouble y_val);
    virtual const tdouble calc_cdf_x(const tdouble& x_val, const bool safeCalc=false);
    virtual const tdouble calc_sf_x(const tdouble& x_val, const bool safeCalc=false);
    virtual const tdouble get_mean_current_config();
    virtual const tdouble get_sd_current_config();
    virtual const bool check_x(const tdouble xV);
    virtual const bool search_circref(FlxFunction* fcr);
};

class PYBIND11_EXPORT RBRV_entry_RV_Binomial : public RBRV_entry_RV_base {
  protected:
    FlxFunction* p;
    FlxFunction* N;
    const bool eval_once;
    tdouble _p;
    tuint _N;
    
  public:
    RBRV_entry_RV_Binomial(const std::string& name, const tuint iID, FlxFunction* p, FlxFunction* N, const bool eval_once);
    ~RBRV_entry_RV_Binomial();

    const std::string get_type() const { return "binomial"; }
    virtual void eval_para();
    virtual const tdouble transform_y2x(const tdouble y_val);
    virtual const tdouble calc_cdf_x(const tdouble& x_val, const bool safeCalc=false);
    virtual const tdouble get_mean_current_config();
    virtual const tdouble get_sd_current_config();
    virtual const bool check_x(const tdouble xV);
    virtual const bool search_circref(FlxFunction* fcr);
};

class PYBIND11_EXPORT RBRV_entry_RV_Cauchy : public RBRV_entry_RV_base {
  protected:
    FlxFunction* loc;
    FlxFunction* scale;
    tdouble locv;
    tdouble scalev;
    
  public:
    RBRV_entry_RV_Cauchy(const std::string& name, const tuint iID, FlxFunction* loc, FlxFunction* scale);
    ~RBRV_entry_RV_Cauchy();

    const std::string get_type() const { return "cauchy"; }
    virtual void eval_para();
    virtual const tdouble transform_y2x(const tdouble y_val);
    virtual const tdouble transform_x2y(const tdouble& x_val);
    virtual const tdouble calc_pdf_x(const tdouble& x_val, const bool safeCalc=false);
    virtual const tdouble calc_cdf_x(const tdouble& x_val, const bool safeCalc=false);
    virtual const tdouble get_mean_current_config();
    virtual const tdouble get_sd_current_config();
    virtual const tdouble get_median_current_config();
    virtual const tdouble get_mode_current_config();
    virtual const bool check_x(const tdouble xV) { return true; }
    virtual const bool search_circref(FlxFunction* fcr);
};

class PYBIND11_EXPORT RBRV_entry_RV_Weibull : public RBRV_entry_RV_base {
  protected:
    const bool is_mean;                // true if p1=mean and p2=sd
    FlxFunction* p1;
    FlxFunction* p2;
    FlxFunction* epsilon;
    const bool eval_once;
    tdouble k;
    tdouble lambda;
    tdouble eps;
    
    const tdouble get_mean_help();
    const tdouble get_sd_help();
    const tdouble get_cov_help();
  public:
    RBRV_entry_RV_Weibull(const std::string& name, const tuint iID, const bool is_mean, FlxFunction* p1, FlxFunction* p2, FlxFunction* epsilon, const bool eval_once);
    ~RBRV_entry_RV_Weibull();

    const std::string get_type() const { return "weibull"; }
    virtual void eval_para();
    virtual const tdouble transform_y2x(const tdouble y_val);
    virtual const tdouble transform_x2y(const tdouble& x_val);
    virtual const tdouble calc_pdf_x(const tdouble& x_val, const bool safeCalc=false);
    virtual const tdouble calc_pdf_x_log(const tdouble& x_val, const bool safeCalc=false);
    virtual const tdouble calc_cdf_x(const tdouble& x_val, const bool safeCalc=false);
    virtual const tdouble calc_sf_x(const tdouble& x_val, const bool safeCalc=false);
    virtual const tdouble calc_entropy();
    virtual const tdouble get_mean_current_config();
    virtual const tdouble get_sd_current_config();
    virtual const tdouble get_median_current_config();
    virtual const tdouble get_mode_current_config();
    virtual const bool check_x(const tdouble xV);
    virtual const bool search_circref(FlxFunction* fcr);
    virtual py::dict info();
    
};

class PYBIND11_EXPORT RBRV_entry_RV_ChiSquared : public RBRV_entry_RV_base {
  protected:
    FlxFunction* p1;
    const bool eval_once;
    tdouble dof;
    
  public:
    RBRV_entry_RV_ChiSquared(const std::string& name, const tuint iID, FlxFunction* p1, const bool eval_once);
    ~RBRV_entry_RV_ChiSquared();

    const std::string get_type() const { return "chisquared"; }
    virtual void eval_para();
    virtual const tdouble transform_y2x(const tdouble y_val);
    virtual const tdouble transform_x2y(const tdouble& x_val);
    virtual const tdouble calc_pdf_x(const tdouble& x_val, const bool safeCalc=false);
    virtual const tdouble calc_pdf_x_log(const tdouble& x_val, const bool safeCalc=false);
    virtual const tdouble calc_cdf_x(const tdouble& x_val, const bool safeCalc=false);
    virtual const tdouble calc_sf_x(const tdouble& x_val, const bool safeCalc=false);
    virtual const tdouble calc_entropy();
    virtual const tdouble get_mean_current_config();
    virtual const tdouble get_sd_current_config();
    virtual const tdouble get_median_current_config();
    virtual const tdouble get_mode_current_config();
    virtual const bool check_x(const tdouble xV);
    virtual const bool search_circref(FlxFunction* fcr);
};

class PYBIND11_EXPORT RBRV_entry_RV_Chi : public RBRV_entry_RV_base {
  protected:
    FlxFunction* p1;
    const bool eval_once;
    tdouble dof;

  public:
    RBRV_entry_RV_Chi(const std::string& name, const tuint iID, FlxFunction* p1, const bool eval_once);
    ~RBRV_entry_RV_Chi();

    const std::string get_type() const { return "chi"; }
    virtual void eval_para();
    virtual const tdouble transform_y2x(const tdouble y_val);
    virtual const tdouble transform_x2y(const tdouble& x_val);
    virtual const tdouble calc_pdf_x(const tdouble& x_val, const bool safeCalc=false);
    virtual const tdouble calc_pdf_x_log(const tdouble& x_val, const bool safeCalc=false);
    virtual const tdouble calc_cdf_x(const tdouble& x_val, const bool safeCalc=false);
    virtual const tdouble calc_sf_x(const tdouble& x_val, const bool safeCalc=false);
    virtual const tdouble calc_entropy();
    virtual const tdouble get_mean_current_config();
    virtual const tdouble get_sd_current_config();
    virtual const tdouble get_median_current_config();
    virtual const tdouble get_mode_current_config();
    virtual const bool check_x(const tdouble xV);
    virtual const bool search_circref(FlxFunction* fcr);
};


class PYBIND11_EXPORT RBRV_entry_RV_StudentsT : public RBRV_entry_RV_base {
  protected:
    FlxFunction* p1;
    bool eval_once;
    tdouble dof;
    
  public:
    RBRV_entry_RV_StudentsT(const std::string& name, const tuint iID, FlxFunction* p1, const bool eval_once);
    RBRV_entry_RV_StudentsT(const std::string& name, const tuint iID, py::dict config);
    ~RBRV_entry_RV_StudentsT();

    const std::string get_type() const { return "studentst"; }
    virtual void eval_para();
    virtual const tdouble transform_y2x(const tdouble y_val);
    virtual const tdouble transform_x2y(const tdouble& x_val);
    virtual const tdouble calc_pdf_x(const tdouble& x_val, const bool safeCalc=false);
    virtual const tdouble calc_cdf_x(const tdouble& x_val, const bool safeCalc=false);
    virtual const tdouble calc_entropy();
    virtual const tdouble get_mean_current_config() { return ZERO; }
    virtual const tdouble get_sd_current_config();
    virtual const tdouble get_median_current_config() { return ZERO; }
    virtual const tdouble get_mode_current_config() { return ZERO; }
    virtual const bool check_x(const tdouble xV) { return true; }
    virtual const bool search_circref(FlxFunction* fcr);
    virtual const tdouble get_HPD(const tdouble p);
    virtual py::dict info();
};

class PYBIND11_EXPORT RBRV_entry_RV_StudentsT_generalized : public RBRV_entry_RV_base {
  protected:
    int pid;        // 0: dof, loc, scale; 1: dof, loc, val_1, pr_1
    FlxFunction* p1;
    FlxFunction* p2;
    FlxFunction* p3;
    FlxFunction* p4;
    tdouble dof;
    tdouble loc;
    tdouble scale;

  public:
    RBRV_entry_RV_StudentsT_generalized(const std::string& name, const tuint iID, FlxFunction* nu, FlxFunction* locf, FlxFunction* scalef);
    RBRV_entry_RV_StudentsT_generalized(const std::string& name, const tuint iID, py::dict config);
    ~RBRV_entry_RV_StudentsT_generalized();

    const std::string get_type() const { return "studentstgen"; }
    virtual void eval_para();
    virtual const tdouble transform_y2x(const tdouble y_val);
    virtual const tdouble transform_x2y(const tdouble& x_val);
    virtual const tdouble calc_pdf_x(const tdouble& x_val, const bool safeCalc=false);
    virtual const tdouble calc_cdf_x(const tdouble& x_val, const bool safeCalc=false);
    virtual const tdouble calc_entropy();
    virtual const tdouble get_mean_current_config();
    virtual const tdouble get_sd_current_config();
    virtual const tdouble get_median_current_config();
    virtual const tdouble get_mode_current_config();
    virtual const bool check_x(const tdouble xV) { return true; }
    virtual const bool search_circref(FlxFunction* fcr);
    virtual const tdouble get_HPD(const tdouble p);
    virtual py::dict info();
};


class PYBIND11_EXPORT RBRV_entry_RV_logt : public RBRV_entry_RV_StudentsT_generalized {
  protected:

  public:
    RBRV_entry_RV_logt(const std::string& name, const tuint iID, FlxFunction* nu, FlxFunction* locf, FlxFunction* scalef);
    RBRV_entry_RV_logt(const std::string& name, const tuint iID, py::dict config);
    ~RBRV_entry_RV_logt();

    const std::string get_type() const { return "logt"; }
    virtual const tdouble transform_y2x(const tdouble y_val);
    virtual const tdouble transform_x2y(const tdouble& x_val);
    virtual const tdouble calc_pdf_x(const tdouble& x_val, const bool safeCalc=false);
    virtual const tdouble calc_cdf_x(const tdouble& x_val, const bool safeCalc=false);
    virtual const tdouble calc_entropy();
    virtual const tdouble get_mean_current_config();
    virtual const tdouble get_sd_current_config();
    virtual const tdouble get_median_current_config();
    virtual const tdouble get_mode_current_config();
    virtual const bool check_x(const tdouble xV) { return true; }
    virtual const tdouble get_HPD(const tdouble p);
    virtual py::dict info();
};


class PYBIND11_EXPORT RBRV_entry_RV_Laplace : public RBRV_entry_RV_base {
  protected:
    FlxFunction* locf;
    FlxFunction* scalef;
    tdouble loc;
    tdouble scale;

  public:
    RBRV_entry_RV_Laplace(const std::string& name, const tuint iID, FlxFunction* locf, FlxFunction* scalef);
    ~RBRV_entry_RV_Laplace();

    const std::string get_type() const { return "laplace"; }
    virtual void eval_para();
    virtual const tdouble transform_y2x(const tdouble y_val);
    virtual const tdouble transform_x2y(const tdouble& x_val);
    virtual const tdouble calc_pdf_x(const tdouble& x_val, const bool safeCalc=false);
    virtual const tdouble calc_pdf_x_log(const tdouble& x_val, const bool safeCalc=false);
    virtual const tdouble calc_cdf_x(const tdouble& x_val, const bool safeCalc=false);
    virtual const tdouble calc_sf_x(const tdouble& x_val, const bool safeCalc=false);
    virtual const tdouble calc_entropy();
    virtual const tdouble get_mean_current_config();
    virtual const tdouble get_sd_current_config();
    virtual const tdouble get_median_current_config();
    virtual const tdouble get_mode_current_config();
    virtual const bool check_x(const tdouble xV) { return true; }
    virtual const bool search_circref(FlxFunction* fcr);
    virtual const tdouble get_HPD(const tdouble p);
};


class PYBIND11_EXPORT RBRV_entry_RV_genpareto : public RBRV_entry_RV_base {
  protected:
    FlxFunction* xif;
    FlxFunction* locf;
    FlxFunction* scalef;
    tdouble xi;
    tdouble loc;
    tdouble scale;
    bool eval_once;

    const tdouble eval_cdf_help(const tdouble x_val);
    void free_mem();
  public:
    RBRV_entry_RV_genpareto(const std::string& name, const tuint iID, py::dict config);
    ~RBRV_entry_RV_genpareto();

    const std::string get_type() const { return "genpareto"; }
    virtual void eval_para();
    virtual const tdouble transform_y2x(const tdouble y_val);
    virtual const tdouble transform_x2y(const tdouble& x_val);
    virtual const tdouble calc_pdf_x(const tdouble& x_val, const bool safeCalc=false);
    virtual const tdouble calc_cdf_x(const tdouble& x_val, const bool safeCalc=false);
    virtual const tdouble calc_sf_x(const tdouble& x_val, const bool safeCalc=false);
    virtual const tdouble calc_entropy();
    virtual const tdouble get_mean_current_config();
    virtual const tdouble get_sd_current_config();
    virtual const tdouble get_median_current_config();
    virtual const tdouble get_mode_current_config();
    virtual const bool check_x(const tdouble xV);
    virtual const bool search_circref(FlxFunction* fcr);
    virtual py::dict info();
};



class PYBIND11_EXPORT RBRV_entry_RV_UserTransform : public RBRV_entry_RV_base {
  protected:
    const bool is_z2x;
    FlxFunction* t1;                   // z2x or y2z
    FlxFunction* t2;                 // x2z or z2y
    FlxFunction* dh;                 // dx2z/dx or dz2y/dz
    FlxFunction* checkXf;        // returns true if x is valid
    RBRV_entry_RV_base* rv_z;
    const bool manage_z;

    tVec tPL;        // constains the values of the parameters
    tdouble* const tPLp;
    
    const tdouble eval_para_fun(FlxFunction* fp, const tdouble value);
    const tdouble eval_cdf_sf(const bool is_cdf, const tdouble& x_val, const bool safeCalc);
  public:
    RBRV_entry_RV_UserTransform(const std::string& name, const tuint iID, const bool is_z2x, FlxFunction* t1, FlxFunction* t2, FlxFunction* dh, FlxFunction* checkXf, RBRV_entry_RV_base* rv_z, const bool manage_z=true);
    ~RBRV_entry_RV_UserTransform();

    const std::string get_type() const { return "usertransform"; }
    virtual const tdouble transform_y2x(const tdouble y_val);
    virtual const tdouble transform_x2y(const tdouble& x_val);
    virtual const tdouble calc_pdf_x(const tdouble& x_val, const bool safeCalc=false);
    virtual const tdouble calc_pdf_x_log(const tdouble& x_val, const bool safeCalc=false);
    virtual const tdouble calc_cdf_x(const tdouble& x_val, const bool safeCalc=false);
    virtual const tdouble calc_sf_x(const tdouble& x_val, const bool safeCalc=false);
    virtual const tdouble get_mean_current_config();
    virtual const tdouble get_sd_current_config();
    virtual const bool check_x(const tdouble xV);
    virtual const bool search_circref(FlxFunction* fcr);
    virtual py::dict info();
    void replace_rv_z(RBRV_entry_RV_base* rv_z_);
};


class PYBIND11_EXPORT RBRV_entry_RV_Truncated : public RBRV_entry_RV_base {
  protected:
    FlxFunction* a;
    FlxFunction* b;
    RBRV_entry_RV_base* rv_z;
    const bool manage_z;
    tdouble aV;
    tdouble bV;
    tdouble q;
    tdouble aV_cdf;

  public:
    RBRV_entry_RV_Truncated(const std::string& name, const tuint iID, FlxFunction* a, FlxFunction* b, RBRV_entry_RV_base* rv_z, const bool manage_z=true);
    ~RBRV_entry_RV_Truncated();

    const std::string get_type() const { return "truncated"; }
    virtual void eval_para();
    virtual const tdouble transform_y2x(const tdouble y_val);
    virtual const tdouble transform_x2y(const tdouble& x_val);
    virtual const tdouble calc_pdf_x(const tdouble& x_val, const bool safeCalc=false);
    virtual const tdouble calc_pdf_x_log(const tdouble& x_val, const bool safeCalc=false);
    virtual const tdouble calc_cdf_x(const tdouble& x_val, const bool safeCalc=false);
    virtual const tdouble get_mean_current_config();
    virtual const tdouble get_sd_current_config();
    virtual const bool check_x(const tdouble xV);
    virtual const bool search_circref(FlxFunction* fcr);
    virtual py::dict info();
};


class PYBIND11_EXPORT RBRV_entry_RV_maxminTransform : public RBRV_entry_RV_base {
  protected:
    const bool is_max;
    FlxFunction* n;                   // number of independent observations
    RBRV_entry_RV_base* rv_z;
    tdouble nV;
    
    const tdouble eval_cdf_sf(const bool is_cdf, const tdouble& x_val, const bool safeCalc);
  public:
    RBRV_entry_RV_maxminTransform(const std::string& name, const tuint iID, const bool is_max, FlxFunction* n, RBRV_entry_RV_base* rv_z);
    ~RBRV_entry_RV_maxminTransform();

    const std::string get_type() const { return "maxmintransform"; }
    virtual void eval_para();
    virtual const tdouble transform_y2x(const tdouble y_val);
    virtual const tdouble transform_x2y(const tdouble& x_val);
    virtual const tdouble calc_pdf_x(const tdouble& x_val, const bool safeCalc=false);
    virtual const tdouble calc_cdf_x(const tdouble& x_val, const bool safeCalc=false);
    virtual const tdouble calc_sf_x(const tdouble& x_val, const bool safeCalc=false);
    virtual const tdouble get_mean_current_config();
    virtual const tdouble get_sd_current_config();
    virtual const bool check_x(const tdouble xV);
    virtual const bool search_circref(FlxFunction* fcr);
};



class PYBIND11_EXPORT RBRV_entry_RV_quantiles : public RBRV_entry_RV_base {
  protected:
    tulong N_bins;
    std::vector<tdouble> pv;
    std::vector<tdouble> qv;
    tdouble* p_vec;   // of size N_bins+1
    tdouble* q_vec;   // of size N_bins+1
    tdouble* bin_rv_params; // of size N_bins
    tuint* N_vec;     // of size N_bins+1  TODO for accounting for prior information
    RBRV_entry_RV_base* tail_up;
    RBRV_entry_RV_base* tail_low;

    enum class interpol_type_t { uniform, pchip, bin_beta, bin_linear, pdf_linear };
    interpol_type_t interpol_type;

    std::optional<boost::math::interpolators::pchip<std::vector<tdouble>>> pchip_cdf;
    std::optional<boost::math::interpolators::pchip<std::vector<tdouble>>> pchip_icdf;

    void free_mem();
  public:
    RBRV_entry_RV_quantiles(const std::string& name, const tuint iID, py::dict config);
    virtual ~RBRV_entry_RV_quantiles();

    const std::string get_type() const { return "quantiles"; }
    virtual void eval_para();
    virtual const tdouble transform_y2x(const tdouble y_val);
    virtual const tdouble transform_x2y(const tdouble& x_val);
    virtual const tdouble calc_pdf_x(const tdouble& x_val, const bool safeCalc=false);
    virtual const tdouble calc_cdf_x(const tdouble& x_val, const bool safeCalc=false);
    virtual const tdouble get_mean_current_config();
    virtual const tdouble get_sd_current_config();
    virtual const tdouble get_median_current_config();
    virtual const tdouble get_mode_current_config();
    virtual const bool check_x(const tdouble xV);
    virtual const bool search_circref(FlxFunction* fcr);
};





//======================= Calculate Expectation numerically 1D ============================================

class FLXLIB_EXPORT calc_expectation_numerical_1D {
  private:
    FunBase* fun;
       
    pdouble calc_Interval(const tdouble& x_o, const tdouble& x_u, const tulong& N, RBRV_entry_RV_base& rnd);
  public:
    calc_expectation_numerical_1D(FunBase* fun) : fun(fun) {}
    /**
    * @brief calculates the expectation of fun
    * @param N_Interv Number of intervals
    * @param N Number of samples to calculate - distributed over intervals
    * if (N==0): just N_Interv are calculated
    */
    tdouble calc_expectation(const tulong N_Interv, const tulong N, const tdouble pRedistribute, RBRV_entry_RV_base& rnd, const tdouble LB, const tdouble UB);
};

//======================= Calculate Expectation numerically MCI ============================================

/**
* @brief calculates the expectation value of a function numerically
* ... for functions depending on more than one variable
*/
class FLXLIB_EXPORT calc_expectation_numerical_MCI : public FlxBoxBaseR {
  private:
    /**
    * @brief The function to calculate the expectation from
    */
    FunBase* funR;
    /**
    * @brief Dimension of the problem
    */
    tuint DIM;
    /**
    * @brief uniform distribution
    */
    RBRV_entry_RV_uniform* uD;
    /**
    * @brief one MCI-Integrationstep
    */
    inline const tdouble Integrationstep(flxVec &tempVec, RBRV_constructor& RndBoxN);
  public:
    /**
    * @brief calculates the expectation of fun
    * @param fun The function to calculate the expectation from
    * @param bv a vector with information about the random variables involved in this function
    * @param N_Interv number of intervals to split the sampling
    * @param N_mciI_1 number of MCI-samples per interval (for the first round)
    * @param N_mciI_2 number of MCI-samples in general (for the second round) - adaptive sampling
    * @return returns the expectation of fun
    */
    tdouble calc_expectation(FunBase* fun, RBRV_constructor& RndBoxN, const tulong N_Interv, const tulong N_mciI_1, const tulong N_mciI_2, const tdouble pRedistribute, const tdouble Bound);
};




