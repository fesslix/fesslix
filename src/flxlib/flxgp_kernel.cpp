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


#include "flxgp_kernel.h"
#include "flxMtx.h"

#if FLX_USE_NLOPT
  #include <nlopt.hpp>
#else
  #include "flxMtx_Eigen.h"
#endif

#include "flxparse.h"

#if FLX_DEBUG
int flxGPBox::Cinst = 0;
#endif


// ------------------------------------------------------------------------------------------------

void flxGP_mean_base::initialize_from_template(const flxVec& vt)
{
  const tuint N = vt.get_N();
  for (tuint i=0;i<N;++i) {
    pVec[i] = vt[i];
  }
}

void flxGP_mean_base::reset_pVec()
{
  if (get_Npara()>0) {
    pVec.set_zero();
    pVec[0] = ONE;
  }
}

flxGP_mean_0::flxGP_mean_0(const tuint Ndim)
: flxGP_mean_base(0,Ndim)
{

}

void flxGP_mean_0::assemble_F(FlxMtx& F, const tdouble* dm_ptr)
{
  throw FlxException_NotImplemented("flxGP_mean_0::assemble_F");
}

void flxGP_mean_0::assemble_f_vec(flxVec& f, const tdouble* dm_ptr)
{
  throw FlxException_NotImplemented("flxGP_mean_0::assemble_f_vec");
}

py::dict flxGP_mean_0::info()
{
  py::dict res;
  res["type"] = "zero";
  return res;
}

flxGP_mean_const::flxGP_mean_const(const tdouble val, const tuint Ndim)
: flxGP_mean_base(0,Ndim), val(val)
{

}

void flxGP_mean_const::assemble_F(FlxMtx& F, const tdouble* dm_ptr)
{
  throw FlxException_NotImplemented("flxGP_mean_const::assemble_F");
}

void flxGP_mean_const::assemble_f_vec(flxVec& f, const tdouble* dm_ptr)
{
  throw FlxException_NotImplemented("flxGP_mean_const::assemble_f_vec");
}

py::dict flxGP_mean_const::info()
{
  py::dict res;
  res["type"] = "const";
  res["value"] = val;
  return res;
}

flxGP_mean_ordinary::flxGP_mean_ordinary(FunReadFunUser* mean_f_, const tuint Ndim)
: flxGP_mean_base(1,Ndim), mean_f(mean_f_->get_fun_ptr())
{
  if (mean_f_->get_numbofpara()!=Ndim) {
    std::ostringstream ssV;
    ssV << "The mean function of the Gaussian process must take the same number of parameters as the dimension of the process (" << Ndim << ") and not " << mean_f_->get_numbofpara() << ".";
    throw FlxException("flxGP_mean_ordinary::flxGP_mean_ordinary", "Wrong number of parameters in mean function.", ssV.str() );
  }
}

flxGP_mean_ordinary::flxGP_mean_ordinary(FlxFunction* mean_f, const tuint Ndim)
: flxGP_mean_base(1,Ndim), mean_f(mean_f)
{

}

void flxGP_mean_ordinary::initialize_pVec(const tdouble mean_output)
{
  pVec[0] = mean_output;
}

const tdouble flxGP_mean_ordinary::eval_mean_f(const tdouble* pV)
{
  const tdouble* const tT = FunPara::ParaList;
  const tuint tTS = FunPara::ParaListSize;
  FunPara::ParaList = pV;
  FunPara::ParaListSize = Ndim;
  const tdouble res = mean_f->calc();
  FunPara::ParaList = tT;
  FunPara::ParaListSize = tTS;
  return res*pVec[0];
}

void flxGP_mean_ordinary::assemble_F(FlxMtx& F, const tdouble* dm_ptr)
{
  const tuint N = F.nrows();
  for (tuint i=0;i<N;++i) {
    F.operator()(i,0) = eval_mean_f(dm_ptr);
    dm_ptr += Ndim;
  }
}

void flxGP_mean_ordinary::assemble_f_vec(flxVec& f, const tdouble* dm_ptr)
{
  f[0] = eval_mean_f(dm_ptr);
}

py::dict flxGP_mean_ordinary::info()
{
  py::dict res;
  res["type"] = "ordinary";
  res["para_vec"] =  py_wrap_array_no_ownership<tdouble>(pVec.get_tmp_vptr(),pVec.get_N());;
  return res;
}

flxGP_mean_universal::flxGP_mean_universal(const tuint type_, const tuint Ndim)
: flxGP_mean_base((type_==0)?1:((type_==1)?(1+Ndim):((type_==2)?(1+Ndim+(Ndim*Ndim+Ndim)/2):(1+Ndim+(Ndim*Ndim+Ndim)/2+(Ndim*Ndim*Ndim+3*Ndim*Ndim+2*Ndim)/6))),Ndim), type(type_), normalizef(ONE)
{

}

void flxGP_mean_universal::initialize_pVec(const tdouble mean_output)
{
  pVec.set_zero();
  pVec[0] = ONE;
  normalizef = mean_output;
}

const tdouble flxGP_mean_universal::eval_mean_f(const tdouble* pV)
{
  tuint c = 0;
  tdouble res = pVec[c++];
  if (type==0) return res*normalizef;
  for (tuint i=0; i<Ndim; ++i) {
    res += pVec[c++]*pV[i];
  }
  if (type==1) return res*normalizef;
  for (tuint i=0; i<Ndim; ++i) {
    for (tuint j=0; j<=i; ++j) {
      res += pVec[c++] * pV[i]*pV[j];
    }
  }
  if (type==2) return res*normalizef;
  for (tuint i=0; i<Ndim; ++i) {
    for (tuint j=0; j<=i; ++j) {
      for (tuint k=0; k<=j; ++k) {
        res += pVec[c++] * pV[i]*pV[j]*pV[k];
      }
    }
  }
  if (c!=pVec.get_N()) throw FlxException_Crude("flxGP_mean_universal::eval_mean_f_01");
  if (type==3) return res*normalizef;
  throw FlxException("flxGP_mean_universal::eval_mean_f_02","Invalid type specifier.");
}

void flxGP_mean_universal::assemble_F(FlxMtx& F, const tdouble* dm_ptr)
{
  const tuint N = F.nrows();
  for (tuint i=0;i<N;++i) {
    tuint c = 0;
    // constant parameter
    F.operator()(i,c++) = normalizef;
    if (type>0) {
      for (tuint n=0; n<Ndim; ++n) {
        F.operator()(i,c++) = dm_ptr[n]*normalizef;
      }
      if (type>1) {
        for (tuint n=0; n<Ndim; ++n) {
          for (tuint j=0; j<=n; ++j) {
            F.operator()(i,c++) = dm_ptr[n]*dm_ptr[j]*normalizef;
          }
        }
        if (type>2) {
          for (tuint n=0; n<Ndim; ++n) {
            for (tuint j=0; j<=n; ++j) {
              for (tuint k=0; k<=j; ++k) {
                F.operator()(i,c++) = dm_ptr[n]*dm_ptr[j]*dm_ptr[k]*normalizef;
              }
            }
          }
        }
      }
    }
    dm_ptr += Ndim;
  }
}

void flxGP_mean_universal::assemble_f_vec(flxVec& f, const tdouble* dm_ptr)
{
  tuint c = 0;
  // constant parameter
  f[c++] = normalizef;
  if (type>0) {
    for (tuint n=0; n<Ndim; ++n) {
      f[c++] = dm_ptr[n]*normalizef;
    }
    if (type>1) {
      for (tuint n=0; n<Ndim; ++n) {
        for (tuint j=0; j<=n; ++j) {
          f[c++] = dm_ptr[n]*dm_ptr[j]*normalizef;
        }
      }
      if (type>2) {
        for (tuint n=0; n<Ndim; ++n) {
          for (tuint j=0; j<=n; ++j) {
            for (tuint k=0; k<=j; ++k) {
              f[c++] = dm_ptr[n]*dm_ptr[j]*dm_ptr[k]*normalizef;
            }
          }
        }
      }
    }
  }
}

py::dict flxGP_mean_universal::info()
{
  py::dict res;
  res["type"] = "universal";
  res["para_vec"] =  py_wrap_array_no_ownership<tdouble>(pVec.get_tmp_vptr(),pVec.get_N());
  res["normalizef"] = normalizef;
  return res;
}


// ------------------------------------------------------------------------------------------------

flxGP_kernel_base::flxGP_kernel_base(const tuint Ndim)
: Ndim(Ndim)
{

}

const tdouble flxGP_kernel_base::eval_kernel_corrl_f(const tdouble* pV)
{
  throw FlxException_NotImplemented("flxGP_kernel_base::eval_kernel_corrl_f");
}

const tdouble flxGP_kernel_base::eval_kernel_sd()
{
  throw FlxException_NotImplemented("flxGP_kernel_base::eval_kernel_sd");
}

void flxGP_kernel_base::set_sd(const tdouble sd)
{
  #if FLX_DEBUG
    if (fabs(ONE-sd)>1e-6) {
      throw FlxException_NotImplemented("flxGP_kernel_base::set_sd");
    }
  #endif
}

flxGP_kernel_user::flxGP_kernel_user(FunReadFunUser* kernel_f_, const tuint Ndim)
: flxGP_kernel_base(Ndim), kernel_f(kernel_f_->get_fun_ptr())
{
  if (kernel_f_->get_numbofpara()!=Ndim*2) {
    std::ostringstream ssV;
    ssV << "The kernel function of the Gaussian process must take the number of parameters twice as large as the dimension of the process (2*" << Ndim << ") and not " << kernel_f_->get_numbofpara() << ".";
    throw FlxException("flxGP_kernel_user::flxGP_kernel_user", "Wrong number of parameters in kernel function.", ssV.str() );
  }
}

flxGP_kernel_user::flxGP_kernel_user(FlxFunction* kernel_f, const tuint Ndim)
: flxGP_kernel_base(Ndim), kernel_f(kernel_f)
{

}

const tdouble flxGP_kernel_user::eval_kernel_f(const tdouble* pV)
{
  const tdouble* const tT = FunPara::ParaList;
  const tuint tTS = FunPara::ParaListSize;
  FunPara::ParaList = pV;
  FunPara::ParaListSize = Ndim*2;
  const tdouble res = kernel_f->calc();
  FunPara::ParaList = tT;
  FunPara::ParaListSize = tTS;
  return res;
}

py::dict flxGP_kernel_user::info()
{
  py::dict res;
  res["type"] = "fun";
  return res;
}


flxGP_kernel_auto::flxGP_kernel_auto(const iVec& tVec)
: flxGP_kernel_base(tVec.size()), pVec(tVec.size()+1), nVec(tVec.size()+1), tVec(tVec)
{
  pVec = ONE;
  nVec = ONE;
}

flxGP_kernel_auto::flxGP_kernel_auto(std::vector<std::string>& kernel_type_lst)
: flxGP_kernel_base(kernel_type_lst.size()), pVec(kernel_type_lst.size()+1), nVec(kernel_type_lst.size()+1), tVec(kernel_type_lst.size())
{
  pVec = ONE;
  nVec = ONE;
  for (tuint i=0;i<Ndim;++i) {
    if (kernel_type_lst[i]=="gauss") {
      tVec[i] = 0;
    } else if (kernel_type_lst[i]=="exp") {
      tVec[i] = 1;
    } else {
      throw FlxException("flxGP_kernel_auto::flxGP_kernel_auto", "unknown kernel type '" + kernel_type_lst[i] + "'");
    }
  }
}

void flxGP_kernel_auto::initialize_pVec(const tdouble sd_output, const flxVec& sd_vec)
{
  pVec = ONE;
  nVec[0] = sd_output;
  for (tuint i=0;i<Ndim;++i) {
    nVec[i+1] = sd_vec[i];
  }
}

void flxGP_kernel_auto::initialize_from_template(const flxVec& vt)
{
  pVec = vt;
}

void flxGP_kernel_auto::reset_pVec()
{
  pVec = ONE;
}

const tdouble flxGP_kernel_auto::eval_kernel_f(const tdouble* pV)
{
  return pow2(eval_kernel_sd()) * eval_kernel_corrl_f(pV);
}

const tdouble flxGP_kernel_auto::eval_kernel_corrl_f(const tdouble* pV)
{
  tdouble res = ONE;
  for (tuint i=0;i<Ndim;++i) {
    const tuint t = tVec[i];
    const tdouble corrl = nVec[i+1]*pVec[i+1];
    switch (t) {
      case 0:   // Gaussian
        res *= exp(-pow2((pV[i]-pV[i+Ndim])/corrl));
        break;
      case 1:   // exponential
        res *= exp(-fabs(pV[i]-pV[i+Ndim])/corrl);
        break;
      default:
        throw FlxException("flxGP_kernel_auto::eval_kernel_f");
    };
  }
  return res;
}

const tdouble flxGP_kernel_auto::eval_kernel_sd()
{
  return nVec[0]*pVec[0];
}

void flxGP_kernel_auto::set_sd(const tdouble sd)
{
  pVec[0] = sd;
  #if FLX_DEBUG
    if (nVec[0]!=ONE && fabs(ONE-sd)>1e-6) {
      std::ostringstream ssV;
      ssV << "nVec[0]=" << nVec[0] << "   sd=" << sd;
      throw FlxException("flxGP_kernel_auto::set_sd", ssV.str());
    }
  #endif
}

py::dict flxGP_kernel_auto::info()
{
  py::dict res;
  py::list tlst;
  for (tuint i=0;i<Ndim;++i) {
    std::string ts;
    switch (tVec[i]) {
      case 0:
        ts = "gauss";
        break;
      case 1:
        ts = "exp";
        break;
      default:
        throw FlxException_Crude("flxGP_kernel_auto::info");
    }
    tlst.append(ts);
  }
  res["type"] = tlst;
  res["para_vec"] =  py_wrap_array_no_ownership<tdouble>(pVec.get_tmp_vptr(),pVec.get_N());
  res["n_vec"] =  py_wrap_array_no_ownership<tdouble>(nVec.get_tmp_vptr(),nVec.get_N());
  res["kernel_sd"] = this->eval_kernel_sd();
  return res;
}


// ------------------------------------------------------------------------------------------------

flxGPProj_base::flxGPProj_base(const std::string& name, const tuint Ndim)
: name(name), Ndim(Ndim)
{

}


// ------------------------------------------------------------------------------------------------

flxGPProj::flxGPProj(const std::string& name, const tuint Ndim, flxGP_mean_base* gp_mean, flxGP_kernel_base* gp_kernel, const bool use_LSE)
: flxGPProj_base(name,Ndim), gp_mean(gp_mean), gp_kernel(gp_kernel), use_LSE(use_LSE), obsv_changed(true), N_obsv(0), dm_ptr(NULL), do_ptr(NULL),
  cov_mtx_obsv(NULL), lt_mtx_obsv(NULL), ltinv_mtx_obsv(NULL), covinv_mtx_obsv(NULL), FRinvF_mtx(NULL), FRinvF_cdc_mtx(NULL),
  noise_val(ZERO), lt_mtx_obsv_ldet(ZERO), lpr_obsv(ZERO), do_0mean(NULL), alpha(NULL), Fmtx(NULL), d_info(NULL)
{
  gp_mean->initialize_pVec(ONE);
}

flxGPProj::~flxGPProj()
{
  if (cov_mtx_obsv) delete cov_mtx_obsv;
  if (covinv_mtx_obsv) delete covinv_mtx_obsv;
  if (lt_mtx_obsv) delete lt_mtx_obsv;
  if (ltinv_mtx_obsv) delete ltinv_mtx_obsv;
  if (FRinvF_mtx) delete FRinvF_mtx;
  if (FRinvF_cdc_mtx) delete FRinvF_cdc_mtx;
  if (do_0mean) delete do_0mean;
  if (alpha) delete alpha;
  if (Fmtx) delete Fmtx;
  if (d_info) delete d_info;
  delete gp_mean;
  delete gp_kernel;
}

const tdouble flxGPProj::assemble_observations_help()
{
  // assemble covariance/correlation matrix
    // Define parameter vector for covariance kernel
      flxVec pV(Ndim*2);
      tdouble* pVp = pV.get_tmp_vptr();
      flxVec pV_1(pVp,Ndim);
      flxVec pV_2(pVp+Ndim,Ndim);
      tdouble* cmop = cov_mtx_obsv->get_VecPointer();
    // loop over data-points
      // for LSE » eval scaling factors
        const tdouble kernel_sd_ini = use_LSE?(gp_kernel->eval_kernel_sd()):ONE;
        const tdouble sigma_Z_2_ini = use_LSE?pow2(kernel_sd_ini):ONE;
        const tdouble sigma_eps_2_ini = pow2(noise_val);
        const tdouble sigma_Y_2_ini = sigma_Z_2_ini + sigma_eps_2_ini;
        const tdouble scale_Z_ini = sigma_Z_2_ini/sigma_Y_2_ini;
        //const tdouble scale_eps = sigma_eps_2_ini/sigma_Y_2_ini;
      for (size_t i=0;i<N_obsv;++i) {
        // update parameter vector
          pV_1 = flxVec(dm_ptr+i*Ndim,Ndim);
        for (size_t j=0;j<=i;++j) {
          // update parameter vector
            pV_2 = flxVec(dm_ptr+j*Ndim,Ndim);
          if (use_LSE) {
            *cmop = scale_Z_ini*gp_kernel->eval_kernel_corrl_f(pVp);
          } else {
            *cmop = gp_kernel->eval_kernel_f(pVp);
          }
          ++cmop;
        }
        if (use_LSE) {
          *(cmop-1) = ONE;   // by definition ( equals » *(cmop-1) += scale_eps;
        } else {
          *(cmop-1) += sigma_eps_2_ini;
        }
      }
  // compute cholesky decomposition
    try {
      lt_mtx_obsv->CholeskyDec(*cov_mtx_obsv,false);
      lt_mtx_obsv_ldet = lt_mtx_obsv->det_log();
    } catch (FlxException &e) {
      return log(ZERO);
    }
  // optimize mean parameters
      if (use_LSE) {
        // fully invert R
          *ltinv_mtx_obsv = *lt_mtx_obsv;
          ltinv_mtx_obsv->Invert();
          covinv_mtx_obsv->assign_LTL(*ltinv_mtx_obsv);
        if (gp_mean->get_Npara()>0) {
          // invert FRinvF
            MtxProd_BTKB_mtx(*Fmtx,*covinv_mtx_obsv,*FRinvF_mtx);
            FRinvF_cdc_mtx->CholeskyDec(*FRinvF_mtx,false);
          flxVec do_vec(do_ptr,N_obsv);
          lt_mtx_obsv->MultInv(do_vec,*do_0mean);
          lt_mtx_obsv->TransMultInv(*do_0mean,*do_0mean);  // do_0mean is used to store temporary results
          // as an alternative to the two lines above (possibly nummerically less stable)
            // covinv_mtx_obsv->MultMv(do_vec,*do_0mean);  // do_0mean is used to store temporary results
          flxVec& tpvm = gp_mean->get_paraVec();
          Fmtx->TransposeMmultVec(*do_0mean,tpvm);
          FRinvF_cdc_mtx->MultInv(tpvm,tpvm);
          FRinvF_cdc_mtx->TransMultInv(tpvm,tpvm);
        }
      }
  // remove mean-trend from assemble_observations
    for (size_t i=0;i<N_obsv;++i) {
      pV_1 = flxVec(dm_ptr+i*Ndim,Ndim);
      do_0mean->operator[](i) = do_ptr[i] - gp_mean->eval_mean_f(pV_1.get_tmp_vptr());
    }
  // compute the log-likelihood of the observation
    pdouble nres(lt_mtx_obsv_ldet);  // = log(det(L)) = 0.5*log(det(sigma_Z^2 R + sigma_eps^2 I))  ; for use_LSE = 0.5*log(det(Q))
    // evaluate the alpha-vector
      alpha->operator=(*do_0mean);
      lt_mtx_obsv->MultInv(*alpha,*alpha);
      lt_mtx_obsv->TransMultInv(*alpha,*alpha);
      // compute remaining terms of log-likelihood
        const tdouble mhd2 = (*alpha) * (*do_0mean);  // we need this when evaluating the log-likelihood
        const tdouble sigma_Y_2 = use_LSE?(mhd2/N_obsv):ONE;
        nres += 0.5*log(2*PI*sigma_Y_2)*N_obsv;
        nres += 0.5*mhd2/sigma_Y_2;
        const tdouble neg_logl = -(nres.cast2double());
      // for use_LSE, evaluate variance estimate
        if (use_LSE) {
          const tdouble sigma_Z_2 = sigma_Y_2 * scale_Z_ini;
          const tdouble sigma_eps_2 = max(sigma_Y_2 - sigma_Z_2,ZERO);
          noise_opt_stream << "»» LSE-results  logl=" << GlobalVar.Double2String(neg_logl) << "  "
            << "»»  sd_obsv [" << GlobalVar.Double2String(sqrt(sigma_Y_2)) << " <- " << GlobalVar.Double2String(sqrt(sigma_Y_2_ini)) << "] "
            << "»»  sd_Z [" << GlobalVar.Double2String(sqrt(sigma_Z_2)) << " <- " << GlobalVar.Double2String(sqrt(sigma_Z_2_ini)) << "] "
            << "»»  sd_noise [" << GlobalVar.Double2String(sqrt(sigma_eps_2)) << " <- " << GlobalVar.Double2String(sqrt(sigma_eps_2_ini)) << "] " << std::endl;
          gp_kernel->set_sd(sqrt(sigma_Z_2));
          noise_val = sqrt(sigma_eps_2);
          alpha->operator/=(sigma_Y_2);
        }
  // GlobalVar.slogcout(1) << "  flxGPProj::assemble_observations_74 " << res << "  " << 0.5*((*alpha) * (*do_0mean))/sigma2 << "  " << (lt_mtx_obsv_ldet+0.5*log(sigma2)*N_obsv) << "  " << 0.5*log(2*PI)*N_obsv << std::endl;
  return neg_logl;
}

double gp_likeli_f_nv (const tdouble lnv, void *params)
{
  flxGPProj *p = (flxGPProj *)params;
  // remember starting values
    const tdouble nv = exp(lnv);
    const tdouble sdZ_ini = p->gp_kernel->eval_kernel_sd();
  // consistency checks
    if (fabs(p->noise_val-nv)/sdZ_ini<=GlobalVar.TOL() && p->obsv_changed==false) return p->get_log_likeli_obsv();
    if (nv<sdZ_ini*1e-6) {
      throw FlxException_math("flxGPProj::likeli_f_b01","Noise value is effectively zero.");
    }
  // assemble/evaluate for specified noise value
    p->unassemble();
    p->noise_val = nv;
    tdouble res;
    try {
      res = p->assemble_observations_help();
    } catch (FlxException& e) {
      // restore original values (we just want to optimize noise conditional on the initial configuration)
        if (p->use_LSE) {
          p->gp_kernel->set_sd(sdZ_ini);
        }
        p->noise_val = nv;
      throw;
    }

  // restore original values (we just want to optimize noise conditional on the initial configuration)
    if (p->use_LSE) {
      p->gp_kernel->set_sd(sdZ_ini);
    }
    p->noise_val = nv;

  // output
    //GlobalVar.slogcout(5) << "      " << GlobalVar.Double2String(res) << ": " << nv << std::endl;

  // error checking
    if (std::isnan(res)) {
        throw FlxException_math("flxGPProj::likeli_f_b02","Negative log-likelihood is NaN.");
    }
    if (std::isinf(res)) {
        throw FlxException_math("flxGPProj::likeli_f_b03","Negative log-likelihood is infinite.");
    }
  // keep track of 'best guess'
    if (res > p->noise_best_guess_logl) {
      p->noise_best_guess_val = nv;
      p->noise_best_guess_sdZ = sdZ_ini;
      p->noise_best_guess_logl = res;
      p->noise_opt_stream << "    *** best guess ***  logl=" << GlobalVar.Double2String(res) << "  noise=" << GlobalVar.Double2String(nv) << "  sd_Z=" << GlobalVar.Double2String(sdZ_ini) << std::endl;
    }
  return -res;
}

void flxGPProj::assemble_observations(const bool initialize_pVec, const bool optimize_noise_val)
{
  if (obsv_changed == false) return;
  // reset logging stream
    noise_opt_stream.str("");
    noise_opt_stream.clear();
  if (d_info==NULL) {
      std::ostringstream ssV;
      ssV << "The Gaussian process '" << name << "' has not yet been conditioned on data.";
      throw FlxException("flxGPProj::assemble_observations_01", "Gaussian process not conditioned on data.", ssV.str() );
  }
  // get the observed model input
    N_obsv = 0;
    dm_ptr = d_info->get_data_ptr_mtx(N_obsv,Ndim);
  // get the observed output vector
    do_ptr = d_info->get_data_ptr_vec(N_obsv);
  // initialize mean function and kernel
      if (initialize_pVec) {
        flxVec dov(do_ptr,N_obsv);
        const tdouble dov_mean = dov.get_Mean();
        // mean
          if (!use_LSE) {
            gp_mean->initialize_pVec(dov_mean);
          } else {
            gp_mean->initialize_pVec(ONE);
          }
        // kernel
          flxVec sd_vec(Ndim);
          flxVec dmv(N_obsv);
          for (tuint i=0;i<Ndim;++i) {
            for (size_t j=0;j<N_obsv;++j) {
              dmv[j] = dm_ptr[j*Ndim+i];
            }
            sd_vec[i] = dmv.get_sd(dmv.get_Mean());
          }
          const tdouble sample_sd_out = dov.get_sd(dov_mean);
          gp_kernel->initialize_pVec(use_LSE?ONE:sample_sd_out,sd_vec);
          if (use_LSE) {
            gp_kernel->set_sd(sample_sd_out);
          }
        // noise_val
          if (optimize_noise_val) {
            noise_val = 0.1*sample_sd_out;
          } else {
            noise_val = 1e-4*sample_sd_out;
          }
      }
  // allocate memory
    // allocate covariance matrix that has correct size
      if (cov_mtx_obsv) {
        if (cov_mtx_obsv->nrows() != N_obsv) {
          delete cov_mtx_obsv;
          cov_mtx_obsv = NULL;
        }
      }
      if (cov_mtx_obsv == NULL) {
        cov_mtx_obsv = new FlxMtxSym(N_obsv);
      }
    // Cholesky decomposition: allocate lower-triangular matrix that has correct size
      if (lt_mtx_obsv) {
        if (lt_mtx_obsv->nrows() != N_obsv) {
          delete lt_mtx_obsv;
          lt_mtx_obsv = NULL;
          if (use_LSE) {
            delete ltinv_mtx_obsv;
            ltinv_mtx_obsv = NULL;
            delete covinv_mtx_obsv;
            covinv_mtx_obsv = NULL;
            delete FRinvF_mtx;
            FRinvF_mtx = NULL;
            delete FRinvF_cdc_mtx;
            FRinvF_cdc_mtx = NULL;
          }
        }
      }
      if (lt_mtx_obsv == NULL) {
        lt_mtx_obsv = new FlxMtxLTri(N_obsv);
        if (use_LSE) {
          ltinv_mtx_obsv = new FlxMtxLTri(N_obsv);
          covinv_mtx_obsv = new FlxMtxSym(N_obsv);
          FRinvF_mtx = new FlxMtxSym(gp_mean->get_Npara());
          FRinvF_cdc_mtx = new FlxMtxLTri(gp_mean->get_Npara());
        }
      }
    // allocate vector to store mean-corrected observations
      if (do_0mean) {
        if (do_0mean->get_N() != N_obsv) {
          delete do_0mean;
          do_0mean = NULL;
        }
      }
      if (do_0mean == NULL) do_0mean = new flxVec(N_obsv);
    // allocate alpha-vector
      if (alpha) {
        if (alpha->get_N() != N_obsv) {
          delete alpha;
          alpha = NULL;
        }
      }
      if (alpha == NULL) alpha = new flxVec(N_obsv);
    // allocate and assemble F-matrix (assembly is done only once!)
      if (use_LSE) {
        if (gp_mean->get_Npara()>0) {
          if (Fmtx) {
            if (Fmtx->nrows()!=N_obsv || Fmtx->ncols()!=gp_mean->get_Npara()) {
              delete Fmtx;
              Fmtx = NULL;
            }
          }
          if (Fmtx==NULL) {
            Fmtx = new FlxMtx(N_obsv, gp_mean->get_Npara());
            gp_mean->assemble_F(*Fmtx,dm_ptr);
          } else {
            if (initialize_pVec) {
              // NOTE this is INCORRECT, if the dm_ptr changes but its size remains the same (unlikely in practice!)
              gp_mean->assemble_F(*Fmtx,dm_ptr);
            }
          }
        }
      }

  // minimize noise parameter
    if (optimize_noise_val) {
      tdouble lnvt = log(noise_val);
      noise_best_guess_val = noise_val;
      noise_best_guess_sdZ = use_LSE?(gp_kernel->eval_kernel_sd()):ONE;
      noise_best_guess_logl = -std::numeric_limits<tdouble>::infinity();
      try {
        flx_optim(log(noise_val/10), log(noise_val*10), lnvt, gp_likeli_f_nv, this, true, true, 100, 20, 1e-4, 1e-4, &noise_opt_stream );
      } catch (FlxException& e) {
        noise_opt_stream << "»»» setting noise to best guess (" << GlobalVar.Double2String(noise_best_guess_val) << ")" << std::endl;
      }
      // ensure that minimum is set
        noise_val = noise_best_guess_val;
        if (use_LSE) {
          gp_kernel->set_sd(noise_best_guess_sdZ);
        }
        lpr_obsv = assemble_observations_help();
        // consistency checks
          if (fabs(lpr_obsv-noise_best_guess_logl)/max(ONE,fabs(noise_best_guess_logl))>1e-6) {
            throw FlxException_Crude("flxGPProj::assemble_observations_02");
          }
    } else {
      lpr_obsv = assemble_observations_help();
    }
  obsv_changed = false;
}

const tdouble flxGPProj::register_observation(const flxGP_data_base& d_info_, const bool initialize_pVec, const bool optimize_noise_val)
{
  if (d_info) delete d_info;
  d_info = d_info_.copy();
  unassemble();
  assemble_observations(initialize_pVec, optimize_noise_val);
  return lpr_obsv;
}

void flxGPProj::register_noise(const tdouble noise_sd)
{
  noise_val = noise_sd;
  unassemble();
}

#if FLX_USE_NLOPT
double gp_likeli_f (unsigned n, const double *v, double *grad, void *params)
#else
double gp_likeli_f(const gsl_vector *v, void *params)
#endif
{
  flxGPProj *p = (flxGPProj *)params;

  #if FLX_USE_NLOPT
  #else
  const size_t n = v->size;
  #endif

  // assign parameters
    bool bchange = false;
    tuint c = 0;
    if (!p->use_LSE) {
      flxVec& pvm = p->gp_mean->get_paraVec();
      for (tuint i=0;i<pvm.get_N();++i) {
        #if FLX_USE_NLOPT
        const tdouble pv = v[c++];
        #else
        const tdouble pv = gsl_vector_get(v, c++);
        #endif
        if (pvm[i] != pv) bchange = true;
        pvm[i] = pv;
      }
    }
    flxVec* pkp = p->gp_kernel->get_paraVec();
    for (tuint i=0;i<p->gp_kernel->get_Npara();++i) {
      if (p->use_LSE && i==0) continue;
      #if FLX_USE_NLOPT
      const tdouble pv = exp(v[c++]);
      #else
      const tdouble pv = exp(gsl_vector_get (v, c++));
      #endif
      // make sure that correlation length is finite
        if (!std::isfinite(pv)) {
          return std::numeric_limits<tdouble>::infinity();
        }
      if (fabs(pkp->operator[](i)-pv)>GlobalVar.TOL()) bchange = true;
      pkp->operator[](i) = pv;
    }
    if (c<n) {
      #if FLX_USE_NLOPT
      const tdouble pv = exp(v[c++]);
      #else
      const tdouble pv = exp(gsl_vector_get (v, c++));
      #endif
      if (fabs(p->noise_val-pv)>GlobalVar.TOL()) bchange = true;
      p->noise_val = pv;
    }

  // store initial parameters (for use_LSE)
    const tdouble noise_sd_ini = (p->use_LSE)?(p->noise_val):ZERO;
    const tdouble kernel_sd_ini = (p->use_LSE)?(p->gp_kernel->eval_kernel_sd()):ZERO;

  // evaluate likelihood
    if (bchange) p->unassemble();
    tdouble res;
    try {
      res = p->get_log_likeli_obsv();
    } catch (FlxException& e) {
      // restore original values (we just want to optimize noise conditional on the initial configuration)
        if (p->use_LSE && p->keep_LSE_results==false) {
          p->gp_kernel->set_sd(kernel_sd_ini);
          p->noise_val = noise_sd_ini;
        }
      throw;
    }

    // restore original values (we just want to optimize noise conditional on the initial configuration)
      if (p->use_LSE && p->keep_LSE_results==false) {
        p->gp_kernel->set_sd(kernel_sd_ini);
        p->noise_val = noise_sd_ini;
      }

  #if FLX_USE_NLOPT
      p->para_opt_stream << "    flxGPProj::likeli_f_89 " << res << "   " << flxVec(v,n) << "   " << (bchange?"yes":"no") << std::endl;
  #endif

  // error checking
    if (std::isnan(res)) {
        throw FlxException_math("flxGPProj::likeli_f_01");
    }
    if (std::isinf(res)) {
        throw FlxException_math("flxGPProj::likeli_f_02");
    }
  return -res;
}


const tdouble flxGPProj::optimize_help(const tdouble step_size, const tuint iterMax, const bool opt_noise, const bool output_ini, std::ostream& ostrm)
{
  // make sure that process is assembeld before log-likelihood is evaluated
    assemble_observations(false,false);
    keep_LSE_results = false;
    // we can only optimize the process, if it is conditioned on some data
  // count the number of uncertain parameters
    const tuint Npara = (use_LSE?(gp_kernel->get_Npara()-1):(gp_mean->get_Npara() + gp_kernel->get_Npara())) + (opt_noise?1:0);
    const tuint NparaMean   = use_LSE?0:(gp_mean->get_Npara());
    const tuint NparaKernel = gp_kernel->get_Npara();

  #if FLX_USE_NLOPT

    nlopt::opt opt(nlopt::LN_NELDERMEAD, Npara);
    // nlopt::opt opt(nlopt::LN_SBPLX, Npara);

    // assign lower bound
      // std::vector<double> lb(Npara, -HUGE_VAL);
      // for (tuint i=0;i<Npara;++i) {
      //   lb[Npara-i-1] = GlobalVar.TOL();
      // }
      // if (!use_LSE) {
      //   lb[Npara-(gp_kernel->get_Npara())] = GlobalVar.TOL();
      // }
      // opt.set_lower_bounds(lb);

    // specify function to minimize
      opt.set_min_objective(gp_likeli_f,this);

    // stopping criteria
      opt.set_xtol_rel(1e-4);
      opt.set_ftol_rel(1e-4);
      opt.set_maxeval(iterMax);

    // starting vector
      std::vector<double> x_ini(Npara);
      {
        tuint c = 0;
        if (!use_LSE) {
          for (tuint i=0;i<NparaMean;++i) {
            const tdouble pv = gp_mean->get_paraVec()[i];
            x_ini[c++] = pv;
          }
        }
        const flxVec* pkp = gp_kernel->get_paraVec();
        for (tuint i=0;i<NparaKernel;++i) {
          if (use_LSE && i==0) continue;
          const tdouble pv = log(pkp->operator[](i));
          if (std::isinf(pv)) {
            throw FlxException("flxGPProj::optimize_help_NLopt_01", GlobalVar.Double2String(pkp->operator[](i)));
          }
          x_ini[c++] = pv;
        }
        if (opt_noise) {
          x_ini[c++] = log(max(noise_val,1e-4*(gp_kernel->eval_kernel_sd())));
        }
      }
      if (output_ini) {
          const tdouble lres_ini = -gp_likeli_f(Npara,&x_ini[0],NULL,this);
          ostrm << "       initial point estimate: " << GlobalVar.Double2String(lres_ini) << " at ( ";
          for (tuint i=0;i<x_ini.size(); ++i) {
              if (i>0) {
                  ostrm << ", ";
              }
              ostrm << GlobalVar.Double2String(x_ini[i]);
          }
          ostrm << " ) dim=" << Npara << std::endl;
      }
      std::vector<double> x(x_ini);

    // perform the minimization
      double minf, lres;

      try {
        nlopt::result result = opt.optimize(x, minf);

        lres = -minf;

          ostrm << "      " << " MLE: " << lres << " at ( ";
          for (tuint i=0;i<Npara; ++i) {
              if (i>0) {
                  ostrm << ", ";
              }
              ostrm << x[i];
          }
          ostrm << " ) :: niter = " << opt.get_numevals() << ", dim = " << Npara << std::endl;

        // make sure MLE-model parameters are actually set
          keep_LSE_results = true;
          const tdouble res_err = fabs(gp_likeli_f(Npara,&x[0],NULL,this)+lres)/max(ONE,fabs(lres));
          if (res_err>=1e-6) {
            std::ostringstream ssV;
            ssV << GlobalVar.Double2String(lpr_obsv) << "  " << GlobalVar.Double2String(lres) << "  " << GlobalVar.Double2String(res_err);
            throw FlxException("flxGPProj::optimize_help_NLopt_02", ssV.str());
          }

      } catch(...) {
        // re-run evaluation for initial model
          keep_LSE_results = true;
          lres = -gp_likeli_f(Npara,&x_ini[0],NULL,this);
        throw;
      }

  #else
    // Initialize optimization
      const gsl_multimin_fminimizer_type *T = gsl_multimin_fminimizer_nmsimplex2;
      gsl_multimin_fminimizer *s = NULL;
      gsl_vector *ss, *x;
      gsl_multimin_function minex_func;

      size_t iter = 0;
      int status;
      double size;
      tdouble lres;

      /* Starting point */
      flxVec pvec_ini(Npara);
      {
        x = gsl_vector_alloc(Npara);
        tuint c = 0;
        if (!use_LSE) {
          for (tuint i=0;i<NparaMean;++i) {
            const tdouble pv = gp_mean->get_paraVec()[i];
            pvec_ini[c] = pv;
            gsl_vector_set (x, c++, pv);
          }
        }
        const flxVec* pkp = gp_kernel->get_paraVec();
        for (tuint i=0;i<NparaKernel;++i) {
          if (use_LSE && i==0) continue;
          const tdouble pv = log(pkp->operator[](i));
          if (std::isinf(pv)) {
            throw FlxException("flxGPProj::optimize_help_GSL_01");
          }
          pvec_ini[c] = pv;
          gsl_vector_set (x, c++, pv);
        }
        gsl_vector_set (x, c++, log(max(noise_val,1e-4*(gp_kernel->eval_kernel_sd()))));
      }

      /* Set initial step sizes to 1 */
      ss = gsl_vector_alloc(Npara);
      gsl_vector_set_all(ss, step_size);

      /* Initialize method and iterate */
      minex_func.n = Npara;
      minex_func.f = gp_likeli_f;
      minex_func.params = this;

    try {

    // output inital point
        lres = -gp_likeli_f(x,this);
        if (output_ini) {
            ostrm << "    initial point estimate: " << GlobalVar.Double2String(lres) << " at ( ";
            for (tuint i=0;i<pvec_ini.get_N(); ++i) {
                if (i>0) {
                    ostrm << ", ";
                }
                ostrm << GlobalVar.Double2String(pvec_ini[i]);
            }
            ostrm << " ) dim=" << Npara << std::endl;
        }

    // perform iteration
      s = gsl_multimin_fminimizer_alloc (T, Npara);
      gsl_multimin_fminimizer_set (s, &minex_func, x, ss);

      do
          {
          iter++;
          status = gsl_multimin_fminimizer_iterate(s);

          if (status)
              break;

          size = gsl_multimin_fminimizer_size (s);
          status = gsl_multimin_test_size (size, 1e-2);

          if (status == GSL_SUCCESS)
          {
              lres = -gp_likeli_f(s->x,this);
              break;
          }

          if (!use_LSE) {
            printf ("%5d %10.3e %10.3e f() = %7.3f size = %.3f\n",
                    (int)iter,
                    gsl_vector_get (s->x, 0),
                    exp(gsl_vector_get (s->x, NparaMean)),
                    s->fval, size);
          }

          }
      while (status == GSL_CONTINUE && iter < iterMax);

      // output results
        if (status!=GSL_SUCCESS) {
            lres = -(s->fval);
        }

        ostrm << "      " << ((status==GSL_SUCCESS)?' ':'~') << "MLE: " << lres << " at ( ";
        for (tuint i=0;i<pvec_ini.get_N(); ++i) {
            if (i>0) {
                ostrm << ", ";
            }
            tdouble pv = gsl_vector_get(s->x, i);
            if (i>=NparaMean) pv = exp(pv);
            ostrm << pv;
        }
        ostrm << " ) :: niter = " << iter << ", dim = " << Npara << std::endl;

      // make sure MLE-model parameters are actually set
        keep_LSE_results = true;
        const tdouble res_err = fabs(gp_likeli_f(s->x,this)+lres)/max(ONE,fabs(lres));
        if (res_err>=1e-6) {
            std::ostringstream ssV;
            ssV << GlobalVar.Double2String(lpr_obsv) << "  " << GlobalVar.Double2String(lres) << "  " << GlobalVar.Double2String(res_err);
          throw FlxException("flxGPProj::optimize_help_GSL_02", ssV.str());
        }

    // if an error occurs, reset parameters to initial values
    } catch (...) {
      // reset parameter values
        for (tuint i=0;i<Npara;++i) {
          gsl_vector_set (x, i, pvec_ini[i]);
        }
      // re-run evaluation for initial model
        keep_LSE_results = true;
        lres = -gp_likeli_f(x,this);
      // free memory
        gsl_vector_free(x);
        gsl_vector_free(ss);
        gsl_multimin_fminimizer_free (s);
      throw;
    }

    gsl_vector_free(x);
    gsl_vector_free(ss);
    gsl_multimin_fminimizer_free (s);

  #endif

  return lres;
}

const tdouble flxGPProj::optimize(const tuint iterMax, const bool opt_noise)
{
  // reset logging-stream
    para_opt_stream.str("");
    para_opt_stream.clear();
  tdouble lres = log(ZERO);
  tdouble step_size = 0.1;
  #if FLX_USE_NLOPT
    lres = optimize_help(step_size,iterMax,opt_noise,true,para_opt_stream);
  #else
    const tuint itermaxMLEstep = 20;
    tuint i = 0;
    bool was_nan = false;
    while (true) {
      // for the final iteration, return log-likelihood of initial configuration
        if (++i>=itermaxMLEstep) {
          lres = get_log_likeli_obsv();
          if (std::isnan(lres)) {  // try to reset the parameters
            if (was_nan) {
              throw FlxException("flxGPProj::optimize","It seems that the log-likelihood cannot be evaluated.");
              return lres;
            }
            i = 0;
            gp_mean->reset_pVec();
            gp_kernel->reset_pVec();
            was_nan = true;
            continue;
          } else {
            break;
          }
        }
      // for all other iterations, perform an optimization
        try {
          lres = optimize_help(step_size,iterMax,opt_noise,i<=1,para_opt_stream);
        } catch (FlxException_math& e) {
            step_size /= 2; // in case of an error, use only half the step size
            continue;
        }
        break;
    }
  #endif
  return lres;
}

const tdouble flxGPProj::get_log_likeli_obsv()
{
  assemble_observations(false,false);
  if (obsv_changed) return log(ZERO);
  return lpr_obsv;
}

void flxGPProj::eval_covar_point(flxVec& K_star, const flxVec& x_vec, tdouble& prior_mean, tdouble& prior_var)
{
  // check dimension of x_vec
    x_vec.check_size(Ndim);
  // make sure that the model is assembled
    assemble_observations(false,false);
    if (obsv_changed) throw FlxException("flxGPProj::eval_covar_point");
  // parameter vector for covariance function
    flxVec pV(Ndim*2);
    tdouble* pVp = pV.get_tmp_vptr();
    flxVec pV_1(pVp,Ndim);
    flxVec pV_2(pVp+Ndim,Ndim);
    pV_1 = x_vec;
  // loop over data-points
    for (size_t i=0;i<N_obsv;++i) {
      pV_2 = flxVec(dm_ptr+i*Ndim,Ndim);
      K_star[i] = gp_kernel->eval_kernel_f(pVp);
    }
  // evaluate mean at point
    prior_mean = gp_mean->eval_mean_f(pVp);
  // evaluate variance at point
    pV_2 = x_vec;
    prior_var = gp_kernel->eval_kernel_f(pVp);
}

void flxGPProj::predict_mean_var(const flxVec& x_vec, const bool predict_noise, tdouble& res_mean, tdouble& res_var)
{
  x_vec.check_size(Ndim);
  // in case there is a problem
    if (obsv_changed) {
      res_mean = std::numeric_limits<tdouble>::infinity();
      res_var = std::numeric_limits<tdouble>::infinity();
      return;
    }
  // compute covariance between x_vec and input-data
    flxVec K_star(N_obsv);
    tdouble prior_mean, prior_var_noise_free;
    eval_covar_point(K_star,x_vec, prior_mean, prior_var_noise_free);
  // evaluate mean
   res_mean = prior_mean + K_star.operator*(*alpha);
  // evaluate variance
    // Obsolete and WRONG???
      // assemble v-vector
        //lt_mtx_obsv->MultInv(K_star,K_star,true);
        //res_var = prior_var - K_star*K_star;
    // evaluate reduction in variance
      flxVec Rinvr(K_star.get_N());
      lt_mtx_obsv->MultInv(K_star,Rinvr);
      lt_mtx_obsv->TransMultInv(Rinvr,Rinvr);
      // as an alternative to the two lines above (maybe faster, but less accurate)
        // covinv_mtx_obsv->MultMv(K_star,Rinvr);
      tdouble var_red = K_star*Rinvr;
      if (use_LSE) {
        tdouble sigma_Y_2 = prior_var_noise_free + pow2(noise_val);    // = sigma_Y_2
        var_red /= sigma_Y_2;
      }
    // evaluate variance
      res_var = prior_var_noise_free - var_red;
      // avoid round-off errors that lead to small negative variances
        // a problem is that due to round-off errors in obtaining Rinvr, K_star*Rinvr>1 !!!
        if (res_var<ZERO && fabs(res_var)<=sqrt(GlobalVar.TOL())) {
          res_var = ZERO;
        }
        if (res_var<ZERO) {
          throw FlxException_Crude("flxGPProj::predict_mean_var");
        }
      if (predict_noise) {
        res_var += pow2(noise_val);
      }

    // if (use_LSE) {
    //   // evaluate u-vector
    //     flxVec Rinvr(K_star.get_N());
    //     lt_mtx_obsv->MultInv(K_star,Rinvr);
    //     lt_mtx_obsv->TransMultInv(Rinvr,Rinvr);
    //     // as an alternative to the two lines above (maybe faster, but less accurate)
    //       // covinv_mtx_obsv->MultMv(K_star,Rinvr);
    //     // evaluate variance -- without impact of gp_mean
    //       res_var = ONE - K_star*Rinvr;
    //     // correct for impact of gp_mean
    //     if (gp_mean->get_Npara()>0) {
    //       flxVec uvec(Fmtx->ncols());
    //       Fmtx->TransposeMmultVec(Rinvr,uvec);
    //       flxVec fvec(Fmtx->ncols());
    //       gp_mean->assemble_f_vec(fvec,x_vec.get_tmp_vptr_const());
    //       uvec -= fvec;
    //       // evaluate last term
    //         flxVec uhelp(uvec);
    //         FRinvF_cdc_mtx->MultInv(uhelp,uhelp);
    //         FRinvF_cdc_mtx->TransMultInv(uhelp,uhelp);
    //       // finally, correct variance
    //         res_var += uvec*uhelp;
    //     }
    //     // avoid round-off errors that lead to small negative variances
    //       // a problem is that due to round-off errors in obtaining Rinvr, K_star*Rinvr>1 !!!
    //       if (res_var<ZERO && fabs(res_var)<=sqrt(GlobalVar.TOL())) {
    //         res_var = ZERO;
    //       }
    //     res_var *= prior_var;
    // } else {
    //   throw FlxException_NotImplemented("flxGPProj::predict_mean_var");
    // }
}

const tdouble flxGPProj::eval_trend(const flxVec& x_vec, const bool predict_noise)
{
  x_vec.check_size(Ndim);
  tdouble prior_mean = gp_mean->eval_mean_f(x_vec.get_tmp_vptr_const());
  return prior_mean;
}

const tdouble flxGPProj::eval_kernel(const flxVec& x_vec, const bool predict_noise)
{
  x_vec.check_size(Ndim*2);
  tdouble prior_cov = gp_kernel->eval_kernel_f(x_vec.get_tmp_vptr_const());
  if (predict_noise) {
    prior_cov += pow2(noise_val);
  }
  return prior_cov;
}

py::dict flxGPProj::info()
{
  py::dict res;
  res["type"] = "singlegp";
  res["name"] = name;
  // parameters of mean
    res["mean"] = gp_mean->info();
  // parameters of kernel
    res["kernel"] = gp_kernel->info();
  res["noise"] = noise_val;
  res["noise_log"] = noise_opt_stream.str();
  res["opt_log"] = para_opt_stream.str();
  if (d_info) {
    res["logl_obsv"] = get_log_likeli_obsv();
  }
  res["N_obsv"] = N_obsv;
  res["obsv_up2date"] = !obsv_changed;
  return res;
}

void flxGPProj::initialize_from_template_trend(const flxGPProj& ref)
{
  const flxVec& vt_trend = ref.gp_mean->get_paraVec();
  gp_mean->initialize_from_template(vt_trend);
}

void flxGPProj::initialize_from_template_kernel(const flxGPProj& ref)
{
  const flxVec* vt_kernel_ptr = ref.gp_kernel->get_paraVec();
  if (vt_kernel_ptr) {
    gp_kernel->initialize_from_template(*vt_kernel_ptr);
  }
}

const tuint flxGPProj::get_AIC_penalty() const
{
  const tdouble Npara = gp_mean->get_Npara() + gp_kernel->get_Npara() + ONE;  // + ONE is for noise
  return Npara + Npara*(Npara+ONE)/(N_obsv-Npara-ONE);
}
// ------------------------------------------------------------------------------------------------

flxGP_avgModel::flxGP_avgModel(const std::string& name, const tuint Ndim, const tuint iterMax, const bool use_LSE)
: flxGPProj_base(name,Ndim), max_poM(3), max_NKernel(2), model_lst(max_poM*tuint(pow(max_NKernel,Ndim)),NULL), modProbVec(model_lst.size()), iterMax(iterMax)
{
  flxGPProj* model_ptr = NULL;
  flxGP_mean_universal* mean_ptr = NULL;
  flxGP_kernel_auto* kernel_ptr = NULL;
  iVec tVec(Ndim);
  try {
    tuint c = 0;
    for (tuint type_=0;type_<max_poM;++type_) {
      tVec = 0; // initialize kernel switch
      do {
        std::string lbl = name + "_" + std::to_string(type_) + "_";
        for (tuint i=0;i<Ndim;++i) {
          lbl += std::to_string(tVec[i]);
        }
        mean_ptr = new flxGP_mean_universal(type_, Ndim);
        kernel_ptr = new flxGP_kernel_auto(tVec);
        model_ptr = new flxGPProj(lbl,Ndim,mean_ptr,kernel_ptr, use_LSE);
          mean_ptr = NULL;
          kernel_ptr = NULL;
        model_lst[c++] = model_ptr;
          model_ptr = NULL;
      } while (increase_kernel_switch(tVec));
    }
    if (c!=model_lst.size()) throw FlxException_Crude("flxGP_avgModel::flxGP_avgModel");
  } catch (...) {
    if (model_ptr) {
      delete model_ptr;
    } else {
      if (mean_ptr) delete mean_ptr;
      if (kernel_ptr) delete kernel_ptr;
    }
    for (auto const& model_ptr : model_lst) {
      if (model_ptr) delete model_ptr;
    }
    throw;
  }
  // initialize vector of model probabilities
    modProbVec = ONE/model_lst.size();
}

flxGP_avgModel::~flxGP_avgModel()
{
  for (auto const& model_ptr : model_lst) {
    delete model_ptr;
  }
}

const bool flxGP_avgModel::increase_kernel_switch(iVec& tVec) const
{
  const tuint Ndim = tVec.size();
  // sum over tVec
    tuint k = 0;
    for (tuint i=0;i<Ndim;++i) {
      k += tVec[i];
    }
    if (k >= Ndim*(max_NKernel-1)) return false;
  // increase first element
    tVec[0] += 1;
  // make sure that all elements are smaller than N
    for (tuint i=0;i<Ndim;++i) {
      if (tVec[i]>=max_NKernel) {
        tVec[i] = 0;
        #if FLX_DEBUG
          if (i+1==Ndim) throw FlxException_Crude("flxGP_avgModel::increase_kernel_switch");
        #endif
        tVec[i+1] += 1;
      } else {
        break;
      }
    }
  return true;
}

const tdouble flxGP_avgModel::register_observation(const flxGP_data_base& d_info_, const bool initialize_pVec, const bool optimize_noise_val)
{
  // reset logging-stream
    para_opt_stream.str("");
    para_opt_stream.clear();
  // register the observation
    for (auto const& model_ptr : model_lst) {
      model_ptr->register_observation(d_info_,initialize_pVec, optimize_noise_val);
    }
  // optimize model parameters
    iVec tVec(Ndim);
    tuint c = 0;
    flxGPProj* ref_model_trend = NULL;
    const tuint N_kernel_variants = tuint(pow(max_NKernel,Ndim));
    for (tuint type_=0;type_<max_poM;++type_) {
      tVec = 0; // initialize kernel switch
      tuint c_inner = 0;
      do {
        flxGPProj* model_ptr = model_lst[c];
        para_opt_stream << "Optimize parameters of model " << model_ptr->get_name() << " ..." << std::endl;
        if (initialize_pVec) {
        // use results from previous optimizations as starting solution
          // trend
            if (ref_model_trend) {
              model_ptr->initialize_from_template_trend(*ref_model_trend);
            }
            if (c_inner++ == 0) {
              ref_model_trend = model_ptr;
            }
          // kernel

        // use results from previous optimizations as starting solution
          flxGPProj* ref_model_kernel = NULL;
          if (type_==0) {
            if (c>0) {
              ref_model_kernel = model_lst[0];
            }
          } else {
            ref_model_kernel = model_lst[c-N_kernel_variants];
          }
          if (ref_model_kernel) {
            model_ptr->initialize_from_template_kernel(*ref_model_kernel);
          }
        }
        // optimize model
          const tdouble lres = model_ptr->optimize(iterMax,false);
        // evaluate model probability
          const tdouble lpenalty = model_ptr->get_AIC_penalty();
          modProbVec[c++] = lres - lpenalty;
        para_opt_stream << "    log-likelihood: " << GlobalVar.Double2String(lres) << " penalty: " << lpenalty  << std::endl;
      } while (increase_kernel_switch(tVec));
    }
  // evaluate model probabilities
    const tdouble AICc_max = modProbVec.get_max();
    modProbVec -= AICc_max;
    for (tuint i=0;i<modProbVec.get_N();++i) {
      modProbVec[i] = exp(modProbVec[i]);
    }
    modProbVec /= modProbVec.get_sum();
    para_opt_stream << "Model probabilities:" << std::endl;
    for (tuint i=0;i<modProbVec.get_N();++i) {
      para_opt_stream << "     " << model_lst[i]->get_name() << ": " << GlobalVar.Double2String(modProbVec[i]) << std::endl;
    }
    return modProbVec.get_Mean();
}

void flxGP_avgModel::register_noise(const tdouble noise_sd)
{
  for (auto const& model_ptr : model_lst) {
    model_ptr->register_noise(noise_sd);
  }
}

void flxGP_avgModel::unassemble()
{
  for (auto const& model_ptr : model_lst) {
    model_ptr->unassemble();
  }
}

const tdouble flxGP_avgModel::optimize(const tuint iterMax, const bool opt_noise)
{
  return log(ZERO);
}

const tdouble flxGP_avgModel::get_log_likeli_obsv()
{
  throw FlxException_NotImplemented("flxGP_avgModel::get_log_likeli_obsv");
}

void flxGP_avgModel::predict_mean_var(const flxVec& x_vec, const bool predict_noise, tdouble& res_mean, tdouble& res_var)
{
  res_mean = ZERO;
  res_var = ZERO;
  for (tuint i=0;i<modProbVec.get_N();++i) {
    const tdouble mp = modProbVec[i];
    if (mp>GlobalVar.TOL()) {
      tdouble tmp_mean, tmp_var;
      model_lst[i]->predict_mean_var(x_vec,predict_noise,tmp_mean,tmp_var);
      res_mean += mp * tmp_mean;
      res_var += pow2(mp) * tmp_var;
    }
  }
}

const tdouble flxGP_avgModel::eval_trend(const flxVec& x_vec, const bool predict_noise)
{
  tdouble res = ZERO;
  for (tuint i=0;i<modProbVec.get_N();++i) {
    const tdouble mp = modProbVec[i];
    if (mp>GlobalVar.TOL()) {
      res += mp * model_lst[i]->eval_trend(x_vec,predict_noise);
    }
  }
  return res;
}

const tdouble flxGP_avgModel::eval_kernel(const flxVec& x_vec, const bool predict_noise)
{
  tdouble res = ZERO;
  for (tuint i=0;i<modProbVec.get_N();++i) {
    const tdouble mp = modProbVec[i];
    if (mp>GlobalVar.TOL()) {
      res += pow2(mp) * model_lst[i]->eval_kernel(x_vec,predict_noise);
    }
  }
  return res;
}

py::dict flxGP_avgModel::info()
{
  py::dict res;
  res["type"] = "mavggp";
  res["name"] = name;
  // TODO return other information
  res["opt_log"] = para_opt_stream.str();
  return res;
}


// ------------------------------------------------------------------------------------------------


flxGPBox::flxGPBox()
{
  #if FLX_DEBUG
    ++Cinst;
    if (Cinst > 1) {
      std::ostringstream ssV;
      ssV << "More than one instance of 'flxGPBox' created ...";
      throw FlxException("flxGPBox::flxGPBox", ssV.str() );
    }
  #endif
}

flxGPBox::~flxGPBox()
{
  for (std::map<std::string, flxGPProj_base*>::iterator pos = box.begin(); pos != box.end(); ++pos) {
    delete pos->second;
  }
}

void flxGPBox::insert(const std::string name, const tuint Ndim, flxGP_mean_base* gp_mean, flxGP_kernel_base* gp_kernel, const bool use_LSE)
{
  std::map<std::string, flxGPProj_base*>::iterator pos = box.find(name);
  if ( pos != box.end() ) {
    std::ostringstream ssV;
    ssV << "A gp-project with the name '" << name << "' does already exist.";
    throw FlxException("flxGPBox::insert_1", "Project already exists", ssV.str() ); 
  } else {
    flxGPProj* valuep = new flxGPProj(name,Ndim,gp_mean,gp_kernel, use_LSE);
    std::pair<std::string, flxGPProj_base*> Element(name, valuep);
    box.insert(Element);
  }
}

void flxGPBox::insert (const std::string name, flxGPProj_base* gp_model)
{
  std::map<std::string, flxGPProj_base*>::iterator pos = box.find(name);
  if ( pos != box.end() ) {
    std::ostringstream ssV;
    ssV << "A gp-project with the name '" << name << "' does already exist.";
    throw FlxException("flxGPBox::insert_2", "Project already exists", ssV.str() );
  } else {
    std::pair<std::string, flxGPProj_base*> Element(name, gp_model);
    box.insert(Element);
  }
}

flxGPProj_base& flxGPBox::get(const std::string& name)
{
  std::map<std::string, flxGPProj_base*>::iterator pos = box.find(name);
  if ( pos != box.end() ) {
    return (*pos->second);
  } else {
    std::ostringstream ssV_2;
    ssV_2 << "The GP-project '" << name << "' does not exist.";
    throw FlxException("flxGPBox::get_1", ssV_2.str(), "You have to define the project before you can use it."); 
  }
}





