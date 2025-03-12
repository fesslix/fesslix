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


#include <cmath>

#include "flxBayDA.h"
#include "flxfunction_ope_calc.h"
#include "flxfunction_fun_calc.h"
#include "flxfunction.h"

#include <gsl/gsl_multimin.h>
FlxBayUpBox* RBRV_entry_read_bayDA::BayUpBox = NULL;

// the number of distribution types supported by bayDA
const tuint flxBayDA_likeli_N_type = 12;


flxBayDA_likeli::flxBayDA_likeli(tuint type, flxVec& data_vec_, const tuint id_transform, FlxRndCreator& RndCreator, tuint Nchain, tuint Npost, tuint N_adapt)
: type(type), M(0), pvec(NULL), rv(NULL), rv_ext(NULL), data_vec(&data_vec_), id_transform(id_transform), id_transform2(0), tcl(ZERO), RndCreator(RndCreator), Nchain(Nchain), Npost(Npost), pchain(NULL), ppost(NULL), pstep(NULL), lcvec(NULL), N(0), N_adapt(N_adapt), N_adapt_(0),
  N_adapt_count(NULL), dim_sample(NULL),chain_stats(NULL)   // ,lsum(Npost,false)
{
    // transformation related
    switch (id_transform) {
        case 0:
            break;
        case 1:
            tcl -= data_vec->get_sum();
            break;
        default:
            throw FlxException_Crude("flxBayDA::flxBayDA_01");
    }
    // type-specific allocation
    switch (type) {
        case 0:     // Normal
        {
            // define internal random variable
            M = 2;
            pvec = new flxVec(M);
            FlxFunction* mu = gen_para_fun(0,0);
            FlxFunction* sd = gen_para_fun(1,1);
            rv = new RBRV_entry_RV_normal("BayDA_Normal",0,0,mu,sd,NULL,NULL,false);
            (*pvec)[0] = data_vec->get_Mean();
            (*pvec)[1] = log(data_vec->get_sd((*pvec)[0]));
            apply_transform2(0);
            break;
        }
        case 1:     // Gumbel
        {
            // define internal random variable
            M = 2;
            pvec = new flxVec(M);
            FlxFunction* mu = gen_para_fun(0,0);
            FlxFunction* sd = gen_para_fun(1,1);
            rv = new RBRV_entry_RV_Gumbel("BayDA_Gumbel",0,1,mu,sd,NULL,NULL,false);
            (*pvec)[0] = data_vec->get_Mean();
            (*pvec)[1] = log(data_vec->get_sd((*pvec)[0]));
            apply_transform2(0);
            break;
        }
        case 2:     // Laplace
        {
            // define internal random variable
            M = 2;
            pvec = new flxVec(M);
            FlxFunction* loc = gen_para_fun(0,0);
            FlxFunction* scale = gen_para_fun(1,1);
            rv = new RBRV_entry_RV_Laplace("BayDA_Laplace",0,loc,scale);
            (*pvec)[0] = data_vec->get_Mean();
            (*pvec)[1] = log(data_vec->get_sd((*pvec)[0])/sqrt(2.));
            apply_transform2(0);
            break;
        }
        case 3:     // Cauchy
        {
            // define internal random variable
            M = 2;
            pvec = new flxVec(M);
            FlxFunction* loc = gen_para_fun(0,0);
            FlxFunction* scale = gen_para_fun(1,1);
            rv = new RBRV_entry_RV_Cauchy("BayDA_Cauchy",0,loc,scale);
            flxVec tmpv(*data_vec);
            tmpv.sort();
            (*pvec)[0] = flx_percentile(tmpv.get_tmp_vptr(),tmpv.get_N(),0.5);
            (*pvec)[1] = log((flx_percentile(tmpv.get_tmp_vptr(),tmpv.get_N(),0.75)-flx_percentile(tmpv.get_tmp_vptr(),tmpv.get_N(),0.25))/2);
            apply_transform2(0);
            break;
        }
        case 4:     // Student's t
        {
            // define internal random variable
            M = 3;
            pvec = new flxVec(M);
            FlxFunction* nu = gen_para_fun(1,0);
            FlxFunction* loc = gen_para_fun(0,1);
            FlxFunction* scale = gen_para_fun(1,2);
            rv = new RBRV_entry_RV_StudentsT_generalized("BayDA_StudentstGen",0,nu,loc,scale);
            (*pvec)[0] = log(2.);
            (*pvec)[1] = data_vec->get_Mean();
            (*pvec)[2] = log(data_vec->get_sd((*pvec)[1]));
            apply_transform2(0);
            break;
        }
        case 5:     // truncated Normal on exp-transformed space
        {
            // define internal random variable
            M = 2;
            pvec = new flxVec(M);
            FlxFunction* mu = gen_para_fun(0,0);
            FlxFunction* sd = gen_para_fun(1,1);
            FlxFunction* a = new FlxFunction(new FunNumber(ZERO));
            rv = new RBRV_entry_RV_normal_trunc("BayDA_exp_truncNormal",0,mu,sd,a,NULL,false);
            apply_transform2(1);
            (*pvec)[0] = data_vec->get_Mean();
            (*pvec)[1] = log(data_vec->get_sd((*pvec)[0]));
            break;
        }
        case 6:     // Gamma on exp-transformed space
        {
            // define internal random variable
            M = 2;
            pvec = new flxVec(M);
            FlxFunction* mu = gen_para_fun(1,0);
            FlxFunction* sd = gen_para_fun(1,1);
            rv = new RBRV_entry_RV_gamma("BayDA_exp_Gamma",0,true,mu,sd,NULL,false);
            apply_transform2(1);
            const tdouble tm = data_vec->get_Mean();
            (*pvec)[0] = log(tm);
            (*pvec)[1] = log(data_vec->get_sd(tm));
            break;
        }
        case 7:     // exponential on exp-transformed space
        {
            // define internal random variable
            M = 1;
            pvec = new flxVec(M);
            FlxFunction* lambda = gen_para_fun(1,0);
            rv = new RBRV_entry_RV_exponential("BayDA_exp_exponential",0,lambda,NULL);
            apply_transform2(1);
            (*pvec)[0] = log(ONE/data_vec->get_Mean());
            break;
        }
        case 8:     // Weibull on exp-transformed space
        {
            // define internal random variable
            M = 2;
            pvec = new flxVec(M);
            FlxFunction* mu = gen_para_fun(1,0);
            FlxFunction* sd = gen_para_fun(1,1);
            rv = new RBRV_entry_RV_Weibull("BayDA_exp_Weibull",0,true,mu,sd,NULL,false);
            apply_transform2(1);
            const tdouble tm = data_vec->get_Mean();
            (*pvec)[0] = log(tm);
            (*pvec)[1] = log(data_vec->get_sd(tm));
            break;
        }
        case 9:     // Gumbel on exp-transformed space
        {
            // define internal random variable
            M = 2;
            pvec = new flxVec(M);
            FlxFunction* mu = gen_para_fun(0,0);
            FlxFunction* sd = gen_para_fun(1,1);
            FlxFunction* a = new FlxFunction(new FunNumber(ZERO));
            rv = new RBRV_entry_RV_Truncated("BayDA_Gumbel_trunc",0, a,NULL, new RBRV_entry_RV_Gumbel("BayDA_Gumbel",0,1,mu,sd,NULL,NULL,false));
            apply_transform2(1);
            (*pvec)[0] = data_vec->get_Mean();
            (*pvec)[1] = log(data_vec->get_sd((*pvec)[0]));
            break;
        }
        case 10:     // Laplace on exp-transformed space
        {
            // define internal random variable
            M = 2;
            pvec = new flxVec(M);
            FlxFunction* loc = gen_para_fun(0,0);
            FlxFunction* scale = gen_para_fun(1,1);
            FlxFunction* a = new FlxFunction(new FunNumber(ZERO));
            rv = new RBRV_entry_RV_Truncated("BayDA_Laplace_trunc",0, a,NULL, new RBRV_entry_RV_Laplace("BayDA_Laplace",0,loc,scale));
            apply_transform2(1);
            (*pvec)[0] = data_vec->get_Mean();
            (*pvec)[1] = log(data_vec->get_sd((*pvec)[0])/sqrt(2.));
            break;
        }
        case 11:     // DUMMY - different variant of Normal distribution
        {
            // define internal random variable
            M = 2;
            pvec = new flxVec(M);
            FlxFunction* mu = gen_para_fun(0,0);
            FlxFunction* sd = gen_para_fun(1,1);
            rv = new RBRV_entry_RV_lognormal("BayDA_LogNormal",0,0,mu,sd,NULL,NULL,NULL,false);
            (*pvec)[0] = data_vec->get_Mean();
            (*pvec)[1] = log(data_vec->get_sd((*pvec)[0]));
            apply_transform2(1);
            break;
        }
        default:
            throw FlxException_Crude("flxBayDA::flxBayDA_02");
    }
    // allocate memory for chains
    pchain = new flxVec(M*Nchain);
    // allocate memory for posterior samples
    ppost = new flxVec(M*Npost);
    // allocate step vector
    pstep = new flxVec(M);
    // allocate dim_sample
    dim_sample = new tuint[M];
    // allocate lcvec
    lcvec = new flxVec(Nchain);
    // allocate N_adapt_count
    N_adapt_count = new tuint[M*2];
    for (tuint i=0;i<2*M;++i) {
        N_adapt_count[i] = 0;
    }
    // allocate ashuffle
    ashuffle = new tuint[Npost];
    RndCreator.shuffle(ashuffle,Npost);
    // allocate chain_stats
    chain_stats = new vdouble[Nchain*M];
}

flxBayDA_likeli::~flxBayDA_likeli() {
    if (pvec) { delete pvec; }
    if (rv) { delete rv; }
    if (pchain) { delete pchain; }
    if (ppost) { delete ppost; }
    if (pstep) { delete pstep; }
    if (dim_sample) { delete [] dim_sample; }
    if (N_adapt_count) { delete [] N_adapt_count; }
    if (lcvec) { delete lcvec; }
    if (ashuffle) { delete [] ashuffle; }
    if (id_transform2!=0) {
        delete data_vec;
        delete rv_ext;
    }
    if (chain_stats) delete [] chain_stats;
}

void flxBayDA_likeli::apply_transform2(const tuint idt)
{
    id_transform2 = idt;
    switch (idt) {
        case 0:
            rv_ext = rv;
            break;
        case 1:
        {
            // define external random variable
            FlxFunction* z2x = NULL;
            FlxFunction* x2z = NULL;
            FlxFunction* dxdz = NULL;
            tdouble s = data_vec->get_Mean();
            if (s <= 10.) {
                s = ONE;
                z2x = new FlxFunction(new FunLn(new FunPara(1)));
                x2z = new FlxFunction(new FunExp(new FunPara(1)));
                dxdz = new FlxFunction(new FunExp(new FunPara(1)));
            } else {
                z2x = new FlxFunction(new FunMult(new FunNumber(s), new FunLn(new FunPara(1))));
                x2z = new FlxFunction(new FunExp(new FunMult_Div(new FunPara(1),new FunNumber(s))));
                dxdz = new FlxFunction(new FunMult(new FunExp(new FunMult_Div(new FunPara(1),new FunNumber(s))),new FunNumber(ONE/s)));
            }
            rv_ext = new RBRV_entry_RV_UserTransform("BayDA_LogNormal_ext",0,true,z2x,x2z,dxdz,NULL,rv,false);
            // transform data-vector
            const tuint N = data_vec->get_N();
            tcl += data_vec->get_sum()/s - log(s)*N;
            data_vec = new flxVec(*data_vec);
            tdouble* pdv = data_vec->get_tmp_vptr();
            for (tuint i=0;i<N;++i) {
                pdv[i] = exp(pdv[i]/s);
            }
            break;
        }
        default:
            throw FlxException_Crude("flxBayDA::apply_transform2");
    }
}

FlxFunction* flxBayDA_likeli::gen_para_fun(const tuint ftype, const tuint pid)
{
    switch (ftype) {
        case 0:
            return new FlxFunction(new FunConst(pvec->get_tmp_vptr()+pid));
        case 1:
            return new FlxFunction(new FunAdd(new FunExp(new FunConst(pvec->get_tmp_vptr()+pid)),new FunNumber(GlobalVar.TOL())));
        default:
            throw FlxException_Crude("flxBayDA_likeli::gen_para_fun");
    }
}

const tdouble flxBayDA_likeli::calc_likeli()
{
    pdouble res;
    const tuint N = data_vec->get_N();
    const tdouble* pdv = data_vec->get_tmp_vptr_const();
    try {
        for (tuint i=0;i<N;++i) {
            const tdouble r1 =  rv->calc_pdf_x_log(pdv[i],true);
            res += r1;
        }
    } catch (FlxException &e) {
        GlobalVar.alert.alert("flxBayDA_likeli::calc_likeli_01",e.what());
        res += log(ZERO);
    }
    res += tcl;
    const tdouble rd = res.cast2double();
    if (!std::isfinite(rd)) {
        if (rd!=log(ZERO)) {
            #if FLX_DEBUG
                for (tuint i=0;i<N;++i) {
                    const tdouble r1 =  rv->calc_pdf_x_log(pdv[i],true);
                    if (!std::isfinite(r1)) {
                        if (r1!=log(ZERO)) {
                            const tdouble r2 =  rv->calc_pdf_x_log(pdv[i],true);
                            std::ostringstream ssV;
                            ssV << "Likelihood is not finite. term=" << N << ", val=" << r2;
                            throw FlxException("flxBayDA_likeli::calc_likeli_02", ssV.str());
                        }
                    }
                }
            #endif
            std::ostringstream ssV;
            ssV << "Likelihood is not finite. (" << GlobalVar.Double2String(rd) << ")";
            throw FlxException("flxBayDA_likeli::calc_likeli_03", ssV.str());
        }
    }
    return rd;
}

void flxBayDA_likeli::initialize_chains()
{
    // step position of chains (to pvec)
    for (tuint i=0;i<Nchain;++i) {
        flxVec ptmp(pchain->get_tmp_vptr()+i*M,M);
        ptmp = *pvec;
    }
    // set initial step size
    for (tuint m=0;m<M;++m) {
        (*pstep)[m] = fabs((*pvec)[m])*0.1;
    }
}

void flxBayDA_likeli::move_chains(const bool is_tuning_phase)
{
    for (tuint i=0;i<Nchain;++i) {
        // position of current sample
        flxVec x(pchain->get_tmp_vptr()+i*M,M);
        *pvec = x;
        // sample dimension vector
        RndCreator.shuffle(dim_sample,M);
        for (tuint m=0;m<M;++m) {
            const tuint m_ = dim_sample[m];
            const tdouble s_ = (*pstep)[m_];
            // find target likelihood value
                tdouble aa;
                do {        // it did happen that aa equals ZERO ... which causes lval equal to -inf
                    aa = RndCreator.gen_smp_uniform();
                } while (aa<=GlobalVar.TOL());
                const tdouble lval = log(aa)+calc_likeli();
            if (!std::isfinite(lval)) {
                throw FlxException_Crude("flxBayDA_likeli::move_chains_01");
            }
            // extraction - lower bound
            tdouble x_low = x[m_];
            tuint c=0;
            tdouble s__ = s_;
            tdouble l_prev = std::numeric_limits<tdouble>::infinity();
            tuint lcc = 0;
            do {
                if (c>10) s__*=2;
                x_low -= s__;
                (*pvec)[m_] = x_low;
                const tdouble l_now = calc_likeli();
                if (fabs(l_now-l_prev)<1e-6) {
                    ++lcc;
                } else {
                    lcc = 0;
                }
                ++c;
                l_prev = l_now;
                // if (c>=10) GlobalVar.slogcout(4) << "     move_chains::x_low=" << x_low << "   c=" << c << "   likeli=" << l_now << "   target=" << lval << "   dim=" << m_ << std::endl;
                if (lcc>4) break;
                if (c>10000) {
                    throw FlxException("flxBayDA_likeli::move_chains_02" );
                }
            } while (l_prev>=lval );
            N_adapt_count[m_*2] += c;
            // extraction - upper bound
            tdouble x_up = x[m_];
            c=0;
            s__ = s_;
            l_prev = std::numeric_limits<tdouble>::infinity();
            lcc = 0;
            do {
                if (c>10) s__*=2;
                x_up += s__;
                (*pvec)[m_] = x_up;
                const tdouble l_now = calc_likeli();
                if (fabs(l_now-l_prev)<1e-6) {
                    ++lcc;
                } else {
                    lcc = 0;
                }
                ++c;
                l_prev = l_now;
                // if (c>=10) GlobalVar.slogcout(4) << "     move_chains::x_up=" << x_up << "   c=" << c << "   likeli=" << calc_likeli() << "   target=" << lval << "   dim=" << m_ << std::endl;
                if (lcc>4) break;
                if (c>10000) {
                    throw FlxException("flxBayDA_likeli::move_chains_03" );
                }
            } while (l_prev>=lval );
            N_adapt_count[m_*2] += c;
            if (std::isnan(x_low) || std::isnan(x_up)) {
                throw FlxException_Crude("flxBayDA_likeli::move_chains_04");
            }
            // sample uniformly and shrink
            c = 0;
            do {
                (*pvec)[m_] = RndCreator.gen_smp_uniform()*(x_up-x_low)+x_low;
                if (std::isnan((*pvec)[m_])) {
                    throw FlxException_Crude("flxBayDA_likeli::move_chains_05");
                }
                ++c;
                (*lcvec)[i] = calc_likeli();
                if ((*lcvec)[i]>=lval) {
                    break;
                }
                if ( (*pvec)[m_]>x[m_] ) {
                    x_up = (*pvec)[m_];
                } else {
                    x_low = (*pvec)[m_];
                }
                if (c>10000) {
                    throw FlxException("flxBayDA_likeli::move_chains_06" );
                }
            } while (true);
            N_adapt_count[m_*2+1] += c;
        }
        // remember current position of chain
        x = *pvec;
    }
    ++N_adapt_;
    if (is_tuning_phase) {
        if (N_adapt_>=N_adapt) {
            for (tuint m=0;m<M;++m) {
                const tdouble Ne = ((double)N_adapt_count[m*2])/2/(N_adapt_*Nchain);
                const tdouble Nc = (double)N_adapt_count[m*2+1]/(N_adapt_*Nchain);
                (*pstep)[m] = 2*((*pstep)[m])*(Ne/(Ne+Nc));
                N_adapt_count[m*2] = 0;
                N_adapt_count[m*2+1] = 0;
            }
            N_adapt_ = 0;
        }
    }
}

void flxBayDA_likeli::fill_post_samples()
{
    if (N>0) { return; }
    while (true) {
        move_chains(false);
        for (tuint i=0;i<Nchain;++i) {
            lsum += (*lcvec)[i];
        }
        for (tuint i=0;i<Nchain;++i) {
            const tuint j = ashuffle[N++];
            flxVec x1(pchain->get_tmp_vptr()+i*M,M);
            flxVec x2(ppost->get_tmp_vptr()+j*M,M);
            x2 = x1;
            // monitor chain convergence
            for (tuint m=0;m<M;++m) {
                chain_stats[i*M+m] += x1[m];
            }
            if (N>=Npost) {
                return;
            }
        }
    }
}

void flxBayDA_likeli::get_post_sample(flxVec* sample_ptr)
{
    // make sure that there are posterior samples left to retrieve
    if (N<=0) {
        fill_post_samples();
    }
    // get the next sample in the shuffled list
    const tuint i = ashuffle[--N];
    *pvec = flxVec(ppost->get_tmp_vptr()+i*M,M);
    if (sample_ptr) {
        *sample_ptr = *pvec;
    }
}

void flxBayDA_likeli::get_posterior_mean(flxVec& mv, flxVec& sv) const {
    const tdouble* pv = ppost->get_tmp_vptr_const();
    for (tuint m=0;m<M;++m) {
        vdouble s;
        for (tuint i=0;i<Npost;++i) {
            s += pv[i*M+m];
        }
        mv[m] = s.get_mean();
        sv[m] = sqrt(s.get_variance());
    }
}

const tdouble flxBayDA_likeli::eval_Gelman_Rubin_convergence(const tuint m) {
    flxVec chain_mean(Nchain);
    flxVec sj(Nchain);
    for (tuint i=0;i<Nchain;++i) {
        chain_mean[i] = chain_stats[i*M+m].get_mean();
        sj[i] = chain_stats[i*M+m].get_variance();
    }
    const tdouble grand_mean = chain_mean.get_Mean();
    const tdouble B = chain_mean.get_Var(grand_mean);
    const tdouble W = sj.get_Mean();
    const tdouble L = chain_stats[m].get_size();
    const tdouble R = ((L-ONE)/L*W+B/L)/W;
    return R;
}

flxBayDA::flxBayDA(const std::string& name, const tuint id_transform, FlxSMtx& smtx, FlxRndCreator& RndCreator, tuint Nchain, tuint Nburn, tuint Ntune, tuint Npost, tuint N_adapt, const tdouble tplaus, std::valarray<int> dtfv, const std::string& pvec_name, const std::string& distid_name)
: name(name), id_transform(id_transform), data_vec(smtx.get_Ncoeff()), RndCreator(RndCreator), Nchain(Nchain), Nburn(Nburn), Ntune(Ntune), Npost(Npost), N_adapt(N_adapt), rv(NULL),
  rv_dummy("bayDA_"+name+"_dummy",0), model_pr(NULL), models(NULL), tplaus(tplaus), dtfv(dtfv), pvec_name(pvec_name), distid_name(distid_name)
{
    if (distid_name!="") {
        data->ConstantBox.declareC(distid_name);
    }
    // assign data_vec
    const tuint N = data_vec.get_N();
    tdouble* dptr = data_vec.get_tmp_vptr();
    tuint c=0;
    for (tuint i=0;i<smtx.get_nrows();++i) {
        for (tuint j=0;j<smtx.get_ncols();++j) {
            dptr[c++] = smtx(i,j);
        }
    }
    if (id_transform==1 && data_vec.get_min()<=ZERO) {
        std::ostringstream ssV;
        ssV << "(" << data_vec.get_min() << "; data: " << data_vec << ")";
        throw FlxException("flxBayDA::gen_samples_02", "Not all entries of the data-vector are positive.", ssV.str() );
    }
    // perform transformation
    FlxFunction* z2x = NULL;
    FlxFunction* x2z = NULL;
    FlxFunction* dxdz = NULL;
    FlxFunction* checkx = NULL;
    try {
        switch (id_transform) {
            case (0):
                z2x = new FlxFunction(new FunPara(1));
                x2z = new FlxFunction(new FunPara(1));
                dxdz = new FlxFunction(new FunNumber(ONE));
                checkx = new FlxFunction(new FunNumber(ONE));
                break;
            case (1):
                for (tuint i=0;i<N;++i) dptr[i] = log(dptr[i]);
                z2x = new FlxFunction(new FunExp(new FunPara(1)));
                x2z = new FlxFunction(new FunLn(new FunPara(1)));
                dxdz = new FlxFunction(new FunMult_Div(new FunNumber(ONE),new FunPara(1)));
                checkx = new FlxFunction(new FunPara(1));
                break;
            default:
                throw FlxException_Crude("flxBayDA::flxBayDA_01");
        };
        // prepare chains
        if (Ntune>=Nburn) {
            throw FlxException("flxBayDA::flxBayDA_02", "Ntune must be smaller than Ntune" );
        }
        // create associated random variable
        rv = new RBRV_entry_RV_UserTransform(name,0,true,z2x,x2z,dxdz,checkx,&rv_dummy,false);
    } catch (FlxException& e) {
        if (z2x) delete z2x;
        if (x2z) delete x2z;
        if (dxdz) delete dxdz;
        if (checkx) delete checkx;
        throw;
    }
}

flxBayDA::~flxBayDA() {
    if (rv) delete rv;
    if (model_pr) delete model_pr;
    free_models();
}

void flxBayDA::free_models()
{
    if (models) {
        for (tuint i=0; i<flxBayDA_likeli_N_type;++i) {
            if (models[i]) delete models[i];
        }
        delete [] models;
        models = NULL;
    }
}

void flxBayDA::gen_samples()
{
    if (model_pr) return;
    if (models) throw FlxException_Crude("flxBayDA::gen_samples_01");
    // posterior sampling conditional on the different distribution types
        GlobalVar.slogcout(4) << "BayDA (" << name << ")" << std::endl;
        try {
            models = new flxBayDA_likeli*[flxBayDA_likeli_N_type];
            for (tuint type=0;type<flxBayDA_likeli_N_type;++type) models[type] = NULL;
            tdouble l_best = log(ZERO);
            for (tuint type=0;type<flxBayDA_likeli_N_type;++type) {
                // check if type is to be processed
                    bool skip = true;
                    for(auto it = std::begin(dtfv); it != std::end(dtfv); ++it) {
                        if (*it<0 || *it==(int)type) {
                            skip = false;
                            break;
                        }
                    }
                    if (skip) continue;
                GlobalVar.slogcout(4) << "  inference for likelihood TYPE " << type << std::endl;
                models[type] = new flxBayDA_likeli(type,data_vec,id_transform,RndCreator,Nchain,Npost,N_adapt);
                flxBayDA_likeli& likeli = *(models[type]);

                // find MLE of current model
                tdouble lres = log(ZERO);
                tdouble step_size_MLE = ONE;
                const tuint itermaxMLEstep = 20;
                for (tuint i=0;i<itermaxMLEstep;++i) {
                    try {
                        lres = find_MLE(likeli,step_size_MLE,i==itermaxMLEstep-1);
                    } catch (FlxException_math& e) {
                        step_size_MLE /= 2;
                        continue;
                    }
                    break;
                }
                if (lres<l_best-tplaus) {
                    GlobalVar.slogcout(4) << "      UNPLAUSIBLE model" << std::endl;
                    delete models[type];
                    models[type] = NULL;
                    continue;
                }
                likeli.initialize_chains();

                // tuning phase
                GlobalVar.slogcout(4) << "      tuning phase ..." << std::endl;
                for (tuint i=0;i<Ntune;++i) {
                    likeli.move_chains(true);
                }
                // finish burn-in phase
                GlobalVar.slogcout(4) << "      burn-in phase ..." << std::endl;
                for (tuint i=0;i<(Nburn-Ntune);++i) {
                    likeli.move_chains(false);
                }
                // generate and store posterior samples
                GlobalVar.slogcout(4) << "      sampling phase ..." << std::endl;
                likeli.fill_post_samples();
                for (tuint m=0;m<likeli.get_M();++m) {
                    GlobalVar.slogcout(4) << "        R("<< (m+1) << ") = " << likeli.eval_Gelman_Rubin_convergence(m) << std::endl;
                }
                lres = likeli.get_fit();
                GlobalVar.slogcout(4) << "    data-fit: " << lres << " [sigma " << GlobalVar.Double2String(likeli.get_fit_sd(),false,2) << "] at ( ";
                flxVec mv(likeli.get_M());
                flxVec sv(likeli.get_M());
                likeli.get_posterior_mean(mv,sv);
                for (tuint m=0;m<likeli.get_M();++m) {
                    if (m>0) GlobalVar.slogcout(4) << ", ";
                    GlobalVar.slogcout(4) << mv[m] << " [" << GlobalVar.Double2String(sv[m],false,2) << "]";
                }
                GlobalVar.slogcout(4) << " )" << std::endl;
                if (lres<l_best-tplaus || std::isnan(lres)) {
                    GlobalVar.slogcout(4) << "      UNPLAUSIBLE model" << std::endl;
                    delete models[type];
                    models[type] = NULL;
                    continue;
                }
                if (lres>l_best) l_best = lres;
            }
        } catch (FlxException& e) {
            free_models();
            throw;
        }
    // evaluate model probabilities
        model_pr = new flxVec(flxBayDA_likeli_N_type);
        for (tuint type=0;type<flxBayDA_likeli_N_type;++type) {
            if (models[type]==NULL) {
                (*model_pr)[type] = log(ZERO);
            } else {
                (*model_pr)[type] = models[type]->get_fit();
            }
        }
        *model_pr -= model_pr->get_max();
        for (tuint type=0;type<flxBayDA_likeli_N_type;++type) {
            (*model_pr)[type] = exp((*model_pr)[type]);
        }
        *model_pr /= model_pr->get_sum();
        GlobalVar.slogcout(4) << "    model probabilities:" << std::endl;
        for (tuint type=0;type<flxBayDA_likeli_N_type;++type) {
            GlobalVar.slogcout(4) << "      TYPE " << type << ": " << (*model_pr)[type] << std::endl;;
        }
        tdouble s = ZERO;
        for (tuint type=0;type<flxBayDA_likeli_N_type;++type) {
            s += (*model_pr)[type];
            (*model_pr)[type] = s;
        }
}

void flxBayDA::sample()
{
    const tuint type = RndCreator.gen_smp_index(*model_pr);
    GlobalVar.slogcout(4) << "BayDA (" << name << ") :: sample :: TYPE " << type << " :: PARA ( ";
    const tuint M = models[type]->get_M();
    flxVec pvec(M);
    models[type]->get_post_sample(&pvec);
    for (tuint m=0;m<M;++m) {
        if (m>0) {
            GlobalVar.slogcout(4) << ", ";
        }
        GlobalVar.slogcout(4) << pvec[m];
    }
    GlobalVar.slogcout(4) << " )" << std::endl;
    rv->replace_rv_z(models[type]->get_rv());
    // store current outcome of parameter vector
        if (pvec_name!="") {
            FlxSMtx* s_pvec = new FlxSMtx(pvec);
            data->ConstMtxBox.insert(pvec_name,s_pvec);
        }
    // store current distribution type
        if (distid_name!="") {
            data->ConstantBox.getRef(distid_name) = (tdouble)type;
        }
}

double likeli_f (const gsl_vector *v, void *params)
{
  flxBayDA_likeli *p = (flxBayDA_likeli *)params;

  const tuint M = p->get_M();
  flxVec& pvec = p->get_pvec();

  for (tuint i=0;i<M;++i) {
      pvec[i] = gsl_vector_get(v, i);
  }
  const tdouble res = p->calc_likeli();
  if (std::isnan(res)) {
      throw FlxException_Crude("flxBayDA::likeli_f_01");
  }
  if (std::isinf(res)) {
      throw FlxException_math("flxBayDA::likeli_f_02");
  }
  return -res;
}

const tdouble flxBayDA::find_MLE(flxBayDA_likeli& likeli, const tdouble step_size, const bool return_ini)
{
  const gsl_multimin_fminimizer_type *T = gsl_multimin_fminimizer_nmsimplex2;
  gsl_multimin_fminimizer *s = NULL;
  gsl_vector *ss, *x;
  gsl_multimin_function minex_func;

  size_t iter = 0;
  int status;
  double size;
  const tuint M = likeli.get_M();
  tdouble lres;

  /* Starting point */
    x = gsl_vector_alloc(M);
    for (tuint i=0;i<M;++i) {
        gsl_vector_set (x, i, likeli.get_pvec()[i]);
    }

  /* Set initial step sizes to 1 */
  ss = gsl_vector_alloc(M);
  gsl_vector_set_all(ss, step_size);

  /* Initialize method and iterate */
  minex_func.n = M;
  minex_func.f = likeli_f;
  minex_func.params = &likeli;

  const flxVec pvec_ini = likeli.get_pvec();
  try {
    // output inital point
        lres = -likeli_f(x,&likeli);
        if (return_ini) {
            gsl_vector_free(x);
            gsl_vector_free(ss);
            return lres;
        }
        if (step_size==ONE) {
            GlobalVar.slogcout(4) << "    initial point estimate: " << lres << " at ( ";
            for (tuint i=0;i<pvec_ini.get_N(); ++i) {
                if (i>0) {
                    GlobalVar.slogcout(4) << ", ";
                }
                GlobalVar.slogcout(4) << pvec_ini[i];
            }
            GlobalVar.slogcout(4) << " ) " << std::endl;
        }

    s = gsl_multimin_fminimizer_alloc (T, M);
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
            lres = -likeli_f(s->x,&likeli);
            break;
        }

    //       printf ("%5d %10.3e %10.3e f() = %7.3f size = %.3f\n",
    //               (int)iter,
    //               gsl_vector_get (s->x, 0),
    //               gsl_vector_get (s->x, 1),
    //               s->fval, size);

        }
    while (status == GSL_CONTINUE && iter < 100);

    // output results
    if (status!=GSL_SUCCESS) {
        lres = -(s->fval);
    }

    GlobalVar.slogcout(4) << "   " << ((status==GSL_SUCCESS)?' ':'~') << "MLE: " << lres << " at ( ";
    const flxVec& pvec = likeli.get_pvec();
    for (tuint i=0;i<pvec.get_N(); ++i) {
        if (i>0) {
            GlobalVar.slogcout(4) << ", ";
        }
        GlobalVar.slogcout(4) << pvec[i];
    }
    GlobalVar.slogcout(4) << " ) " << std::endl;

  } catch (...) {
    gsl_vector_free(x);
    gsl_vector_free(ss);
    gsl_multimin_fminimizer_free (s);
    likeli.get_pvec() = pvec_ini;
    throw;
  }

  gsl_vector_free(x);
  gsl_vector_free(ss);
  gsl_multimin_fminimizer_free (s);

  return lres;

}

const tdouble FunRVcheckX::calc()
{
    return rv->check_x(child_1->calc());
}

const std::string FunRVcheckX::write_v()
{
    throw FlxException_Crude("FunRVcheckX::write_v");
}

RBRV_entry_read_bayDA::RBRV_entry_read_bayDA(const bool readName, const bool readBrakets)
: RBRV_entry_read_base(readName, readBrakets)
{
    name_bayDA = new FlxString(false,false);
}

RBRV_entry_read_bayDA::~RBRV_entry_read_bayDA()
{
    delete name_bayDA;
}

RBRV_entry* RBRV_entry_read_bayDA::generate_entry(const std::string& family, tuint& running_iID)
{
  const std::string name = family + nameF->eval_word(true,false,true);
  flxBayDA& da_obj = RBRV_entry_read_bayDA::BayUpBox->get_DA(name_bayDA->eval_word(true));
  FlxFunction* z2x = new FlxFunction(new FunPara(1));
  FlxFunction* x2z = new FlxFunction(new FunPara(1));
  FlxFunction* dxdz = new FlxFunction(new FunNumber(ONE));
  FlxFunction* checkx = new FlxFunction(new FunRVcheckX(da_obj.get_rv_ptr(),new FunPara(1)));
  return new RBRV_entry_RV_UserTransform(name,running_iID++, true, z2x, x2z, dxdz, checkx, da_obj.get_rv_ptr(), false );
}




