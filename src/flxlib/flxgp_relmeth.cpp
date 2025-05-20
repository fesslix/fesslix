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

#include "flxgp_relmeth.h"
#include "flxdefault.h"

#include <algorithm>
#if FLX_PARALLEL
    #include <thread>
    #include <mutex>
#endif

rng_type randgen_local;

const tuint N_MCS_tqi = 1000;


flxGP_MCI::flxGP_MCI(flxGPProj_base& gp, const tuint N_init, const tuint Nreserve, const tuint user_seed_int, const tuint user_init_calls, std::ostream& logStream)
: gp(gp), Ndim(gp.get_Ndim()), user_seed_int(user_seed_int), user_init_calls(user_init_calls), tqi_vec(N_MCS_tqi), tqi_vec_rv_u(N_MCS_tqi), id_next_point(0), logStream(logStream),
  static_sum(ZERO),last_m(ZERO), last_n(ZERO), N_init(N_init), lh_samples(N_init*Ndim)
{
    dmV.reserve(Ndim*Nreserve);
    doV.reserve(Nreserve);
    mcs_pi.reserve(1000000);

    // generate pseudo-random samples for tqi_vec_rv_u
    rv_initialize(false,true,13560607,0,&randgen_local,false);
    for (tuint i=0;i<N_MCS_tqi;++i) {
        tqi_vec_rv_u[i] = rv_uniform(randgen_local);
    }

    init_RNG();
    FlxRndCreator rndc(&randgen_local);
    rndc.latin_hypercube(lh_samples.get_tmp_vptr(),N_init,Ndim);
    lh_samples *= 6.;
    lh_samples -= 3.;
}

void flxGP_MCI::init_RNG()
{
    // initialize the pseudo random number generator
        // so that it generates the same set of random numbers in each simulation run
        rv_initialize(false,true,user_seed_int,user_init_calls,&randgen_local,false);
}

const bool flxGP_MCI::is_point_unique(const flxVec& uvec_) const
{
    const size_t N = doV.size();
    size_t c = 0;
    for (size_t i=0;i<N;++i) {
        const flxVec svec(&(dmV[c]),Ndim,false);
        const tdouble t = svec.get_Norm2_NOroot_of_distance(uvec_)/(2*Ndim);
        if (t<=1e-6) return false;
        c += Ndim;
    }
    return true;
}

void flxGP_MCI::generate_sample(flxVec &uvec_)
{
    rv_normal(uvec_,randgen_local);
}

void flxGP_MCI::get_next_point(flxVec& uvec_)
{
    if (N_init>0) {
        --N_init;
        flxVec tv(lh_samples.get_tmp_vptr()+N_init*Ndim,Ndim,false,false);
        uvec_ = tv;
        if ( !is_point_unique(uvec_) ) {
            get_next_point(uvec_);
        }
        return;
    }
    init_RNG();
    for (tulong i=0;i<=id_next_point;++i) {
        generate_sample(uvec_);
    }
    ++id_next_point;
    while ( !is_point_unique(uvec_) ) {
        generate_sample(uvec_);
        ++id_next_point;
    };
}

void flxGP_MCI::register_sample(const tdouble lsfval, const flxVec& uvec_, const bool condition_GP, const bool register_pvec, const bool learn_noise)
{
    doV.push_back(lsfval);
    for (tuint i=0;i<Ndim;++i) {
        dmV.push_back(uvec_[i]);
    }
    if (condition_GP) {
        gp.register_observation(flxGP_data_ptr(this),register_pvec,learn_noise);
    }
}

void flxGP_MCI::optimize_gp_para(const tuint iterMax)
{
    gp.optimize(iterMax,false);
}

const tdouble flxGP_MCI::tqi_eval(const tdouble m, const tdouble n) const
{
    //return m/n;
    //return (m+ONE)/(n+2*ONE);
    // return iBeta_reg_inv(m+ONE,n-m+ONE,0.95);
    return iBeta_reg_inv(m+ONE,n-m+ONE,0.99);
}

struct tqi_struct {
    flxVec* tqi_vec;
    tulong n;
    tdouble pf_ref;
    tdouble q;
};

double tqi_rsfun(const tdouble x, void *params)
{
    tqi_struct *p = (tqi_struct *)params;
    flxVec& tqi_vec = *(p->tqi_vec);
    const tulong n = p->n;
    const tdouble pf_ref = p->pf_ref;
    const tdouble q = p->q;
    tdouble s = ZERO;
    const tdouble pf = exp(x)*pf_ref;
    if (pf>=ONE) {
        throw FlxException_Crude("flxgp::tqi_rsfun");
    }
    for (tuint i=0;i<N_MCS_tqi;++i) {
        s += iBeta_reg(tqi_vec[i]+ONE,n-tqi_vec[i]+ONE,pf);
    }
    return s/N_MCS_tqi-q;
}

const tdouble flxGP_MCI::get_mean_tqi(const tdouble ref_m, const tulong n, const tulong* skip_id, const tuint nrep)
{
    // Variant: in case of skip_id, assume we can fully reduce Kiriging-uncertainty
        if (skip_id) {
            return tqi_eval(ref_m*nrep,n*nrep);
        }
    // check whether a specific id needs to be skipped
        tdouble skip_pr = ZERO;
        if (skip_id) {
            if (*skip_id < mcs_pi.size()) {
                skip_pr = mcs_pi[*skip_id];
                mcs_pi[*skip_id] = ZERO;
            }
        }
    // perform MCS to sample realizations of m
        const size_t N_ = mcs_pi.size();
        tqi_vec = tqi_vec_rv_u;
        std::for_each(tqi_vec.begin(),tqi_vec.end(),[&](auto&& val){
            tdouble m = static_sum;
            // assume maximum dependence
                for (tulong j=0;j<N_;++j) {
                    if (val<=mcs_pi[j]) {
                        m += ONE;
                    }
                }
                m *= nrep;
            val = m + skip_pr*nrep;
        });
        // for (tuint i=0;i<N_MCS_tqi;++i) {
        //     tdouble m = static_sum;
        //     // assume full independence
        //         // for (tuint k=0;k<nrep;++k) {
        //         //     for (tulong j=0;j<n;++j) {
        //         //         const tdouble r = rv_uniform(randgen_local);
        //         //         if (r<=mcs_pi[j]) {
        //         //             ++m;
        //         //         }
        //         //     }
        //         // }
        //     // assume maximum dependence
        //         const tdouble r = tqi_vec[i];
        //         for (tulong j=0;j<N_;++j) {
        //             if (r<=mcs_pi[j]) {
        //                 m += ONE;
        //             }
        //         }
        //         m *= nrep;
        //     tqi_vec[i] = m + skip_pr*nrep;
        // }
    // scale m-values
        const tdouble tqvm = tqi_vec.get_Mean();
        if (tqvm>GlobalVar.TOL()) {
            tqi_vec *= (ref_m*nrep)/tqvm;
        }
    // root search of 99%-quantile
        tqi_struct tqis;
        tqis.tqi_vec = &tqi_vec;
        tqis.n = n*nrep;
        tqis.pf_ref = (ref_m+ONE)/(n*nrep+2*ONE);
        tqis.q = 0.99;
        const tdouble x_rs_res = flx_RootSearch_RegulaFalsi(tqi_rsfun,&tqis,log(1.),log(std::min(100.,(ONE-1e-7)/tqis.pf_ref)),1e-4,1e-4*tqis.pf_ref,NULL);  // &GlobalVar.slogcout(1)
        return exp(x_rs_res)*tqis.pf_ref;

    // perform MCS to sample realizations of tqi
    //     tdouble s = ZERO;
    //     for (tuint i=0;i<N_MCS_tqi;++i) {
    //         tqi_vec[i] = tqi_eval(tqi_vec[i],n*nrep);
    //         s += tqi_vec[i];
    //     }
    //
    //
    // // revert the skipped id
    //     if (skip_pr>ZERO) {
    //         mcs_pi[*skip_id] = skip_pr;
    //     }
    // // evaluate quantity of interest
    //     // in this case, the 95% quantile of the tqi-sample-vector
    //     tqi_vec.sort();
    //     return tqi_vec[tuint(0.90*N_MCS_tqi)];
    //     // return tqi_vec[tuint(0.95*N_MCS_tqi)];
    //     //return tqi_vec.get_Var(s/N_MCS_tqi);
}

const tdouble flxGP_MCI::tqi_eval_covar() const
{
    pdouble s1, s2, s3;
    for (tuint i=0;i<N_MCS_tqi;++i) {
        const tdouble m = tqi_vec[i];
        s1 += (m+ONE)*(last_n-m+ONE)/(pow2(last_n+2*ONE)*(last_n+3*ONE));
        const tdouble mu_i = (m+ONE)/(last_n+2*ONE);
        s2 += pow2(mu_i);
        s3 += mu_i;
    }
    s1 /= N_MCS_tqi;
    s2 /= N_MCS_tqi;
    s3 /= N_MCS_tqi;
    const tdouble mu = s3.cast2double();
    s3 *= s3;
    s1 += s2;
    s1 -= s3;
    return sqrt(s1.cast2double())/mu;
}

const tdouble flxGP_MCI::tqi_eval_pr(const tdouble p) const
{
    pdouble s;
    for (tuint i=0;i<N_MCS_tqi;++i) {
        s += iBeta_reg(tqi_vec[i]+ONE,last_n-tqi_vec[i]+ONE,p);
    }
    s /= N_MCS_tqi;
    return s.cast2double();
}

const tdouble flxGP_MCI::simulate_GP_mci(const tulong Nsmpls, tdouble& err, int& proposed_action_id)
{
    proposed_action_id = 0;
    init_RNG();
    tulong id_worst_point = 0;
    tulong id_worst_mcspi = 0;
    const tdouble pi_TOL = 1e-8;
    tdouble Uval_worst_point = std::numeric_limits<tdouble>::infinity();
    // ensure that 'mcs_pi' is large engough
        mcs_pi.clear();
    // conduct Monte Carlo simulation
        tdouble sum = ZERO;
        static_sum = ZERO;
        tdouble vsum = ZERO;
        #if FLX_PARALLEL
            std::mutex mutex_smpling;
            std::mutex mutex_postp;
            tulong i=0;
            auto thread_fun = [&](const tuint pid) {
                flxVec uvec(Ndim);
                while (true) {
                    tuint j;
                    // generate sample
                    {
                        std::lock_guard<std::mutex> guard(mutex_smpling);
                        if (i<Nsmpls) {
                            j = i++;
                            generate_sample(uvec);
                        } else {
                            break;
                        }
                    }
                    // evaluate surrogate model
                        tdouble smpl_mean, smpl_var;
                        gp.predict_mean_var(uvec,true,smpl_mean,smpl_var);
                    // probability of value being smaller than zero
                        const tdouble y = smpl_mean/sqrt(smpl_var);
                        const tdouble p_i = rv_Phi(-y);
                    {
                        std::lock_guard<std::mutex> guard(mutex_postp);
                        sum += p_i;
                        vsum += p_i*(ONE-p_i);
                    // learning function
                        if (fabs(y)<Uval_worst_point) {
                            id_worst_point = j;
                            id_worst_mcspi = mcs_pi.size();
                            Uval_worst_point = fabs(y);
                        }
                    // append to mcs_pi
                        if (p_i<pi_TOL || p_i>ONE-pi_TOL) {
                            static_sum += p_i;
                        } else {
                            mcs_pi.push_back(p_i);
                        }
                    }
                }
            };
            std::vector<std::thread> threadVec;
            const tuint Nt = GlobalVar.max_parallel_threads;
            // start parallel threads
                for (tuint i=0;i<Nt;++i) {
                    threadVec.push_back(std::thread(thread_fun,i));
                }
            // wait for threads to finish
                for (tuint i=0;i<Nt;++i) {
                    threadVec[i].join();
                }
        #else
            flxVec uvec(Ndim);
            for (tulong i=0;i<Nsmpls;++i) {
                generate_sample(uvec);
                // evaluate surrogate model
                    tdouble smpl_mean, smpl_var;
                    gp.predict_mean_var(uvec,true,smpl_mean,smpl_var);
                // probability of value being smaller than zero
                    const tdouble y = smpl_mean/sqrt(smpl_var);
                    const tdouble p_i = rv_Phi(-y);
                    sum += p_i;
                    vsum += p_i*(ONE-p_i);
                // learning function
                    if (fabs(y)<Uval_worst_point) {
                        id_worst_point = i;
                        id_worst_mcspi = mcs_pi.size();
                        Uval_worst_point = fabs(y);
                    }
                // append to mcs_pi
                    if (p_i<pi_TOL || p_i>ONE-pi_TOL) {
                        static_sum += p_i;
                    } else {
                        mcs_pi.push_back(p_i);
                    }
            }
        #endif
    // set id_next_point
        id_next_point = id_worst_point;
        N_init = 0;  // no more latin-hypercube samples are used
    // quantify uncertainty about estimate
        last_m = sum;
        last_n = Nsmpls;
        const tdouble pf_mle = sum/Nsmpls;      // MLE estimate of Pf
        const tdouble mean_m = sum;             // sample mean of hits
        const tdouble mean_pf_bayesian = (sum+ONE)/(Nsmpls+2);
        // E[tqi]: expected target quantity of interest
            const tdouble tqi = get_mean_tqi(mean_m,Nsmpls)/mean_pf_bayesian;
        // estimated E[tqi] if one more LSF call is invested
            const tdouble tqi_1LSF = get_mean_tqi(mean_m,Nsmpls,&id_worst_mcspi)/mean_pf_bayesian;
        // estimated E[tqi] if twice as much surrogate samples are unsed
            const tdouble tqi_2N = get_mean_tqi(mean_m,Nsmpls,NULL,2)/mean_pf_bayesian;
    // recommend an action
        // 0: continue
        // 1: increase surrogate samples
        // 2: stop
        if (tqi_2N<=tqi_1LSF) {
            proposed_action_id = 1;
        }
        const tdouble af = (tqi_2N-tqi_1LSF)/tqi;
        err = tqi-ONE;
    // output
        logStream << "  MCS on surrogate: E[pf]=" << GlobalVar.Double2String(mean_pf_bayesian,false,2,8)
            << "  Pr[pf<" << GlobalVar.Double2String(tqi*mean_pf_bayesian,false,1,6) << "]≈99%"
            << "  err=" << GlobalVar.Double2String(err,false,2,4)
            << "  minU=" << GlobalVar.Double2String(Uval_worst_point,false,2,5);
        GlobalVar.Double2String_setType(3);
        logStream << "  af=" << GlobalVar.Double2String(af,false,2,5);
        GlobalVar.Double2String_setType(DEFAULT_FPC_TYPE);
        logStream << "  N=" << GlobalVar.Double2String(Nsmpls) << std::endl;
    return pf_mle;
}

void flxGP_MCI::output_summary()
{
    tdouble t;
    logStream << "  Model calls taken into account: " << doV.size() << std::endl;
    logStream << "  Samples in surrogate MCS:       " << last_n << std::endl;

    logStream << "  Unbiased estimate of P_f (from maximum likelihood estimation, MLE):" << std::endl;
    logStream << "      Expectation of P_f                       = " << GlobalVar.Double2String(last_m/last_n,false,2,8) << std::endl;
    GlobalVar.Double2String_setType(3);
    t = sqrt((ONE-last_m/last_n)/last_m);
    logStream << "      Coefficient of Variation (C.o.V.)        = " << GlobalVar.Double2String(t*100,false,2,5) << "%" << std::endl;
    if (last_m>GlobalVar.TOL()) {
        t = -rv_InvPhi(last_m/last_n);
        logStream << "      Corresponding reliability index (Beta)   = " << GlobalVar.Double2String(t,false,2,5) << std::endl;
    }
    GlobalVar.Double2String_setType(DEFAULT_FPC_TYPE);

    logStream << "  Bayesian statistics for P_f:" << std::endl;
    const tdouble pf_bayes = (last_m+ONE)/(last_n+2*ONE);
    logStream << "      Expectation for P_f                      = " << GlobalVar.Double2String(pf_bayes,false,2,8) << std::endl;
    GlobalVar.Double2String_setType(3);
    get_mean_tqi(last_m,last_n);
    const tdouble CoV_bayes = tqi_eval_covar();
    logStream << "      C.o.V. for P_f                           = " << GlobalVar.Double2String(CoV_bayes*100,false,2,5) << "%" << std::endl;
    // reliability index :: 1st order approximation
        const tdouble Var_bayes = pow2(CoV_bayes*pf_bayes);
        const tdouble u_ = rv_InvPhi_noAlert(pf_bayes);
        const tdouble phi_u_ = rv_phi(u_);
        const tdouble a1 = -ONE/phi_u_;
        const tdouble a2 = -u_/pow2(phi_u_);
        // second-order approximation of mean
        const tdouble SMean = -u_+a2/2*Var_bayes;
        // first-order approximation of variance
        const tdouble SVar = pow2(a1)*Var_bayes;
        // compute CoV of samples
        const tdouble SCoV = sqrt(SVar)/SMean;
    logStream << "      Expectation for reliability index (beta) = " << GlobalVar.Double2String(SMean,false,2,5) << std::endl;
    logStream << "      C.o.V. for beta                          = " << GlobalVar.Double2String(SCoV*100,false,2,5) << "%" << std::endl;
    GlobalVar.Double2String_setType(DEFAULT_FPC_TYPE);
    logStream << "      upper credible intervals:" << std::endl;
    tdouble parr [5] = { 0.5, 0.75, 0.9, 0.95, 0.99 };
    for (tuint i=0;i<5;++i) {
        // root search of quantile
            tqi_struct tqis;
            tqis.tqi_vec = &tqi_vec;
            tqis.n = last_n;
            tqis.pf_ref = pf_bayes;
            tqis.q = parr[i];
            const tdouble x_rs_res = flx_RootSearch_RegulaFalsi(tqi_rsfun,&tqis,log(0.5),log(std::min(100.,(ONE-1e-7)/tqis.pf_ref)),1e-4,1e-4*tqis.pf_ref,NULL);  // &GlobalVar.slogcout(1)
            t = exp(x_rs_res)*tqis.pf_ref;
        logStream << "          Pr[ P_f < " << GlobalVar.Double2String(t,false,2,8) << " ] = ";
        GlobalVar.Double2String_setType(3);
        logStream << GlobalVar.Double2String(parr[i],false,2,5);
        logStream << "       Pr[beta > " << GlobalVar.Double2String(-rv_InvPhi(t),false,3,6) << " ] = " << GlobalVar.Double2String(parr[i],false,2,5) << std::endl;
        GlobalVar.Double2String_setType(DEFAULT_FPC_TYPE);
    }
    logStream << std::endl;
    tdouble qarr [8] = { 1e-2, 1e-3, 1e-4, 1e-5, 1e-6, 1e-7, 1e-8, 1e-9 };
    for (tuint i=0;i<5;++i) {
        t = tqi_eval_pr(qarr[i]);
        logStream << "          Pr[ P_f < " << GlobalVar.Double2String(qarr[i],false,2,8) << " ] = " << GlobalVar.Double2String(t,false,2,9);
        GlobalVar.Double2String_setType(3);
        logStream << "   Pr[beta > " << GlobalVar.Double2String(-rv_InvPhi(qarr[i]),false,3,6) << " ] = ";
        GlobalVar.Double2String_setType(DEFAULT_FPC_TYPE);
        logStream << GlobalVar.Double2String(t,false,2,8) << std::endl;
    }

}

// ------------------------------------------------------------------------------------------------

flxGP_data_base* flxGP_data_ptr::copy() const
{
    return new flxGP_data_ptr(akmcs);
}

const tdouble* flxGP_data_ptr::get_data_ptr_mtx(tuint &N_obsv, const tuint N_dim)
{
    #if FLX_DEBUG
        if (N_dim != akmcs->get_Ndim()) {
            throw FlxException_Crude("flxGP_data_ptr::get_data_ptr_mtx_01");
        }
    #endif
    N_obsv = akmcs->get_N_model_calls();
    if (N_obsv==0) {
        throw FlxException("flxGP_data_ptr::get_data_ptr_mtx_02", "No data registered to condition on.");
    }
    return &(akmcs->dmV[0]);
}

const tdouble* flxGP_data_ptr::get_data_ptr_vec(const tuint N_obsv)
{
    #if FLX_DEBUG
        if (N_obsv!=akmcs->get_N_model_calls()) {
            throw FlxException_Crude("flxGP_data_ptr::get_data_ptr_vec");
        }
    #endif
    return &(akmcs->doV[0]);
}

