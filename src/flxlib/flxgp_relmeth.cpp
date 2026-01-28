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

#include "flxgp_relmeth.h"
#include "flxdefault.h"
#include "flxparse.h"

#include <algorithm>
#if FLX_PARALLEL
    #include <thread>
    #include <mutex>
    #ifdef HAS_STD_EXECUTION_PAR
        #include <execution>
    #endif
    #include <ranges>
#endif

rng_type randgen_local;

const tuint N_MCS_tqi = 1000;


flxGP_MCI::flxGP_MCI(flxGPProj_base& gp, const tuint Nreserve, const tuint user_seed_int, const tuint user_init_calls, const tdouble tqi_val, const bool allow_decrease_of_N, const bool account4noise)
: gp(gp), Ndim(gp.get_Ndim()), user_seed_int(user_seed_int), user_init_calls(user_init_calls), tqi_val(tqi_val), allow_decrease_of_N(allow_decrease_of_N),account4noise(account4noise),
  tqi_vec(N_MCS_tqi), tqi_vec_rv_u(N_MCS_tqi), id_next_point(0),
  static_sum(ZERO),last_m(ZERO), last_n(ZERO)
{
    dmV.reserve(Ndim*Nreserve);
    doV.reserve(Nreserve);
    mcs_pi.reserve(1000000);

    // generate pseudo-random samples for tqi_vec_rv_u
    rv_initialize(false,true,13560607,0,&randgen_local,false);
    for (tuint i=0;i<N_MCS_tqi;++i) {
        tqi_vec_rv_u[i] = rv_uniform(randgen_local);
    }
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
        if (t<=1e-8) {
            return false;
        }
        c += Ndim;
    }
    return true;
}

void flxGP_MCI::generate_sample(flxVec &uvec_)
{
    rv_normal(uvec_,randgen_local);
}

void flxGP_MCI::assemble_lh_samples(flxVec& lh_samples, const tdouble box_bounds)
{
    init_RNG();
    FlxRndCreator rndc(&randgen_local);
    const tuint N_init = lh_samples.get_N()/Ndim;
    if (N_init*Ndim != lh_samples.get_N()) {
        throw FlxException_Crude("flxGP_MCI::assemble_lh_samples");
    }
    rndc.latin_hypercube(lh_samples.get_tmp_vptr(),N_init,Ndim);
    lh_samples *= 2*box_bounds;
    lh_samples -= box_bounds;
}

void flxGP_MCI::get_next_point(flxVec& uvec_)
{
    init_RNG();
    for (tulong i=0;i<=id_next_point;++i) {
        generate_sample(uvec_);
    }
    ++id_next_point;
    while ( !is_point_unique(uvec_) ) {
        logStream << "WARNING [flxGP_MCI::get_next_point]: selected point is not unique" << std::endl;
        generate_sample(uvec_);
        ++id_next_point;
    };
}

void flxGP_MCI::register_sample(const tdouble lsfval, const flxVec& uvec_)
{
    doV.push_back(lsfval);
    for (tuint i=0;i<Ndim;++i) {
        dmV.push_back(uvec_[i]);
    }
}

void flxGP_MCI::condition_on_data(const bool register_pvec)
{
    gp.register_observation(flxGP_data_ptr(this),register_pvec,register_pvec&&account4noise);
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
    return iBeta_reg_inv(m+ONE,n-m+ONE,tqi_val);
}

struct tqi_struct {
    flxVec* tqi_vec;
    tulong n;
    tdouble pf_ref;
    tdouble q;
};

/**
* @brief root-search-function for getting a specified quantile of the expected posterior probability of failure
*/
tdouble tqi_rsfun(const tdouble x, void *params)
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
        const tulong n_hits = min(tqi_vec[i],n);   // ensure that 'tqi_vec[i]' is not larger than 'n'
        s += iBeta_reg(n_hits+ONE,n-n_hits+ONE,pf);
    }
    return s/N_MCS_tqi-q;
}

const tdouble flxGP_MCI::get_mean_tqi(const tdouble ref_m, const tulong n, const tulong* skip_id, const tdouble nrep)
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
        #ifdef HAS_STD_EXECUTION_PAR
            std::vector<size_t> indices(tqi_vec.get_N());
            std::iota(indices.begin(), indices.end(), 0);
            std::for_each(std::execution::par, indices.begin(), indices.end(), [&](size_t i) {
                auto& val = tqi_vec[i];
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
        #else
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
        #endif
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
    // root search of tqi_val=99%-quantile
        tqi_struct tqis;
        tqis.tqi_vec = &tqi_vec;
        tqis.n = n*nrep;
        tqis.pf_ref = (ref_m+ONE)/(n*nrep+2*ONE);
        tqis.q = tqi_val;
        try {
            const tdouble x_rs_res = flx_RootSearch_RegulaFalsi(tqi_rsfun,&tqis,log(1.),log(std::min(100.,(ONE-1e-7)/tqis.pf_ref)),1e-4,1e-4*tqis.pf_ref,nullptr);  // &GlobalVar.slogcout(1)
            return exp(x_rs_res)*tqis.pf_ref;
        } catch (FlxException& e) {
            return ONE;
        }

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

py::dict flxGP_MCI::simulate_GP_mci(const tulong Nsmpls, tdouble& err, int& proposed_action_id)
{
    py::dict res;
    proposed_action_id = 0;
    init_RNG();
    tulong id_worst_point = 0;
    tulong id_worst_mcspi = 0;
    const tdouble pi_TOL = 1e-8;
    tdouble Uval_worst_point = std::numeric_limits<tdouble>::infinity();
    tdouble sd_worst_point = std::numeric_limits<tdouble>::infinity();
    // ensure that 'mcs_pi' is large engough
        mcs_pi.clear();
    // conduct Monte Carlo simulation
        pdouble sum;
        static_sum = ZERO;
        tulong sum_no_Krig_unc = 0;
        pdouble sd_avg;
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
                        gp.predict_mean_var(uvec,false,smpl_mean,smpl_var);
                    // probability of value being smaller than zero
                        const tdouble smpl_sd = sqrt(smpl_var);
                        const tdouble y = smpl_mean/smpl_sd;
                        const tdouble p_i = rv_Phi(-y);
                    {
                        std::lock_guard<std::mutex> guard(mutex_postp);
                        sum += p_i;
                        vsum += p_i*(ONE-p_i);
                        if (smpl_mean<=ZERO) {
                            ++sum_no_Krig_unc;
                        }
                        sd_avg += smpl_sd;
                    // learning function
                        if (fabs(y)<Uval_worst_point) {
                            id_worst_point = j;
                            id_worst_mcspi = mcs_pi.size();
                            Uval_worst_point = fabs(y);
                            sd_worst_point = smpl_sd;
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
                    gp.predict_mean_var(uvec,false,smpl_mean,smpl_var);
                // probability of value being smaller than zero
                    const tdouble smpl_sd = sqrt(smpl_var);
                    const tdouble y = smpl_mean/smpl_sd;
                    const tdouble p_i = rv_Phi(-y);
                    sum += p_i;
                    vsum += p_i*(ONE-p_i);
                    if (smpl_mean<=ZERO) {
                        ++sum_no_Krig_unc;
                    }
                    sd_avg += smpl_sd;
                // learning function
                    if (fabs(y)<Uval_worst_point) {
                        id_worst_point = i;
                        id_worst_mcspi = mcs_pi.size();
                        Uval_worst_point = fabs(y);
                        sd_worst_point = smpl_sd;
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
    // quantify uncertainty about estimate
        last_m = sum.cast2double();
        last_n = Nsmpls;
        const tdouble pf_mle = sum.cast2double()/Nsmpls;      // MLE estimate of Pf
        const tdouble mean_m = sum.cast2double();             // sample mean of hits
        const tdouble mean_pf_bayesian = (sum.cast2double()+ONE)/(Nsmpls+2);
        // E[tqi]: expected target quantity of interest
            const tdouble tqi = get_mean_tqi(mean_m,Nsmpls)/mean_pf_bayesian;
        // estimated E[tqi] if one more LSF call is invested
            const tdouble tqi_1LSF = get_mean_tqi(mean_m,Nsmpls,&id_worst_mcspi)/mean_pf_bayesian;
        // estimated E[tqi] if twice as much surrogate samples are unsed
            const tdouble tqi_2N = get_mean_tqi(mean_m,Nsmpls,nullptr,2*ONE)/mean_pf_bayesian;
        // estimated E[tqi] if one more LSF call is invested AND the number of surrogate samples is halfed
            const tdouble tqi_1LSF_N_half = get_mean_tqi(mean_m,Nsmpls,&id_worst_mcspi,ONE/2)/mean_pf_bayesian;
    // recommend an action
        // 0: call actual LSF/model
        // 1: increase surrogate samples
        // 2: stop
        // 3: call actual LSF/model AND reduce surrogate samples
        if (allow_decrease_of_N && tqi_1LSF_N_half < tqi_2N) {
            proposed_action_id = 3;
        } else if (tqi_2N<=tqi_1LSF) {
            proposed_action_id = 1;
        }
        const tdouble af = (tqi_2N-tqi_1LSF)/tqi;
        err = tqi-ONE;
    // output
        res["mean_pf_bayesian"] = mean_pf_bayesian;
        res["pf_mle"] = pf_mle;
        res["pf_no_Kriging_uncertainty"] = tdouble(sum_no_Krig_unc)/Nsmpls;
        res["Pr_q_tqi"] = tqi*mean_pf_bayesian;  // Pr[pf<p]≈tqi=99%
        res["err"] = err;
        res["af"] = af;
        res["r"] = tqi;
        res["r_increase_N_surrogate"] = tqi_2N;
        res["r_no_Kriging_uncertainty"] = tqi_eval(mean_m,Nsmpls)/mean_pf_bayesian;
        res["r_no_Kriging_uncertainty_AND_N_half"] = tqi_1LSF_N_half;
        res["propose_to_increase_N_smpls_surrogate"] = proposed_action_id;
        res["N"] = Nsmpls;
        res["N_model_calls"] = get_N_model_calls();
        res["Uval_worst_point"] = Uval_worst_point;
        res["sd_worst_point"] = sd_worst_point;
        res["sd_avg"] = sd_avg.cast2double()/Nsmpls;
        // try to extract information about the kernel
        {
            py::dict gp_info = gp.info();
            if (gp_info.contains("kernel")) {
                py::dict kernel_info = gp_info["kernel"];
                if (kernel_info.contains("kernel_sd")) {
                    const tdouble kernel_sd = parse_py_para_as_float("kernel_sd",kernel_info,true);
                    res["kernel_sd"] = kernel_sd;
                    if (kernel_info.contains("para_vec") && kernel_info.contains("n_vec")) {
                        flxVec kernel_para_vec = parse_py_obj_as_flxVec(kernel_info["para_vec"],"kernel::para_vec");
                        flxVec kernel_n_vec = parse_py_obj_as_flxVec(kernel_info["n_vec"],"kernel::n_vec");
                        const tuint kernel_N = kernel_para_vec.get_N();
                        if (kernel_N<=2) {
                            throw FlxException_Crude("flxGP_MCI::simulate_GP_mci_80");
                        }
                        // create Python-array
                            // Allocate memory for the return array
                            auto corrl_res_buf = py::array_t<tdouble>(kernel_N-1);
                            // Get the buffer info to access the underlying return data
                            py::buffer_info res_buf_info = corrl_res_buf.request();
                            tdouble* corrl_res_ptr = static_cast<tdouble*>(res_buf_info.ptr);
                        for (tuint i=1;i<kernel_N;++i) {
                            corrl_res_ptr[i-1] = kernel_para_vec[i]*kernel_n_vec[i];
                        }
                        res["kernel_corrl"] = corrl_res_buf;
                    }
                }
            }
            if (gp_info.contains("noise")) {
                res["kernel_noise"] = gp_info["noise"];
            }
        }
    return res;
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
    get_mean_tqi(last_m,last_n);  // to assemble tqi_vec for the call of tqi_eval_covar()
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
            try {
                const tdouble x_rs_res = flx_RootSearch_RegulaFalsi(tqi_rsfun,&tqis,log(0.5),log(std::min(100.,(ONE-1e-7)/tqis.pf_ref)),1e-4,1e-4*tqis.pf_ref,nullptr);  // &GlobalVar.slogcout(1)
                t = exp(x_rs_res)*tqis.pf_ref;
            } catch (FlxException& e) {
                t = ONE;
            }
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

