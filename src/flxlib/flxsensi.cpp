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


#include <cmath>

#include "flxsensi.h"
#include "flxdata.h"
#include "flxobjrbrv.h"


flx_sensi_s1o::flx_sensi_s1o(const std::string name, const size_t N_learn, const tuint x_dim)
: name(name), N_learn(N_learn), x_dim(x_dim), x_learn(x_dim), y_learn(new flxVec(N_learn)), Nbrv(0), brecord_vec(NULL), x_record(x_dim), y_record(true), eval_last(ZERO), has_eval(false)
{
    for (tuint i=0;i<x_dim;++i) {
        x_learn[i] = new flxVec(N_learn);
    }
}

flx_sensi_s1o::~flx_sensi_s1o()
{
    if (brecord_vec) {
        for (tuint i=0;i<Nbrv;++i) {
            if (brecord_vec[i]) delete brecord_vec[i];
        }
        delete [] brecord_vec;
    }
    if (y_learn) delete y_learn;
    for (tuint i=0;i<x_dim;++i) {
        if (x_learn[i]) delete x_learn[i];
    }
}

void flx_sensi_s1o::allocate_brecord()
{
    if (brecord_vec) return;
    const tuint Nsplit_vec [] = { 1,2,4,8,16,32,64,100,200,400,1000 };
    Nbrv = sizeof(Nsplit_vec)/sizeof(*Nsplit_vec);
    brecord_vec = new flx_sensi_splitter*[Nbrv];
    for (tuint i=0;i<Nbrv;++i) brecord_vec[i] = NULL;
    const size_t n = y_record.get_size();
    for (tuint i=0;i<Nbrv;++i) {
        brecord_vec[i]  = new flx_sensi_splitter(Nsplit_vec[i],x_dim,x_learn,n);
    }
    flxVec v(x_dim);
    for (tuint i=0;i<Nbrv;++i) {
        flx_sensi_splitter& brecord = *(brecord_vec[i]);
        for (size_t j=0;j<n;++j) {
            for (tuint k=0;k<x_dim;++k) {
                v[k] = (*x_learn[k])[j];
            }
            brecord.record_value(v,(*y_learn)[j]);
        }
    }
    if (y_learn==NULL) throw FlxException_Crude("flx_sensi_s1o::allocate_brecord");
    delete y_learn; y_learn = NULL;
    for (tuint i=0;i<x_dim;++i) {
        delete x_learn[i]; x_learn[i] = NULL;
    }
}

void flx_sensi_s1o::record_value(const flxVec& x_vec, const tdouble y)
{
    // consistency check
        if (x_vec.get_N() != x_dim) {
            throw FlxException("flx_sensi_s1o::record_value", "dimension mismatch");
        }
    // update global statistics
        for (tuint i=0;i<x_dim;++i) {
        x_record[i] += x_vec[i];
        }
        y_record += y;
    // record values in corresponding batch
        if (brecord_vec==NULL) {
            const size_t n = y_record.get_size();
            for (tuint i=0;i<x_dim;++i) {
                (*x_learn[i])[n-1] = x_vec[i];
            }
            (*y_learn)[n-1] = y;
            // initiate the splitting (so that we do not need to store the samples anymore)
            if (n==N_learn) {
                allocate_brecord();
            }
        } else {
            for (tuint i=0;i<Nbrv;++i) {
                flx_sensi_splitter& brecord = *(brecord_vec[i]);
                brecord.record_value(x_vec,y);
            }
        }
}

const tdouble flx_sensi_s1o::eval()
{
    const pdouble y_mean = y_record.get_mean_p();
    const pdouble y_var = y_record.get_variance_p();
    const size_t N = y_record.get_size();
    allocate_brecord();
    flxVec Si_vec(Nbrv);
    flxVec diff_vec(Nbrv-1);
    diff_vec.set_nan();
    Si_vec.set_nan();
    size_t Nbatch_prev = 0;
    for (tuint i=0;i<Nbrv;++i) {
        flx_sensi_splitter& brecord = *(brecord_vec[i]);
        // make sure to stop if not enough batches are available (otherwise the sensitivity estimate will converge to ONE)
            if (brecord.get_Nbatches()*10>N) break;
            if (brecord.get_Nbatches()==Nbatch_prev) break;
            Nbatch_prev = brecord.get_Nbatches();
        const tdouble Si = brecord.eval(y_mean,y_var);
        Si_vec[i] = Si;
        if (i>0) diff_vec[i-1] = fabs(Si-Si_vec[i-1])/max(Si,Si_vec[i-1]);
        GlobalVar.slogcout(5) << "   batch-set " << i << "   " << Si << "   " << brecord.get_Nbatches();
        if (i>0) GlobalVar.slogcout(5) << "   " << diff_vec[i-1];
        GlobalVar.slogcout(5) << std::endl;
    }
    const size_t mid = diff_vec.get_minID();
    if (diff_vec[mid]>0.1) {
        std::ostringstream ssV;
        ssV << "Sensitivity estimate is likely inaccurate. It is recommended to increase the number of recorded values. (";
        ssV << "min_diff=" << diff_vec[mid] << ")";
        GlobalVar.alert.alert("flx_sensi_s1o::eval",ssV.str());
    }
    eval_last = Si_vec[mid+1];
    has_eval = true;
    return Si_vec[mid+1];
}

void flx_sensi_s1o::eval_dist(flxVec& svec, FlxRndCreator &rndCreator)
{
    if (!has_eval) eval();
    const size_t N = y_record.get_size();
    size_t Nbatch_prev = 0;
    flxVec tvec(svec.get_N());
    tdouble dmin = ZERO;
    for (tuint i=0;i<Nbrv;++i) {
        flx_sensi_splitter& brecord = *(brecord_vec[i]);
        // make sure to stop if not enough batches are available (otherwise the sensitivity estimate will converge to ONE)
            if (brecord.get_Nbatches()*100>N) break;
            if (brecord.get_Nbatches()==Nbatch_prev) break;
            Nbatch_prev = brecord.get_Nbatches();
        brecord.eval_dist(tvec,rndCreator,y_record);
        const tdouble dmin_ = fabs(tvec.get_Mean()-eval_last);
        GlobalVar.slogcout(5) << "   batch-set " << i << "   " << brecord.get_Nbatches() << "   mean=" << tvec.get_Mean() << "   sd=" << tvec.get_sd(tvec.get_Mean()) << std::endl;
        // stopping criterion
            if (i==0 || dmin_<dmin) {
                dmin = dmin_;
                svec = tvec;
            }
    }
}

flx_sensi_splitter_el::flx_sensi_splitter_el(const size_t Nsplit, const tdouble* vvecp, const size_t n)
: Nbatches(0)
{
    splits.reserve(Nsplit);
    // make sure to get a vector without NaN
        flxVec vvec(vvecp,n,false);
        const size_t nanc = vvec.count_nan();
        const size_t N = vvec.get_N()-nanc;
        if (N==0) {
            Nbatches = 1;
            return;
        }
        flxVec svec(N);
        vvec.copy_vals_without_nan(svec);
    // sort the vector
        svec.sort();
    // split into batches
        tdouble last_split = ZERO;  // dummy initialization
        const tdouble dq = ONE/(Nsplit+1);
        pdouble q = dq;
        const tdouble Nscale = ONE*1.5;
        for (tuint i=0;i<Nsplit;++i) {
            // find splitting value
                tdouble p, sv;
                if (q.cast2double()<=ONE/2) {
                    // map to percentile
                    p = iBeta_reg(Nscale,Nscale,q.cast2double());
                    // map percentile to original space
                    sv = flx_percentile(svec.get_tmp_vptr_const(), n, p);
                } else {
                    // map to percentile
                    p = iBeta_reg(Nscale,Nscale,(pdouble(ONE)-q).cast2double());
                    // map percentile to original space
                    sv = flx_percentile(svec.get_tmp_vptr_const(), n, p, true);
                }
            // add split
                if (i==0 || last_split!=sv) {
                    splits.push_back(sv);
                    last_split = sv;
                }
            // increase counter
                q += dq;
        }
    if (splits.size()>0) {
        Nbatches = splits.size() + 2;
    } else {
        Nbatches = 1;
    }
}

const size_t flx_sensi_splitter_el::get_batch_index(const tdouble val) const
{
    if (std::isnan(val)) {
        return Nbatches - 1;
    }
    if (Nbatches<=1) return 0;
    const size_t Nsplit = Nbatches - 2;   // as the last batch is reserved for NaN
    #if FLX_DEBUG
        if (Nsplit != splits.size()) {
            throw FlxException_Crude("flx_sensi_splitter_el::get_batch_index_01");
        }
    #endif
    if (Nsplit==0) {
        #if FLX_DEBUG
           throw FlxException_Crude("flx_sensi_splitter_el::get_batch_index_02");
        #endif
        return 0;
    }
    const tdouble* vp = &(splits[0]);
    return flx_find_pos(vp, Nsplit, val);

}

flx_sensi_splitter::flx_sensi_splitter(const size_t Nsplit, const tuint x_dim, const std::valarray< flxVec*>& vvecp, const size_t n)
: Nbatches(1), x_dim(x_dim), splits(x_dim), batches(NULL)
{
    for (tuint i=0;i<x_dim;++i) splits[i] = NULL;
    try {
        for (tuint i=0;i<x_dim;++i) {
            splits[i] = new flx_sensi_splitter_el(Nsplit,vvecp[i]->get_tmp_vptr_const(),n);
            Nbatches *= splits[i]->get_Nbatches();
        }
        batches = new std::valarray<flx_sensi_batch>(flx_sensi_batch(x_dim),Nbatches);
    } catch (...) {
        for (tuint i=0;i<x_dim;++i) {
            if (splits[i]) {
                delete splits[i]; splits[i] = NULL;
            }
        }
    }
}

flx_sensi_splitter::~flx_sensi_splitter()
{
    for (tuint i=0;i<x_dim;++i) {
        if (splits[i]) delete splits[i];
    }
    if (batches) delete batches;
}


void flx_sensi_splitter::record_value(const flxVec& x_vec, const tdouble y)
{
    size_t c = 0;
    size_t n = 1;
    for (tuint i=0;i<x_dim;++i) {
        const size_t si = splits[i]->get_batch_index(x_vec[i]);
        c += si*n;
        n *= splits[i]->get_Nbatches();
    }
    (*batches)[c].record_value(x_vec,y);
}

const tdouble flx_sensi_splitter::eval(const pdouble &y_mean, const pdouble &y_var)
{
    // evaluate variance of mean
        pdouble s;
        vdouble w;
        for (auto &&batch: *batches) {
            vdouble& y = batch.get_y_ref();
            if (y.get_size()==0) continue;
            const tdouble z = (tdouble)y.get_size();
            s += pow2(y.get_mean_p()-y_mean)*z;
            w += z;
        }
        const pdouble td = w.get_sum_p() - w.get_sum_of_squares_p()/w.get_sum_p();
        const tdouble p = (s/td).cast2double();
    // normalize by variance
        return p/y_var.cast2double();
}

void flx_sensi_splitter::eval_dist(flxVec& svec, FlxRndCreator &rndCreator, vdouble& y_record)
{
    // assemble Dirichlet distribution
        flxVec avec(Nbatches);
        for (tuint i=0;i<Nbatches;++i) {
            vdouble& y = (*batches)[i].get_y_ref();
            avec[i] = ONE + y.get_size();
        }
        RBRV_dirichlet rvs = RBRV_dirichlet(true, "flx_sensi_splitter::eval_dist", true, Nbatches, NULL, 0, NULL, &avec);
    // draw samples
        flxVec yvec(Nbatches);  // y-transform of Dirichlet distribution for batch weights (and y-transform to sample batch mean)
        flxVec wvec(Nbatches);   // batch probabilities
        flxVec mvec(Nbatches);   // the batch mean vector
        const tuint N = svec.get_N();
        for (tuint i=0;i<N;++i) {
            // sample the batch probabilities
                rndCreator.gen_smp(yvec);
                rvs.set_y(yvec.get_tmp_vptr());
                rvs.transform_y2x();
                rvs.get_x(wvec.get_tmp_vptr());
            // sample the batch mean
                rndCreator.gen_smp(yvec);
                pdouble mu;
                for (tuint j=0;j<Nbatches;++j) {
                    vdouble& y = (*batches)[j].get_y_ref();
                    if (y.get_size()>1) {
                        mvec[j] = y.get_mean_sample(yvec[j]);
                    } else {
                        mvec[j] = y_record.get_mean();
                    }
                    mu += mvec[j]*wvec[j];
                }
            // get the variance of the batch mean
                pdouble bv;
                for (tuint j=0;j<Nbatches;++j) {
                    pdouble tmp(mvec[j]);
                    tmp -= mu;
                    tmp *= tmp;
                    tmp *= wvec[j];
                    bv += tmp;
                }
            // sample total variance
                const tdouble y_var = y_record.get_var_sample(rndCreator.gen_smp());
            // evaluate Sobol coefficient
                svec[i] = bv.cast2double()/y_var;
        }
}

flx_sensi_batch::flx_sensi_batch(const tuint x_dim)
: x_record(x_dim)
{

}

void flx_sensi_batch::record_value(const flxVec& x_vec, const tdouble y)
{
    const tuint N = x_vec.get_N();
    if (x_record.size() != N) {
        throw FlxException("flx_sensi_batch::record_value", "Array dimensions do not match.");
    }
    for (tuint i=0;i<N;++i) {
        x_record[i] += x_vec[i];
    }
    y_record += y;
}
