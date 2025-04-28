"""Tools module for Fesslix.

"""

# Fesslix - Stochastic Analysis
# Copyright (C) 2010-2025 Wolfgang Betz
#
# Fesslix is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# Fesslix is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with Fesslix.  If not, see <http://www.gnu.org/licenses/>. 

import fesslix as flx

import numpy as np
from scipy import stats as scipy_stats
from scipy import optimize as scipy_opt
import scipy.interpolate


##################################################
# discretization                                 #
##################################################

def detect_bounds_x(rv, config_dict, q_low=1e-3, q_up=None, mode='ignore'):
    """Makes sure that x_low and x_up are assigned in config_dict.

    """
    if q_up is None:
        q_up = 1.-q_low
    if 'x_low' not in config_dict:
        config_dict['x_low'] = rv.icdf(q_low)
    else:
        if (config_dict['x_low'] is None):
            config_dict['x_low'] = rv.icdf(q_low)
        else:
            if mode=='overwrite':
                config_dict['x_low'] = rv.icdf(q_low)
            elif mode=='minmax':
                config_dict['x_low'] = min(config_dict['x_low'], rv.icdf(q_low))
    if 'x_up' not in config_dict:
        config_dict['x_up'] = rv.icdf(q_up)
    else:
        if (config_dict['x_up'] is None):
            config_dict['x_up'] = rv.icdf(q_up)
        else:
            if mode=='overwrite':
                config_dict['x_up'] = rv.icdf(q_up)
            elif mode=='minmax':
                config_dict['x_up'] = max(config_dict['x_up'], rv.icdf(q_up))


def discretize_x(x_low, x_up, x_disc_N=int(1e3), x_disc_shift=False, x_disc_on_log=False):
    """Returns an array with discretized values for the x-axis."""
    if x_disc_on_log:
        if x_low<0. or x_up<0.:
            raise NameError(f'ERROR 202202071550: {x_low} {x_up}')
        x_low = np.log(x_low)
        x_up = np.log(x_up)
    x, dx = np.linspace(x_low,x_up,num=x_disc_N,endpoint=(not x_disc_shift),retstep=True)
    if x_disc_shift:
        x += dx/2
    if x_disc_on_log:
        x = np.exp(x)
    return x

def discretize_x_get_diff(x_low, x_up, x_disc_N=int(1e3), x_disc_on_log=False):
    """Returns an array with discretized values for the x-axis and additionally returns an vector with the size of the elements."""
    x = discretize_x(x_low=x_low, x_up=x_up, x_disc_N=x_disc_N, x_disc_shift=True, x_disc_on_log=x_disc_on_log)
    N = len(x)
    dx = np.empty(N)
    x_prev = x_low
    for i in range(N):
        if i+1<N:
            x_next = (x[i]+x[i+1])/2
        else:
            x_next = x_up
        dx[i] = (x_next-x_prev)
        x_prev = x_next
    return x, dx

def discretize_stdNormal_space(q_low=1e-3, q_up=None, x_disc_N=int(1e3)):
    """Returns an array with discretized values on U-space (standard Normal space).

    """
    if q_up is None:
        q_up = 1. - q_low
    xl = flx.cdfn_inv(q_low)
    xu = flx.cdfn_inv(q_up)
    return discretize_x(x_low=xl, x_up=xu, x_disc_N=x_disc_N, x_disc_shift=False, x_disc_on_log=False)



##################################################
# Working with float arrays                      #
##################################################

def fit_tail_to_data(tail_data_transformed, bound=None):
    res = { 'models':{} }
    def neg_log_likelihood(dist, data, params):
        return -np.sum(dist.logpdf(data, *params))
    ## ======================
    ## Fit generalized pareto
    ## ======================
    ## Fit GPD to the tail data
    gpd_params = scipy_stats.genpareto.fit(tail_data_transformed,floc=0.)
    res_ = { 'type':'genpareto', 'xi':gpd_params[0], 'scale':gpd_params[2] }
    ## Kolmogorov–Smirnov Test
    D_gpd, p_gpd = scipy_stats.kstest(tail_data_transformed, 'genpareto', args=gpd_params)
    res_['kstest_D'] = D_gpd
    res_['kstest_p'] = p_gpd
    ## Log-Likelihood
    res_['nll'] = neg_log_likelihood(scipy_stats.genpareto, tail_data_transformed, gpd_params)
    ## store results
    res['models']['genpareto'] = res_
    ## ==========================
    ## Fit lognormal distribution
    ## ==========================
    ## Fit log-Normal distribution to the tail data
    logn_params = scipy_stats.lognorm.fit(tail_data_transformed,floc=0.)
    res_ = { 'type':'logn', 'lambda':logn_params[2], 'zeta':logn_params[0] }
    ## Kolmogorov–Smirnov Test
    D_logn, p_logn = scipy_stats.kstest(tail_data_transformed, 'lognorm', args=logn_params)
    res_['kstest_D'] = D_logn
    res_['kstest_p'] = p_logn
    ## Log-Likelihood
    res_['nll'] = neg_log_likelihood(scipy_stats.lognorm, tail_data_transformed, logn_params)
    ## store results
    res['models']['logn'] = res_
    ## ==========================
    ## fit beta distribution in case there is a bound
    ## ==========================
    if bound is not None:
        beta_params = scipy_stats.beta.fit(tail_data_transformed,floc=0.,fscale=bound)
        res_ = { 'type':'beta', 'alpha':beta_params[0], 'beta':beta_params[1], 'a':beta_params[2], 'b':beta_params[3] }
        ## Kolmogorov–Smirnov Test
        D_beta, p_beta = scipy_stats.kstest(tail_data_transformed, 'beta', args=beta_params)
        res_['kstest_D'] = D_beta
        res_['kstest_p'] = p_beta
        ## Log-Likelihood
        res_['nll'] = neg_log_likelihood(scipy_stats.beta, tail_data_transformed, beta_params)
        ## store results
        res['models']['beta'] = res_
    ## ==========================
    ## Select model with best fit
    ## ==========================
    best_fit = None
    best_model = None
    for model_type, model in res['models'].items():
        if best_fit is None:
            best_fit = model['nll']
            best_model = model_type
        else:
            if best_fit > model['nll']:
                best_fit = model['nll']
                best_model = model_type
    res['best_model'] = best_model  ## model with best fit
    res['use_model'] = best_model   ## use the model with best fit
    return res


def _fit_linear_inclined(x_data):
    def neg_log_likelihood(m):
        # Keep inside domain to avoid log of negative numbers
        if not (-1 <= m <= 1):
            return np.inf
        pdf_values = 1 + m * (2 * x_data - 1)
        if np.any(pdf_values <= 0):
            return np.inf
        return -np.sum(np.log(pdf_values))

    result = scipy_opt.minimize_scalar(neg_log_likelihood, bounds=(-1., 1.), method='bounded')
    return result.x



def get_quantiles_from_data(data, p_vec=None, N_points_per_bin=100, data_is_sorted=False, lower_bound=None, upper_bound=None):
    """Extract quantiles from data."""
    res = {}
    ## ===============
    ## Sort data array
    ## ===============
    if data_is_sorted:
        sdata = data
    else:
        sdata = np.sort(data, axis=None)
    N_total = sdata.size
    res['N_total'] = N_total    ## total number of samples considered
    ## ==============
    ## Assemble p_vec
    ## ==============
    if p_vec is None:
        ## ---------------------
        ## select number of bins
        ## ---------------------
        if N_total < N_points_per_bin*10:
            raise NameError(f'ERROR 202504241439: Not enough data. {N_total = }, {N_points_per_bin = }')
        N_bins = int(N_total / N_points_per_bin)
        ## --------------------
        ## assign probabilities
        ## --------------------
        p_vec = discretize_x(x_low=0., x_up=1., x_disc_N=N_bins+1, x_disc_shift=False, x_disc_on_log=False)
    else:
        N_bins = p_vec.size()-1
    if p_vec[0]!=0. or p_vec[-1]!=1.:
        raise NameError(f"ERROR 202504241513: First and last value of 'p_vec' must be 0.0 and 1.0, respectively.")
    if N_bins < 5:
        raise NameError(f'ERROR 202504250716: Not enough bins. {N_bins = }')
    res['N_bins'] = N_bins
    res['p_vec'] = p_vec    ## vector of probabilities (for quantile evaluation)
    ## ==============
    ## Assemble q_vec
    ## ==============
    res['nf_points_per_bin'] = N_total/float(N_bins)  ## (average) number of points per bin
    q_vec = np.empty(N_bins+1)
    N_vec = np.empty(N_bins,dtype=int)
    j_first = None
    j_last  = None
    j_prev = 0
    for i in range(N_bins+1):
        p = p_vec[i]
        if p<=0.:
            if lower_bound is None:
                q_vec[i] = sdata[0] - (sdata[1]-sdata[0])/2
            else:
                q_vec[i] = lower_bound
                if lower_bound>sdata[0]:
                    raise NameError(f"ERROR 202504250831: lower bound is larger than minimum in data. {lower_bound = }, {sdata[0] = }")
            j = 0
        elif p>=1.:
            if upper_bound is None:
                q_vec[i] = sdata[-1] + (sdata[-1]-sdata[-2])/2
            else:
                q_vec[i] = upper_bound
                if upper_bound<sdata[-1]:
                    raise NameError(f"ERROR 202504250832: upper bound is smaller than maimum in data. {upper_bound = }, {sdata[-1] = }")
            j = N_total
        else:
            ## linear interpolation
            nf = p*N_total - 0.5
            j = int(nf)
            nf -= j
            q_vec[i] = sdata[j] + (sdata[j+1]-sdata[j])*nf
            ## remember relevant indices for tail-fitting
            if i==1:
                j_first = j+1
                if nf<1e-12:  ## if quantile matches 'exactly' the relevant data-point it is not part of the tail
                    j_first -= 1
            if i==N_bins-1:
                j_last = j+1
            ## increase j by 1 for assembling N_vec
            j += 1
        ## store number of samples per bin
        if i>0:
            N_vec[i-1] = j-j_prev
            j_prev = j
    ## make sure that only unique values are in q_vec
    for i in range(N_bins):
        if q_vec[i]==q_vec[i+1]:
            raise NameError(f"ERROR 202504250837: quantiles are not unique. {i = }")
    res['q_vec'] = q_vec    ## vector of quantiles (associated with p_vec)
    ## make sure total number of samples is consistent
    if sum(N_vec)!=N_total:
        raise NameError(f"ERROR 202504250841: {sum(N_vec)}, {N_vec = }")
    res['N_vec'] = N_vec
    ## ==============================================
    ## fit beta distribution to individual bins
    ## ==============================================
    bin_rvbeta_params = np.ones(N_bins*2)*-1.
    bin_rvlinear_params = np.zeros(N_bins)
    for i in range(N_bins):
        ## transform bin-data to [0.,1.] values
        x_low = q_vec[i]
        x_up  = q_vec[i+1]
        if (x_up-x_low)<1e-12: ## avoid fit if values are 'almost' equivalent
            continue
        data_bin = data[np.logical_and( (x_low<data), (x_up>data) )]
        data_bin -= x_low
        data_bin /= (x_up-x_low)
        ## fit distribution
        beta_params = scipy_stats.beta.fit(data_bin,floc=0.,fscale=1.)
        bin_rvbeta_params[i*2] = beta_params[0]
        bin_rvbeta_params[i*2+1] = beta_params[1]
        ## fit linear distribution
        bin_rvlinear_params[i] = _fit_linear_inclined(data_bin)
    res['bin_rvbeta_params'] = bin_rvbeta_params
    res['bin_rvlinear_params'] = bin_rvlinear_params
    ## ==============================================
    ## fit upper tail
    ## ==============================================
    Q_tail = q_vec[-2]
    tail_data_transformed = sdata[j_last:] - Q_tail
    if upper_bound is None:
        bound_transformed = upper_bound
    else:
        bound_transformed = upper_bound - Q_tail
    res['tail_upper'] = fit_tail_to_data(tail_data_transformed,bound_transformed)
    ## ==============================================
    ## fit lower tail
    ## ==============================================
    Q_tail = q_vec[1]
    tail_data_transformed = Q_tail - sdata[:j_first]
    if lower_bound is None:
        bound_transformed = lower_bound
    else:
        bound_transformed = Q_tail - lower_bound
    res['tail_lower'] = fit_tail_to_data(tail_data_transformed,bound_transformed)
    ## ==============================================
    ## return
    ## ==============================================
    res['type'] = 'quantiles'
    res['interpol'] = "uniform"
    res['use_tail_fit'] = True
    return res




