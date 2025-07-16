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


#include "flxgp_kernel.h"
#include "flxgp_relmeth.h"
#include "flxrbrv.h"
#include "flxobjrandom.h"

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>  // For NumPy support

#define PYBIND11_DETAILED_ERROR_MESSAGES

namespace py = pybind11;



// #################################################################################
// Gaussian processes
// #################################################################################

class flxGP_data_from_Py : public flxGP_data_base {
  protected:
    py::array_t<tdouble> dm_py;
    py::array_t<tdouble> do_py;
  public:
    flxGP_data_from_Py(py::array_t<tdouble> dm_py, py::array_t<tdouble> do_py) : dm_py(dm_py), do_py(do_py) {}

    virtual flxGP_data_base* copy() const;
    virtual const tdouble* get_data_ptr_mtx(tuint &N_obsv, const tuint N_dim);
    virtual const tdouble* get_data_ptr_vec(const tuint N_obsv);
};

class PYBIND11_EXPORT flxPyGP {
  private:
    flxGPProj_base* gp_ptr;
    bool mem_managed;       // true, if memory is managed by this class
    std::string descr;

  public:
    flxPyGP(py::dict config);
    flxPyGP(flxGPProj_base* gp_ptr, const bool mem_managed=false);
    flxPyGP() = delete;
    flxPyGP(flxPyGP& rhs);
    flxPyGP(flxPyGP&& rhs);
    ~flxPyGP();

    flxPyGP& operator=(const flxPyGP& rhs) = delete;
    flxGPProj_base* get_gp_ptr() { return gp_ptr; }

    const std::string get_name() const;
    const std::string& get_descr() const;

    const tdouble condition_on(py::array_t<tdouble> dm_in, py::array_t<tdouble> dv_out, const bool init_pvec, const bool opt_noise);
    void noise_white(const tdouble noise_sd);
    const tdouble optimize(const tuint itermax, const bool opt_noise);
    py::object predict(py::array_t<tdouble> arr, const std::string& type, const bool predict_noise);
    void unassemble();
    py::dict info();

};


// #################################################################################
// AK-MCS
// #################################################################################

enum class akmcs_status {
    undefined,            // initial state -> gp is unconditioned
    defined,              // internal state, nothing to do
    evalLSF,              // requires a model call to worst point
    increase_N_surrogate, // requires an increase of surrogate samples
    stop_success,         // stop is recommended
    stop_iterLimit        // maximum number of surrogate samples is exceeded
};

class PYBIND11_EXPORT flxGP_AKMCS {
  private:
    // sampler
      RBRV_constructor* RndBox;
      py::object Sampler_obj;
    // dataBox
      flxDataBox* dBox_ptr;
      py::object dataBox_obj;
    // model
      FlxFunction* lsf;
    // Gaussian process
      flxGPProj_base* gp_ptr;   // the Gaussian process used to approximate the limit-state function
      py::object gp_obj;    // if gp_obj is not None, memory of gp is not managed by this class!
    // AK-MCS handler
      flxGP_MCI* gp_mci;    // the ak-mcs handler
      tuint iterMax;
    // internal attributes
      // maximum number of samples to use in surrogate model (if 0: infinite)
        tulong NmaxSur;
      // current number of surrogate samples
        tulong Nsmpls;
      // last status after calling simulate()
        akmcs_status last_state;
      // threshold for stopping criterion
        tdouble err_thresh;
      // number of LSF-calls in the current instance
        tuint N_model_calls;

    void free_mem();
    const bool eval_model(flxVec& y_vec);
  public:
    // results of run of surrogate model
      py::dict res;

    flxGP_AKMCS(py::dict config);
    flxGP_AKMCS() = delete;
    flxGP_AKMCS(flxGP_AKMCS& rhs);
    flxGP_AKMCS(flxGP_AKMCS&& rhs);
    ~flxGP_AKMCS();

    flxGP_AKMCS& operator=(const flxGP_AKMCS& rhs) = delete;

    void initialize_with_LHS(tuint N);
    akmcs_status simulate();
    flxPyGP get_GP();

    const tuint get_N_model_calls(const bool only_from_current_run) const;

};




