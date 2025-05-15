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
    bool mem_managed;
    std::string descr;

  public:
    flxPyGP(py::dict config);
    flxPyGP(flxGPProj_base* gp_ptr, const bool mem_managed=false);
    flxPyGP() = delete;
    flxPyGP(flxPyGP& rhs);
    flxPyGP(flxPyGP&& rhs);
    ~flxPyGP();

    flxPyGP& operator=(const flxPyGP& rhs) = delete;

    const std::string get_name() const;
    const std::string& get_descr() const;

    void condition_on(py::array_t<tdouble> dm_in, py::array_t<tdouble> dv_out, const bool init_pvec, const bool opt_noise);
    void noise_white(const tdouble noise_sd);
    const tdouble optimize(const tuint itermax, const bool log_opt);
    py::object predict(py::array_t<tdouble> arr, const std::string& type, const bool predict_noise);
    void unassemble();
    py::dict info();

};



