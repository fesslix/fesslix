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

#include "python_interface_gpr.h"

#include "flxgp_kernel.h"

#include "flxparse.h"


// #################################################################################
// Gaussian processes
// #################################################################################

flxGP_data_base* flxGP_data_from_Py::copy() const
{
    return new flxGP_data_from_Py(dm_py, do_py);
}

const tdouble* flxGP_data_from_Py::get_data_ptr_mtx(tuint& N_obsv, const tuint N_dim)
{
    // Request a buffer info from the numpy array (zero-copy, safe)
    py::buffer_info buf = dm_py.request();

    if (buf.ndim != 2) {
        throw FlxException("flxGP_data_from_Py::get_data_ptr_mtx_01", "Expected a 2D array");
    }

    tdouble* ptr = static_cast<tdouble*>(buf.ptr); // raw pointer to the data

    size_t rows = buf.shape[0];
    size_t cols = buf.shape[1];

    if (cols!=N_dim) {
        throw FlxException("flxGP_data_from_Py::get_data_ptr_mtx_02", "Dimension of matrix (number of columns) does not match dimension of gaussian process.");
    }

    N_obsv = rows;


    return ptr;
}

const tdouble* flxGP_data_from_Py::get_data_ptr_vec(const tuint N_obsv)
{
    // Access the input data as a raw pointer
    py::buffer_info buf_info = do_py.request();
    tdouble* ptr = static_cast<tdouble*>(buf_info.ptr);

    // Get the size of the input array
    size_t size = buf_info.size;

    if (size!=N_obsv) {
        throw FlxException("flxGP_data_from_Py::get_data_ptr_vec", "Mismatch in number of observations.");
    }
    return ptr;
}


flxPyGP::flxPyGP(py::dict config)
: gp_ptr(nullptr), mem_managed(false)
{
    flxGP_mean_base* gp_mean = nullptr;
        flxGP_kernel_base* gp_kernel = nullptr;
    // retrieve type
        const std::string gp_type = parse_str_as_word(parse_py_para_as_string("type",config,false,"singleGP"),true);
    // retrieve name
        std::string gp_name = parse_str_as_word(parse_py_para_as_string("name",config,false,"name_unspecified"),true);
    // retrieve description
        descr = parse_py_para_as_string("descr",config,false);
    // retrieve dimension
        const tuint Ndim = parse_py_para_as_tuintNo0("Ndim",config,true);
    // other parameters
        const bool useLSE = parse_py_para_as_bool("useLSE", config, false, true);
    if (gp_type=="singlegp") {
        try {
        // retrieve mean
            const std::string mean_type = parse_py_para_as_word("mean_type",config,false,true,false,false,"zero");
            if (mean_type=="zero") {
                gp_mean = new flxGP_mean_0(Ndim);
            } else if (mean_type=="const") {
                const tdouble val = parse_py_para_as_float("mean_value",config,true);
                gp_mean = new flxGP_mean_const(val, Ndim);
            } else if (mean_type=="universal") {
                const tuint tid = parse_py_para_as_tuint("mean_polyo",config,true);
                gp_mean = new flxGP_mean_universal(tid,Ndim);
            } else if (mean_type=="ordinary") {
                FlxFunction* mean_f = parse_py_para("mean_fun", config, true, Ndim);
                gp_mean = new flxGP_mean_ordinary(mean_f,Ndim);
            } else {
              throw FlxException("flxPyGP::flxPyGP_19", "Unknown type specified in 'mean_type': " + mean_type);
            }
        // retrieve kernel
            if (config.contains("kernel_lst")) {
                std::vector<std::string> kernel_type_lst;
                parse_py_para_as_word_lst(kernel_type_lst,"kernel_lst",config,true,true);
                if (kernel_type_lst.size()!=Ndim) {
                    throw FlxException("flxPyGP::flxPyGP_20","Size of the array of kernel types must match the number of specified dimensions.");
                }
                gp_kernel = new flxGP_kernel_auto(kernel_type_lst);
            } else {
                FlxFunction* kernel_f = parse_py_para("kernel_fun", config, true, 2*Ndim);
                gp_kernel = new flxGP_kernel_user(kernel_f,Ndim);
            }
        // create Gaussian process
            gp_ptr = new flxGPProj(gp_name,Ndim,gp_mean,gp_kernel,useLSE);
        } catch (...) {
            if (gp_mean) delete gp_mean;
            if (gp_kernel) delete gp_kernel;
            throw;
        }
    } else if (gp_type=="mavggp") {
        // retrieve maximum number of iterations
            const tuint itermax = parse_py_para_as_tuintNo0("itermax",config,false,500);
        // create Gaussian process
            gp_ptr = new flxGP_avgModel(gp_name,Ndim,itermax,useLSE);
    } else {
        throw FlxException("flxPyGP::flxPyGP", "Unknown type '" + gp_type + "' specified in type.");
    }
}

flxPyGP::flxPyGP(flxGPProj_base* gp_ptr, const bool mem_managed)
: gp_ptr(gp_ptr), mem_managed(mem_managed)
{
}

flxPyGP::flxPyGP(flxPyGP& rhs)
: gp_ptr(rhs.gp_ptr), mem_managed(rhs.mem_managed)
{
    if (mem_managed) {
        throw FlxException_NotImplemented("flxPyGP::flxPyGP(flxPyGP& rhs)");
    }
}

flxPyGP::flxPyGP(flxPyGP&& rhs)
: gp_ptr(rhs.gp_ptr), mem_managed(rhs.mem_managed)
{
    rhs.gp_ptr = nullptr;
    rhs.mem_managed = false;
}

flxPyGP::~flxPyGP()
{
    if (gp_ptr) {
        if (mem_managed) {
            delete gp_ptr;
        }
    }
}

const std::string flxPyGP::get_name() const
{
    return gp_ptr->get_name();
}

const std::string& flxPyGP::get_descr() const
{
    return descr;
}

void flxPyGP::condition_on(py::array_t<tdouble> dm_in, py::array_t<tdouble> dv_out, const bool init_pvec, const bool opt_noise)
{
    gp_ptr->register_observation(flxGP_data_from_Py(dm_in,dv_out),init_pvec,opt_noise,GlobalVar.slogcout(3));
}

void flxPyGP::noise_white(const tdouble noise_sd)
{
    if (noise_sd<ZERO) {
        throw FlxException("flxPyGP::noise_white", "Standard deviation of noise must not be negative.");
    }
    gp_ptr->register_noise(noise_sd);
}

const tdouble flxPyGP::optimize(const tuint itermax, const bool log_opt)
{
    std::ostream& ostrm = log_opt?(GlobalVar.slogcout(3)):(GlobalVar.slog_dummy());
    return gp_ptr->optimize(itermax,ostrm);
}

py::object flxPyGP::predict(py::array_t<tdouble> arr, const std::string& type, const bool predict_noise)
{
    // Access the input data as a raw pointer
    py::buffer_info buf_info = arr.request();
    tdouble* input_ptr = static_cast<tdouble*>(buf_info.ptr);
    // Get the size of the input array
    size_t size = buf_info.size;
    // define as flxVec
    flxVec x_vec(input_ptr,size,false);

    if (type=="mean_sd") {
      tdouble res_mean, res_var;
      gp_ptr->predict_mean_var(x_vec,predict_noise,res_mean,res_var);
      return py::make_tuple(res_mean, sqrt(res_var));
    }
    if (type=="mean") {
      tdouble res_mean, res_var;
      gp_ptr->predict_mean_var(x_vec,predict_noise,res_mean,res_var);
      return py::float_(res_mean);
    }
    if (type=="sd") {
      tdouble res_mean, res_var;
      gp_ptr->predict_mean_var(x_vec,predict_noise,res_mean,res_var);
      return py::float_(sqrt(res_var));
    }
    if (type=="var") {
      tdouble res_mean, res_var;
      gp_ptr->predict_mean_var(x_vec,predict_noise,res_mean,res_var);
      return py::float_(res_var);
    }
    if (type=="trend") {
      return py::float_(gp_ptr->eval_trend(x_vec,predict_noise));
    }
    if (type=="kernel") {
      return py::float_(gp_ptr->eval_kernel(x_vec,predict_noise));
    }
    throw FlxException_NeglectInInteractive("flxPyGP::predict", "Unknown type-keyword '" + type + "'.");
}

void flxPyGP::unassemble()
{
    gp_ptr->unassemble();
}

py::dict flxPyGP::info()
{
    return gp_ptr->info();
}


// #################################################################################
// only for debugging purposes
// #################################################################################

double add(double a, double b) {
    return a + b;
}


// #################################################################################
// Expose interface to Python
// #################################################################################

PYBIND11_MODULE(gpr, m) {
    // ====================================================
    // Gaussian processes
    // ====================================================
        py::class_<flxPyGP>(m, "gp")
            .def(py::init<py::dict>())
            .def("get_name", &flxPyGP::get_name, "get name of gaussian process")
            .def("get_descr", &flxPyGP::get_descr, "get description of gaussian process")
            .def("condition_on", &flxPyGP::condition_on, pybind11::arg("dm_in"), pybind11::arg("dv_out"), pybind11::arg("init_pvec") = true, pybind11::arg("opt_noise") = true, "conditions the tp-model on an observation")
            .def("noise_white", &flxPyGP::noise_white, "model error of the gp-model")
            .def("optimize", &flxPyGP::optimize, pybind11::arg("itermax") = 500, pybind11::arg("log_opt") = false, "optimizes the parameters of the gp-model")
            .def("predict", &flxPyGP::predict, pybind11::arg("x_vec"), pybind11::arg("type"), pybind11::arg("predict_noise") = false, "evaluate/predict quantities of the gp-model conditioned on the observations")
            .def("unassemble", &flxPyGP::unassemble, "forces the gp-model to re-assemble the covariance matrix on the next call")
            .def("info", &flxPyGP::info, "return a dict with properties of the gp-model")
            ;
    // ====================================================
    // only for debugging purposes (TODO remove at some point)
    // ====================================================
        m.def("add", &add, "A function that adds two numbers");
        m.attr("the_answer") = 42;
}


