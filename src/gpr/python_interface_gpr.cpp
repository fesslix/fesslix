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
#include "flxobjrandom.h"

void register_dataBox_post_processors_gpr()
{
  flxDataBox::postProc_map["akmcs"] = [](const py::dict& cfg, flxDataBox& dBox) { return new post_proc_akmcs(cfg, dBox); };
}

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
: gp_ptr(nullptr), mem_managed(true)
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

const tdouble flxPyGP::condition_on(py::array_t<tdouble> dm_in, py::array_t<tdouble> dv_out, const bool init_pvec, const bool opt_noise)
{
    return gp_ptr->register_observation(flxGP_data_from_Py(dm_in,dv_out),init_pvec,opt_noise);
}

void flxPyGP::noise_white(const tdouble noise_sd)
{
    if (noise_sd<ZERO) {
        throw FlxException("flxPyGP::noise_white", "Standard deviation of noise must not be negative.");
    }
    gp_ptr->register_noise(noise_sd);
}

const tdouble flxPyGP::optimize(const tuint itermax, const bool opt_noise)
{
    return gp_ptr->optimize(itermax, opt_noise);
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
// AK-MCS
// #################################################################################

flxGP_AKMCS::flxGP_AKMCS(py::dict config)
: RndBox(nullptr), dBox_ptr(nullptr), lsf(nullptr), gp_ptr(nullptr), gp_mci(nullptr),
    iterMax(500), NmaxSur(10000000), Nsmpls(1000000), last_state(akmcs_status::undefined), err_thresh(0.3), N_model_calls(0)
{
    try {
        // ====================================================
        // initialize the sampler
        // ====================================================
        if (config.contains("sampler")) {
            py::object py_sampler_class = py::module_::import("fesslix.core").attr("sampler");  // Get the actual class from the original module
            py::object tmp_obj = config["sampler"];
            if (py::isinstance(tmp_obj, py_sampler_class)) {
                gp_obj = tmp_obj;
                flxPySampler &sampler_obj_ref = gp_obj.cast<flxPySampler &>();
                RndBox = sampler_obj_ref.get_ptr_RndBox();
            } else {
                throw FlxException_NeglectInInteractive("flxGP_AKMCS::flxGP_AKMCS_s01", "'sampler' in 'config' is not of type <flx.sampler>.");
            }
        } else {
            throw FlxException_NeglectInInteractive("flxGP_AKMCS::flxGP_AKMCS_s02", "'sampler' not found in 'config'.");
        }
        const tuint M = RndBox->get_NRV();
        // ====================================================
        // initialize the dataBox
        // ====================================================
        if (config.contains("data_box")) {
            py::object py_dataBox_class = py::module_::import("fesslix.core").attr("dataBox");  // Get the actual class from the original module
            py::object tmp_obj = config["data_box"];
            if (py::isinstance(tmp_obj, py_dataBox_class)) {
                dataBox_obj = tmp_obj;
                flxDataBox &dataBox_obj_ref = dataBox_obj.cast<flxDataBox &>();
                // consistency checks
                    dataBox_obj_ref.ensure_M_in(M);
                    dataBox_obj_ref.ensure_M_out(1);
                dBox_ptr = &dataBox_obj_ref;
            } else {
                throw FlxException_NeglectInInteractive("flxGP_AKMCS::flxGP_AKMCS_d01", "'data_box' in 'config' is not of type <flx.dataBox>.");
            }
        }
        // ====================================================
        // limit-state function
        // ====================================================
        lsf = parse_py_para("lsf",config,true);
        // ====================================================
        // initialize the Gaussian process
        // ====================================================
        if (config.contains("gp")) {
            py::object tmp_obj = config["gp"];
            if (py::isinstance<flxPyGP>(tmp_obj)) {
                gp_obj = tmp_obj;
                flxPyGP &gp_obj_ref = gp_obj.cast<flxPyGP &>();
                gp_ptr = gp_obj_ref.get_gp_ptr();
            } else {
                throw FlxException_NeglectInInteractive("flxGP_AKMCS::flxGP_AKMCS_g01", "'gp' in 'config' is not of type <flx.gpr.gp>.");
            }
        } else {
            // mean function
                flxGP_mean_base* gp_mean = new flxGP_mean_0(M);
            // kernel
                iVec tVec(M);
                tVec = 0;
                flxGP_kernel_base* gp_kernel = new flxGP_kernel_auto(tVec);
            gp_ptr = new flxGPProj("gp4akmcs",M,gp_mean,gp_kernel,true);
        }
        if (gp_ptr->get_Ndim() != M) {
            throw FlxException("flxGP_AKMCS::flxGP_AKMCS_g02", "Specified Gaussian process has wrong dimension.");
        }
        // ====================================================
        // initialize AK-MCS
        // ====================================================
        const tuint seed_id = parse_py_para_as_tuint("seed",config,false,0);
        const tuint N_RNG_init = parse_py_para_as_tuint("N_RNG_init",config,false,0);
        const tuint N_reserve = parse_py_para_as_tuintNo0("N_reserve",config,false,10000);
        gp_mci = new flxGP_MCI(*gp_ptr,N_reserve,seed_id,N_RNG_init);
        iterMax = parse_py_para_as_tuintNo0("itermax",config,false,500);
        NmaxSur = parse_py_para_as_tulong("NmaxSur",config,false,10000000);
        Nsmpls = parse_py_para_as_tulong("Nsmpls",config,false,1000000);
        if (Nsmpls<1000) Nsmpls = 1000;
        const tdouble err_thresh_ = parse_py_para_as_float("err_thresh",config,false,0.3);
        if (err_thresh_>ZERO) {
            err_thresh = err_thresh_;
        }

    } catch (...) {
        free_mem();
        throw;
    }

}

flxGP_AKMCS::~flxGP_AKMCS()
{
    free_mem();
}

void flxGP_AKMCS::free_mem()
{
    if (gp_obj.is_none()) {     // delete pointer only if GP-process is not handled by Python
        if (gp_ptr) {
            delete gp_ptr;
            gp_ptr = nullptr;
        }
    }
    if (lsf) {
        delete lsf;
        lsf = nullptr;
    }
}

const bool flxGP_AKMCS::eval_model(flxVec& y_vec)
{
    // make sure point has not been considered before
        if (gp_mci->is_point_unique(y_vec)==false) {
            return false;
        }
    // evaluate the model
        RndBox->set_smp(y_vec);
        const tdouble lsf_val = lsf->calc();
        ++N_model_calls;
    // register the call
        gp_mci->register_sample(lsf_val,y_vec);
    // register call in dataBox
        if (dBox_ptr) {
            dBox_ptr->vec_in = y_vec;
            dBox_ptr->vec_out[0] = lsf_val;
            dBox_ptr->append_data();
        }
    return true;
}

void flxGP_AKMCS::initialize_with_LHS(tuint N)
{
    // generate the set of samples
        const tuint M = RndBox->get_NRV();
        if (N==0) N = 2*M;
        flxVec lh_samples(N*M);
        gp_mci->assemble_lh_samples(lh_samples);
    // run the actual model for each sample
        for (tuint i=0;i<N;++i) {
            flxVec tv(lh_samples.get_tmp_vptr()+i*M,M,false,false);
            eval_model(tv);
        }
    last_state = akmcs_status::undefined;
}

void flxGP_AKMCS::initialize_with_sample(const flxVec& y_vec, const tdouble lsf_val)
{
    // make sure point has not been considered before
        if (gp_mci->is_point_unique(y_vec)==false) {
            throw FlxException("flxGP_AKMCS::initialize_with_sample_01", "Sample to add is not unique.");
        }
    gp_mci->register_sample(lsf_val,y_vec);
    last_state = akmcs_status::undefined;
}

akmcs_status flxGP_AKMCS::simulate()
{
    // perform action based on recommendation form last call of simulate()
        switch (last_state) {
            case akmcs_status::undefined:
            {
                // condition GP on data and set initial parameter values
                    gp_mci->condition_on_data(true,false);
                    gp_mci->optimize_gp_para(iterMax);
                last_state = akmcs_status::defined;
                break;
            }
            case akmcs_status::defined:
                break;
            case akmcs_status::evalLSF:
            {
                const tuint M = RndBox->get_NRV();
                flxVec uvec(M);
                gp_mci->get_next_point(uvec);
                if (!eval_model(uvec)) {
                    throw FlxException_Crude("flxGP_AKMCS::simulate_01");
                }
                gp_mci->condition_on_data(false,false);
                gp_mci->optimize_gp_para(iterMax);
                last_state = akmcs_status::defined;
                break;
            }
            case akmcs_status::increase_N_surrogate:
            {
                Nsmpls *= 2;
                if (NmaxSur>0 && Nsmpls>NmaxSur) Nsmpls = NmaxSur;
                last_state = akmcs_status::defined;
                break;
            }
            case akmcs_status::stop_success:
            case akmcs_status::stop_iterLimit:
                return last_state;
            default:
                throw FlxException_Crude("flxGP_AKMCS::simulate_02");
        }
        if (last_state!=akmcs_status::defined) {
            throw FlxException_Crude("flxGP_AKMCS::simulate_03");
        }
    // perform the sampling
        tdouble err;
        int proposed_action_id;
        res = gp_mci->simulate_GP_mci(Nsmpls,err,proposed_action_id);
    // decide on what to do
        if (err<=err_thresh) {
            last_state = akmcs_status::stop_success;
            return last_state;
        }
        if (NmaxSur>0 && Nsmpls>=NmaxSur) {
            last_state = akmcs_status::stop_iterLimit;
            return last_state;
        }
        if (proposed_action_id==1) {
            last_state = akmcs_status::increase_N_surrogate;
            return last_state;
        }
        last_state = akmcs_status::evalLSF;
        return last_state;
}

flxPyGP flxGP_AKMCS::get_GP()
{
    return flxPyGP(gp_ptr);
}

const tuint flxGP_AKMCS::get_N_model_calls(const bool only_from_current_run) const
{
    if (only_from_current_run) {
        return N_model_calls;
    } else {
        return gp_mci->get_N_model_calls();
    }
}

post_proc_akmcs::post_proc_akmcs(py::dict config, flxDataBox& dBox)
: akmcs_ptr(nullptr)
{
    if (config.contains("akmcs")) {
        py::object tmp_obj = config["akmcs"];
        if (py::isinstance<flxGP_AKMCS>(tmp_obj)) {
            akmcs_obj = tmp_obj;
            flxGP_AKMCS &akmcs_obj_ref = akmcs_obj.cast<flxGP_AKMCS &>();
            akmcs_ptr = &akmcs_obj_ref;
        } else {
            throw FlxException_NeglectInInteractive("post_proc_akmcs::post_proc_akmcs_01", "'akmcs' in 'config' is not of type <flx.gpr.akmcs>.");
        }
    } else {
        throw FlxException_NeglectInInteractive("post_proc_akmcs::post_proc_akmcs_99", "'akmcs' not found in 'config'.");
    }
}

void post_proc_akmcs::append_data(const flxVec& vec_full)
{
    if (akmcs_ptr==nullptr) {
        throw FlxException_Crude("post_proc_akmcs::append_data");
    }
    const tdouble lsf_val = vec_full[0];
    flxVec y_vec(vec_full.get_tmp_vptr_const()+1,vec_full.get_N()-1);
    akmcs_ptr->initialize_with_sample(y_vec,lsf_val);
}

py::dict post_proc_akmcs::eval()
{
    py::dict res;
    res["akmcs"] = akmcs_obj;
    return res;
}

// #################################################################################
// Expose interface to Python
// #################################################################################

PYBIND11_MODULE(gpr, m) {
    register_dataBox_post_processors_gpr();
    // ====================================================
    // Gaussian processes
    // ====================================================
        py::class_<flxPyGP>(m, "gp")
            .def(py::init<py::dict>())
            .def("get_name", &flxPyGP::get_name, "get name of gaussian process")
            .def("get_descr", &flxPyGP::get_descr, "get description of gaussian process")
            .def("condition_on", &flxPyGP::condition_on, pybind11::arg("dm_in"), pybind11::arg("dv_out"), pybind11::arg("init_pvec") = true, pybind11::arg("opt_noise") = true, "conditions the tp-model on an observation")
            .def("noise_white", &flxPyGP::noise_white, "model error of the gp-model")
            .def("optimize", &flxPyGP::optimize, pybind11::arg("itermax") = 500, pybind11::arg("opt_noise") = false, "optimizes the parameters of the gp-model")
            .def("predict", &flxPyGP::predict, pybind11::arg("x_vec"), pybind11::arg("type"), pybind11::arg("predict_noise") = false, "evaluate/predict quantities of the gp-model conditioned on the observations")
            .def("unassemble", &flxPyGP::unassemble, "forces the gp-model to re-assemble the covariance matrix on the next call")
            .def("info", &flxPyGP::info, "return a dict with properties of the gp-model")
            ;
    // ====================================================
    // AK-MCS
    // ====================================================
        py::enum_<akmcs_status>(m, "akmcs_status")
            .value("undefined", akmcs_status::undefined)
            .value("defined", akmcs_status::defined)
            .value("evalLSF", akmcs_status::evalLSF)
            .value("increase_N_surrogate", akmcs_status::increase_N_surrogate)
            .value("stop_success", akmcs_status::stop_success)
            .value("stop_iterLimit", akmcs_status::stop_iterLimit)
            .export_values();
        py::class_<flxGP_AKMCS>(m, "akmcs")
            .def(py::init<py::dict>())
            .def("initialize_with_LHS", &flxGP_AKMCS::initialize_with_LHS, pybind11::arg("N")=0, "Initialize AK-MCS using N Latin-hypercube samples.")
            .def("simulate", &flxGP_AKMCS::simulate, "Perform a single simulation step using the surrogate model.")
            .def("get_GP", &flxGP_AKMCS::get_GP, "Retrieve a reference to the internal Gaussian process.")
            .def("get_N_model_calls", &flxGP_AKMCS::get_N_model_calls, pybind11::arg("only_from_current_run") = true, "Retrieve total number of calls of the actual limit-state function.")
            .def_readwrite("res", &flxGP_AKMCS::res, "The result dictionary of 'simulate()'.");
            ;
}


