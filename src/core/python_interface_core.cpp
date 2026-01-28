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

#include "python_interface_core.h"
#include "fesslix.h"
#include "flxfunction_data.h"
#include "flxrbrv_rvs_read.h"
#include "flxstringfun_fun.h"
#include "flxobjrbrv.h"
#include "flxobjrandom.h"
#include "flxphys.h"


#include <iostream>


void finalize_call()
{
    std::cout << std::flush;
}

// #################################################################################
// 'global' functions
// #################################################################################

std::string Double2String(tdouble a) {
    return GlobalVar.Double2String(a);
}

void print_info()
{
    std::ostream& lout = GlobalVar.slogcout(3);
    // check state of Engine
        lout << "Engine: ";
        if (get_FlxEngine_ptr()) {
            lout << "up and running";
        } else {
            lout << "NOT running";
        }
        lout << std::endl;
    // log information
    fesslix_logInfo(lout);
    lout.flush();
}


// #################################################################################
// load configuration
// #################################################################################

void set_exe_dir(const std::string& exe_dir)
{
    GlobalVar.set_exe_dir(exe_dir);
}

void process_config(py::dict config) {
    // ======================================================
    // General options
    // ======================================================
    if (config.contains("flx")) {
        if (py::isinstance<py::dict>(config["flx"])) {
            py::dict pdl1 = config["flx"].cast<py::dict>();
            // ----------------------------------------------
            // prgbar
            // ----------------------------------------------
            if (pdl1.contains("prgbar")) {
                if (py::isinstance<py::bool_>(pdl1["prgbar"])) {
                    GlobalVar.prgBar = pdl1["prgbar"].cast<bool>();
                } else {
                    GlobalVar.slogcout(2) << "ERROR [process_config::flx::prgbar]: expected an entry of type <bool>" << std::endl;
                }
            }
            // ----------------------------------------------
            // leak-check
            // ----------------------------------------------
            if (pdl1.contains("leak-check")) {
                if (py::isinstance<py::bool_>(pdl1["leak-check"])) {
                    const bool lcheck = pdl1["leak-check"].cast<bool>();
                    if (lcheck) {
                        GlobalVar.set_leak_check_mode();
                    }
                } else {
                    GlobalVar.slogcout(2) << "ERROR [process_config::flx::leak-check]: expected an entry of type <bool>" << std::endl;
                }
            }
        } else {
            GlobalVar.slogcout(2) << "ERROR [process_config::flx]: expected an entry of type <dict>" << std::endl;
        }
    }
    // ======================================================
    // TOL
    // ======================================================
    if (config.contains("TOL")) {
        if (py::isinstance<py::float_>(config["TOL"])) {
            GlobalVar.set_TOL( config["TOL"].cast<tdouble>() );
        } else {
            GlobalVar.slogcout(2) << "ERROR [process_config::TOL]: expected an entry of type <float>, got <" << py::str(py::type::of(config["TOL"])) << ">" << std::endl;
        }
    }
    // ======================================================
    // Floating-point conversion
    // ======================================================
    if (config.contains("fpc")) {
        if (py::isinstance<py::dict>(config["fpc"])) {
            py::dict pdl1 = config["fpc"].cast<py::dict>();
            // ----------------------------------------------
            // prec
            // ----------------------------------------------
            if (pdl1.contains("prec")) {
                if (py::isinstance<py::int_>(pdl1["prec"])) {
                    GlobalVar.Double2String_setPrec(pdl1["prec"].cast<tuint>());
                } else {
                    GlobalVar.slogcout(2) << "ERROR [process_config::fpc::prec]: expected an entry of type <unsigned int>" << std::endl;
                }
            }
            // ----------------------------------------------
            // type
            // ----------------------------------------------
            if (pdl1.contains("type")) {
                if (py::isinstance<py::int_>(pdl1["type"])) {
                    GlobalVar.Double2String_setType(pdl1["type"].cast<tuint>());
                } else {
                    GlobalVar.slogcout(2) << "ERROR [process_config::fpc::type]: expected an entry of type <unsigned int>" << std::endl;
                }
            }
            // ----------------------------------------------
            // bvalu
            // ----------------------------------------------
            if (pdl1.contains("bvalu")) {
                if (py::isinstance<py::float_>(pdl1["bvalu"])) {
                    GlobalVar.Double2String_setBValU( pdl1["bvalu"].cast<tdouble>() );
                } else {
                    GlobalVar.slogcout(2) << "ERROR [process_config::fpc::bvalu]: expected an entry of type <float>" << std::endl;
                }
            }
            // ----------------------------------------------
            // bvall
            // ----------------------------------------------
            if (pdl1.contains("bvall")) {
                if (py::isinstance<py::float_>(pdl1["bvall"])) {
                    GlobalVar.Double2String_setBValL( pdl1["bvall"].cast<tdouble>() );
                } else {
                    GlobalVar.slogcout(2) << "ERROR [process_config::fpc::bvall]: expected an entry of type <float>" << std::endl;
                }
            }
            // ----------------------------------------------
            // del0
            // ----------------------------------------------
            if (pdl1.contains("del0")) {
                if (py::isinstance<py::bool_>(pdl1["del0"])) {
                    GlobalVar.Double2String_setDel0( pdl1["del0"].cast<bool>() );
                } else {
                    GlobalVar.slogcout(2) << "ERROR [process_config::fpc::del0]: expected an entry of type <bool>" << std::endl;
                }
            }
            // ----------------------------------------------
            // delp
            // ----------------------------------------------
            if (pdl1.contains("delp")) {
                if (py::isinstance<py::bool_>(pdl1["delp"])) {
                    GlobalVar.Double2String_setDel0( pdl1["delp"].cast<bool>() );
                } else {
                    GlobalVar.slogcout(2) << "ERROR [process_config::fpc::delp]: expected an entry of type <bool>" << std::endl;
                }
            }
        } else {
            GlobalVar.slogcout(2) << "ERROR [process_config::fpc]: expected an entry of type <dict>" << std::endl;
        }
    }
    // ======================================================
    // Gauss points
    // ======================================================
    if (config.contains("gauss")) {
        if (py::isinstance<py::dict>(config["gauss"])) {
            py::dict pdl1 = config["gauss"].cast<py::dict>();
            // ----------------------------------------------
            // file
            // ----------------------------------------------
            if (pdl1.contains("file")) {
                if (py::isinstance<py::str>(pdl1["file"])) {
                    fesslix_gaussFile = pdl1["file"].cast<std::string>();
                } else {
                    GlobalVar.slogcout(2) << "ERROR [process_config::gauss::file]: expected an entry of type <string>" << std::endl;
                }
            }
            // ----------------------------------------------
            // maxnumb
            // ----------------------------------------------
            if (pdl1.contains("maxnumb")) {
                if (py::isinstance<py::int_>(pdl1["maxnumb"])) {
                    GaussIntegration::GaussPointMaxArraySize = pdl1["maxnumb"].cast<tuint>();
                } else {
                    GlobalVar.slogcout(2) << "ERROR [process_config::gauss::maxnumb]: expected an entry of type <unsigned int>" << std::endl;
                }
            }
            // ----------------------------------------------
            // numb
            // ----------------------------------------------
            if (pdl1.contains("numb")) {
                if (py::isinstance<py::int_>(pdl1["numb"])) {
                    GaussIntegration::GaussPointArraySize = pdl1["numb"].cast<tuint>();
                } else {
                    GlobalVar.slogcout(2) << "ERROR [process_config::gauss::numb]: expected an entry of type <unsigned int>" << std::endl;
                }
            }
        } else {
            GlobalVar.slogcout(2) << "ERROR [process_config::gauss]: expected an entry of type <dict>" << std::endl;
        }
    }
    // ======================================================
    // Control the logging behavior
    // ======================================================
    if (config.contains("log")) {
        if (py::isinstance<py::dict>(config["log"])) {
            py::dict pdl1 = config["log"].cast<py::dict>();
            // ----------------------------------------------
            // input
            // ----------------------------------------------
            if (pdl1.contains("input")) {
                if (py::isinstance<py::bool_>(pdl1["input"])) {
                    GlobalVar.prelog_activated( pdl1["input"].cast<bool>() );
                } else {
                    GlobalVar.slogcout(2) << "ERROR [process_config::log::input]: expected an entry of type <bool>" << std::endl;
                }
            }
            // ----------------------------------------------
            // file
            // ----------------------------------------------
            if (pdl1.contains("file")) {
                if (py::isinstance<py::str>(pdl1["file"])) {
                    fesslix_logFile = pdl1["file"].cast<std::string>();
                } else {
                    GlobalVar.slogcout(2) << "ERROR [process_config::log::file]: expected an entry of type <string>" << std::endl;
                }
            }
            // ----------------------------------------------
            // level
            // ----------------------------------------------
            if (pdl1.contains("level")) {
                if (py::isinstance<py::int_>(pdl1["level"])) {
                    GlobalVar.logLevel = pdl1["level"].cast<tuint>();
                } else {
                    GlobalVar.slogcout(2) << "ERROR [process_config::log::level]: expected an entry of type <integer>" << std::endl;
                }
            }
            // ----------------------------------------------
            // output
            // ----------------------------------------------
            if (pdl1.contains("output")) {
                if (py::isinstance<py::bool_>(pdl1["output"])) {
                    fesslix_logOutput = pdl1["output"].cast<bool>();
                } else {
                    GlobalVar.slogcout(2) << "ERROR [process_config::log::output]: expected an entry of type <bool>" << std::endl;
                }
            }
            // ----------------------------------------------
            // logTrunc
            // ----------------------------------------------
            if (pdl1.contains("logTrunc")) {
                if (py::isinstance<py::bool_>(pdl1["logTrunc"])) {
                    fesslix_logTrunc = pdl1["logTrunc"].cast<bool>();
                } else {
                    GlobalVar.slogcout(2) << "ERROR [process_config::log::logTrunc]: expected an entry of type <bool>" << std::endl;
                }
            }
        } else {
            GlobalVar.slogcout(2) << "ERROR [process_config::log]: expected an entry of type <dict>" << std::endl;
        }
    }
    // ======================================================
    // MT19937
    // ======================================================
    if (config.contains("MT19937")) {
        if (py::isinstance<py::dict>(config["MT19937"])) {
            py::dict mt19937 = config["MT19937"].cast<py::dict>();

            // ----------------------------------------------
            // init
            // ----------------------------------------------
            if (mt19937.contains("init")) {
                if (py::isinstance<py::dict>(mt19937["init"])) {
                    py::dict init = mt19937["init"].cast<py::dict>();
                    // --------------------------------------
                    // calls
                    // --------------------------------------
                    if (init.contains("calls")) {
                        if (py::isinstance<py::int_>(init["calls"])) {
                            GlobalVar.MT19937_init_calls = init["calls"].cast<tuint>();
                        } else {
                            GlobalVar.slogcout(2) << "ERROR [process_config::MT19937::init::calls]: expected an entry of type <unsigned int>" << std::endl;
                        }
                    }
                    // --------------------------------------
                    // rand
                    // --------------------------------------
                    if (init.contains("rand")) {
                        if (py::isinstance<py::bool_>(init["rand"])) {
                            GlobalVar.MT19937_init_RAND = init["rand"].cast<bool>();
                        } else {
                            GlobalVar.slogcout(2) << "ERROR [process_config::MT19937::init::rand]: expected an entry of type <bool>" << std::endl;
                        }
                    }
                    // --------------------------------------
                    // seed
                    // --------------------------------------
                    if (init.contains("seed")) {
                        if (py::isinstance<py::bool_>(init["seed"])) {
                            GlobalVar.MT19937_init_seed = init["seed"].cast<bool>();
                        } else {
                            GlobalVar.slogcout(2) << "ERROR [process_config::MT19937::init::seed]: expected an entry of type <bool>" << std::endl;
                        }
                    }
                    // --------------------------------------
                    // seedvalue
                    // --------------------------------------
                    if (init.contains("seedvalue")) {
                        if (py::isinstance<py::int_>(init["seedvalue"])) {
                            GlobalVar.MT19937_init_seedvalue = init["seedvalue"].cast<tuint>();
                        } else {
                            GlobalVar.slogcout(2) << "ERROR [process_config::MT19937::init::seedvalue]: expected an entry of type <unsigned int>" << std::endl;
                        }
                    }
                } else {
                    GlobalVar.slogcout(2) << "ERROR [process_config::MT19937::init]: expected an entry of type <dict>" << std::endl;
                }
            }
        } else {
            GlobalVar.slogcout(2) << "ERROR [process_config::MT19937]: expected an entry of type <dict>" << std::endl;
        }
    }

    // ======================================================
    // Legendre polynomials
    // ======================================================
    if (config.contains("Legendre")) {
        if (py::isinstance<py::dict>(config["Legendre"])) {
            py::dict legendre = config["Legendre"].cast<py::dict>();
            // ----------------------------------------------
            // numb
            // ----------------------------------------------
            if (legendre.contains("numb")) {
                if (py::isinstance<py::int_>(legendre["numb"])) {
                    GlobalVar.LegendrePolyH_init_numb = legendre["numb"].cast<tuint>();
                } else {
                    GlobalVar.slogcout(2) << "ERROR [process_config::Legendre::numb]: expected an entry of type <unsigned int>" << std::endl;
                }
            }
        } else {
            GlobalVar.slogcout(2) << "ERROR [process_config::Legendre]: expected an entry of type <dict>" << std::endl;
        }
    }
}


// #################################################################################
// logging
// #################################################################################
flxPyLogger pyLogger;

const char * PythonLoggerBuffer::logLevelToString(int logLevel)
{
    switch (logLevel) {
        //case 0: return "critical";
        case 1: return "critical";
        case 2: return "warning";
        case 3: return "info";
        case 4: return "info";
        case 5: return "debug";
        default: return "";
    }
}

int PythonLoggerBuffer::sync()
{
    std::string msg = str();
    if (!msg.empty()) {
        if (py_logger.is_none()) {
            throw FlxException_Crude("PythonLoggerBuffer::sync");
        }
        // remove newline at the end
            if (msg.back()=='\n') {
                msg.pop_back();
            }
        // Forward the message to the Python logger when buffer is flushed
            // Note: However, we need a special treatment of line breaks
            std::stringstream oss(msg);
            std::string line;
            while (std::getline(oss, line)) {
                py_logger.attr(logLevelToString(log_Level))(line);
            }

        // const int i = msg[msg.size()-2];
        // py_logger.attr(logLevelToString(log_Level))(msg + "[" + GlobalVar.Double2String(i) + "]");
    }
    str(""); // Clear the buffer
    return 0;
}

void PythonLoggerBuffer::set_logLevel(const int logLevel_)
{
    log_Level = logLevel_;
}

std::ostream& flxPyLogger::slog(const int logLevel_)
{
   streamPy.set_logLevel(logLevel_);
   return streamPy;
}

void flxPyLogger::set_logger(py::object logger_obj)
{
    streamPy.set_logger(logger_obj);
}

void set_logger(py::object logger_obj)
{
    pyLogger.set_logger(logger_obj);
    GlobalVar.set_logger(pyLogger);
}

void slog(const std::string& message, int logLevel) {
    check_engine_state();
    GlobalVar.slogcout(logLevel) << message << std::endl;
    GlobalVar.slogcout(logLevel).flush();
    finalize_call();
}


// #################################################################################
// Fesslix Engine
// #################################################################################

void error_out(std::string errMsg)
{
  // write to the actual error stream
    *(GlobalVar.get_true_cerr()) << std::endl << errMsg << std::endl;
  // write to the default error stream of Fesslix
    if (GlobalVar.get_true_cerr() != GlobalVar.get_cerr() ) {
      *(GlobalVar.get_cerr()) << std::endl << errMsg << std::endl;
    }
  // write to the log-file
    GlobalVar.logLevel_strong_reset();
    if (GlobalVar.check_logNOTcout()) {
      GlobalVar.slog(1) << std::endl << errMsg << std::endl;
      flx_print_RETURN_ERROR();
    }
  flxengine_unload();
  flx_delete_FLXcout();
}

void unload_engine()
{
    flxengine_unload();
    flx_delete_FLXcout();
}

int load_engine()
{
    // make sure it is not already running
        if (get_FlxEngine_ptr()) {
            unload_engine();
        }

    // start the Fesslix-engine
        const int itmp = flxengine_init();
        finalize_call();
        if (itmp>=0) return itmp;

    // return success
        return RETURN_SUCCESS;
}

class ModuleCleanup {
  public:
    ~ModuleCleanup() {
        unload_engine();
    }
};
static ModuleCleanup module_cleanup;

// #################################################################################
// Standard functions
// #################################################################################

void set_const(const std::string const_name, const tdouble value)
{
    FlxEngine().dataBox.ConstantBox.insert(const_name,value);
}

void set_var(const std::string var_name, py::object fun)
{
    FlxFunction* var_fun = parse_function(fun, "<flxPara> in 'set_var'");
    try {
        FlxEngine().dataBox.VarBox.insert(var_name,var_fun);
    } catch (FlxException& e) {
        delete var_fun;
        throw;
    }
}


// #################################################################################
// random variables
// #################################################################################




// #################################################################################
// sets of random variables
// #################################################################################

flxPyRVset::flxPyRVset(flxPyRVset& rhs)
: rvset_ptr(rhs.rvset_ptr), name_of_set(rhs.name_of_set)
{
}

flxPyRVset::flxPyRVset(flxPyRVset&& rhs)
: rvset_ptr(rhs.rvset_ptr), name_of_set(rhs.name_of_set)
{
    rhs.rvset_ptr = nullptr;
    rhs.name_of_set = "";
}

flxPyRVset& flxPyRVset::operator=(const flxPyRVset& rhs)
{
    rvset_ptr = rhs.rvset_ptr;
    name_of_set = rhs.name_of_set;
    return *this;
}

const std::string flxPyRVset::get_name() const
{
    return name_of_set;
}

const tuint flxPyRVset::get_NRV() const
{
    return rvset_ptr->get_NRV_only_this();
}

const tuint flxPyRVset::get_NOX() const
{
    return rvset_ptr->get_NOX_only_this();
}

py::array_t<tdouble> flxPyRVset::get_values(const std::string mode)
{
    // process mode
        enum RBRVvecGetType { x, u, mean, sd };
        RBRVvecGetType gType;
        if (mode=="x") {
            gType = RBRVvecGetType::x;
        } else if (mode=="u") {
            gType = RBRVvecGetType::u;
        } else if (mode=="mean") {
            gType = RBRVvecGetType::mean;
        } else if (mode=="sd") {
            gType = RBRVvecGetType::sd;
        } else {
            throw FlxException_NeglectInInteractive("flxPyRVset::get_values_01","Unkown mode '" + mode + "' for mode.");
        }
    // prepare result array
        const tuint NOX = rvset_ptr->get_NOX_only_this();
        const tuint NRV = rvset_ptr->get_NRV_only_this();
        if ( (gType==u&&NRV==0) || NOX==0 ) {
          std::ostringstream ssV;
          ssV << "The set '" << name_of_set << "' does not contain any random variables.";
          throw FlxException("FlxObjRBRV_vec_get::task_2", ssV.str() );
        }
        const tuint N = (gType==u)?NRV:NOX;
        // Allocate memory for the return array
        auto res_buf = py::array_t<tdouble>(N);
        // Get the buffer info to access the underlying return data
        py::buffer_info res_buf_info = res_buf.request();
        tdouble* res_ptr = static_cast<tdouble*>(res_buf_info.ptr);
    // assign to the array
      switch(gType) {
        case x:
          rvset_ptr->get_x_only_this(res_ptr);
          break;
        case u:
          rvset_ptr->get_y_only_this(res_ptr);
          break;
        case mean:
          rvset_ptr->get_mean_only_this(res_ptr);
          break;
        case sd:
          rvset_ptr->get_sd_only_this(res_ptr);
          break;
      }
    // Return the array
    return res_buf;
}

py::list flxPyRVset::get_rvs()
{
    check_engine_state();
    py::list res;
    const tuint N = rvset_ptr->get_NOX_only_this();
    for (tuint i=0;i<N;++i) {
        const std::string rv_name = rvset_ptr->get_rv_name(i);
        res.append( get_rv_from_set(rv_name) );
    }
    return res;
}

void flxPyRVset::set_y_vec(py::array_t<tdouble> arr) {
    const tuint NRV = rvset_ptr->get_NRV_only_this();
    // prepare input array
        // Access the input data as a raw pointer
        py::buffer_info buf_info = arr.request();
        tdouble* input_ptr = static_cast<tdouble*>(buf_info.ptr);
        // Get the size of the input array
        size_t size = buf_info.size;
        if (size!=NRV) {
            std::ostringstream ssV;
            ssV << "Input array has size " << size << ", whereas an input array of size " << NRV << " is expected.";
            throw FlxException_NeglectInInteractive("flxPyRVset::set_y_vec", ssV.str() );
        }
    rvset_ptr->set_is_valid(false);
    rvset_ptr->set_y_only_this(input_ptr);
    rvset_ptr->transform_y2x();
}

void flxPyRVset::set_x_vec(py::array_t<tdouble> arr) {
    const tuint NOX = rvset_ptr->get_NOX_only_this();
    // prepare input array
        // Access the input data as a raw pointer
        py::buffer_info buf_info = arr.request();
        tdouble* input_ptr = static_cast<tdouble*>(buf_info.ptr);
        // Get the size of the input array
        size_t size = buf_info.size;
        if (size!=NOX) {
            std::ostringstream ssV;
            ssV << "Input array has size " << size << ", whereas an input array of size " << NOX << " is expected.";
            throw FlxException_NeglectInInteractive("flxPyRVset::set_x_vec", ssV.str() );
        }
    rvset_ptr->set_is_valid(false);
    rvset_ptr->set_x_only_this(input_ptr);
    rvset_ptr->transform_x2y();
}

const tdouble flxPyRVset::pdf_log(py::array_t<tdouble> arr)
{
    this->set_x_vec(arr);
    return rvset_ptr->get_pdf_x_eval_log();
}

py::array_t<tdouble> flxPyRVset::eval_rp_psd(py::array_t<tdouble> arr)
{
    // ensure that process is based on "power spectral density"
        RBRV_set_psd* sb_psd = dynamic_cast<RBRV_set_psd*>(rvset_ptr);
        if (sb_psd==nullptr) {
            std::ostringstream ssV;
            ssV << "The rbrv-set '" << name_of_set << "' is not a random process (with specified power spectral density function).";
            throw FlxException_NeglectInInteractive("flxPyRVset::eval_rp_psd_01", ssV.str() );
        }
    // prepare input array
        // Access the input data as a raw pointer
        py::buffer_info buf_info = arr.request();
        tdouble* input_ptr = static_cast<tdouble*>(buf_info.ptr);
        // Get the size of the input array
        size_t size = buf_info.size;
    // prepare the output array
        // Allocate memory for the return array
        auto res_buf = py::array_t<tdouble>(size);
        // Get the buffer info to access the underlying return data
        py::buffer_info res_buf_info = res_buf.request();
        tdouble* res_ptr = static_cast<tdouble*>(res_buf_info.ptr);
    // Fill the array with values
        for (size_t i = 0; i < size; ++i) {
            res_ptr[i] = sb_psd->eval_realization(input_ptr[i]);    // TODO avoid re-evaluating the parameters of the random variable
        }
    // Return the array
        return res_buf;
}

flxPyRVset rbrv_set(py::dict config, py::list rv_list)
{
    check_engine_state();
    // eval and check name
        std::string set_name = parse_py_para_as_word("name",config,true,true);
        RBRV_entry_read_base::generate_set_base_check_name(FlxEngine().dataBox.rbrv_box, set_name);
    // prepare rv-set-creator
        const bool is_Nataf = parse_py_para_as_bool("is_Nataf", config, false, false);
        std::optional<FlxObjRBRV_set_creator> crtr;
        if (is_Nataf==false) {  // Rosenblatt transformation
            RBRV_set_baseDPtr parents = nullptr;
            try {
                // identify parents
                    std::vector<std::string> set_parents;
                    parse_py_para_as_word_lst(set_parents,"parents",config,false,true);
                    const tuint Nparents = set_parents.size();
                    RBRV_entry_read_base::generate_set_base(FlxEngine().dataBox.rbrv_box, set_parents,parents);
                const bool allow_x2y = parse_py_para_as_bool("allow_x2y", config, false, false);
                crtr.emplace(set_name,parents,Nparents,allow_x2y);
            } catch (FlxException& e) {
                if (parents) delete [] parents;
                throw;
            }
        } else {    // Nataf transformation
            if (config.contains("parents")) {
                throw FlxException_NeglectInInteractive("rbrv_set_01", "Python <dict> of 'rbrv_set' for Nataf transformation MUST NOT contain 'parents'.");
            }
            const bool is_Nataf_only_once = parse_py_para_as_bool("is_Nataf_only_once", config, false, true);
            crtr.emplace(set_name,is_Nataf_only_once);
        }
    // create the random variables
        const tuint Nrv = rv_list.size();
        RBRV_entry* entry = nullptr;
        FlxFunction* csVal = nullptr;
        std::string csNam;
        bool csFix = false;
        const std::string family = set_name + "::";
        try {
            tuint iID = 0;
            for (tuint i=0;i<Nrv;++i) {
                // cast config of RV into a dict
                    std::ostringstream ssV;
                    ssV << "configuration for random variable " << (i+1);
                    py::dict rv_config = parse_py_obj_as_dict(rv_list[i], ssV.str());
                // generate RV from config
                    if (crtr->get_Nataf_evalOnce()) {
                        rv_config["eval_once"] = true;
                    }
                    entry = parse_py_obj_as_grv(rv_config, iID, family, ssV.str());
                // correlation (for Rosenblatt)
                    if (rv_config.contains("corr")) {
                        py::dict corr_config = parse_py_obj_as_dict(rv_config["corr"], "'corr' in " + ssV.str());
                        csNam = family + parse_py_para_as_string("rv_name", corr_config, true);
                        csVal = parse_py_para("value", corr_config, true);
                        csFix = parse_py_para_as_bool("fix", corr_config, false, false);
                    }
                // register the random variable in the set
                    try {
                        crtr->add_entry(FlxEngine().dataBox.rbrv_box,entry, csVal, csNam, csFix);
                        csVal = nullptr;
                        entry = nullptr;
                    } catch (FlxException& e) {
                        // memory management is taken care of by crtr->add_entry(...)
                            entry = nullptr;
                            csVal = nullptr;
                        throw;
                    }
            }
        } catch (FlxException& e) {
            if (entry) delete entry;
            if (csVal) delete csVal;
            throw;
        }
    // set up the correlation matrix (in case of Nataf transformation)
        if (is_Nataf) {
            if (config.contains("corr")) {
                py::list corr_lst = parse_py_obj_as_list(config["corr"],"key 'corr'");
                const size_t Nc = corr_lst.size();
                for (size_t i=0; i<Nc;++i) {
                    std::ostringstream ssV;
                    ssV << "entry '" << (i+1) << "'";
                    py::dict corr_entry = parse_py_obj_as_dict( corr_lst[i],ssV.str() );
                    const std::string rv1 = family + parse_py_para_as_word("rv_1",corr_entry, true, true);
                    const std::string rv2 = family + parse_py_para_as_word("rv_2",corr_entry, true, true);
                    const tdouble rho = parse_py_para_as_float("value",corr_entry,true);
                    const bool corr_approx = parse_py_para_as_bool("corr_approx",corr_entry,false,true);
                    const bool rhogauss = parse_py_para_as_bool("rhogauss",corr_entry,false,false);
                    crtr->add_corr(rv1, rv2, rho, corr_approx, rhogauss, false);
                }
            }
        }
    // register the set in the engine
        RBRV_set_base* set_ptr = crtr->register_set(FlxEngine().dataBox.rbrv_box,true);
        flxPyRVset res(set_ptr, set_name);
    finalize_call();
    return res;
}

flxPyRVset rbrv_set_noise(py::dict config, py::dict rv_config)
{
    check_engine_state();
    // eval and check name
        std::string set_name = parse_py_para_as_word("name",config,true,true);
        RBRV_entry_read_base::generate_set_base_check_name(FlxEngine().dataBox.rbrv_box, set_name);
        const std::string family = set_name + "::";
    // number of random variables in the set
        const tuint Ndim = parse_py_para_as_tuintNo0("N", config, true);
    RBRV_set_baseDPtr parents = nullptr;
    RBRV_entry_RV_base* entry = nullptr;
    RBRV_set_noise* ts = NULL;
    try {
        // identify parents
            std::vector<std::string> set_parents;
            parse_py_para_as_word_lst(set_parents,"parents",config,false,true);
            const tuint Nparents = set_parents.size();
            RBRV_entry_read_base::generate_set_base(FlxEngine().dataBox.rbrv_box, set_parents,parents);
        // retrieve base type / distribution of random variable
            entry = parse_py_obj_as_rv(rv_config, false, 0, family, "rv_config");
        // generate set
            ts = new RBRV_set_noise(false,Ndim,set_name,false,entry,Nparents,parents);
            parents = nullptr;
            entry = nullptr;
            FlxEngine().dataBox.rbrv_box.register_set(ts);
    } catch (FlxException& e) {
        FLXMSG("rbrv_set_noise_01",1);
        if (parents) delete [] parents;
        if (entry) delete entry;
        if (ts) delete ts;
        throw;
    }
    // create and return Python-object of generated set
        flxPyRVset res(ts, set_name);
        finalize_call();
        return res;
}

flxPyRVset rbrv_set_proc(py::dict config, py::dict rv_config)
{
    check_engine_state();
    // eval and check name
        std::string set_name = parse_py_para_as_word("name",config,true,true);
        RBRV_entry_read_base::generate_set_base_check_name(FlxEngine().dataBox.rbrv_box, set_name);
        const std::string family = set_name + "::";
    // number of random variables in the set
        const tuint Ndim = parse_py_para_as_tuintNo0("N", config, true);
    // "distance" between two points
        const tdouble dx = parse_py_para_as_floatPosNo0("dx",config,false,ONE);
    const tuint M = parse_py_para_as_tuint("M", config, false, 0);
    const tuint evtype = parse_py_para_as_tuint("evtype", config, false, 2);
    const bool only_once = parse_py_para_as_bool("only_once", config, false, true );
    const bool rhoGauss = parse_py_para_as_bool("rhogauss", config, false, false );
    RBRV_set_baseDPtr parents = nullptr;
    RBRV_entry_RV_base* entry = nullptr;
    RBRV_set_proc* ts = nullptr;
    FlxFunction* rho = nullptr;
    try {
        // auto-correlation coefficient function
            rho = parse_py_para("rho",config,true);
        // identify parents
            std::vector<std::string> set_parents;
            parse_py_para_as_word_lst(set_parents,"parents",config,false,true);
            const tuint Nparents = set_parents.size();
            RBRV_entry_read_base::generate_set_base(FlxEngine().dataBox.rbrv_box, set_parents,parents);
        // retrieve base type / distribution of random variable
            entry = parse_py_obj_as_rv(rv_config, false, 0, family, "rv_config");
        // generate set
            ts = new RBRV_set_proc(false,Ndim,M,set_name,false,entry,rho,dx,Nparents,parents,evtype,only_once,rhoGauss);
            parents = nullptr;
            entry = nullptr;
            rho = nullptr;
            FlxEngine().dataBox.rbrv_box.register_set(ts);
    } catch (FlxException& e) {
        FLXMSG("rbrv_set_proc",1);
        if (rho) delete rho;
        if (parents) delete [] parents;
        if (entry) delete entry;
        if (ts) delete ts;
        throw;
    }
    // create and return Python-object of generated set
        flxPyRVset res(ts, set_name);
        finalize_call();
        return res;
}

flxPyRVset rbrv_set_psd(py::dict config)
{
    check_engine_state();
    // eval and check name
        std::string set_name = parse_py_para_as_word("name",config,true,true);
        RBRV_entry_read_base::generate_set_base_check_name(FlxEngine().dataBox.rbrv_box, set_name);
        const std::string family = set_name + "::";
    // number of random variables in the set
        const tuint N = parse_py_para_as_tuintNo0("N", config, true);
    // bounds
        const tdouble lb = parse_py_para_as_float("lb",config,true);
        const tdouble ub = parse_py_para_as_float("ub",config,true);
        if (ub<=lb) {
            throw FlxException("rbrv_set_psd_01", "Lower bound must be smaller than upper bound.");
        }
    RBRV_set_baseDPtr parents = nullptr;
    RBRV_set_psd* ts = NULL;
    FlxFunction* psd_fun = nullptr;
    try {
        // power spectral density function
            psd_fun = parse_py_para("psd",config,true);
        // identify parents
            std::vector<std::string> set_parents;
            parse_py_para_as_word_lst(set_parents,"parents",config,false,true);
            const tuint Nparents = set_parents.size();
            RBRV_entry_read_base::generate_set_base(FlxEngine().dataBox.rbrv_box, set_parents,parents);
        // generate set
            ts = new RBRV_set_psd(false,set_name,N,psd_fun,lb,ub,Nparents,parents,FlxEngine().dataBox.ConstantBox.getRef("gx"));
            parents = nullptr;
            psd_fun = nullptr;
            FlxEngine().dataBox.rbrv_box.register_set(ts);
    } catch (FlxException& e) {
        FLXMSG("rbrv_set_psd_55",1);
        if (psd_fun) delete psd_fun;
        if (parents) delete [] parents;
        if (ts) delete ts;
        throw;
    }
    // create and return Python-object of generated set
        flxPyRVset res(ts, set_name);
        finalize_call();
        return res;
}

flxPyRVset rbrv_set_sphere(py::dict config)
{
    check_engine_state();
    // eval and check name
        std::string set_name = parse_py_para_as_word("name",config,true,true);
        RBRV_entry_read_base::generate_set_base_check_name(FlxEngine().dataBox.rbrv_box, set_name);
        const std::string family = set_name + "::";
    // number of random variables in the set
        const tuint N = parse_py_para_as_tuintNo0("N", config, true);
    // radius of sphere
        FlxFunction* radius = parse_py_para("radius",config,true);
    RBRV_set_baseDPtr parents = nullptr;
    RBRV_set_sphere* ts = NULL;
    try {
        // identify parents
            std::vector<std::string> set_parents;
            parse_py_para_as_word_lst(set_parents,"parents",config,false,true);
            const tuint Nparents = set_parents.size();
            RBRV_entry_read_base::generate_set_base(FlxEngine().dataBox.rbrv_box, set_parents,parents);
        // generate set
            ts = new RBRV_set_sphere(false,N,set_name,false,Nparents,parents,radius);
            parents = nullptr;
            radius = nullptr;
            FlxEngine().dataBox.rbrv_box.register_set(ts);
    } catch (FlxException& e) {
        FLXMSG("rbrv_set_sphere_55",1);
        if (radius) delete radius;
        if (parents) delete [] parents;
        if (ts) delete ts;
        throw;
    }
    // create and return Python-object of generated set
        flxPyRVset res(ts, set_name);
        finalize_call();
        return res;
}

flxPyRVset rbrv_set_vfun(py::dict config)
{
    check_engine_state();
    // eval and check name
        std::string set_name = parse_py_para_as_word("name",config,true,true);
        RBRV_entry_read_base::generate_set_base_check_name(FlxEngine().dataBox.rbrv_box, set_name);
        const std::string family = set_name + "::";
    // number of random variables in the set
        const tuint N = parse_py_para_as_tuintNo0("N", config, true);
    // flxMtxFun
        FlxMtxFun_base* vecfun = parse_py_para_as_flxMtxFun(N,"vecfun",config);
    RBRV_set_baseDPtr parents = nullptr;
    RBRV_vfset* ts = NULL;
    try {
        // identify parents
            std::vector<std::string> set_parents;
            parse_py_para_as_word_lst(set_parents,"parents",config,false,true);
            const tuint Nparents = set_parents.size();
            RBRV_entry_read_base::generate_set_base(FlxEngine().dataBox.rbrv_box, set_parents,parents);
        // generate set
            ts = new RBRV_vfset(false,set_name,false,N,vecfun,Nparents,parents);
            vecfun = nullptr;
            parents = nullptr;
            FlxEngine().dataBox.rbrv_box.register_set(ts);
    } catch (FlxException& e) {
        FLXMSG("rbrv_set_vfun_55",1);
        if (vecfun) delete vecfun;
        if (parents) delete [] parents;
        if (ts) delete ts;
        throw;
    }
    // create and return Python-object of generated set
        flxPyRVset res(ts, set_name);
        finalize_call();
        return res;
}

flxPyRVset rbrv_set_dirichlet(py::dict config)
{
    check_engine_state();
    // eval and check name
        std::string set_name = parse_py_para_as_word("name",config,true,true);
        RBRV_entry_read_base::generate_set_base_check_name(FlxEngine().dataBox.rbrv_box, set_name);
        const std::string family = set_name + "::";
    // number of random variables in the set
        flxVec alpha_vec = parse_py_para_as_flxVec("alpha",config,true);
    RBRV_set_baseDPtr parents = nullptr;
    RBRV_dirichlet* ts = NULL;
    try {
        // identify parents
            std::vector<std::string> set_parents;
            parse_py_para_as_word_lst(set_parents,"parents",config,false,true);
            const tuint Nparents = set_parents.size();
            RBRV_entry_read_base::generate_set_base(FlxEngine().dataBox.rbrv_box, set_parents,parents);
        // generate set
            ts = new RBRV_dirichlet(false,set_name,false,alpha_vec.get_N(),nullptr,Nparents,parents,&alpha_vec);
            parents = nullptr;
            FlxEngine().dataBox.rbrv_box.register_set(ts);
    } catch (FlxException& e) {
        FLXMSG("rbrv_set_vfun_55",1);
        if (parents) delete [] parents;
        if (ts) delete ts;
        throw;
    }
    // create and return Python-object of generated set
        flxPyRVset res(ts, set_name);
        finalize_call();
        return res;
}

flxPyRVset rbrv_set_multinomial(py::dict config)
{
    check_engine_state();
    // eval and check name
        std::string set_name = parse_py_para_as_word("name",config,true,true);
        RBRV_entry_read_base::generate_set_base_check_name(FlxEngine().dataBox.rbrv_box, set_name);
        const std::string family = set_name + "::";
    // number of random variables in the set
        const tuint N = parse_py_para_as_tuintNo0("N", config, true);
    // number of random variables in the set
        const tuint Ntrial = parse_py_para_as_tuintNo0("Ntrial", config, true);
    // flxMtxFun
        FlxMtxFun_base* vecfun = parse_py_para_as_flxMtxFun(N,"pvec",config);
    RBRV_set_baseDPtr parents = nullptr;
    RBRV_multinomial* ts = NULL;
    try {
        // identify parents
            std::vector<std::string> set_parents;
            parse_py_para_as_word_lst(set_parents,"parents",config,false,true);
            const tuint Nparents = set_parents.size();
            RBRV_entry_read_base::generate_set_base(FlxEngine().dataBox.rbrv_box, set_parents,parents);
        // generate set
            ts = new RBRV_multinomial(false,set_name,false,N,vecfun,Nparents,parents,Ntrial);
            vecfun = nullptr;
            parents = nullptr;
            FlxEngine().dataBox.rbrv_box.register_set(ts);
    } catch (FlxException& e) {
        FLXMSG("rbrv_set_multinomial_55",1);
        if (vecfun) delete vecfun;
        if (parents) delete [] parents;
        if (ts) delete ts;
        throw;
    }
    // create and return Python-object of generated set
        flxPyRVset res(ts, set_name);
        finalize_call();
        return res;
}

flxPyRV get_rv_from_set(const std::string& rv_name)
{
    check_engine_state();
    RBRV_entry* rv_ptr = FlxEngine().dataBox.rbrv_box.get_entry(rv_name,true);
    flxPyRV res(rv_ptr);
    finalize_call();
    return res;
}


// #################################################################################
// advanced features
// #################################################################################

tdouble eval_fun(py::object expr)
{
    check_engine_state();
    FlxFunction* fun = parse_function(expr,"parameter 'expr' of 'flx.eval_fun'");
    tdouble res;
    try {
        res = fun->calc();
        delete fun;
    } catch (FlxException& e) {
        delete fun;
        GlobalVar.slogcout(1) << std::endl << e.what() << std::endl;
        res = std::numeric_limits<tdouble>::quiet_NaN();
    }
    return res;
}

void eval_code(py::object expr)
{
    check_engine_state();
    FlxObjBase* code = parse_code(expr,"parameter 'expr' of 'flx.eval_code'");
    try {
        code->exec();
        delete code;
    } catch (FlxException& e) {
        delete code;
        throw;
    }
}

py::array_t<tdouble> eval_vfun(const tuint N, py::object expr)
{
    check_engine_state();
    FlxMtxFun_base* vfun = parse_FlxMtxFun(N,expr,"parameter 'expr' of 'flx.eval_vfun'");
    try {
        vfun->eval();
        const flxVec& res_ = vfun->get_res_vec();
        // Allocate memory for the return array
        auto res_buf = py::array_t<tdouble>(N);

        // Get the buffer info to access the underlying return data
        py::buffer_info res_buf_info = res_buf.request();
        tdouble* res_ptr = static_cast<tdouble*>(res_buf_info.ptr);
        flxVec res(res_ptr,N);

        // copy values into result vector
        res = res_;
        delete vfun;
        return res_buf;
    } catch (FlxException& e) {
        delete vfun;
        throw;
    }
}

// #################################################################################
// only for debugging purposes
// #################################################################################

double add(double a, double b) {
    check_engine_state();
    return a + b;
}


// #################################################################################
// Expose interface to Python
// #################################################################################

PYBIND11_MODULE(core, m) {
    register_dataBox_post_processors();
    // ====================================================
    // 'global' functions
    // ====================================================
        m.def("Double2String", &Double2String, "Convert a double into a string");   // TODO docu
        m.def("print_info", &print_info, "output configuration options of Fesslix");

    // ====================================================
    // load configuration
    // ====================================================
        m.def("set_exe_dir", &set_exe_dir, "Set the directory which contains the Fesslix module");
        m.def("process_config", &process_config, "Process a dictionary containing configuration options");

    // ====================================================
    // logging
    // ====================================================
        m.def("set_logger", &set_logger, "Set the logger object");
        m.def("slog", &slog, pybind11::arg("message"), pybind11::arg("logLevel") = 4, "Log a message at a specified level");

    // ====================================================
    // Fesslix Engine
    // ====================================================
        m.def("load_engine", &load_engine, "Load the Fesslix Engine");

    // ====================================================
    // Standard functions
    // ====================================================
        m.def("set_const", &set_const, "assigns a value to a const-variable");
        m.def("set_var", &set_var, "assigns a value to a var-variable");

        m.def("cdfn", &rv_Phi, "the CDF of the standard Normal distribution");
        m.def("cdfn_inv", &rv_InvPhi, "inverse of the CDF of the standard Normal distribution");

        m.def("randStr", &generate_randStr, "returns a random string of length N");
        m.def("get_time_since_start", []() {return GlobalVar.get_time_since_start(); }, "returns the time since start-up in seconds");

    // ====================================================
    // random variables
    // ====================================================
        py::class_<flxPyRV>(m, "rv")
            .def(py::init<py::dict>())
            .def("get_name", &flxPyRV::get_name, "get name of random variable")
            .def("get_descr", &flxPyRV::get_descr, "get description of random variable")
            .def("get_type", &flxPyRV::get_type, "get type of random variable")
            .def("get_value", &flxPyRV::get_value, "get the value currently assigned to the random variable")
            .def("x2y", &flxPyRV::x2y, "transformation from 'original space' to standard normal space")
            .def("y2x", &flxPyRV::y2x, "transformation from standard normal space into 'original space'")
            .def("sample", &flxPyRV::sample, "generate a random realization of the random variable")
            .def("sample_array", &flxPyRV::sample_array, "generate a vector of random realizations of the random variable")
            .def("pdf", &flxPyRV::pdf, pybind11::arg("x_val"), pybind11::arg("safeCalc") = true, "evaluates the pdf of the random variable at x_val")
            .def("pdf_array", &flxPyRV::pdf_array, pybind11::arg("x_vec"), pybind11::arg("safeCalc") = true, "evaluates the pdf of the random variable for array x_vec")
            .def("pdf_log", &flxPyRV::pdf_log, pybind11::arg("x_val"), pybind11::arg("safeCalc") = true, "evaluates the log-pdf of the random variable at x_val")
            .def("cdf", &flxPyRV::cdf, pybind11::arg("x_val"), pybind11::arg("safeCalc") = true, "evaluates the cdf of the random variable at x_val")
            .def("cdf_array", &flxPyRV::cdf_array, pybind11::arg("x_vec"), pybind11::arg("safeCalc") = true, "evaluates the cdf of the random variable for array x_vec")
            .def("icdf", &flxPyRV::icdf, "evaluates the inverse of the cdf of the random variable for probability p")
            .def("icdf_array", &flxPyRV::icdf_array, "evaluates the inverse of the cdf of the random variable for array p_vec")
            .def("sf", &flxPyRV::sf, pybind11::arg("x_val"), pybind11::arg("safeCalc") = true, "returns the survival function of the random variable; i.e., 1-cdf(x_val)")
            .def("sf_array", &flxPyRV::sf_array, pybind11::arg("x_vec"), pybind11::arg("safeCalc") = true, "evaluates the survival function of the random variable for array x_vec")
            .def("entropy", &flxPyRV::entropy, "returns the entropy of the random variable")
            .def("mean", &flxPyRV::mean, "returns the mean of the random variable")
            .def("sd", &flxPyRV::sd, "returns the standard deviation of the random variable")
            .def("median", &flxPyRV::median, "returns the median of the random variable")
            .def("mode", &flxPyRV::mode, "returns the mode of the random variable")
            .def("check_x", &flxPyRV::check_x, "check if x_val is inside of the valid domain of the random variable")
            .def("get_HPD", &flxPyRV::get_HPD, "returns the lower quantile value of the HPD (highest probability density) interval of the distribution")
            .def("info", &flxPyRV::info, "return information about random variable");

    // ====================================================
    // sets of random variables
    // ====================================================
        py::class_<flxPyRVset>(m, "rvset")
            .def("get_name", &flxPyRVset::get_name, "get name of set of random variables")
            .def("get_NRV", &flxPyRVset::get_NRV, "return number of random variables (in standard Normal space) in the set of random variables")
            .def("get_NOX", &flxPyRVset::get_NOX, "return number of random variables (in original space) in the set of random variables")
            .def("get_values", &flxPyRVset::get_values, pybind11::arg("mode")="x", "returns a vector of quantities of all entries contained in the set of random variables")
            .def("get_rvs", &flxPyRVset::get_rvs, "retrieve a list of all random variables in the set")
            .def("set_y_vec", &flxPyRVset::set_y_vec, "assign y_vec as standard Normal values of the random variables in the set")
            .def("set_x_vec", &flxPyRVset::set_x_vec, "assign x_vec as values of the random variables in the set")
            .def("pdf_log", &flxPyRVset::pdf_log, "evaluates the log-pdf of the set of random variable at x_vec")
            .def("eval_rp_psd", &flxPyRVset::eval_rp_psd, "evaluates the realization of a random process at time t with given power spectral density function");

        m.def("rv_set", &rbrv_set, "creates a set of general random variables");
        m.def("rv_set_noise", &rbrv_set_noise, "creates a set of independent random variables that have all the same distribution");
        m.def("rv_set_proc", &rbrv_set_proc, "creates a discrete random process with given correlation structure");
        m.def("rv_set_psd", &rbrv_set_psd, "creates a Gaussian process with given power spectral density function");
        m.def("rv_set_sphere", &rbrv_set_sphere, "creates a random point that is uniformly distribted in a hyper-sphere");
        m.def("rv_set_vfun", &rbrv_set_vfun, "creates a vector that is an explicit function of random variables");
        m.def("rv_set_dirichlet", &rbrv_set_dirichlet, "creates a vector that is the outcome of a Dirichlet distributed random variable");
        m.def("rv_set_multinomial", &rbrv_set_multinomial, "creates a vector that is the outcome of a multinomial distributed random variable");

        m.def("get_rv_from_set", &get_rv_from_set, "retrieve a random variable from a set of random variables");

        py::class_<flxPySampler>(m, "sampler")
            .def(py::init<py::list>())
            .def("sample", &flxPySampler::sample, "generate a random sample for a collection of sets of random variables")
            .def("get_NRV", &flxPySampler::get_NRV, "return number of random variables (in standard Normal space) in the collection of random variables")
            .def("get_NOX", &flxPySampler::get_NOX, "return number of random variables (in original space) in the collection of random variables")
            .def("get_values", &flxPySampler::get_values, pybind11::arg("mode")="x", "returns a vector of quantities of all entries contained in the sampler")
            .def("assign_u", &flxPySampler::assign_u, "assign y_vec as standard Normal values of the random variables in the sampler and perform the transformation to original space")
            .def("assign_x", &flxPySampler::assign_x, pybind11::arg("x_vec"), pybind11::arg("transform") = true, "assign x_vec as values of the random variables in the sampler")
            .def("perform_MCS", &flxPySampler::perform_MCS, "perform a Monte Carlo simulation")
            ;

    // ====================================================
    // dataBox
    // ====================================================
        py::class_<post_proc_base>(m, "postProc")
            .def("eval", &post_proc_base::eval,"retrieve the state of the post-processor")
            ;
        py::class_<flxDataBox>(m, "dataBox")
            .def(py::init<tuint, tuint>())
            .def("write2mem", &flxDataBox::write2mem, "allocate memory for storing data")
            .def("extract_col_from_mem", &flxDataBox::extract_col_from_mem, py::return_value_policy::reference_internal, "return a numpy-array that points to the memory of 'col'")
            .def("free_mem", &flxDataBox::free_mem, "free the allocated memory for storing data")
            .def("write2file", &flxDataBox::write2file, "send samples to a file")
            .def("read_from_file", &flxDataBox::read_from_file, "import data from a file")
            .def("close_file", &flxDataBox::close_file, "close the file stream")
            .def("register_post_processor", &flxDataBox::register_post_processor, py::return_value_policy::reference_internal, "register a new post-processor")
            ;

    // ====================================================
    // Structural reliability analysis
    // ====================================================
            m.def("perform_Line_Sampling", &perform_Line_Sampling, "Perform Line Sampling (structural reliability analysis)");
            m.def("perform_FORM", &perform_FORM, "Perform FORM (First Order Reliability Method)");

    // ====================================================
    // physical analysis
    // ====================================================
        m.def("phys_temp2svp", &phys_temp2svp, "Evaluates the saturation vapor pressure as a function of temperature");
        m.def("phys_dewpoint", &flxPhys_dewpoint, "evaluates the dew point as a function of temperature and humidity");
        m.def("phys_tauphi2temp", &flxPhys_tauphi2temp, "evaluates the temperature associated with a given dew point and humidity");
        m.def("phys_tauptemp2phi", &flxPhys_tautemp2phi, "evaluates the humidity associated with a given dew point and temperature");
        m.def("phys_abs_humidity", &flxPhys_abs_humidity, "evaluates the absolute humidity as a function of temperature and humidity");

    // ====================================================
    // Advanced features
    // ====================================================
        m.def("eval_fun", &eval_fun, "Evaluates an expression and returns the result.");
        m.def("eval_code", &eval_code, "Evaluates code in the Fesslix engine.");
        m.def("eval_vfun", &eval_vfun, "Evaluates an expression and returns a vector.");

    // ====================================================
    // only for debugging purposes (TODO remove at some point)
    // ====================================================
        m.def("add", &add, "A function that adds two numbers");
        m.attr("the_answer") = 42;
}


