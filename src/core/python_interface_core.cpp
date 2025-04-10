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

#include "python_interface_core.h"
#include "fesslix.h"
#include "flxfunction_data.h"
#include "flxrbrv_rvs_read.h"
#include <iostream>

void check_engine_state()
{
    if (FlxEngine==nullptr) {
        throw FlxException_NeglectInInteractive("check_engine_state","Fesslix Engine is not running", "Please start the engine using load_engine()");
    }
}

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
        if (FlxEngine) {
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
            GlobalVar.slogcout(2) << "ERROR [process_config::TOL]: expected an entry of type <float>, got <" << py::str(config["TOL"].get_type()) << ">" << std::endl;
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
        if (FlxEngine) {
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
// random variables
// #################################################################################

RBRV_entry_RV_base * parse_py_obj_as_rv(py::dict config, const bool name_required, const tuint iID, const std::string family, std::string descr)
{
    RBRV_entry_RV_base* rv_ptr = nullptr;
    // retrieve type
        const std::string rv_type = parse_str_as_word(parse_py_para_as_string("type",config,true),true);
    // retrieve name
        std::string rv_name = family + parse_str_as_word(parse_py_para_as_string("name",config,name_required,"name_unspecified"),true);
    // select rbrv-class based on type
        if (rv_type=="stdn") {
            rv_ptr = new RBRV_entry_RV_stdN(rv_name,iID,config);
        } else if (rv_type=="normal") {
            rv_ptr = new RBRV_entry_RV_normal(rv_name,iID,config);
        } else if (rv_type=="logn") {
            rv_ptr = new RBRV_entry_RV_lognormal(rv_name,iID,config);
        } else if (rv_type=="uniform") {
            rv_ptr = new RBRV_entry_RV_uniform(rv_name,iID,config);
        // } else if (rv_type=="gumbel") {
        //     rv_ptr = new RBRV_entry_read_Gumbel(rv_name,iID,config);
        // } else if (rv_type=="normal_trunc") {
        //     rv_ptr = new RBRV_entry_read_normal_trunc(rv_name,iID,config);
        } else if (rv_type=="beta") {
            rv_ptr = new RBRV_entry_RV_beta(rv_name,iID,config);
        // } else if (rv_type=="exponential") {
        //     rv_ptr = new RBRV_entry_read_exponential(rv_name,iID,config);
        // } else if (rv_type=="gamma") {
        //     rv_ptr = new RBRV_entry_read_gamma(rv_name,iID,config);
        // } else if (rv_type=="poisson") {
        //     rv_ptr = new RBRV_entry_read_Poisson(rv_name,iID,config);
        // } else if (rv_type=="binomial") {
        //     rv_ptr = new RBRV_entry_read_Binomial(rv_name,iID,config);
        // } else if (rv_type=="cauchy") {
        //     rv_ptr = new RBRV_entry_read_Cauchy(rv_name,iID,config);
        // } else if (rv_type=="weibull") {
        //     rv_ptr = new RBRV_entry_read_Weibull(rv_name,iID,config);
        // } else if (rv_type=="chisquared") {
        //     rv_ptr = new RBRV_entry_read_ChiSquared(true,rv_name,iID,config);
        // } else if (rv_type=="chi") {
        //     rv_ptr = new RBRV_entry_read_ChiSquared(false,rv_name,iID,config);
        } else if (rv_type=="studentst") {
            rv_ptr = new RBRV_entry_RV_StudentsT(rv_name,iID,config);
        } else if (rv_type=="studentstgen") {
            rv_ptr = new RBRV_entry_RV_StudentsT_generalized(rv_name,iID,config);
        } else if (rv_type=="logt") {
            rv_ptr = new RBRV_entry_RV_logt(rv_name,iID,config);
        // } else if (rv_type=="laplace") {
        //     rv_ptr = new RBRV_entry_read_Laplace(rv_name,iID,config);
        // } else if (rv_type=="usertransform") {
        //     rv_ptr = new RBRV_entry_read_UserTransform(rv_name,iID,config);
        // } else if (rv_type=="truncated") {
        //     rv_ptr = new RBRV_entry_read_Truncated(rv_name,iID,config);
        // } else if (rv_type=="maxmintransform") {
        //     rv_ptr = new RBRV_entry_read_maxminTransform(rv_name,iID,config);
        // } else if (rv_type=="bayda") {
        //     rv_ptr = new RBRV_entry_read_bayDA(rv_name,iID,config);
        } else {
            std::ostringstream ssV;
            ssV << "Unknown random variable type '" << rv_type << "'.";
            throw FlxException("flxPyRV::flxPyRV_50", ssV.str() );
        }
    return rv_ptr;
}

flxPyRV::flxPyRV(py::dict config)
: rv_ptr(nullptr), mem_managed(true)
{
    rv_ptr_ = parse_py_obj_as_rv(config, false, 0, "", "flx.rv");
    rv_ptr = rv_ptr_;
    finalize_call();
}

flxPyRV::flxPyRV(RBRV_entry* rv_ptr)
: rv_ptr(rv_ptr), rv_ptr_(dynamic_cast<RBRV_entry_RV_base*>(rv_ptr)), mem_managed(false)
{

}

flxPyRV::flxPyRV(flxPyRV& rhs)
: rv_ptr(rhs.rv_ptr), rv_ptr_(rhs.rv_ptr_), mem_managed(rhs.mem_managed)
{
    if (mem_managed) {
        throw FlxException_NotImplemented("flxPyRV::flxPyRV(flxPyRV& rhs)");
    }
}

flxPyRV::flxPyRV(flxPyRV && rhs)
: rv_ptr(rhs.rv_ptr), rv_ptr_(rhs.rv_ptr_), mem_managed(rhs.mem_managed)
{
    rhs.rv_ptr = nullptr;
    rhs.rv_ptr_ = nullptr;
    rhs.mem_managed = false;
}

flxPyRV::~flxPyRV()
{
    if (rv_ptr) {
        if (mem_managed) {
            delete rv_ptr;
        }
    }
}

void flxPyRV::ensure_is_a_basic_rv()
{
    if (rv_ptr_==nullptr) {
        throw FlxException_NeglectInInteractive("flxPyRV::ensure_is_a_basic_rv", "'" + rv_ptr->name + "' is not a basic random variable.");
    }
}

const std::string flxPyRV::get_name() const
{
    return rv_ptr->name;
}

const std::string flxPyRV::get_type() const
{
    return rv_ptr->get_type();
}

const tdouble flxPyRV::get_value() const
{
    return rv_ptr->get_value();
}

const tdouble flxPyRV::x2y(const tdouble x_val)
{
    rv_ptr->eval_para();
    return rv_ptr->transform_x2y(x_val);
}

const tdouble flxPyRV::y2x(const tdouble y_val)
{
    ensure_is_a_basic_rv();
    rv_ptr_->eval_para();
    return rv_ptr_->transform_y2x(y_val);
}

const tdouble flxPyRV::sample()
{
    ensure_is_a_basic_rv();
    rv_ptr_->eval_para();
    const tdouble y = FlxEngine->dataBox.RndCreator.gen_smp();
    return rv_ptr_->transform_y2x(y);
}

void flxPyRV::sample_array(py::array_t<tdouble> arr)
{
    ensure_is_a_basic_rv();
    rv_ptr_->eval_para();

    // Access the data as a raw pointer
    py::buffer_info buf_info = arr.request();
    tdouble* res_ptr = static_cast<tdouble*>(buf_info.ptr);

    // Get the size of the input array
    size_t size = buf_info.size;

    // Generate the samples (in standard Normal space)
    flxVec res_vec(res_ptr,size,false,false);
    FlxEngine->dataBox.RndCreator.gen_smp(res_vec);

    // transform the samples to original space
    for (size_t i = 0; i < size; ++i) {
        res_ptr[i] = rv_ptr_->transform_y2x(res_ptr[i]);    // TODO avoid re-evaluating the parameters of the random variable
    }
}

const tdouble flxPyRV::pdf(const tdouble x_val, const bool safeCalc)
{
    rv_ptr->eval_para();
    return rv_ptr->calc_pdf_x(x_val,safeCalc);
}

py::array_t<tdouble> flxPyRV::pdf_array(py::array_t<tdouble> arr, const bool safeCalc)
{
    rv_ptr->eval_para();
    // Access the input data as a raw pointer
    py::buffer_info buf_info = arr.request();
    tdouble* input_ptr = static_cast<tdouble*>(buf_info.ptr);

    // Get the size of the input array
    size_t size = buf_info.size;

    // Allocate memory for the return array
    auto res_buf = py::array_t<tdouble>(size);

    // Get the buffer info to access the underlying return data
    py::buffer_info res_buf_info = res_buf.request();
    tdouble* res_ptr = static_cast<tdouble*>(res_buf_info.ptr);

    // Fill the array with values
    for (size_t i = 0; i < size; ++i) {
        res_ptr[i] = rv_ptr->calc_pdf_x(input_ptr[i],safeCalc);    // TODO avoid re-evaluating the parameters of the random variable
    }

    // Return the array
    return res_buf;


}

const tdouble flxPyRV::pdf_log(const tdouble x_val, const bool safeCalc)
{
    rv_ptr->eval_para();
    return rv_ptr->calc_pdf_x_log(x_val,safeCalc);
}

const tdouble flxPyRV::cdf(const tdouble x_val, const bool safeCalc)
{
    rv_ptr->eval_para();
    return rv_ptr->calc_cdf_x(x_val,safeCalc);
}

py::array_t<tdouble> flxPyRV::cdf_array(py::array_t<tdouble> arr, const bool safeCalc)
{
    rv_ptr->eval_para();

    // Access the input data as a raw pointer
    py::buffer_info buf_info = arr.request();
    tdouble* input_ptr = static_cast<tdouble*>(buf_info.ptr);

    // Get the size of the input array
    size_t size = buf_info.size;

    // Allocate memory for the return array
    auto res_buf = py::array_t<tdouble>(size);

    // Get the buffer info to access the underlying return data
    py::buffer_info res_buf_info = res_buf.request();
    tdouble* res_ptr = static_cast<tdouble*>(res_buf_info.ptr);

    // Fill the array with values
    for (size_t i = 0; i < size; ++i) {
        res_ptr[i] = rv_ptr->calc_cdf_x(input_ptr[i],safeCalc);    // TODO avoid re-evaluating the parameters of the random variable
    }

    // Return the array
    return res_buf;


}

const tdouble flxPyRV::icdf(const tdouble p)
{
    ensure_is_a_basic_rv();
    rv_ptr_->eval_para();
    const tdouble y = rv_InvPhi_noAlert( p );
    return rv_ptr_->transform_y2x(y);
}

const tdouble flxPyRV::sf(const tdouble x_val, const bool safeCalc)
{
    rv_ptr->eval_para();
    return rv_ptr->calc_sf_x(x_val,safeCalc);
}

const tdouble flxPyRV::entropy()
{
    rv_ptr->eval_para();
    return rv_ptr->calc_entropy();
}

const tdouble flxPyRV::mean()
{
    rv_ptr->eval_para();
    return rv_ptr->get_mean_current_config();
}

const tdouble flxPyRV::sd()
{
    rv_ptr->eval_para();
    return rv_ptr->get_sd_current_config();
}

const tdouble flxPyRV::median()
{
    rv_ptr->eval_para();
    return rv_ptr->get_median_current_config();
}

const tdouble flxPyRV::mode()
{
    rv_ptr->eval_para();
    return rv_ptr->get_mode_current_config();
}

const bool flxPyRV::check_x(const tdouble xV)
{
    rv_ptr->eval_para();
    return rv_ptr->check_x(xV);
}

const tdouble flxPyRV::get_HPD(const tdouble p)
{
    ensure_is_a_basic_rv();
    rv_ptr_->eval_para();
    return rv_ptr_->get_HPD(p);
}

py::dict flxPyRV::info()
{
    rv_ptr->eval_para();
    return rv_ptr->info();
}




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

py::array_t<tdouble> flxPyRVset::get_values(const std::string mode)
{
    // process mode
        enum RBRVvecGetType { x, y, mean, sd };
        RBRVvecGetType gType;
        if (mode=="x") {
            gType = RBRVvecGetType::x;
        } else if (mode=="y") {
            gType = RBRVvecGetType::y;
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
        if ( (gType==y&&NRV==0) || NOX==0 ) {
          std::ostringstream ssV;
          ssV << "The set '" << name_of_set << "' does not contain any random variables.";
          throw FlxException("FlxObjRBRV_vec_get::task_2", ssV.str() );
        }
        const tuint N = (gType==y)?NRV:NOX;
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
        case y:
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

flxPyRVset rbrv_set(py::dict config, py::list rv_list)
{
    // eval and check name
        std::string set_name = parse_py_para_as_word("name",config,true,true);
        RBRV_entry_read_base::generate_set_base_check_name(FlxEngine->dataBox.rbrv_box, set_name);
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
                    RBRV_entry_read_base::generate_set_base(FlxEngine->dataBox.rbrv_box, set_parents,parents);
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
        RBRV_entry_RV_base* entry = nullptr;
        FlxFunction* csVal = nullptr;
        std::string csNam;
        bool csFix = false;
        const std::string family = set_name + "::";
        try {
            for (tuint i=0;i<Nrv;++i) {
                // cast config of RV into a dict
                    std::ostringstream ssV;
                    ssV << "configuration for random variable " << (i+1);
                    py::dict rv_config = parse_py_obj_as_dict(rv_list[i], ssV.str());
                // generate RV from config
                    if (crtr->get_Nataf_evalOnce()) {
                        rv_config["eval_once"] = true;
                    }
                    entry = parse_py_obj_as_rv(rv_config, true, i, family, ssV.str());
                // correlation (for Rosenblatt)
                    if (rv_config.contains("corr")) {
                        py::dict corr_config = parse_py_obj_as_dict(rv_config["corr"], "'corr' in " + ssV.str());
                        csNam = family + parse_py_para_as_string("rv_name", corr_config, true);
                        csVal = parse_py_para("value", corr_config, true);
                        csFix = parse_py_para_as_bool("fix", corr_config, false, false);
                    }
                // register the random variable in the set
                    try {
                        crtr->add_entry(FlxEngine->dataBox.rbrv_box,entry, csVal, csNam, csFix);
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
        RBRV_set_base* set_ptr = crtr->register_set(FlxEngine->dataBox.rbrv_box,true);
        flxPyRVset res(set_ptr, set_name);
    finalize_call();
    return res;
}

flxPyRV get_rv_from_set(const std::string& rv_name)
{
    RBRV_entry* rv_ptr = FlxEngine->dataBox.rbrv_box.get_entry(rv_name,true);
    flxPyRV res(rv_ptr);
    finalize_call();
    return res;
}

flxPySampler::flxPySampler(py::list rvsets)
: RndBox(nullptr)
{
    std::vector<std::string> set_str_vec;
    set_str_vec.reserve(rvsets.size());
    for (ssize_t i = 0; i < rvsets.size(); ++i) {
        py::object obj = rvsets[i];
        const std::string entry = parse_str_as_word(parse_py_obj_as_string(obj, "list entry"), true);
        set_str_vec.push_back(entry);
    }
    RndBox = new RBRV_constructor(set_str_vec,FlxEngine->dataBox.rbrv_box);
}

flxPySampler::~flxPySampler()
{
    if (RndBox) delete RndBox;
}

void flxPySampler::sample()
{
    RndBox->gen_smp();
}


// #################################################################################
// advanced features
// #################################################################################

tdouble eval_fun(py::object expr)
{
    FlxFunction* fun = parse_function(expr,"parameter 'expr' of 'flx.eval_fun'");
    tdouble res;
    try {
        res = fun->calc();
    } catch (FlxException& e) {
        GlobalVar.slogcout(1) << std::endl << e.what() << std::endl;
        res = std::numeric_limits<tdouble>::quiet_NaN();
    }
    return res;
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
        m.def("cdfn_inv", &rv_InvPhi, "inverse of the CDF of the standard Normal distribution");

    // ====================================================
    // random variables
    // ====================================================
        py::class_<flxPyRV>(m, "rv")
            .def(py::init<py::dict>())
            .def("get_name", &flxPyRV::get_name, "get name of random variable")
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
            .def("sf", &flxPyRV::sf, pybind11::arg("x_val"), pybind11::arg("safeCalc") = true, "returns the survival function of the random variable; i.e., 1-cdf(x_val)")
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
            .def("get_values", &flxPyRVset::get_values, pybind11::arg("mode")="x", "returns a vector of quantities of all entries contained in the set of random variables");

        m.def("rv_set", &rbrv_set, "creates a set of general random variables");
        m.def("get_rv_from_set", &get_rv_from_set, "retrieve a random variable from a set of random variables");

        py::class_<flxPySampler>(m, "sampler")
            .def(py::init<py::list>())
            .def("sample", &flxPySampler::sample, "generate a random sample for a collection of sets of random variables");

    // ====================================================
    // Advanced features
    // ====================================================
        m.def("eval_fun", &eval_fun, "Evaluates an expression and returns the result.");

    // ====================================================
    // only for debugging purposes (TODO remove at some point)
    // ====================================================
        m.def("add", &add, "A function that adds two numbers");
        m.attr("the_answer") = 42;
}


