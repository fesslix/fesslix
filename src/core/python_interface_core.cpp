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

flxPyRV::flxPyRV(py::dict config)
: rv_ptr(nullptr)
{
    // retrieve type
        if (config.contains("type")==false) {
            throw FlxException_NeglectInInteractive("flxPyRV::flxPyRV_01", "Python <dict> does not contain key 'type'.");
        }
        if (py::isinstance<py::str>(config["type"])==false) {
            throw FlxException_NeglectInInteractive("flxPyRV::flxPyRV_02", "Key 'type' in Python <dict> is not of type 'str'.");
        }
        const std::string rv_type = config["type"].cast<std::string>();
    // retrieve name
        std::string rv_name = "name_unspecified";
        if (config.contains("name")) {
            if (py::isinstance<py::str>(config["name"])) {
                rv_name = config["name"].cast<std::string>();
            } else {
                throw FlxException_NeglectInInteractive("flxPyRV::flxPyRV_03", "Key 'name' in Python <dict> is not of type 'str'.");
            }
        }
    // select rbrv-class based on type
        if (rv_type=="stdn") {
            rv_ptr = new RBRV_entry_RV_stdN(rv_name,0,config);
        } else if (rv_type=="normal") {
            rv_ptr = new RBRV_entry_RV_normal(rv_name,0,config);
        // } else if (rv_type=="logn") {
        //     rv_ptr = new RBRV_entry_read_logn(readName,readBrakets);
        // } else if (rv_type=="uniform") {
        //     rv_ptr = new RBRV_entry_read_uniform(readName,readBrakets);
        // } else if (rv_type=="gumbel") {
        //     rv_ptr = new RBRV_entry_read_Gumbel(readName,readBrakets);
        // } else if (rv_type=="normal_trunc") {
        //     rv_ptr = new RBRV_entry_read_normal_trunc(readName,readBrakets);
        // } else if (rv_type=="beta") {
        //     rv_ptr = new RBRV_entry_read_beta(readName,readBrakets);
        // } else if (rv_type=="exponential") {
        //     rv_ptr = new RBRV_entry_read_exponential(readName,readBrakets);
        // } else if (rv_type=="gamma") {
        //     rv_ptr = new RBRV_entry_read_gamma(readName,readBrakets);
        // } else if (rv_type=="poisson") {
        //     rv_ptr = new RBRV_entry_read_Poisson(readName,readBrakets);
        // } else if (rv_type=="binomial") {
        //     rv_ptr = new RBRV_entry_read_Binomial(readName,readBrakets);
        // } else if (rv_type=="cauchy") {
        //     rv_ptr = new RBRV_entry_read_Cauchy(readName,readBrakets);
        // } else if (rv_type=="weibull") {
        //     rv_ptr = new RBRV_entry_read_Weibull(readName,readBrakets);
        // } else if (rv_type=="chisquared") {
        //     rv_ptr = new RBRV_entry_read_ChiSquared(true,readName,readBrakets);
        // } else if (rv_type=="chi") {
        //     rv_ptr = new RBRV_entry_read_ChiSquared(false,readName,readBrakets);
        } else if (rv_type=="studentst") {
            rv_ptr = new RBRV_entry_RV_StudentsT(rv_name,0,config);
        } else if (rv_type=="studentstgen") {
            rv_ptr = new RBRV_entry_RV_StudentsT_generalized(rv_name,0,config);
        // } else if (rv_type=="laplace") {
        //     rv_ptr = new RBRV_entry_read_Laplace(readName,readBrakets);
        // } else if (rv_type=="usertransform") {
        //     rv_ptr = new RBRV_entry_read_UserTransform(readName,readBrakets);
        // } else if (rv_type=="truncated") {
        //     rv_ptr = new RBRV_entry_read_Truncated(readName,readBrakets);
        // } else if (rv_type=="maxmintransform") {
        //     rv_ptr = new RBRV_entry_read_maxminTransform(readName,readBrakets);
        // } else if (rv_type=="bayda") {
        //     rv_ptr = new RBRV_entry_read_bayDA(readName,readBrakets);
        } else {
            std::ostringstream ssV;
            ssV << "Unknown random variable type '" << rv_type << "'.";
            throw FlxException("flxPyRV::flxPyRV_50", ssV.str() );
        }
}

flxPyRV::~flxPyRV()
{
    delete rv_ptr;
}

const std::string flxPyRV::get_type() const
{
    return rv_ptr->get_type();
}

const tdouble flxPyRV::x2y(const tdouble x_val)
{
    return rv_ptr->transform_x2y(x_val);
}

const tdouble flxPyRV::y2x(const tdouble y_val)
{
    return rv_ptr->transform_y2x(y_val);
}

const tdouble flxPyRV::pdf(const tdouble& x_val, const bool safeCalc)
{
    return rv_ptr->calc_pdf_x(x_val,safeCalc);
}

const tdouble flxPyRV::pdf_log(const tdouble& x_val, const bool safeCalc)
{
    return rv_ptr->calc_pdf_x_log(x_val,safeCalc);
}

const tdouble flxPyRV::cdf(const tdouble& x_val, const bool safeCalc)
{
    return rv_ptr->calc_cdf_x(x_val,safeCalc);
}

const tdouble flxPyRV::sf(const tdouble& x_val, const bool safeCalc)
{
    return rv_ptr->calc_sf_x(x_val,safeCalc);
}

const tdouble flxPyRV::entropy()
{
    return rv_ptr->calc_entropy();
}

const tdouble flxPyRV::mean()
{
    return rv_ptr->get_mean_current_config();
}

const tdouble flxPyRV::sd()
{
    return rv_ptr->get_sd_current_config();
}

const tdouble flxPyRV::median()
{
    return rv_ptr->get_median_current_config();
}

const tdouble flxPyRV::mode()
{
    return rv_ptr->get_mode_current_config();
}

const bool flxPyRV::check_x(const tdouble xV)
{
    return rv_ptr->check_x(xV);
}

const tdouble flxPyRV::get_HPD(const tdouble p)
{
    return rv_ptr->get_HPD(p);
}

const std::string flxPyRV::info()
{
    std::ostringstream ssV;
    rv_ptr->info(ssV);
    return ssV.str();
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
    // random variables
    // ====================================================
        py::class_<flxPyRV>(m, "rv")
            .def(py::init<py::dict>())
            .def("get_type", &flxPyRV::get_type, "get type of random variable")
            .def("x2y", &flxPyRV::x2y, "transformation from 'original space' to standard normal space")
            .def("y2x", &flxPyRV::y2x, "transformation from standard normal space into 'original space'")
            .def("pdf", &flxPyRV::pdf, pybind11::arg("x_val"), pybind11::arg("safeCalc") = true, "evaluates the pdf of the random variable at x_val")
            .def("pdf_log", &flxPyRV::pdf_log, pybind11::arg("x_val"), pybind11::arg("safeCalc") = true, "evaluates the log-pdf of the random variable at x_val")
            .def("cdf", &flxPyRV::cdf, pybind11::arg("x_val"), pybind11::arg("safeCalc") = true, "evaluates the cdf of the random variable at x_val")
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
    // only for debugging purposes (TODO remove at some point)
    // ====================================================
        m.def("add", &add, "A function that adds two numbers");
        m.attr("the_answer") = 42;
}


