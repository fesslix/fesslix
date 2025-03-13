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


// #################################################################################
// 'global' functions
// #################################################################################

std::string Double2String(tdouble a) {
    return GlobalVar.Double2String(a);
}

void print_info()
{
    std::ostream& lout = GlobalVar.slogcout(3);
    fesslix_logInfo(lout);
    lout.flush();
}


// #################################################################################
// load configuration
// #################################################################################

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

void slog(int logLevel, const std::string& message) {
    GlobalVar.slogcout(logLevel) << message;
    GlobalVar.slogcout(logLevel).flush();
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

PYBIND11_MODULE(core, m) {
    // ====================================================
    // 'global' functions
    // ====================================================
        m.def("Double2String", &Double2String, "Convert a double into a string");   // TODO docu
        m.def("print_info", &print_info, "output configuration options of Fesslix");

    // ====================================================
    // load configuration
    // ====================================================
        m.def("process_config", &process_config, "Process a dictionary containing configuration options");

    // ====================================================
    // logging
    // ====================================================
        m.def("set_logger", &set_logger, "Set the logger object");
        m.def("slog", &slog, "Log a message at a specified level");

    // ====================================================
    // only for debugging purposes (TODO remove at some point)
    // ====================================================
        m.def("add", &add, "A function that adds two numbers");
        m.attr("the_answer") = 42;
}


