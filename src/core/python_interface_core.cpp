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


std::string Double2String(tdouble a) {
    return GlobalVar.Double2String(a);
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
    const std::string msg = str();
    if (!msg.empty()) {
        if (py_logger.is_none()) {
            throw FlxException_Crude("PythonLoggerBuffer::sync");
        }
        // Forward the message to the Python logger when buffer is flushed
        py_logger.attr(logLevelToString(log_Level))(msg);
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
    GlobalVar.slogcout(logLevel) << message << std::endl;
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
    GlobalVar.slogcout(1) << "Fesslix::core » loaded" << std::endl;
    // ====================================================
    // 'global' functions
    // ====================================================
        m.def("Double2String", &Double2String, "Convert a double into a string");

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


