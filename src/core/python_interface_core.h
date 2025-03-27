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


#include "flxmath.h"
#include "flxrbrv_rvs.h"

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;


/**
* @brief ensure that the Fesslix engine is up and running; otherwise, an error is thrown
*/
void check_engine_state();
/**
* @brief run this after a call that uses the engine
*/
void finalize_call();


// #################################################################################
// 'global' functions
// #################################################################################

/**
* @brief Convert a double into a string
*/
std::string Double2String(tdouble a);

/**
* @brief output configuration options of Fesslix
*/
void print_info();


// #################################################################################
// load configuration
// #################################################################################

void set_exe_dir(const std::string& exe_dir);
void process_config(py::dict config);


// #################################################################################
// logging
// #################################################################################

class PythonLoggerBuffer : public std::stringbuf {
  private:
    py::object py_logger;
    int log_Level;

    PythonLoggerBuffer( const PythonLoggerBuffer& );
    PythonLoggerBuffer& operator=( const PythonLoggerBuffer& ); // Kopieren verhindern

    const char* logLevelToString(int logLevel);
  public:
    PythonLoggerBuffer ( ) : py_logger(py::none()), log_Level(0) { }

    int sync() override;

    void set_logLevel(const int logLevel_);
    void set_logger(py::object py_logger_) { py_logger = std::move(py_logger_); }
};

class flxStreamPy : public std::ostream {
  private:
    PythonLoggerBuffer theBuf;

  public:
    flxStreamPy ( ) : std::ostream(&theBuf) { }

    void set_logLevel(const int logLevel_) { theBuf.set_logLevel(logLevel_); }
    void set_logger(py::object py_logger_) { theBuf.set_logger(std::move(py_logger_)); }
};

class flxPyLogger : public flxLoggerBase {
  private:
    flxStreamPy streamPy;

  public:
    flxPyLogger() {};
    virtual ~flxPyLogger() {};

    /**
    * @brief returns the logging-stream
    * (1) alerts and errors are logged (ALERT)
    * (2) warnings are logged (WARNING)
    * (3) normal but significant information is logged (NOTICE)
    * (4) informational (INFO)
    * (5) debug-level messages (DEBUG)
    */
    virtual std::ostream& slog(const int logLevel_);

    void set_logger(py::object logger_obj);
};

/**
* @brief set Python-object as logger
*/
void set_logger(py::object logger_obj);

/**
* @brief log some <message> with <logLevel>
*/
void slog(const std::string& message, int logLevel=4);


// #################################################################################
// Fesslix Engine
// #################################################################################

int load_engine();


// #################################################################################
// random variables
// #################################################################################

class flxPyRV {
  private:
    RBRV_entry_RV_base* rv_ptr;
  public:
    flxPyRV(py::dict config);
    ~flxPyRV();

    const std::string get_type() const;
    const tdouble x2y(const tdouble x_val);
    const tdouble y2x(const tdouble y_val);

    const tdouble pdf(const tdouble x_val, const bool safeCalc);
    const tdouble pdf_log(const tdouble x_val, const bool safeCalc);
    const tdouble cdf(const tdouble x_val, const bool safeCalc);
    const tdouble icdf(const tdouble p);
    const tdouble sf(const tdouble x_val, const bool safeCalc);
    const tdouble entropy();
    const tdouble mean();
    const tdouble sd();
    const tdouble median();
    const tdouble mode();
    const bool check_x(const tdouble xV);
    const tdouble get_HPD(const tdouble p);
    py::dict info();
};




