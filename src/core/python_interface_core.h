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
#include "flxobjrandom.h"



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
// Standard functions
// #################################################################################

void set_const(const std::string const_name, const tdouble value);
void set_var(const std::string var_name, py::object fun);


// #################################################################################
// random variables
// #################################################################################


// #################################################################################
// sets of random variables
// #################################################################################

class flxPyRVset {
  private:
    RBRV_set_base* rvset_ptr;   // memory needs to be managed externally
    std::string name_of_set;
  public:
    flxPyRVset(RBRV_set_base* rvset_ptr, const std::string& name_of_set) : rvset_ptr(rvset_ptr), name_of_set(name_of_set) {}
    flxPyRVset() = delete;
    flxPyRVset(flxPyRVset& rhs);
    flxPyRVset(flxPyRVset&& rhs);
    ~flxPyRVset() {}

    flxPyRVset& operator=(const flxPyRVset& rhs);

    const std::string get_name() const;
    const tuint get_NRV() const;
    const tuint get_NOX() const;

    py::array_t<tdouble> get_values(const std::string mode="x");
    py::list get_rvs();
    void set_y_vec(py::array_t<tdouble> arr);
    void set_x_vec(py::array_t<tdouble> arr);
    const tdouble pdf_log(py::array_t<tdouble> arr);

    /**
    * @brief evaluates the realization of a random process at time t with given power spectral density function
    */
    py::array_t<tdouble> eval_rp_psd(py::array_t<tdouble> arr);

};

flxPyRVset rbrv_set(py::dict config, py::list rv_list);
flxPyRVset rbrv_set_noise(py::dict config, py::dict rv_config);
flxPyRVset rbrv_set_proc(py::dict config, py::dict rv_config);
flxPyRVset rbrv_set_psd(py::dict config);
flxPyRVset rbrv_set_sphere(py::dict config);
flxPyRVset rbrv_set_vfun(py::dict config);
flxPyRVset rbrv_set_dirichlet(py::dict config);
flxPyRVset rbrv_set_multinomial(py::dict config);

flxPyRV get_rv_from_set(const std::string& rv_name);



// #################################################################################
// advanced features
// #################################################################################

tdouble eval_fun(py::object expr);
void eval_code(py::object expr);
py::array_t<tdouble> eval_vfun(const tuint N, py::object expr);
