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

#include "pdouble.h"
#include "flxplatform.h"
#include "flxexception.h"
#include "flxMemoryManager.h"

#include <boost/format.hpp>

#include <valarray>
typedef std::valarray<tdouble> tVec;
// typedef std::valarray<pdouble> pVec;
typedef std::valarray<bool> bVec;
typedef std::valarray<tuint> iVec;
typedef std::valarray<tnlong> lVec;

#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>  // For NumPy support
namespace py = pybind11;

FLXLIB_EXPORT size_t bVec_NumbTrue(const bVec& bv);

FLXLIB_EXPORT void FlxError(bool errSerious, std::string errnumber, std::string Titel = "", std::string Msg = "");

#if FLX_DEBUG
#define FLXMSG(M,I) GlobalVar.slogcout(I) << M << std::endl;
#else
  #define FLXMSG(M,I)
#endif


class FLXLIB_EXPORT flxLoggerBase {
  public:
    virtual ~flxLoggerBase() {};

    virtual std::ostream& slog(const int logLevel_) = 0;
};


class FLXLIB_EXPORT strGlobalVar {
  private:
    const ostreamp true_stdcout;              // this pointer will always point to std::cout
    const ostreamp true_cerr;                 // this pointer will always point to std::cerr
    ostreamp sdummy;
    ostreamp slogP;
    /**
    * @brief stores a pointer to the standard output stream
    */
    ostreamp stdcout;
    ostreamp stdcerr;
    /**
    * @brief distributes a stream to cout and log!
    */
    ostreamp slogcoutP;
    
    tuint D2S_Prec;
    tuint D2S_Type;
    tdouble D2S_BValU;
    tdouble D2S_BValL; 
    bool D2S_Del0;
    bool D2S_DelP;
    tdouble TOLval;
    /**
    * @brief stores if the pre-log is active
    */
    bool prelog_active;
    /**
    * @brief stores if the pre-log is active
    */
    std::string prelog_stream;
    /**
    * @brief this is the default loglevel
    */
    int logLevel_strong;
    /**
    * @brief stores the directory of the binary (executable)
    */
    std::string exe_dir;
    /**
    * @brief truncates new_line at the end of prelog
    */
    void prelog_prepare_write();
    /**
    * @brief true, if fesslix is run in leak-check-mode
    */
    bool leak_check;
    flxLoggerBase* Logger_ptr;
    
  public:
    tuint logLev_counter;
    FlxMemoryManager MemMngr;
    
    /**
    * @brief this is the temporary loglevel - can be reset by objects
    */
    int logLevel;
    tuint MT19937_init_calls;
    bool MT19937_init_RAND;
    bool MT19937_init_seed;
    tuint MT19937_init_seedvalue;
    size_t LegendrePolyH_init_numb;
    bool prgBar;                // true if progress-bar is activated for timeconsuming jobs
    FlxAlert alert;
    const tdouble sqrtEps;                // square-root of precision of epsilon
    tuint max_parallel_threads;
    void set_TOL(const tdouble tolV);
    /**
    * @brief if a value is smaller or equal than TOL, it can be considered to be zero
    */
    const tdouble& TOL() const { return TOLval; }
    
    strGlobalVar();
    ~strGlobalVar();

    /**
    * @brief updates Logger_ptr
    * @note the memory management has to be taken care of externally
    */
    void set_logger(flxLoggerBase& Logger) { Logger_ptr = &Logger; }
    const bool has_logger() const { return Logger_ptr!=nullptr; }
    /**
    * @brief updates the log and the cout stream
    * @note the memory management has to be taken care of externally
    */
    void set_slogcout(ostreamp slogV, ostreamp stdcoutV);
//     void set_slog(ostreamp slogV);
//     void set_stdcout(ostreamp stdcoutV) { stdcout = stdcoutV; }
    void set_stdcerr(ostreamp stdcerrV) { stdcerr= stdcerrV; }
    const ostreamp& get_cout() { return stdcout; }                     // TODO obsolete? delete?
    const ostreamp& get_cerr() { return stdcerr; }                     // TODO obsolete? delete?
    const ostreamp& get_log() { return slogP; }
    const ostreamp& get_true_cout() { return true_stdcout; }           // TODO obsolete? delete?
    const ostreamp& get_true_cerr() { return true_cerr; }              // TODO obsolete? delete?
    const std::string& get_exe_dir() { return exe_dir; }               // TODO obsolete? delete?
    /**
    * @returns the directory of the executable - must end with an '/' !!!
    */
    void set_exe_dir(const std::string& exe_dirV) { exe_dir = exe_dirV; }
    /**
    * @brief returns the logging-stream
    * (1) alerts and errors are logged (ALERT)
    * (2) warnings are logged (WARNING)
    * (3) normal but significant information is logged (NOTICE)
    * (4) informational (INFO)
    * (5) debug-level messages (DEBUG)
    */
    std::ostream& slog(const int logLevel_);
    std::ostream& slogcout(const int logLevel_);
    void slog_flush() { slogP->flush(); }                             // TODO obsolete? delete?                      // TODO obsolete? delete?
    std::ostream& slog_dummy() { return *sdummy; }
    /**
    * @returns true if log-stream is not cout
    */
    const bool check_logNOTcout() { return (slogP != stdcout); };     // TODO obsolete? delete?
    const std::string get_flxPrompt() { return "fesslix:> "; }

    /**
    * @returns send double to string in maximum precision
    */
    const std::string Double2String_maxPrec(tdouble dv);
    /**
    * @brief transforms a double to a string
    * @param dv the double value to transfrom
    * @param checkTOL true, if values close to zero should be considered as zero
    * @param prec the precision of floating-point numbers to use
    *                 -1: use the default setting
    *                 all other (positive) values: the precision of floating-point numbers
    * @param fixW output a string of fixed length
    *                 -1: no, output the string as short as possible
    *                 0: fixW is calculated depending on the precision value
    *                         fixW = prec+7
    *                 all other (positive) values: that's the fixed width
    *                 Note: if the string is longer than fixW, it is NOT truncated!
    */
    const std::string Double2String(tdouble dv, const bool checkTOL = false, const int prec=-1, int fixW=-1);
    void Double2String_setPrec(const tuint prec);
    void Double2String_setType(const tuint type);
    void Double2String_setBValU(const tdouble bv);
    void Double2String_setBValL(const tdouble bv);
    void Double2String_setDel0(const bool d0);
    void Double2String_setDelP(const bool dP);
    void Double2String_log();
    const int D2S_get_fixW(const int prec, const int fixW) const;
    const std::string Double2String_sci(tdouble dv, const int prec=-1, int fixW=-1);
    
    const boost::basic_format<char> D2S_totalPrec(const tdouble dv);
    /**
    * @brief set the strong loglevel
    */
    void logLevel_strong_set(const int logL);                                             // TODO obsolete? delete?
    /**
    * @brief deactivates/activates logging
    */
    void logLevel_log_deactivate(const bool deact);                                       // TODO obsolete? delete?
    /**
    * @brief set the temporary loglevel to the strong loglevel
    */
    void logLevel_strong_reset() { logLevel = logLevel_strong; logLev_counter=0; };       // TODO obsolete? delete?
    /**
    * @brief activates/deactivates the prelog
    */
    void prelog_activated(const bool b);
    /**
    * @brief returns true if the prelog is active
    */
    bool prelog_isActive() { return prelog_active; }
    /**
    * @brief write the pre-log to the log - and clear it afterwards
    */
    void prelog_write();
    /**
    * @brief clear the pre-log
    */
    void prelog_clear();
    /**
    * @brief force writing the pre-log to the log - do not clear it afterwards
    */
    std::string prelog_force_write();
    /**
    * @brief check if the pre-log is not empty
    */
    bool prelog_isNOTempty() { if(prelog_stream!="") return true; else return false; }
    /**
    * @brief append a character to the pre-log
    */
    void prelog_append(char c);
    /**
    * @brief append a string to the pre-log
    */
    void prelog_append(std::string str);
    /**
    * @brief sets leak_check to true
    */
    void set_leak_check_mode() { leak_check = true; }
    /**
    * @brief returns whether fesslix is run in leak-check mode
    */
    const bool is_leak_check() const { return leak_check; }
    
    friend class FlxOstreamBox;
};
#ifdef fesslix_flxglobal_CPP
  FLXLIB_EXPORT strGlobalVar GlobalVar;
#else
  FLXLIB_EXPORT extern strGlobalVar GlobalVar;
#endif


class FLXLIB_EXPORT FlxProgress {
  private:
    const bool active;
    std::ostream& op;
    tuint N;        // total size
    clock_t last_t;        // last time of action
    tuint last_p;
    bool is_running;
    tuint size_e;        // number of chars written during last step
    
    void tick_t(const tuint idx);
  public:
    FlxProgress(std::ostream& op, const bool doLog) : active(GlobalVar.prgBar&&doLog),op(op),last_p(0),is_running(false),size_e(0) { }
    
    void start(const tuint N_v);
    inline void tick(const tuint idx) { if (!active) return;
      if (tdouble(clock()-last_t)<0.2*CLOCKS_PER_SEC) return;
      tick_t(idx);
    };
    void stop();
};


FLXLIB_EXPORT const std::string bool2string(const bool b);
  
FLXLIB_EXPORT void fesslix_logInfo(std::ostream& lout);



  
