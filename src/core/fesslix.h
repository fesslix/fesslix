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

#include "flxobjcommon.h"
#include "flxdefault.h"

#ifdef fesslix_fesslix_CPP
  bool fesslix_logOutput=false;        // true, if output is logged and log-stream is not cout
  bool fesslix_useScreen=DEFAULT_FLX_USESCREEN;
  std::string fesslix_gaussFile=DEFAULT_GAUSS_FILE;
  std::string fesslix_logFile=DEFAULT_LOG_FILE;
  bool fesslix_logTrunc=DEFAULT_LOG_TRUNC;
#else
  extern bool fesslix_logOutput;
  extern bool fesslix_useScreen;
  extern std::string fesslix_gaussFile;
  extern std::string fesslix_logFile;
  extern bool fesslix_logTrunc;
#endif

/**
* @brief The Main-Class
*/
class FesslixMain : public FlxReaderBase {
  public:
    // FlxAddOnManager AddOnManager;
    FlxData dataBox;
    FlxFunctionReader *funReader;
    FlxEvaluateCmd *EvaluateCmd;
    FlxDefParaBox AllDefParaBox;
  private:
    static int Cinst;
    
    void load_FlxReaders(FlxCreateObjReaders* cor);
  public:
    FesslixMain ();
    ~FesslixMain ();
    
    /**
    * @brief Initialize some streams in Fesslix
    * @param coutV (cout)
    * @param cerrV (cerr)
    * @param gaussI (the input file stream to read in the gauss points - will be autom. deleted after input!)
    */
    void initialize(const ostreamp coutV, const ostreamp cerrV, const std::string& gaussFile);
    
    /**
    * @brief Evaluate the File (this will take some time - depending on the size of the file).
    * @Note this is a top-level command; don't use this recursively!
    */
    void evaluate ();
    void evaluate (ReadStream *readerV); 
    /**
    * @brief Read - but not evaluate - the File
    * @Note this is a top-level command; don't use this recursively!
    */
    FlxObjBase* read_objects (ReadStream *readerV);
    
    /**
    * @brief Evaluate an expression - when 'evaluate' is already active at a higher level
    */
    void evaluate_expr(const std::string exprStr);
    
    #if FLX_DEBUG
      void test ();
    #endif
};
#ifdef fesslix_fesslix_CPP
  FesslixMain* FlxEngine = NULL;
#else
  extern FesslixMain* FlxEngine;
#endif



std::string flx_get_version();
void flx_print_current_time(std::ostream& os);

/**
* @brief initialize Fesslix
* @returns -1 if end of function was reached
*/
const int flxengine_init();
/**
* @brief deallocate the FesslixEngine
*/
void flxengine_unload();
/**
* @brief must be executed after flxengine_unload
*/
void flx_delete_FLXcout();

/**
* @brief evaluates the parameter file FileName
* @param is_file: true, if 'Expr' is FileName; false, if 'Expr' is command
* @param catch_all: true, if all errors should be catched
* @return only relevant if catch_all=true
*           if RETURN_RETURN_SUCCESS, the program can continue
*         otherwise, the execution should be stopped
*/
const int flx_parse_file(const std::string& Expr, const bool is_file, const bool catch_all, bool do_log=false, int tabWidth=8 );

/**
* @brief same as 'flx_parse_file', however, the file is only read but not executed
*/
FlxObjBase* flx_read_file(const std::string& Expr, const bool is_file, const bool catch_all, bool do_log=false, int tabWidth=8 );

void flx_print_RETURN_SUCCESS();
void flx_print_RETURN_ERROR();



  
  
  


