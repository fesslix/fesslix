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

#define fesslix_fesslix_CPP

#include "fesslix.h"
#include "flxobjrandom.h"
#include "flxobjrbrv.h"
#include "flxobjmtx.h"
#include "flxdefault.h"
#include "flxobjects.h"
#include "flxstringfun_fun.h"
#include "flxMtx_Eigen.h"
#include "flxphys.h"

#include <fstream>
#include <limits>
#include <iostream>


using namespace std;

int FesslixMain::Cinst = 0;      // true, cout is std::cout (otherwise cout is redirected to the log-file)
std::ostream* flxcout = nullptr;
int flx_engine_init_count = 0;        

FesslixMain* FlxEngine_ptr = nullptr;

/**
* @brief execute this command at the very beginning
*/
void flxengine_load();
/**
* @brief opens the log-file and outputs the header of the log
*/
void flx_create_log(const std::string& logFile, const bool logTrunc);



  
  
FesslixMain::FesslixMain() : initialized(false)
{
  // check for multiple instances
    ++Cinst;
    if (Cinst > 1) {
      std::ostringstream ssV;
      ssV << "More than one instance of 'FlxData' created ...";
      throw FlxException("FesslixMain::FesslixMain", ssV.str() );
    }
    
  // create readers ...
    funReader = new FlxFunctionReader();
      FlxReadManager::set_funReader(funReader);
    set_ReadManager(&(dataBox.ReadManager));
  set_RndCreator_ptr( &(dataBox.RndCreator) );
  // FlxObjReadLoadLib::AddOnManager = &AddOnManager;

  // create the evaluation class
    EvaluateCmd = new FlxEvaluateCmd();

  // create ObjReaders
    // The LoadLib-Objects
      // EvaluateCmd->get_ObjReadBox()->insert("loadlib", new FlxObjReadLoadLib(EvaluateCmd->get_ObjReadBox()));
      // EvaluateCmd->get_ObjReadBox()->insert("loadlibfun", new FlxObjReadLoadLibFun(EvaluateCmd->get_ObjReadBox()));
    // Common Objects
      load_FlxReaders( new FlxCreateObjReaders_Common() );
    // Mtx Objects
      load_FlxReaders( new FlxCreateObjReaders_Mtx() );
    // RBRV Objects
      load_FlxReaders( new FlxCreateObjReaders_RBRV() );
    // FlxString Objects
      load_FlxReaders( new FlxCreateObjReaders_FlxString() );
    // FlxPhys Objects
      load_FlxReaders( new FlxCreateObjReaders_FlxPhys() );

}

void FesslixMain::initialize(const ostreamp coutV, const ostreamp cerrV, const string& gaussFile)
{
  if (initialized) return;
  // initialise output-streams
    dataBox.OstreamBox.insert("cout", coutV);
    dataBox.OstreamBox.insert("cerr", cerrV);
    dataBox.OstreamBox.insert("dummy",new flxDummyOstream());
    dataBox.OstreamBox.insert("log", &(GlobalVar.slog(0)) );
  // Read Gauss-Point Information
    dataBox.GaussInt.ReadGP(0,gaussFile);
  initialized = true;
}

void FesslixMain::load_FlxReaders(FlxCreateObjReaders* cor)
{
  cor->createFunReaders( &dataBox);    
  cor->createObjReaders( EvaluateCmd->get_ObjReadBox() );
  delete cor; 
}

FesslixMain::~FesslixMain () {
  // AddOnManager.delete_AddOns();        // free the loaded AddOns
  delete EvaluateCmd;
  delete funReader;
  --Cinst;
}

void FesslixMain::evaluate () {
  FlxObjBase* ob=NULL;
  if (reader == NULL) {        // Check if a reader is available
    ostringstream ssV_2;
    ssV_2 << "No reader available.";
    throw FlxException("FesslixMain::evaluate_1", ssV_2.str() ); 
  }
  try {
    while ( reader->getNextType() != ReadStream::ENDOFFILE ) {
      // print the prompt for output of input to log
        GlobalVar.prelog_append(GlobalVar.get_flxPrompt());
      ob = EvaluateCmd->evaluateCmd();
      ob->exec();
      GlobalVar.logLevel_strong_reset();
      GlobalVar.slog_flush();
      delete ob; ob=NULL;
    }
  } catch ( FlxEndE &e) {
    FLXMSG("FesslixMain::evaluate_2",1);
    GlobalVar.logLevel_strong_reset();
    if (ob) delete ob;
    throw;
  } catch ( FlxExitE &e) {
    FLXMSG("FesslixMain::evaluate_3",1);
    GlobalVar.logLevel_strong_reset();
    GlobalVar.slog(4) << "exit: 'exit;' was called" << std::endl;
    if (ob) delete ob;
  } catch ( FlxException &e ) {
    FLXMSG("FesslixMain::evaluate_4",1);
    GlobalVar.logLevel_strong_reset();
    if (ob) delete ob;
    throw;
  }
}

void FesslixMain::evaluate(ReadStream* readerV)
{
  dataBox.ReadManager.push(readerV);
  try {
    evaluate();
    dataBox.ReadManager.pop();
  } catch ( FlxEndE &e) {
    FLXMSG("FesslixMain::evaluate_101",1);
    dataBox.ReadManager.pop();
    throw;
  } catch ( FlxException &e ) {
    FLXMSG("FesslixMain::evaluate_102",1);
    dataBox.ReadManager.pop();
    throw;
  }
}

FlxObjBase* FesslixMain::read_objects(ReadStream* readerV)
{
  dataBox.ReadManager.push(readerV);
  FlxCodeBlock* obR = new FlxCodeBlock(true);
  try {
    if (reader == NULL) {        // Check if a reader is available
      ostringstream ssV_2;
      ssV_2 << "No reader available.";
      throw FlxException("FesslixMain::read_objects_1", ssV_2.str() ); 
    }
    while ( reader->getNextType() != ReadStream::ENDOFFILE ) {
      FlxObjBase* ob = EvaluateCmd->evaluateCmd();
      obR->attach_obj(ob);
    }
    dataBox.ReadManager.pop();
    return obR;
  } catch ( FlxException &e) {
    FLXMSG("FesslixMain::read_objects_101",1);
    delete obR;
    dataBox.ReadManager.pop();
    throw;
  }
}

void FesslixMain::evaluate_expr(const string exprStr)
{
  ReadStream Freader(exprStr);
  dataBox.ReadManager.push(&Freader);
  FlxObjBase* ob=NULL;
  try {
    while ( Freader.getNextType() != ReadStream::ENDOFFILE ) {
      // print the prompt for output of input to log
        GlobalVar.prelog_append(GlobalVar.get_flxPrompt());
      ob = EvaluateCmd->evaluateCmd();
      ob->exec();
      delete ob; ob=NULL;
    }
  } catch ( FlxExitE &e) {
    FLXMSG("FesslixMain::evaluate_expr_1",1);
    if (ob) delete ob;
  } catch ( FlxEndE &e) {
    FLXMSG("FesslixMain::evaluate_expr_2",1);
    if (ob) delete ob;
    dataBox.ReadManager.pop();
    throw;
  } catch (...) {
    if (ob) delete ob;
    dataBox.ReadManager.pop();
    throw;
  };
  dataBox.ReadManager.pop();
}


#if FLX_DEBUG
void FesslixMain::test()
{  
  GlobalVar.slog(5) << "Fesslix: test-section ------------ START -----------------------" << std::endl;
  
  
  GlobalVar.slog(5) << "Fesslix: test-section ------------ END -------------------------" << std::endl;
}
#endif


void flx_create_log(const std::string& logFile, const bool logTrunc) {
  // check whether a 'special' logger is already registered
    if (GlobalVar.has_logger()) return;
  if (logFile != "{screen}") {
    std::ofstream* ofLog;
    if ( logTrunc ) {
      ofLog = new std::ofstream(logFile.c_str());
    } else {
      ofLog = new std::ofstream(logFile.c_str(), std::ios_base::out | std::ios_base::app); 
    }
    if ( ofLog == NULL || !ofLog->is_open() ) {
      std::ostringstream ssV;
      ssV << "File (" << logFile << ") could not be opened.";
      throw FlxException("create_log_1", ssV.str() );
    } else {
      if (fesslix_useScreen) {
        GlobalVar.set_slogcout(ofLog,GlobalVar.get_cout());
      } else {
        GlobalVar.set_slogcout(ofLog,ofLog);
        GlobalVar.set_stdcerr(ofLog);
      }
    }
  }
  
  if (logFile == "{screen}" && fesslix_useScreen==false) {
    throw FlxException("flx_create_log_2","Configuration options 'flx.usescreen=false' and 'log.file=\"{screen}\" are not allowed simultaneously");
  }
}

void flx_print_current_time_help(tuint i, ostream& os) {
  if (i<10) {
    os << "0";
  } 
  os << i;
}

void flx_print_current_time(ostream& os)
{
  time_t Zeitstempel = time(0);
  tm *nun = localtime(&Zeitstempel);
  os << nun->tm_year+1900 << "-";
  flx_print_current_time_help(nun->tm_mon+1,os);
  os << "-";
  flx_print_current_time_help(nun->tm_mday,os);
  os << " - ";
  flx_print_current_time_help(nun->tm_hour,os);
  os << ':';
  flx_print_current_time_help(nun->tm_min,os);
}

void flx_delete_FLXcout() {
  if (flx_engine_init_count>0) return;
  if (flxcout==nullptr) return;
  flxStreamAlloc* fsa = dynamic_cast<flxStreamAlloc*>(flxcout);
  if (fsa) {
    delete fsa;
  }
  flxcout = nullptr;
}

string flx_get_version()
{
  std::ostringstream ssV;
  ssV << FLX_VERSION;
  return ssV.str();
}

void flxengine_load()
{
  FlxEngine_ptr = new FesslixMain;
}

void flxengine_unload()
{
  if (flx_engine_init_count==0) return;
  flx_engine_init_count--;
  if (flx_engine_init_count>0) return;
  if (FlxEngine_ptr) {
    delete FlxEngine_ptr;
    FlxEngine_ptr = nullptr;
  }
}

void check_engine_state()
{
    if (FlxEngine_ptr==nullptr) {
        throw FlxException_NeglectInInteractive("check_engine_state","Fesslix Engine is not running", "Please start the engine using load_engine()");
    }
}

FesslixMain& FlxEngine()
{
  check_engine_state();
  return *FlxEngine_ptr;
}

FesslixMain* get_FlxEngine_ptr()
{
  return FlxEngine_ptr;
}


const int flxengine_init()
{
  flx_engine_init_count++;
  if (flx_engine_init_count>1) return RETURN_SUCCESS;

  // run pre-checks
    if (flxcout) {
      throw FlxException_Crude("flxengine_init_01");
    }
  // set some global options
    GlobalVar.logLevel_strong_set(GlobalVar.logLevel);

  // Deal with the log-files
    flx_create_log(fesslix_logFile, fesslix_logTrunc);

    bool logout = fesslix_logOutput;
    if (!GlobalVar.check_logNOTcout()) logout=false;
    if (logout) {
      flxcout = new flxStreamAlloc(GlobalVar.get_cout(), GlobalVar.get_log());
    } else {
      flxcout = GlobalVar.get_cout();
    }

  // create main object
    flxengine_load();
    FlxEngine_ptr->initialize(flxcout, GlobalVar.get_cerr(), fesslix_gaussFile);

    #if FLX_DEBUG
    // testing
      FlxEngine_ptr->test();
    #endif

  return -1;
}

const int flx_parse_file(const string& Expr, const bool is_file, const bool catch_all, bool do_log, int tabWidth)
{
  ReadStream *reader = NULL;
  try {
    if (is_file) {
      GlobalVar.slog(3) << "Fesslix: start parsing parameter file '" << Expr << "'" << std::endl;
      reader = new ReadStream( Expr.c_str(), do_log, tabWidth );
    } else {
      GlobalVar.slog(3) << "Fesslix: start evaluating the command from the command line: " << Expr << std::endl;
      reader = new ReadStream(Expr, do_log, tabWidth);  
    }
    FlxEngine_ptr->evaluate(reader);
    delete reader; reader = NULL;
    GlobalVar.logLevel_strong_reset();
  } 
  catch ( FlxEndE const &e) {
    FLXMSG("flx_parse_file_1",1);
    if (reader) { delete reader; reader=NULL; }
    GlobalVar.logLevel_strong_reset();
    GlobalVar.slog(4) << "end: 'end' was called" << std::endl;
    return RETURN_SUCCESS;
  }
  catch ( FlxException_NeglectInInteractive &e) {
    FLXMSG("flx_parse_file_2",1);
    if (reader) { delete reader; reader=NULL; }
    *(GlobalVar.get_true_cerr()) << std::endl << e.what() << std::endl;
    GlobalVar.logLevel_strong_reset();
    if (GlobalVar.check_logNOTcout()) {
      GlobalVar.slog(1) << std::endl << e.what() << std::endl;
    }
    return RETURN_ERROR;
  }
  catch ( FlxException &e) {
    FLXMSG("flx_parse_file_3",1);
    if (reader) { delete reader; reader=NULL; }
    if (!catch_all) throw;
    *(GlobalVar.get_true_cerr()) << std::endl << e.what() << std::endl;
    if (GlobalVar.check_logNOTcout()) {
      GlobalVar.slog(1) << std::endl << e.what() << std::endl;
    }
    GlobalVar.logLevel_strong_reset();
    return RETURN_ERROR;
  }
  catch (std::exception const &e) {
    FLXMSG("flx_parse_file_4",1);
    if (reader) { delete reader; reader=NULL; }
    if (!catch_all) throw;
    *(GlobalVar.get_true_cerr()) << std::endl << "ERROR: Whoops, something went wrong: " << e.what() << std::endl;
    GlobalVar.logLevel_strong_reset();
    if (GlobalVar.check_logNOTcout()) {
      GlobalVar.slog(1) << std::endl << "ERROR: Whoops, something went wrong: " << e.what() << std::endl;
    }
    return RETURN_ERROR;
  }
  catch(...) { 
    FLXMSG("flx_parse_file_5",1);
    if (reader) { delete reader; reader=NULL; }
    if (!catch_all) throw;
    *(GlobalVar.get_true_cerr()) << std::endl << "ERROR: whoops!" << std::endl; 
    GlobalVar.logLevel_strong_reset();
    if (GlobalVar.check_logNOTcout()) {
      GlobalVar.slog(1) << std::endl << "ERROR: whoops!" << std::endl; 
    }
    return RETURN_ERROR;
  }
  return RETURN_SUCCESS;
}

FlxObjBase* flx_read_file(const string& Expr, const bool is_file, const bool catch_all, bool do_log, int tabWidth)
{
  ReadStream *reader = NULL;
  try {
    if (is_file) {
      GlobalVar.slog(3) << "Fesslix: start reading parameter file '" << Expr << "'" << std::endl;
      reader = new ReadStream( Expr.c_str(), do_log, tabWidth );
    } else {
      GlobalVar.slog(3) << "Fesslix: start reading the command from the command line: " << Expr << std::endl;
      reader = new ReadStream(Expr, do_log, tabWidth);  
    }
    FlxObjBase* ob = FlxEngine_ptr->read_objects(reader);
    delete reader; reader = NULL;
    GlobalVar.logLevel_strong_reset();
    return ob;
  } 
  catch ( FlxException &e) {
    FLXMSG("flx_read_file_3",1);
    if (reader) { delete reader; reader=NULL; }
    if (!catch_all) throw;
    *(GlobalVar.get_true_cerr()) << std::endl << e.what() << std::endl;
    if (GlobalVar.check_logNOTcout()) {
      GlobalVar.slog(1) << std::endl << e.what() << std::endl;
    }
    GlobalVar.logLevel_strong_reset();
    return NULL;
  }
  return NULL;
}

void flx_print_RETURN_SUCCESS() {
  GlobalVar.slogcout(4) << std::endl 
            << "Fesslix: Calculation completed. Fesslix(RETURN_SUCCESS)" << std::endl
            << "-------------------------------------------------------" << std::endl << std::endl;
}

void flx_print_RETURN_ERROR() {
  GlobalVar.slogcout(3) << "Fesslix - Version: "  << flx_get_version() << std::endl;
  GlobalVar.slogcout(4) << "Current Time: ";
  flx_print_current_time(GlobalVar.slogcout(4));
  GlobalVar.slogcout(4) << std::endl << std::endl;
  GlobalVar.slogcout(4) << "Fesslix: Termination due to an error. Fesslix(RETURN_ERROR)" << std::endl
                    << "-----------------------------------------------------------" << std::endl << std::endl;
}




