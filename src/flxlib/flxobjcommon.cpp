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

#include "flxrandom.h"
#include "flxobjcommon.h"
#include "flxobjmtx.h"

#include <fstream>
#include <time.h>

#include "flxostools_files.h"
#include "flxstringfun_fun.h"
#include "flxfunction_fun_calc.h"


void FlxCreateObjReaders_Common::createObjReaders(FlxObjectReadBox* objReadBox) {
  objReadBox->insert("var", new FlxObjReadVar());
  objReadBox->insert("var_", new FlxObjReadVar(false));
  objReadBox->insert("fun", new FlxObjReadFun());
  objReadBox->insert("interpolate", new FlxObjReadInterpolate());
  objReadBox->insert("calc", new FlxObjReadCalc());
  objReadBox->insert("echo", new FlxObjReadEcho());
  objReadBox->insert("run_external", new FlxObjReadRunExternal());
  objReadBox->insert("run_files", new FlxObjReadRunExternal_Files());
  objReadBox->insert("file_filter_cv", new FlxObjReadFileFilterCV());
  objReadBox->insert("file_filter_sofistik", new FlxObjReadFileFilterSOFiSTiK());
  objReadBox->insert("filestream", new FlxObjReadFileStream());
  objReadBox->insert("stringstream", new FlxObjReadStringStream());
  objReadBox->insert("ostream_close", new FlxObjReadOStream_close());
  objReadBox->insert("ifstream", new FlxObjReadInputFileStream());
  objReadBox->insert("ifstream_combine", new FlxObjReadInputFileStreamCombine());
  objReadBox->insert("ivstream", new FlxObjReadInputVectorStream());
  objReadBox->insert("ivstream_append", new FlxObjReadivstream_append());
  objReadBox->insert("ivstream_clear", new FlxObjReadivstream_clear());
  objReadBox->insert("ivstream_reset", new FlxObjReadivstream_clear(true));
  objReadBox->insert("istream_write", new FlxObjReadistream_write());
  objReadBox->insert("dstream", new FlxObjReadDistributorStream());
  objReadBox->insert("default", new FlxObjReadDefault());
  objReadBox->insert("if", new FlxObjReadIf());
  objReadBox->insert("if_no_read", new FlxObjReadIf_no_read());
  objReadBox->insert("while", new FlxObjReadWhile());
  objReadBox->insert("for", new FlxObjReadFor());
  objReadBox->insert("sfor", new FlxObjReadSFor());
  objReadBox->insert("for_each", new FlxObjReadForEach());
  objReadBox->insert("catch_error", new FlxObjReadCatch_Error());
  objReadBox->insert("timer", new FlxObjReadTimer());
  objReadBox->insert("time", new FlxObjReadTime());
  objReadBox->insert("funplot", new FlxObjReadFunPlot());
  objReadBox->insert("funplot_header", new FlxObjReadFunPlot_header());
  objReadBox->insert("end", new FlxObjReadEnd());
  objReadBox->insert("exit", new FlxObjReadExit());
  objReadBox->insert("return", new FlxObjReadReturn());
  objReadBox->insert("continue", new FlxObjReadContinue());
  objReadBox->insert("break", new FlxObjReadBreak());
  objReadBox->insert("throw", new FlxObjReadThrow());
  objReadBox->insert("read", new FlxObjReadReadFile());
  objReadBox->insert("intervalcount", new FlxObjReadIntervalCount());
  objReadBox->insert("filter", new FlxObjReadFilter());
  objReadBox->insert("warranty", new FlxObjReadWarranty());
  objReadBox->insert("isread_vec", new FlxObjReadISread_vec());
  objReadBox->insert("sleep", new FlxObjReadSleep());
  objReadBox->insert("rnd_smp", new FlxObjReadRndSmp());
  objReadBox->insert("rnd_seed", new FlxObjReadRndSeed());
  objReadBox->insert("rnd_track", new FlxObjReadRndTrack());
}

void FlxCreateObjReaders_Common::createFunReaders(FlxData* dataBox)
{
  FlxDataBase::set_data(dataBox);
  dataBox->FunBox.insert("ivstream_size", new FunReadFunIvStream_size() );
  dataBox->FunBox.insert("isread", new FunReadFunISread() );
  dataBox->FunBox.insert("objexec", new FunReadObjExec() );
  dataBox->FunBox.insert("catch_error", new FunReadFunCatchError() );
  dataBox->FunBox.insert("root", new FunReadFunRoot() );
  dataBox->FunBox.insert("lines_in_file", new FunReadFunLinesInFile() );
  dataBox->FunBox.insert("rnd_vec_id", new FunReadFunRndVecID() );
  
  flxString_fun_insert(dataBox->StrFunBox);
  
  FlxBoxBaseR::set_GI(&(dataBox->GaussInt));
}

void FlxObjFun::task()
{
  #if FLX_DEBUG
    if (fun == NULL) {
      std::ostringstream ssV;
      ssV << "ERROR";
      throw FlxException("FlxObjFun::task", ssV.str() );
    }
  #endif
  fun->initialize();
  data->FunBox.insert(cname, fun);
  fun = NULL;
  GlobalVar.slog(4) << "fun: Function '" << cname << "' declared." << std::endl;
}

void FlxObjVar::task()
{
  if (put_on_IL) {
    data->VarBox.insert(cname, fun); 
    fun=NULL;
  } else {
    data->VarBox.insert(cname, new FlxFunction(*fun) ); 
  }
}

FlxObjVar::~FlxObjVar()
{
  if (put_on_IL) {
    if (fun) delete fun;
  } else {
    delete fun;
  }
}


FlxObjFun::~FlxObjFun()
{
  if ( fun ) delete fun;
}

FlxObjBase* FlxObjReadVar::read() {   
  const std::string cname = reader->getWord(true,false);
  isdefined(cname, 'V',false);
  reader->getChar('=',false);
  FlxFunction* fi = NULL;
  FlxFunction* fV = NULL;
  try {
    fi = new FlxFunction(funReader,false);
    fV = data->VarBox.get(cname);
    if ( fi->search_circref(fV) ) {
      std::ostringstream ssV;
      ssV << "Circular reference in '" << cname << "'.";
      throw FlxException("FlxObjReadVar::read_1", ssV.str(), reader->getCurrentPos() );
    }
  } catch ( FlxException &e) {
    FLXMSG("FlxObjReadVar::read_2",1);
    if (fi) delete fi;
    throw;
  }
  read_optionalPara(false);
  FlxObjBase* obj1 = new FlxObjVar( get_doLog(), cname, fi, put_on_IL );
  data->VarBox.declareV(cname);
  return obj1;
}

FlxObjBase* FlxObjReadCalc::read() {
  FlxFunction* f1 = new FlxFunction(funReader,false);
  try {
    read_optionalPara(false);
    return new FlxObjCalc(get_doLog(), f1, get_stream(), get_checkTOL() );
  } catch (FlxException &e) {
    FLXMSG("FlxObjReadCalc::read",1);
    delete f1;
    throw;
  }
}


FlxObjReadEcho::FlxObjReadEcho(): FlxObjReadOutputBase()
{
  // newline
    AllDefParaBox->insert(new FlxOptionalParaBool(true,"echo::newline"));
    ParaBox.insert("newline", "echo::newline" );
}

FlxObjBase* FlxObjReadEcho::read()
{
  FlxString* strV = NULL;
  try {
    strV = new FlxString(true,false);
    read_optionalPara(false);
    return new FlxObjEcho(get_doLog(), strV, get_stream(), get_optPara_bool("newline") );
  } catch (FlxException &e) {
    FLXMSG("FlxObjReadEcho::read",1);
    if (strV) delete strV;
    throw;
  }
}

FlxObjBase* FlxObjReadWarranty::read()
{        
  read_optionalPara(false);
  return new FlxObjWarranty(get_doLog(), get_stream() );
}


FlxObjReadFileFilterCV::FlxObjReadFileFilterCV()
{
  // s_init
    AllDefParaBox->insert(new FlxOptionalParaText("@{","file_filter_cv::s_init"));
    ParaBox.insert("s_init", "file_filter_cv::s_init" );
  // s_end
    AllDefParaBox->insert(new FlxOptionalParaText("}","file_filter_cv::s_end"));
    ParaBox.insert("s_end", "file_filter_cv::s_end" );   
  // totalprec
    AllDefParaBox->insert(new FlxOptionalParaBool(true,"file_filter_cv::totalprec"));
    ParaBox.insert("totalprec", "file_filter_cv::totalprec" ); 
}

FlxObjBase* FlxObjReadFileFilterCV::read()
{
  FlxString* fn = NULL;
  try {
    reader->getChar('(',false);
    fn = new FlxString(false,false);
    reader->getChar(')',false);
    read_optionalPara(false);
    return new FlxObjFileFilterCV(get_doLog(),fn,get_stream(),get_optPara_string("s_init"),get_optPara_string("s_end"),get_optPara_bool("totalprec"));
  } catch (FlxException &e) {
    FLXMSG("FlxObjReadFileFilterCV::read",1);
    if (fn) delete fn;
    throw;
  }
}

FlxObjBase* FlxObjReadFileFilterSOFiSTiK::read()
{
  FlxString* fn = NULL;
  FlxObjBase* block = NULL;
  // for defining the materials
    FlxMtxConstFun* mid = NULL;
    FlxFunction* mid_start = NULL;
  // for setting up a correlation matrix
    FlxFunction* rho = NULL;
  try {
    reader->getChar('(',false);
    fn = new FlxString(false,false);        // input file (the system)
    reader->getChar(')',false);
    const std::string syst_stream = reader->getWord(true,false);
    reader->getChar(',',false);
    const std::string mat_stream = reader->getWord(true,false);
    reader->getChar(',',false);
    tdouble& cvar = *(data->ConstantBox.get(reader->getWord(true,false),true));
    reader->getChar(',',false);
    tdouble& cvar2 = *(data->ConstantBox.get(reader->getWord(true,false),true));
    reader->getChar(',',false);
    const std::string mat_string = reader->getText(true,false);
    reader->getChar(',',false);
    mid = new FlxMtxConstFun(true);
    reader->getChar(',',false);
    mid_start = new FlxFunction(funReader);
    block = FlxObjReadCodeBlock::read_block(true,false);
    read_optionalPara(false);
    return new FlxObjFileFilterSOFiSTiK(get_doLog(),fn,syst_stream,mat_stream,cvar,cvar2,mat_string,block,mid,mid_start);
  } catch (FlxException &e) {
    FLXMSG("FlxObjReadFileFilterSOFiSTiK::read",1);
    if (fn) delete fn;
    if (block) delete block;
    if (mid) delete mid;
    if (mid_start) delete mid_start;
    if (rho) delete rho;
    throw;
  }
}

FlxObjReadRunExternal::FlxObjReadRunExternal()
{
  AllDefParaBox->insert(new FlxOptionalParaBool(true,"runext::throw"));
  ParaBox.insert("throw", "runext::throw" );
}

FlxObjBase* FlxObjReadRunExternal::read()
{
  FlxString* fn = new FlxString(false,false);
  try {
    read_optionalPara(false);
    return new FlxObjRunExternal(get_doLog(),fn,get_stream(),get_optPara_bool("throw"));
  } catch (FlxException &e) {
    FLXMSG("FlxObjReadRunExternal::read",1);
    if (fn) delete fn;
    throw;
  }
}

FlxObjBase* FlxObjReadRunExternal_Files::read()
{
  const std::string jid = reader->getWord(true,false);
  FlxString* fn = new FlxString(false,false);
  FlxString* dest = NULL;
  if ( jid=="delete" || jid=="mkdir" ) {
  } else if ( jid=="copy" || jid=="move" ) {
    dest = new FlxString(false,false);
  }
  try {
    read_optionalPara(false);
    return new FlxObjRunExternal_Files(get_doLog(),jid,fn,dest,get_stream());
  } catch (FlxException &e) {
    FLXMSG("FlxObjReadRunExternal_Files::read",1);
    delete fn;
    if (dest) delete dest;
    throw;
  }
}


FlxObjBase* FlxObjReadTimer::read() {
  const std::string task = reader->getWord(true,false);
  const std::string timer = reader->getWord(true,false);
  if ( task == "start" ) {
    read_optionalPara(false);
    return new FlxObjTimerStart(get_doLog(),timer);
  } else if ( task == "stop" ) {
    read_optionalPara(false);
    return new FlxObjTimerStop(get_doLog(),timer);
  } else if ( task == "print" ) {
    read_optionalPara(false);
    return new FlxObjTimerPrint(get_doLog(),timer, get_stream());
  } else if ( task == "define" ) {
    read_optionalPara(false);
    return new FlxObjTimerDefine(get_doLog(),timer);
  } else if ( task == "delete" ) {
    read_optionalPara(false);
    return new FlxObjTimerDelete(get_doLog(),timer);
  } else {
    std::ostringstream ssV;
    ssV << "Unknown action '" << task << "'.";
    throw FlxException_NeglectInInteractive("FlxObjReadTimer::read_1", ssV.str(), reader->getCurrentPos() );
  }
}

FlxObjReadFunPlot::FlxObjReadFunPlot(): FlxObjReadOutputBase()
{
  // sep_str
    AllDefParaBox->insert(new FlxOptionalParaText("","flxoutput::sep_str"));
    ParaBox.insert("sep_str", "flxoutput::sep_str" );
  // init_str
    AllDefParaBox->insert(new FlxOptionalParaText("","flxoutput::init_str"));
    ParaBox.insert("init_str", "flxoutput::init_str" );
  // end_str
    AllDefParaBox->insert(new FlxOptionalParaText("","flxoutput::end_str"));
    ParaBox.insert("end_str", "flxoutput::end_str" );
  
  AllDefParaBox->insert(new FlxOptionalParaBool(false,"flxoutput::binary"));
  ParaBox.insert("binary", "flxoutput::binary" );
  
  AllDefParaBox->insert(new FlxOptionalParaBool(false,"flxoutput::binaryfloat"));
  ParaBox.insert("binaryfloat", "flxoutput::binaryfloat" );
}

FlxObjBase* FlxObjReadFunPlot::read()
{
  std::vector<FlxFunction*> V;
  std::vector<FlxMtxConstFun*> M;
  std::vector<FlxString*> S;
  std::vector<int> B;
  bool b=false;
  try {
    do {
      if (b) reader->getChar(',',false);
      else b=true;
      if (reader->whatIsNextChar()=='{') {
        reader->getChar('{',false);
        FlxMtxConstFun* f = new FlxMtxConstFun(true);
        M.push_back(f);
        B.push_back(2);
        reader->getChar('}',false);
      } else if (reader->whatIsNextChar()=='[') {
        reader->getChar('[',false);
        FlxString* f = new FlxString(false,false);
        S.push_back(f);
        B.push_back(3);
        reader->getChar(']',false);
      } else {
        FlxFunction* f = new FlxFunction(funReader,false);
        V.push_back(f);
        B.push_back(1);
      }
    } while (reader->whatIsNextChar() == ',');
    read_optionalPara(false);
    return new FlxObjFunPlot(get_doLog(),B,V,M,S,get_stream(),get_checkTOL(),get_prec(),get_fixW(),get_boost_str(),get_optPara_string("sep_str"),get_optPara_string("init_str"),get_optPara_string("end_str"),get_optPara_bool("binary"),get_optPara_bool("binaryfloat"));
  } catch ( FlxException &e) {
    FLXMSG("FlxObjReadFunPlot::read",1);
    for (size_t i = 0; i < V.size(); ++i) {
      delete V[i];
    }
    for (size_t i = 0; i < M.size(); ++i) {
      delete M[i];
    }
    for (size_t i = 0; i < S.size(); ++i) {
      delete S[i];
    }
    throw;
  }
}

FlxObjReadFunPlot_header::FlxObjReadFunPlot_header(): FlxObjReadOutputBase()
{
  // only_once
    AllDefParaBox->insert(new FlxOptionalParaBool(false,"funplot_header::only_once"));
    ParaBox.insert("only_once", "funplot_header::only_once" );
}

FlxObjBase* FlxObjReadFunPlot_header::read()
{
  std::vector<std::string> hdr;
  bool b=false;
  do {
    if (b) reader->getChar(',',false);
    else b=true;
    const std::string hestr = reader->getText(false,false);
    hdr.push_back(hestr);
  } while (reader->whatIsNextChar() == ',');
  read_optionalPara(false);
  return new FlxObjFunPlot_header(get_doLog(),hdr,get_stream(),get_prec(),get_fixW(),get_optPara_bool("only_once"));
}


FlxObjBase* FlxObjReadIntervalCount::read()
{
  FlxObjBase* Ivblock = NULL;
  reader->getChar('(',false);
  FlxFunction* Inumb = new FlxFunction(funReader,false);
  try {
    reader->getChar(')',false);
    Ivblock = FlxObjReadCodeBlock::read_block(true,false);
    read_optionalPara(false);
    return new FlxObjIntervalCount(get_doLog(),Inumb, Ivblock);
  } catch (FlxException &e) {
    FLXMSG("FlxObjReadIntervalCount::read",1);
    delete Inumb;
    if (Ivblock!=NULL) delete Ivblock;
    throw;
  }
}

FlxObjBase* FlxObjReadFilter::read()
{
  reader->getChar('(',false);
  
  if (reader->getNextType() != ReadStream::STRING) {
    std::ostringstream ssV_2;
    ssV_2 << "Name of the 'const' variable to use expected.";
    throw FlxException_NeglectInInteractive("FunReadMtxFunFilter::read_1", ssV_2.str(), reader->getCurrentPos());
  }
  tdouble* cv = data->ConstantBox.get(reader->getWord(true),true);
  reader->getChar(';',false);
  FlxMtxConstFun* seqMtx = new FlxMtxConstFun(true);
  FlxCodeBlock* block = NULL;
  try {
    reader->getChar(')',false);
    block = FlxObjReadCodeBlock::read_block(true,false);
    block->activate_continue_catch();
    read_optionalPara(false);
    return new FlxObjFilter(get_doLog(),cv, seqMtx, block);
  } catch (FlxException &e) {
    FLXMSG("FlxObjReadFilter::read",1);
    delete seqMtx;
    if (block) delete block;
    throw;
  }
}

void FlxObjTimerPrint::task()
{
  if ( data->TimerBox.get(tname)->IsRunning() ) {
    std::ostringstream ssV;
    ssV << "Timer '" << tname << "' is running.";
    throw FlxException("FlxObjTimerPrint::task_1", ssV.str(), "To output information out of a timer you have to stop it first.");
  }
  const tdouble ct = data->TimerBox.get(tname)->get_time();
  std::string tstr = GlobalVar.Double2String(ct);
  sout() << "timer(" << tname << ") = " << tstr << "sec. " << std::endl; //[" << ct << "ticks]" << endl;
  GlobalVar.slog(4) << "timer: timer '" << tname << "' has a value of t=" << tstr << "s." << std::endl;
  *(data->ConstantBox.get("ans",true)) = ct;
}

void FlxObjWhile::task()
{
  try {
    if ( maxCycles <= 0 ) {
        while ( funWhile->calc() > ZERO ) { InternListLoop->exec(); }
      } else {
        unsigned int curCycles = 0;
        while ( funWhile->calc() > ZERO && maxCycles > curCycles) { InternListLoop->exec(); curCycles++; } 
        if ( funWhile->calc() > ZERO ) {
          std::ostringstream ssV;
          ssV << "While-Loop: maximum number of loop-passes exceeded (" << maxCycles << ").";
          throw FlxException("FlxObjWhile::task_1", ssV.str(), "This might be an infinite loop ..." );
        }
      }
  } catch (FlxBreakE &e) { }
}

void FlxObjFor::task()
{
  const tdouble i_old = *theConst;
  ConstDef->exec();
  try {
    if ( maxCycles <= 0 ) {
      while ( funCond->calc() > ZERO ) {
        InternListLoop->exec();
        *theConst = funConst->calc();
      }
    } else {
      unsigned int curCycles = 0;
      while ( funCond->calc() > ZERO && maxCycles > curCycles ) {
        InternListLoop->exec();
        *theConst = funConst->calc();
        curCycles++;
      }
      if ( funCond->calc() > ZERO ) {
        std::ostringstream ssV;
        ssV << "For-Loop: maximum number of loop-passes exceeded (" << maxCycles << ").";
        throw FlxException("FlxObjFor::task_1", ssV.str(), "This might be an infinite loop ..." );
      }
    }
  } catch (FlxBreakE &e) { }
  *theConst = i_old;
}

void FlxObjSFor::task()
{
  const size_t n = funTo->cast2tulongW0(false);
  const tdouble i_old = *theConst;
  try {
    for (size_t i = (start0?0:1); i <= n; ++i) {
      *theConst = i;
      InternListLoop->exec();
    }
  } catch (FlxBreakE &e) { }
  *theConst = i_old;
}

void FlxObjForEach::task()
{
  const std::string s_old = loop_var;
  const size_t sep_len = sep_char.size();
  size_t cpos = 0;
  size_t npos;
  const std::string iter = iterstr->eval(false);
  try {
    do {
      npos = iter.find(sep_char,cpos);
      loop_var = iter.substr(cpos, (npos==std::string::npos)?npos:(npos-cpos) );
      cpos = npos + sep_len;
      trim(loop_var);
      InternListLoop->exec();
    } while (npos!=std::string::npos);
  } catch (FlxBreakE &e) { }
  loop_var = s_old;
}

void FlxObjCatch_Error::task()
{
  bool catched = false;
  if (errSerious) {
    try {
      ill->exec();
    } catch (FlxException& e) {
      FLXMSG("FlxObjCatch_Error::task_1",1);
      if (!NOTdolog) GlobalVar.slogcout(1) << std::endl << e.what() << std::endl;
      catched = true;
    }
  } else {
    try {
      ill->exec();
    } catch (FlxException_NeglectInInteractive& e) {
      FLXMSG("FlxObjCatch_Error::task_2",1);
      if (!NOTdolog) GlobalVar.slogcout(1) << std::endl << e.what() << std::endl;
      catched = true;
    }
  }
  if (catched) {
    if (catchblock) catchblock->exec();
  }
}

void FlxObjTime::task()
{
  // initial difinitions
    FlxTimer t;
    time_t now1;
    time_t now2;
  // perform the measurement  
    time(&now1);
    t.start();
    InternListTime->exec();
    t.stop();
    time(&now2);
    const tdouble mt = t.get_time();
    
  // output the results
    sout() << "time = " << GlobalVar.Double2String(mt) << "sec. " << std::endl;
    if (GlobalVar.check_logNOTcout()) {
      GlobalVar.slog(4) << "time: measured time: " << GlobalVar.Double2String(mt) << "sec. (ticks: " << t.get_ticks() << ")" << std::endl;
    }
    if (rt>=ZERO) {
      sout() << "time to read = " << GlobalVar.Double2String(rt) << "sec. " << std::endl;
      if (GlobalVar.check_logNOTcout()) {
        GlobalVar.slog(4) << "time to read: measured time: " << GlobalVar.Double2String(rt) << "sec." << std::endl;
      }
      rt = -ONE;
    }
    const tdouble pds = difftime(now2,now1);
    if (pds>ZERO) {
      sout() << "physical time = " << GlobalVar.Double2String(pds) << "sec. " << std::endl;
      if (GlobalVar.check_logNOTcout()) {
        GlobalVar.slog(4) << "physical time = " << GlobalVar.Double2String(pds) << "sec. " << std::endl;
      }
    }
    *(data->ConstantBox.get("ans",true)) = (spt)?pds:mt;
}

FlxObjTime::~FlxObjTime()
{
  delete InternListTime;
}

void FlxObjFunPlot::check_first(std::ostream& sout_, bool& not_first)
{
  if (not_first) {
    if (sep_str.empty()) {
      sout_ << ((fixW<0)?'\t':' ');
    } else {
      sout_ << sep_str;
    }
  } else {
    not_first = true;
  }
}

void FlxObjFunPlot::task()
{
  std::vector<FlxFunction*>::iterator j = V.begin();
  std::vector<FlxMtxConstFun*>::iterator k = M.begin();
  std::vector<FlxString*>::iterator s = S.begin();
  std::ostream& sout_ = sout();
  
  if (!binary) {
    sout_ << init_str;
  }
  bool not_first = false;
  for (std::vector<int>::iterator i = B.begin(); i != B.end(); ++i) {
    if (*i == 1) {
      if (binary) {
        if (binaryfloat) {
          const float d = (*j)->calc();
          sout_.write((char *)&d,sizeof(float));
        } else {
          const tdouble d = (*j)->calc();
          sout_.write((char *)&d,sizeof(tdouble));
        }
      } else {
        check_first(sout_,not_first);
        write((*j)->calc(),sout_);
      }
      ++j;
    } else if (*i == 2) {
      FlxSMtx* mb = data->ConstMtxBox.get((*k)->eval(),true);
      for (tuint l = 0; l < mb->get_nrows(); ++l) {
        for (tuint m = 0; m < mb->get_ncols(); ++m) {
          if (binary) {
            if (binaryfloat) {
              const float d = mb->operator()(l,m);
              sout_.write((char *)&d,sizeof(float));
            } else {
              const tdouble d = mb->operator()(l,m);
              sout_.write((char *)&d,sizeof(tdouble));
            }
          } else {
            check_first(sout_,not_first);
            write(mb->operator()(l,m),sout_);
          }
        }
      }
      ++k;
    } else if (*i == 3) {
      if (binary) {
        throw FlxException("FlxObjFunPlot::task_01","The output of string expressions is not available, if parameter 'binary' is set to 'true'.");
      } else {
        check_first(sout_,not_first);
        (*s)->eval(sout_);
      }
      ++s;
    } else {
      throw FlxException_Crude("FlxObjFunPlot::task_02");
    }
  }
  if (!binary) {
    if (end_str.empty()) {
      sout_ << std::endl;
    } else {
      sout_ << end_str;
    }
  }
}

FlxObjFunPlot::~FlxObjFunPlot() {
  for (std::vector<FlxFunction*>::iterator i = V.begin(); i != V.end(); ++i) {
    delete (*i);
  }
  for (std::vector<FlxMtxConstFun*>::iterator i = M.begin(); i != M.end(); ++i) {
    delete (*i);
  }
  for (std::vector<FlxString*>::iterator i = S.begin(); i != S.end(); ++i) {
    delete (*i);
  }
}

void FlxObjFunPlot_header::task()
{
  if (only_once) {
    if (executed) return;
    executed = true;
  }
  const tuint N = hdr.size();
  for (tuint i=0;i<N;++i) {
    write_entry(hdr[i],sout(),prec,fixW,(i==0));
  }
  sout() << std::endl;
}

void FlxObjFunPlot_header::write_entry(std::string hestr, std::ostream& os, const int prec, const int fixW, const bool start)
{
  if (start) os << '#';
  int Wi = GlobalVar.D2S_get_fixW(prec,fixW);
  if (Wi<0) {
    os << hestr << '\t';
    return;
  }
  if (start && Wi>0) Wi -= 1;        // at this point, W should always be >0!
  const tuint W(Wi);
  if (hestr.length()<=W) {
    if (hestr.length()+2<=W && start==false) {
      hestr.insert(0,1,' ');
    }
    if (hestr.length()<W) {
      hestr.insert(hestr.length(),W-hestr.length(),' ');
    }
    os << hestr;
    os << ' ';
  } else {
    os << hestr.substr(0,W);
    os << '.';
  }
}


FlxObjCalc::FlxObjCalc( const bool dolog, FlxFunction* funV, const std::string& ostreamV, const bool checkTOL )
: FlxObjOutputBase(dolog,ostreamV,false,checkTOL), fun(funV), ansptr(data->ConstantBox.get("ans",true))
{

}

void FlxObjCalc::task()
{
  sout() << fun->write() << " = ";
  *ansptr=fun->calc(); 
  sout() << GlobalVar.Double2String(*ansptr,checkTOL) << std::endl;
}

void FlxObjWarranty::task()
{
  sout() << std::endl;
  sout() << "Fesslix  * Copyright (C) 2010-2022 Wolfgang Betz" << std::endl << std::endl;
  sout() << "This program is free software; you can redistribute it and/or modify" << std::endl;
  sout() << "it under the terms of the GNU General Public License as published by" << std::endl;
  sout() << "the Free Software Foundation; either version 3 of the License, or" << std::endl;
  sout() << "(at your option) any later version." << std::endl << std::endl;
  sout() << "This program is distributed in the hope that it will be useful," << std::endl;
  sout() << "but WITHOUT ANY WARRANTY; without even the implied warranty of" << std::endl;
  sout() << "MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the" << std::endl;
  sout() << "GNU General Public License for more details." << std::endl  << std::endl;
  sout() << "You should have received a copy of the GNU General Public License" << std::endl;
  sout() << "along with this program.  If not, see <http://www.gnu.org/licenses/>." << std::endl  << std::endl; 
}

void FlxObjRunExternal::task()
{
  const std::string cmd_str = ext_cmd->eval(false);
  const tuint res = system(cmd_str.c_str());
  if (bthrow && res>0) {
    std::ostringstream ssV;
    ssV << "The command \"" << cmd_str << "\" was not executed successfully. "
    << "The returned error-code is " << res << ".";
    throw FlxException_NeglectInInteractive("FlxObjRunExternal::task_1", "'run_external' caused an error" , ssV.str());
  }
  if (!NOTdolog) {
    GlobalVar.slog(4) << "run_external: \"" << cmd_str << "\" returned " << res << std::endl;
  }
}

FlxObjRunExternal::~FlxObjRunExternal()
{
  delete ext_cmd;
}

FlxObjRunExternal_Files::~FlxObjRunExternal_Files()
{
  delete filename;
  if (destination) delete destination;
}

void FlxObjRunExternal_Files::task()
{
  const std::string strfn = filename->eval(false);
  const std::string strdest = (destination)?(destination->eval(false)):"";
  if (jid=="delete") {
    const tuint i = delDir(strfn);
    if (i>0) {
      GlobalVar.slogcout(4) << "run_files: deleted '" << strfn << "'" << std::endl;
    } else {
      GlobalVar.slogcout(4) << "run_files: file '" << strfn << "' not deleted, because it does not exist." << std::endl;
    }
  } else if (jid=="mkdir") {
    createDir(strfn);
  } else if (jid=="copy") {
    copyFile(strfn,strdest,true);
  } else if (jid=="move") {
    moveFile(strfn,strdest);
  } else {
    std::ostringstream ssV;
    ssV << "ID '" << strfn << "' not recognized.";
    throw FlxException_NeglectInInteractive("FlxObjRunExternal_Files::task_1", ssV.str());
  }
}

void FlxObjIntervalCount::task()
{
  if (funNp != NULL) {
    Np = funNp->cast2tuint(false);
    delete funNp;
    funNp = NULL;
  }
  countI++;
  if (countI >= Np) {
    InternList->exec();
    countI = 0;
  }
}

FlxObjIntervalCount::~FlxObjIntervalCount()
{
  if (funNp != NULL) delete funNp;
  delete InternList;
}

FlxObjFilter::~FlxObjFilter()
{
  delete seqMtx;
  delete block;
}

void FlxObjFilter::task()
{
  const tdouble cold = *cv;
  FlxSMtx* mtx = data->ConstMtxBox.get(seqMtx->eval(),true);
  const tuint nr = mtx->get_nrows();
  const tuint nc = mtx->get_ncols();
  try {
    for (tuint j = 0; j < nc; ++j) {
      for (tuint i = 0; i < nr; ++i) {
        *cv = mtx->operator()(i,j);
        block->exec();
      }
    }
  } catch (FlxBreakE &e) { }
  *cv = cold;
}


void FlxObjFileFilterCV::task() 
{
  const tuint LENGTH = 256;
  // open the file
    const std::string fn = filename->eval(false);
    std::ifstream inpf(fn.c_str());
    if ( ! inpf.is_open() ) {
      std::ostringstream ssV;
      ssV << "File (" << fn << ") could not be opened.";
      throw FlxException("FlxObjFileFilterSOFiSTiK::task_1", ssV.str() );
    }
  // read line for line
    char line[LENGTH];
    std::string str; 
      str.reserve(LENGTH);
    std::ostream& sot = sout();
    while(!inpf.eof()) {
      inpf.getline(line,LENGTH);
      str = line;
      parse_str(str,sot);
    }
}

void FlxObjFileFilterCV::parse_str(const std::string& str, std::ostream& sot)
{
  const tuint LENGTH = 256;
  std::size_t cp=0, c;
  std::string str_num;
  while ((c=str.find(s_init,cp))!=std::string::npos) {
    const size_t c1 = c + s_init.length();
    const size_t c2 = str.find(s_end,c+1);
    sot << str.substr(cp,c-cp);
    if (str[c1]=='!') {
      ReadStream rs(str.substr(c1+1,c2-c1-1), false);  
      data->ReadManager.push(&rs);
      try {
        const tdouble d = data->rbrv_box.get_entry(rs.getWord(true,true,false,true),true)->get_value();
        if (tprec) {
        sot << GlobalVar.D2S_totalPrec(d);
        } else {
        sot << GlobalVar.Double2String(d);
        }
        data->ReadManager.pop();
      } catch ( FlxException &e ) {
        FLXMSG("FlxObjFileFilterCV::parse_str_1",1);
        data->ReadManager.pop();
        throw;
      }
    } else if (str[c1]=='#') {
      const std::string includeFN = str.substr(c1+1,c2-c1-1);
      std::ifstream incl(includeFN.c_str());
      if ( ! incl.is_open() ) {
        std::ostringstream ssV;
        ssV << "File (" << includeFN << ") could not be opened.";
        throw FlxException("FlxObjFileFilterCV::parse_str_2", ssV.str() );
      }
      char line2[LENGTH];
      while(!incl.eof()) {
        incl.getline(line2,LENGTH);
        sot << line2;
        if (!incl.eof()) sot << std::endl;
      }
    } else if (str[c1]=='$') {
      std::string str_cv = str.substr(c1+1,c2-(c1+1));
      std::transform(str_cv.begin(), str_cv.end(), str_cv.begin(), (int(*)(int)) std::tolower);
      FlxSMtx* t_mtx = data->ConstMtxBox.get(str_cv,true);
      SMtxBase_write_fullPrec(sot,*t_mtx);
    } else if (str[c1]=='?') {
      std::string str_cv = str.substr(c1+1,c2-(c1+1));
      std::transform(str_cv.begin(), str_cv.end(), str_cv.begin(), (int(*)(int)) std::tolower);
      sot << data->strConstBox.get(str_cv);
    } else {
      size_t c1m = c1;
      bool asInt=false;
      if (str[c1]=='%') {
        ++c1m;
        asInt = true;
      }
      std::string str_cv = str.substr(c1m,c2-c1m);
      std::transform(str_cv.begin(), str_cv.end(), str_cv.begin(), (int(*)(int)) std::tolower);
      const tdouble* dp = data->ConstantBox.get(str_cv);
      tdouble d;
      if (dp) {
        d = (*dp);
      } else {
        FlxFunction* fp = data->VarBox.get(str_cv);
        if (fp) {
        d = fp->calc();
        } else {
        std::ostringstream ssV;
        ssV << "'" << str_cv << "' is not defined.";
        throw FlxException("FlxObjFileFilterCV::parse_str_3", ssV.str() );
        }
      }
      if (asInt) {
        sot << tulong_from(d,"Number");
      } else {
        if (tprec) {
        sot << GlobalVar.D2S_totalPrec(d);
        } else {
        sot << GlobalVar.Double2String(d);
        }
      }
    }
    cp = c2+1;
  }
  sot << str.substr(cp) << std::endl;
}

FlxObjFileFilterSOFiSTiK::FlxObjFileFilterSOFiSTiK(bool dolog, FlxString* filename, const std::string& syst_stream, const std::string& mat_stream, tdouble& cvar, tdouble& cvar2, const std::string& mat_string, FlxObjBase* block, FlxMtxConstFun* midF, FlxFunction* mid_startF)
: FlxObjBase(dolog), filename(filename), syst_stream(syst_stream), 
  mat_stream(mat_stream), cvar(cvar), cvar2(cvar2), mat_string(mat_string), block(block),
  fcv(new FlxObjFileFilterCV(false,NULL,mat_stream,"@{","}",true)), midF(midF), mid_startF(mid_startF)
{
  
}


void FlxObjFileFilterSOFiSTiK::task() 
{
  const tuint LENGTH = 256;
  const iVec midV = data->ConstMtxBox.get_Vec_cast2tuint(midF->eval(),false);
  tulong mid_start = mid_startF->cast2tulongW0();
  // open the file
    const std::string fn = filename->eval(false);
    std::ifstream inpf(fn.c_str());
    if ( ! inpf.is_open() ) {
      std::ostringstream ssV;
      ssV << "File (" << fn << ") could not be opened.";
      throw FlxException("FlxObjFileFilterCV::task_1", ssV.str() );
    }
  // read line for line
    char line[LENGTH];
    std::string str; 
      str.reserve(LENGTH);
    std::ostream& sstrm = *(data->OstreamBox.get(syst_stream));
    std::ostream& mstrm = *(data->OstreamBox.get(mat_stream));
    while(!inpf.eof()) {
      inpf.getline(line,LENGTH);
      str = line;
      if (str.substr(0,4)=="QUAD") {  // check for key-word (QUAD)
      // get index
        tulong id;
        {
          ReadStream rs(str.substr(4), false);  
          data->ReadManager.push(&rs);
          try {
            id = rs.get_UInt<tulong>(true);
            data->ReadManager.pop();
          } catch ( FlxException &e ) {
            FLXMSG("FlxObjFileFilterCV::task_2",1);
            data->ReadManager.pop();
            throw;
          }
        }
      // detect material-id (mid)
        const size_t c = str.find("MNR",0);
        const size_t c2 = str.find("MBW", 0);
        if (c==std::string::npos || c2==std::string::npos) {
          throw FlxException_Crude("FlxObjFileFilterCV::task_3");
        }
        tulong idm;
        {
          ReadStream rs(str.substr(c+3,c2-c-3), false);  
          data->ReadManager.push(&rs);
          try {
            idm = rs.get_UInt<tulong>(true);
            data->ReadManager.pop();
          } catch ( FlxException &e ) {
            FLXMSG("FlxObjFileFilterCV::task_4",1);
            data->ReadManager.pop();
            throw;
          }
        }
        bool b_mid=false;
        for (tuint i=0;i<midV.size();++i) {
          if (midV[i]==idm) {
            b_mid = true;
            break;
          }
        }
        if (b_mid)  {
          // publish index
            cvar = id;
            cvar2 = mid_start;
          // rewrite the line
            sstrm << str.substr(0,c+3) << "  " << mid_start <<  "  " << str.substr(c2) << std::endl;
          // parse the material
            block->exec();
            fcv->parse_str(mat_string,mstrm);
            ++mid_start;
        } else {
          sstrm << str << std::endl;
        }
      } else {
        sstrm << str << std::endl;
      }
    }
}

void FlxObjFileStream::task() {
  std::ofstream *theStream = NULL;
  const std::string sn = streamname->eval_word(true);
  const std::string fn = filename->eval(false);
  data->OstreamBox.close(sn,false);        // close existing file-stream first
  if ( trunc ) {
    theStream = new std::ofstream(fn.c_str(), std::ios_base::out | std::ios_base::binary );
  } else {
    theStream = new std::ofstream(fn.c_str(), std::ios_base::out | std::ios_base::binary | std::ios_base::app ); 
  }
  try {
    if ( theStream == NULL || ! theStream->is_open() ) {
      std::ostringstream ssV;
      ssV << "File (" << fn << ") could not be opened.";
      throw FlxException("FlxObjFileStream::task_1", ssV.str() );
    }
    data->OstreamBox.insert(sn, theStream);
  } catch (FlxException &e) {
    FLXMSG("FlxObjFileStream::task_2",1);
    if (theStream) delete theStream;
    throw;
  }
  GlobalVar.slog(4) << "filestream: stream '" << sn << "' is directed into file '" << fn << "'." << std::endl;
}

void FlxObjStringStream::task() {
  const std::string sn = streamname->eval_word(true);
  data->OstreamBox.close(sn,false);        // close existing file-stream first
  std::ostringstream *theStream = new std::ostringstream;
  try {
    data->OstreamBox.insert(sn, theStream);
  } catch (FlxException &e) {
    FLXMSG("FlxObjStringStream::task_2",1);
    if (theStream) delete theStream;
    throw;
  }
  GlobalVar.slog(4) << "stringstream: stream '" << sn << "' created." << std::endl;
}

void FlxObjOStream_close::task() {
  const std::string sn = streamname->eval_word(true);
  try {
    data->OstreamBox.close(sn);
    GlobalVar.slog(4) << "ostream_close: output-stream '" << sn << "' closed." << std::endl;
  } catch (FlxException &e) {
    FLXMSG("FlxObjOStream_close::task_1",1);
    std::ostringstream ssV;
    ssV << "Output-stream '" << sn << "' does not exist.";
    GlobalVar.alert.alert("FlxObjOStream_close::task_2",ssV.str());
  }
}

FlxObjInputFileStream::~FlxObjInputFileStream()
{
  delete streamname; 
  if (filename) delete filename; 
  delete bs; 
  delete CnumbF; 
  delete CvecF;
}

const tuint FlxObjInputFileStream::get_Columns(std::vector< tuint >& Cvec, const bool errSerious)
{
  // get the number of columns
    tuint Cnumb = CnumbF->cast2tuint(errSerious);
  // get the columns to read
    std::string strT = CvecF->eval(true);
    if (strT!="") {
      ReadStream rs(strT);
      while (rs.getNextType()!=ReadStream::ENDOFFILE) {
        Cvec.push_back(rs.get_UInt<tuint>(false));
        if (rs.getChar()!=',') break;
      }
    }
  // check the columns to read
    if (Cvec.empty() && Cnumb==1) {
      Cvec.push_back(1);
    } else if (Cvec.empty() && Cnumb > 1) {
      std::ostringstream ssV;
      ssV << "Error reading parameter 'pcol'";
      FlxError(errSerious,"FlxObjInputFileStream::get_Columns_1", ssV.str());
    }
    tuint istart=0;
    for (size_t i=0;i<Cvec.size();++i) {
      if (Cvec[i]<=istart && Cvec[i]>=Cnumb) {
        std::ostringstream ssV;
        ssV << "Error reading parameter 'pcol'";
        FlxError(errSerious,"FlxObjInputFileStream::get_Columns_2", ssV.str());
      } else {
        istart = Cvec[i];
      }
    }
    if (Cvec.size() == Cnumb) {
      Cvec.clear();
      Cvec.push_back(1);
      Cnumb = 1;
    }
  return Cnumb;
}

void FlxObjInputFileStream::task()
{
  std::vector<tuint> Cvec;
  const tuint Cnumb = get_Columns(Cvec,false);  
  const std::string sn = streamname->eval_word(true);
  const std::string fns = filename->eval(false);
  ReadStream* rs=NULL;
  FlxIstream_file* istrm=NULL;
  try {
    if (binary) {
      FlxIstream_file_binary* istrmf = new FlxIstream_file_binary(sn,fns,erreof,bs->cast2tuint(false),Cnumb,Cvec,binaryfloat);
      istrm = istrmf;
      tdouble* ifstream_binary_size = data->ConstantBox.get("ifstream_binary_size",true);
      *ifstream_binary_size = (tdouble)istrmf->get_N_numbers();
    } else {
      const char* test = fns.c_str();
      rs = new ReadStream(test,false,8,false);
      istrm = new FlxIstream_file(sn,rs,erreof,bs->cast2tuint(false),Cnumb,Cvec);
      rs = NULL;
    } 
    data->IstreamBox.insert(sn, istrm, false);
    GlobalVar.slog(3) << "ifstream: file '" << fns << "' directed to input stream '" << sn << "'." << std::endl;
  } catch (FlxException &e) {
    FLXMSG("FlxObjInputFileStream::task",1);
    if (rs) delete rs;
    if (istrm) delete istrm;
    throw;
  }
}

FlxObjInputFileStreamCombine::FlxObjInputFileStreamCombine(const bool dolog, FlxString* streamname, std::vector< FlxString* >& filename_vec, std::vector< FlxFunction* >& weightfun_vec, FlxFunction* bs, FlxFunction* CnumbF, FlxString* CvecF, const bool erreofV)
: FlxObjInputFileStream(dolog, streamname, NULL, bs, CnumbF, CvecF, erreofV), filename_vec(filename_vec), weightfun_vec(weightfun_vec)
{
  
}

FlxObjInputFileStreamCombine::~FlxObjInputFileStreamCombine()
{
  for (size_t i=0;i<filename_vec.size();++i) {
    delete filename_vec[i];
  }
  for (size_t i=0;i<weightfun_vec.size();++i) {
    delete weightfun_vec[i];
  }
}

void FlxObjInputFileStreamCombine::task()
{
  std::vector<tuint> Cvec;
  const tuint Cnumb = get_Columns(Cvec,false);  
  const std::string sn = streamname->eval_word(true);
  // prepare the stream-vector
    if (filename_vec.size()!=weightfun_vec.size()) throw FlxException_Crude("FlxObjInputFileStreamCombine::task_01");
    std::vector<ReadStream*> iReader_vec;
    FlxIstream_file_combine* istrm=NULL;
    flxVec weight_vec(filename_vec.size());
    try {
      for (size_t i=0;i<filename_vec.size();++i) {
        FlxFunction* fp = weightfun_vec[i];
        weight_vec[i] = fp->cast2positive(false);
        const std::string fns = filename_vec[i]->eval(false);
        const char* test = fns.c_str();
        iReader_vec.push_back(new ReadStream(test,false,8,false));
      }
      istrm = new FlxIstream_file_combine(sn,iReader_vec,weight_vec,erreof,bs->cast2tuint(false),Cnumb,Cvec);
      iReader_vec.clear();
      data->IstreamBox.insert(sn, istrm, false);
    } catch (FlxException &e) {
      FLXMSG("FlxObjInputFileStreamCombine::task_02",1);
      for (size_t i=0;i<iReader_vec.size();++i) delete iReader_vec[i];
      if (istrm) delete istrm;
      throw;
    } 
}

void FlxObjInputVectorStream::task()
{
  FlxIstream *istrm_input = NULL;
  std::string istrm_input_str = "";
  const std::string sn = streamname->eval_word(true);
  if (inputStreamName!=NULL) {
    istrm_input_str = inputStreamName->eval_word(true);
  }
  if (istrm_input_str!="") {
    istrm_input = &(data->IstreamBox.get(istrm_input_str));
  }
  FlxIstream_vector* istrm=NULL;
  try {
    istrm = new FlxIstream_vector(sn,istrm_input,erreof,Nreserve->cast2tulong(false));
    data->IstreamBox.insert(sn, istrm);
    GlobalVar.slog(3) << "ivstream: created vector stream '" << sn << "'." << std::endl;
  } catch (FlxException &e) {
    FLXMSG("FlxObjInputVectorStream::task",1);
    if (istrm) delete istrm;
    throw;
  }
}

FlxObjInputVectorStream::~FlxObjInputVectorStream()
{
  delete streamname;
  if (inputStreamName) delete inputStreamName; 
  delete Nreserve;
}

void FlxObjivstream_append::task()
{
  #if FLX_DEBUG
    if (fun && mcf) throw FlxException_Crude("FlxObjivstream_append::task_1");
  #endif
  if (istream==NULL) {
    istream = dynamic_cast<FlxIstream_vector*>(&(data->IstreamBox.get(sn->eval_word(true))));
    if (istream==NULL) {
      std::ostringstream ssV;
      ssV << "Input-stream '" << sn << "' is not a vector-input stream!";
      throw FlxException_NeglectInInteractive("FlxObjivstream_append::task_2", ssV.str() );
    }
  }
  if (fun) {
    istream->appendNumber(fun->calc());
  } else {
    tuint Nr = 0;
    tuint Nc = 0;
    const tdouble * mtxp = data->ConstMtxBox.get_Mtx(mcf->eval(),Nr,Nc);
    for (size_t i=0;i<Nr*Nc;++i) {
      istream->appendNumber(mtxp[i]);
    }
  }
}

void FlxObjivstream_clear::task()
{
  const std::string strmnam = sn->eval_word(true);
  FlxIstream_vector* isv = dynamic_cast<FlxIstream_vector*>(&(data->IstreamBox.get(strmnam)));
  if (isv==NULL) {
    std::ostringstream ssV;
    ssV << "Input-stream '" << sn << "' ins not a vector-input stream!";
    throw FlxException_NeglectInInteractive("FlxObjivstream_clear::task", ssV.str() );
  }
  if (is_reset) {
    isv->reset_stream();
  } else {
    isv->clear();
  }
}

void FlxObjistream_write::task()
{
  FlxIstream &is = data->IstreamBox.get(isname->eval_word(true));
  tdouble d;
  while (is.get_value(d,true)) {
    sout() << GlobalVar.Double2String(d) << std::endl;
  };
}


void FlxObjReadFile::task() {
  // open the file
    std::string fn = filename->eval(false);
    ReadStream *Freader = new ReadStream(fn.c_str(),true,8,false); 
    data->ReadManager.push(Freader);
    GlobalVar.slog(3) << "read: start parsing parameter file '" << fn << "'." << std::endl;
  // read from the file
    FlxObjBase* ob=NULL;
    try {
      while ( Freader->getNextType() != ReadStream::ENDOFFILE ) {
        // print the prompt for output of input to log
          GlobalVar.prelog_append(GlobalVar.get_flxPrompt());
        ob = EvaluateCmd->evaluateCmd();
        ob->exec();
        delete ob; ob=NULL;
      }
    } catch ( FlxExitE &e) {
      FLXMSG("FlxObjReadFile::task_1",1);
      GlobalVar.slog(4) << "exit: 'exit;' was called - file '" << fn << "' will be closed." << std::endl;
      if (ob) delete ob;
    } catch ( FlxEndE &e) {
      FLXMSG("FlxObjReadFile::task_2",1);
      if (ob) delete ob;
      data->ReadManager.pop();
      delete Freader;
      GlobalVar.slog(3) << "read: stop parsing parameter file '" << fn << "'." << std::endl;
      throw;
    } catch (...) {
      if (ob) delete ob;
      data->ReadManager.pop();
      delete Freader;
      throw;
    };
  // close the file
    data->ReadManager.pop();
    delete Freader;
    GlobalVar.slog(3) << "read: stop parsing parameter file '" << fn << "'." << std::endl;
}

void FlxObjDistributorStream::task() {
  const std::string sn = streamname->eval_word(true);
  const std::string s1 = stream1->eval_word(true);
  const std::string s2 = stream2->eval_word(true);
  if (sn==s1 || sn==s2) {
    throw FlxException_NeglectInInteractive("FlxObjDistributorStream::task","Operation not allowed.");
  }
  const ostreamp& stream_1 = data->OstreamBox.get(s1);
  const ostreamp& stream_2 = data->OstreamBox.get(s2);
  std::ostream *theStream = new flxStreamAlloc(stream_1, stream_2);
  data->OstreamBox.insert(sn, theStream);
  GlobalVar.slog(4) << "dstream: stream '" << sn << "' is directed into streams '" << s1 << "' and '" << s2 << "'." << std::endl;
}

void FlxObjDefault::task() {
  para->set_default(value);
}

FlxObjDefault::~FlxObjDefault() {
  para->free_value(value);
}


FlxObjReadFileStream::FlxObjReadFileStream() {
  AllDefParaBox->insert(new FlxOptionalParaBool(true,"filestream::truncate"));
  ParaBox.insert("truncate", "filestream::truncate" );
}

FlxObjBase* FlxObjReadFileStream::read() {
  FlxString* fn = NULL;
  FlxString* sn = new FlxString(false,false);
  try {
    reader->getChar('(',false);
    fn = new FlxString(false,false);
    reader->getChar(')',false);
    read_optionalPara(false);
    const bool bt = *(static_cast<bool*>( ParaBox.get("truncate")->get() ));
    return new FlxObjFileStream(get_doLog(),sn, fn, bt);
  } catch (FlxException &e) {
    FLXMSG("FlxObjReadFileStream::read",1);
    delete sn;
    if (fn) delete fn;
    throw;
  }
}

FlxObjBase* FlxObjReadStringStream::read() {
  FlxString* sn = new FlxString(false,false);
  try {
    read_optionalPara(false);
    return new FlxObjStringStream(get_doLog(),sn);
  } catch (FlxException &e) {
    FLXMSG("FlxObjReadFileStream::read",1);
    delete sn;
    throw;
  }
}

FlxObjBase* FlxObjReadOStream_close::read() {
  FlxString* sn = new FlxString(false,false);
  try {
  read_optionalPara(false);
  return new FlxObjOStream_close(get_doLog(),sn);
  } catch (FlxException &e) {
    delete sn;
    throw;
  }
}

FlxObjReadInputFileStream::FlxObjReadInputFileStream()
{
  AllDefParaBox->insert(new FlxOptionalParaFun(1e3,"ifstream::blocksize"));
  ParaBox.insert("blocksize", "ifstream::blocksize" );
  
  AllDefParaBox->insert(new FlxOptionalParaBool(true,"istream::erreof"));
  ParaBox.insert("erreof", "istream::erreof" );
  
  AllDefParaBox->insert(new FlxOptionalParaFun(1,"ifstream::colnumb"));
  ParaBox.insert("colnumb", "ifstream::colnumb" );
  
  AllDefParaBox->insert(new FlxOptionalParaFlxString("","ifstream::pcol",false));
  ParaBox.insert("pcol", "ifstream::pcol" );
  
  AllDefParaBox->insert(new FlxOptionalParaBool(false,"ifstream::binary"));
  ParaBox.insert("binary", "ifstream::binary" );
  
  AllDefParaBox->insert(new FlxOptionalParaBool(false,"ifstream::binaryfloat"));
  ParaBox.insert("binaryfloat", "ifstream::binaryfloat" );
  
  data->ConstantBox.declareC("ifstream_binary_size");
}

FlxObjBase* FlxObjReadInputFileStream::read()
{
  FlxString* fn = NULL;
  FlxString* sn = new FlxString(false,false);
  try {
    reader->getChar('(',false);
    fn = new FlxString(false,false);
    reader->getChar(')',false);
    read_optionalPara(false);
    return new FlxObjInputFileStream(get_doLog(),sn, fn, get_optPara_FlxFunction("blocksize"), get_optPara_FlxFunction("colnumb"), get_optPara_FlxString("pcol"), get_optPara_bool("erreof"), get_optPara_bool("binary"), get_optPara_bool("binaryfloat"));
  } catch (FlxException &e) {
    FLXMSG("FlxObjReadInputFileStream::read",1);
    delete sn;
    if (fn) delete fn;
    throw;
  }
}

FlxObjBase* FlxObjReadInputFileStreamCombine::read()
{
  FlxString* sn = new FlxString(false,false);
  std::vector<FlxString*> filename_vec;
  std::vector<FlxFunction*> weight_vec;
  try {
    reader->getChar(':');
    bool b = true;
    while (b) {
      filename_vec.push_back(new FlxString(false,false));
      reader->getChar('(',false);
      weight_vec.push_back(new FlxFunction(funReader,false));
      reader->getChar(')',false);
      b = (reader->whatIsNextChar()==',');
      if (b) reader->getChar(',');
    };
    return new FlxObjInputFileStreamCombine(get_doLog(),sn,filename_vec,weight_vec, get_optPara_FlxFunction("blocksize"), get_optPara_FlxFunction("colnumb"), get_optPara_FlxString("pcol"), get_optPara_bool("erreof"));
  } catch (FlxException &e) {
    FLXMSG("FlxObjReadInputFileStreamCombine::read",1);
    delete sn;
    for (size_t i=0;i<filename_vec.size();++i) delete filename_vec[i];
    for (size_t i=0;i<weight_vec.size();++i) delete weight_vec[i];
    throw;
  }
}

FlxObjReadInputVectorStream::FlxObjReadInputVectorStream()
{
  AllDefParaBox->insert(new FlxOptionalParaBool(true,"istream::erreof"));
  ParaBox.insert("erreof", "istream::erreof" );
  
  AllDefParaBox->insert(new FlxOptionalParaFun(1e5,"ivstream::nreserve"));
  ParaBox.insert("nreserve", "ivstream::nreserve" );
}

FlxObjBase* FlxObjReadInputVectorStream::read()
{
  FlxString* fn = NULL;
  FlxString* sn = new FlxString(false,false);
  try {
  reader->getChar('(',false);
  if (reader->whatIsNextChar()!=')') { 
    fn = new FlxString(false,false);
  }
  reader->getChar(')',false);
  read_optionalPara(false);
  return new FlxObjInputVectorStream(get_doLog(),sn, fn, get_optPara_FlxFunction("nreserve"), get_optPara_bool("erreof"));
  } catch (FlxException &e) {
    FLXMSG("FlxObjReadInputVectorStream::read",1);
    if (fn) delete fn;
    throw;
  }
}

FlxObjBase* FlxObjReadivstream_append::read()
{
  FlxString* sn = new FlxString(false,false);
  FlxFunction *fun = NULL;
  FlxMtxConstFun* mcf = NULL;
  try {
    reader->getChar('+');
    reader->getChar('=');
    if (reader->whatIsNextChar()=='{') {
      reader->getChar('{');
      mcf = new FlxMtxConstFun(true);
      reader->getChar('}');
    } else {
      fun = new FlxFunction(funReader,false);
    }
    return new FlxObjivstream_append(get_doLog(),sn,fun,mcf);
  } catch (FlxException &e) {
    FLXMSG("FlxObjReadivstream_append::read",1);
    delete sn;
    if (fun) delete fun;
    if (mcf) delete mcf;
    throw;
  }
}

FlxObjBase* FlxObjReadivstream_clear::read()
{
  FlxString* sn = new FlxString(false,false);
  try {
    read_optionalPara(false);
    return new FlxObjivstream_clear(get_doLog(),sn,is_reset);
  } catch (FlxException& e) {
    delete sn;
    throw;
  }
}

FlxObjBase* FlxObjReadistream_write::read()
{
  FlxString *isname = new FlxString(false,false);
  try {
    read_optionalPara(false);
    return new FlxObjistream_write(get_doLog(),isname,get_stream());
  } catch (FlxException &e) {
    FLXMSG("FlxObjReadistream_write::read",1);
    delete isname;
    throw;
  }
}

FlxObjBase* FlxObjReadReadFile::read() {
  FlxString* fn = new FlxString(false,false);
  read_optionalPara(false);
  return new FlxObjReadFile(get_doLog(),fn);
}

FlxObjBase* FlxObjReadDistributorStream::read() {
  FlxString* stream1=NULL; FlxString* stream2=NULL;
  FlxString* dstream = new FlxString(false,false);
  try {
    reader->getChar('(',false);
    FlxString* stream1 = new FlxString(false,false);
    reader->getChar(',',false);
    FlxString* stream2 = new FlxString(false,false);
    reader->getChar(')',false);
    read_optionalPara(false);
    return new FlxObjDistributorStream(get_doLog(),dstream, stream1, stream2);
  } catch (FlxException &e) {
    delete dstream;
    if (stream1) delete stream1;
    if (stream2) delete stream2;
    throw;
  }
}

void FlxObjFloatingPointConversion::task()
{
  switch (id) {
    case 0:
      GlobalVar.Double2String_setPrec(fun->cast2tuintW0(false));
      break;
    case 1:
      GlobalVar.Double2String_setType(fun->cast2tuintW0(false));
      break;
    case 2:
      GlobalVar.Double2String_setBValU(fun->calc());
      break;
    case 3:
      GlobalVar.Double2String_setBValL(fun->calc());
      break;
    case 4:
      if ( fun->calc() == ZERO ) {
        GlobalVar.Double2String_setDel0(false);
      } else {
        GlobalVar.Double2String_setDel0(true);
      }
      break;
    case 5:
      if ( fun->calc() == ZERO ) {
        GlobalVar.Double2String_setDelP(false);
      } else {
        GlobalVar.Double2String_setDelP(true);
      }
      break;
    case 6:
      GlobalVar.set_TOL(fun->calc());
  }
}

void FlxObjSetVariousDefault::task()
{
  switch (ID) {
    case 1:
      GlobalVar.logLevel_strong_set(fun->cast2int());
      break;
    case 2:
      FlxFunDeg::set_deg(fun->cast2tuintW0());
      break;
    default:
      #if FLX_DEBUG
        throw FlxException_Crude("FlxObjSetVariousDefault::task");
      #endif
      return;
  }
}

FlxObjBase* FlxObjReadDefault::read_special(std::string& key)
{
  // floating point conversion
    if (key == "flxoutput::float::prec") {
      FlxFunction* f1 = new FlxFunction(funReader,false);
      read_optionalPara(false);
      return new FlxObjFloatingPointConversion(get_doLog(),f1,0);
    } else if (key == "flxoutput::float::type") {
      FlxFunction* f1 = new FlxFunction(funReader,false);
      read_optionalPara(false);
      return new FlxObjFloatingPointConversion(get_doLog(),f1,1);
    } else if (key == "flxoutput::float::bvalu") {
      FlxFunction* f1 = new FlxFunction(funReader,false);
      read_optionalPara(false);
      return new FlxObjFloatingPointConversion(get_doLog(),f1,2);
    } else if (key == "flxoutput::float::bvall") {
      FlxFunction* f1 = new FlxFunction(funReader,false);
      read_optionalPara(false);
      return new FlxObjFloatingPointConversion(get_doLog(),f1,3);
    } else if (key == "flxoutput::float::del0") {
      FlxFunction* f1 = new FlxFunction(funReader,false);
      read_optionalPara(false);
      return new FlxObjFloatingPointConversion(get_doLog(),f1,4);
    } else if (key == "flxoutput::float::delp") {
      FlxFunction* f1 = new FlxFunction(funReader,false);
      read_optionalPara(false);
      return new FlxObjFloatingPointConversion(get_doLog(),f1,5);
    } else if (key == "flxoutput::float::tol") {
      FlxFunction* f1 = new FlxFunction(funReader,false);
      read_optionalPara(false);
      return new FlxObjFloatingPointConversion(get_doLog(),f1,6);
    }
  // loglevel
    if (key == "log::level") {
      FlxFunction* f1 = new FlxFunction(funReader);
      read_optionalPara(false);
      return new FlxObjSetVariousDefault(get_doLog(),1,f1);
    }
  // default polynomial degree of FlxFunDeg
    if (key == "flxfundeg::degree") {
      FlxFunction* f1 = new FlxFunction(funReader);
      read_optionalPara(false);
      return new FlxObjSetVariousDefault(get_doLog(),2,f1);
    }
    
    return NULL;
}

FlxObjBase* FlxObjReadDefault::read() {
  // extract the DefSetName
    std::string sn = reader->getWord(true,false);
    while ( reader->whatIsNextChar() == ':' ) {
      reader->getChar(':',false);
      reader->getChar(':',false);
      sn += "::";
      sn += reader->getWord(true,false);
    }
    std::transform(sn.begin(), sn.end(), sn.begin(), (int(*)(int)) std::tolower);
    reader->getChar('=',false);
  // check for special treatment
    FlxObjBase* obj1 = read_special(sn);
    if (obj1) return obj1;
  // get the Parameter-Class
    FlxOptionalParaBase* pb = AllDefParaBox->get(sn);
    if ( pb == NULL ) {
      std::ostringstream ssV;
      ssV << "Unknown parameter '" << sn << "'.";
      throw FlxException_NeglectInInteractive("FlxObjReadDefault::read_1", ssV.str(), reader->getCurrentPos());
    }
    void *vp = pb->read_value(false);
    read_optionalPara(false);
    obj1 = new FlxObjDefault(get_doLog(),vp, pb);
    return obj1;
}

FlxObjIf::~FlxObjIf() {
  delete funIf;
  delete InternListThen;
  delete InternListElse;
}

FlxObjWhile::~FlxObjWhile() {
  delete funWhile;
  delete InternListLoop;
}

FlxObjFor::~FlxObjFor() {
  delete funCond;
  delete funConst;
  delete InternListLoop;
  delete ConstDef;
}

FlxObjBase* FlxObjReadWhile::read() {
  FlxFunction* f1 = NULL;
  FlxCodeBlock* loopblock = NULL;
  try {
    reader->getChar('(',false);
    f1 = new FlxFunction(funReader,false);
    reader->getChar(')',false);
    loopblock = FlxObjReadCodeBlock::read_block(true,false);
    loopblock->activate_continue_catch();
    read_optionalPara(false);
    return new FlxObjWhile(get_doLog(),f1, loopblock, get_maxpasses());
  } catch ( FlxException& e) {
    FLXMSG("FlxObjReadWhile::read",1);
    if (f1) delete f1;
    if (loopblock) delete loopblock;
    throw;
  }
}

FlxObjBase* FlxObjReadFor::read() {
  FlxObjReadConst* objRC = NULL;
  FlxObjConst* obj1 = NULL;
  FlxFunction* funCond = NULL; FlxFunction* funConst = NULL;
  FlxCodeBlock* loopblock = NULL;
  try {
    reader->getChar('(',false);
    const std::string cname = reader->getWord(true,false);
    objRC = new FlxObjReadConst();
    obj1 = dynamic_cast<FlxObjConst*> ( objRC->read(cname) );
    delete objRC; objRC=NULL;
    tdouble* d1 = data->ConstantBox.get(cname);
    reader->getChar(';',false);
    funCond = new FlxFunction(funReader,false);
    reader->getChar(';',false);
    funConst = new FlxFunction(funReader,false);
    reader->getChar(')',false);
    loopblock = FlxObjReadCodeBlock::read_block(true,false);
    loopblock->activate_continue_catch();
    read_optionalPara(false);
    return new FlxObjFor(get_doLog(),d1, obj1, funCond, funConst, loopblock, get_maxpasses());
  } catch (FlxException &e) {
    FLXMSG("FlxObjReadFor::read",1);
    if (objRC) delete objRC;
    if (obj1) delete obj1;
    if (funCond) delete funCond;
    if (funConst) delete funConst;
    if (loopblock) delete loopblock;
    throw;
  }
}

FlxObjBase* FlxObjReadSFor::read() {
  FlxFunction* funTo = NULL;
  FlxCodeBlock* loopblock = NULL;
  try {
    reader->getChar('(',false);
    const std::string cname = reader->getWord(true,false);
    reader->getChar(';',false);
    funTo = new FlxFunction(funReader,false);
    bool start0 = false;
    if (reader->whatIsNextChar()==';') {
      reader->getChar(';',false);
      start0 = reader->getBool(false);
    }
    reader->getChar(')',false);
    tdouble* d1 = data->ConstantBox.get(cname,true);
    loopblock = FlxObjReadCodeBlock::read_block(true,false);
    loopblock->activate_continue_catch();
    return new FlxObjSFor(get_doLog(),d1, funTo, start0, loopblock);
  } catch (FlxException &e) {
    FLXMSG("FlxObjReadSFor::read",1);
    if (funTo) delete funTo;
    if (loopblock) delete loopblock;
    throw;
  }
}

FlxObjReadForEach::FlxObjReadForEach(): FlxObjReadBase()
{
  AllDefParaBox->insert(new FlxOptionalParaBool(true,"for_each::trim"));
    ParaBox.insert("trim", "for_each::trim" );
}

FlxObjBase* FlxObjReadForEach::read()
{
  std::string& loop_var = data->strConstBox.get_ref( reader->getWord(true,false) );
  reader->getWord("in",false);
  reader->getChar('(',false);
  FlxString* iterstr = new FlxString(true,false);
  std::string sep_char = ";";  
  FlxCodeBlock* loopblock = NULL;
  try {
    if (reader->whatIsNextChar()==';') {
      reader->getChar();
      sep_char = reader->getText(true,false);
    }
    reader->getChar(')');
    loopblock = FlxObjReadCodeBlock::read_block(true,false);
    loopblock->activate_continue_catch();
    return new FlxObjForEach(get_doLog(),loop_var, iterstr, sep_char, loopblock, get_optPara_bool("trim"));
  } catch (FlxException &e) {
    FLXMSG("FlxObjReadForEach::read",1);
    delete iterstr;
    if (loopblock) delete loopblock;
    throw;
  }
}

FlxObjReadCatch_Error::FlxObjReadCatch_Error()
{
  AllDefParaBox->insert(new FlxOptionalParaBool(false,"catch_error::errserious"));
    ParaBox.insert("errserious", "catch_error::errserious" );
}

FlxObjBase* FlxObjReadCatch_Error::read() {
  FlxObjBase* loopblock = FlxObjReadCodeBlock::read_block(true,false);
  FlxObjBase* catchblock = NULL;
  try {
    if (reader->whatIsNextString(6,true)=="handle") {
      reader->getWord("handle");
      catchblock = FlxObjReadCodeBlock::read_block(true,false);
    }
    read_optionalPara(false);
  } catch (FlxException& e) {
    delete loopblock;
    if (catchblock) delete catchblock;
    throw;
  }
  return new FlxObjCatch_Error(get_doLog(),loopblock,catchblock,get_optPara_bool("errserious"));
}

FlxObjReadTime::FlxObjReadTime()
{
  // Store CPU-time or physical time?
    AllDefParaBox->insert(new FlxOptionalParaBool(false,"time::store_physical"));
    ParaBox.insert("store_physical", "time::store_physical" );
}

FlxObjBase* FlxObjReadTime::read()
{ 
  FlxTimer t;
  t.start();  
  FlxObjBase* timeblock = FlxObjReadCodeBlock::read_block(false,false);
  read_optionalPara(false);
  t.stop();
  const tdouble rt = t.get_time();
  return new FlxObjTime(get_doLog(),timeblock,get_stream(),rt,get_optPara_bool("store_physical"));
}

FlxObjBase* FlxObjReadIf::read() {
  FlxFunction* f1 = NULL;
  FlxCodeBlock* thenblock = NULL; FlxCodeBlock* elseblock = NULL;
  try {
    reader->getChar('(',false);
    f1 = new FlxFunction(funReader,false);
    reader->getChar(')',false);
    thenblock = FlxObjReadCodeBlock::read_block(false,false);
    if ( reader->getNextType() == ReadStream::STRING ) {
      const std::string strV = reader->getWord(true,false);
      if ( strV != "else" ) {
        std::ostringstream ssV;
        ssV << "Expected 'else' or ';' (not '" << strV << "').";
        throw FlxException_NeglectInInteractive("FlxObjReadIf::read_1", ssV.str(), reader->getCurrentPos());
      }
      elseblock = FlxObjReadCodeBlock::read_block(false,false);
    }
    thenblock->deactivate_return_catch();
    if (elseblock) elseblock->deactivate_return_catch();
    read_optionalPara(false);
    return new FlxObjIf(get_doLog(),f1, thenblock, elseblock);
  } catch ( FlxException &e) {
    FLXMSG("FlxObjReadTime::read",1);
    if (f1) delete f1;
    if (thenblock) delete thenblock;
    if (elseblock) delete elseblock;
    throw;
  }
}

FlxObjBase* FlxObjReadIf_no_read::read() {
  FlxFunction* f1 = NULL;
  FlxCodeBlock* block = NULL;
  try {
    reader->getChar('(',false);
    f1 = new FlxFunction(funReader,false);
    reader->getChar(')',false);
    const bool b = (f1->calc() > ZERO);
    delete f1; f1=NULL;
    // if-block
      if (b) {
        block = FlxObjReadCodeBlock::read_block(false,false);
      } else {
        reader->getChar('{');
        reader->ignore_bracket('}');
      }
    // else-block
      if ( reader->getNextType() == ReadStream::STRING ) {
        const std::string strV = reader->getWord(true,false);
        if ( strV != "else" ) {
          std::ostringstream ssV;
          ssV << "Expected 'else' or ';' (not '" << strV << "').";
          throw FlxException_NeglectInInteractive("FlxObjReadIf::FlxObjReadIf_no_read_1", ssV.str(), reader->getCurrentPos());
        }
        if (b) {
          reader->getChar('{');
          reader->ignore_bracket('}');
        } else {
          block = FlxObjReadCodeBlock::read_block(false,false);
        }
      }
    if (block) {
      block->deactivate_return_catch();
      return block;
    } else {
      return new FlxObjDummy();
    }
  } catch ( FlxException &e) {
    FLXMSG("FlxObjReadTime::read",1);
    if (f1) delete f1;
    if (block) delete block;
    throw;
  }
}

FlxObjBase* FlxObjReadFun::read() {
  const std::string cname = get_name();
  // number of parameters
    reader->getChar('(',false);
    
    FlxFunction* i1f = NULL;
    tuint i1=0;
    if (reader->whatIsNextChar()!=')') {
      try {
        i1f = new FlxFunction(funReader,false);
        i1 = i1f->cast2tuintW0(false);
        delete i1f;
            } catch (FlxException &e) {
        FLXMSG("FlxObjReadFun::read_2",1);
        if (i1f) delete i1f; 
        throw;
      }
    }
    reader->getChar(')',false);
  reader->getChar('=',false);
  
  FunReadPara::set_NumbOfPara(i1);
  FlxFunction *fp = NULL;
  try {
    fp = new FlxFunction(funReader,false);
    read_optionalPara(false);
    FlxObjBase* obj1 = new FlxObjFun( get_doLog(), cname, new FunReadFunUser(cname, fp, i1) );
    FunReadPara::set_NumbOfPara(0);
    data->FunBox.declareF(cname);
    return obj1;
  } catch (FlxException &e) {
    FLXMSG("FlxObjReadFun::read_3",1);
    if (fp) delete fp;
    FunReadPara::set_NumbOfPara(0);
    throw;
  }
}

const std::string FlxObjReadFun::get_name()
{
  const std::string cname = reader->getWord(true,false);
  // check if function is already defined
    if ( data->FunBox.get(cname) != NULL ) {
      std::ostringstream ssV_2;
      ssV_2 << "Function '" << cname << "' was already defined - redefinition of functions is not allowed.";
      throw FlxException_NeglectInInteractive("FlxObjReadFun::read_1", ssV_2.str(), reader->getCurrentPos() );
    }
  // check if a variable with that name exists already
    isdefined(cname, 'F',false);
  return cname;
}

FlxObjBase* FlxObjReadInterpolate::read()
{
  const std::string cname = get_name();
  reader->getChar('(',false);
  FlxMtxConstFun* mtxfun = new FlxMtxConstFun(false);
  reader->getChar(')',false);
  try {
    read_optionalPara(false);
    FlxObjBase* obj1 = new FlxObjFun( get_doLog(), cname, new FunReadFunInterpolate(cname, mtxfun) );
    data->FunBox.declareF(cname);
    return obj1;
  } catch (FlxException &e) {
    FLXMSG("FlxObjReadInterpolate::read",1);
    FunReadPara::set_NumbOfPara(0);
    throw;
  }
}

void FunReadFunInterpolate::initialize()
{
  data.initialize();
}

FunBase* FunReadFunInterpolate::read(bool errSerious)
{
  return new FunInterpolate( read_parameters( 1, errSerious ), fname, &data );
}

const tdouble FunInterpolate::calc()
{
  const tdouble xval = child_1->calc();
  const flxVec& xvec = data->get_xvec();
  const flxVec& yvec = data->get_yvec();
  const size_t Nel = xvec.get_N();
  return flx_interpolate_linear(xval,xvec.get_tmp_vptr_const(),yvec.get_tmp_vptr_const(),Nel);
}

Interpolate_help::~Interpolate_help()
{
  if (mtxfun) delete mtxfun;
  if (xvec) delete xvec;
  if (yvec) delete yvec;
}

void Interpolate_help::initialize()
{
  if (xvec||yvec) throw FlxException_Crude("Interpolate_help::initialize_01");
  const std::string mtxname = mtxfun->eval();
  FlxSMtx* mtx = data->ConstMtxBox.get(mtxname,true);
  // check if matrix has two columns
    if (mtx->get_ncols()!=2) {
      std::ostringstream ssV;
      ssV << "Matrix '" << mtxname << "' must have two columns, but has " << mtx->get_ncols() << " columns.";
      throw FlxException("Interpolate_help::initialize_02", ssV.str() );
    }
  const size_t Nr = mtx->get_nrows();
  xvec = new flxVec(Nr);
  yvec = new flxVec(Nr);
  for (size_t i=0;i<Nr;++i) {
    xvec->operator[](i) = mtx->operator()(i,0);
    yvec->operator[](i) = mtx->operator()(i,1);
  }
  delete mtxfun; mtxfun=NULL;
}

flxVec& Interpolate_help::get_xvec()
{
  if (xvec==NULL) throw FlxException_Crude("Interpolate_help::get_xvec");
  return *xvec;
}

flxVec& Interpolate_help::get_yvec()
{
  if (yvec==NULL) throw FlxException_Crude("Interpolate_help::get_yvec");
  return *yvec;
}

FlxObjBase* FlxObjReadEnd::read()
{
  read_optionalPara(false);
  return new FlxObjEnd(get_doLog());
}

FlxObjBase* FlxObjReadExit::read()
{
  read_optionalPara(false);
  return new FlxObjExit(get_doLog());
}

FlxObjBase* FlxObjReadReturn::read()
{
  read_optionalPara(false);
  return new FlxObjReturn(get_doLog());
}

FlxObjBase* FlxObjReadContinue::read()
{
  read_optionalPara(false);
  return new FlxObjContinue(get_doLog());
}

FlxObjBase* FlxObjReadBreak::read()
{
  read_optionalPara(false);
  return new FlxObjBreak(get_doLog());
}

FlxObjBase* FlxObjReadThrow::read()
{
  read_optionalPara(false);
  return new FlxObjThrow(get_doLog());
}

void FlxObjISread_vec::set_istrm()
{
  #if FLX_DEBUG
  if (strV==NULL) throw FlxException_Crude("FlxObjISread_vec::set_istrm");
  #endif
  strS = strV->eval_word(true);
  istrm = &(data->IstreamBox.get(strS));
  delete strV;
  strV = NULL;
}

void FlxObjISread_vec::task()
{
  if (istrm==NULL) {
    set_istrm();
  }
  tuint dim = (dimF?dimF->cast2tuintW0(false):0);
  const std::string vname = VecConstStr->eval_word(true);
  tdouble* vp = data->ConstMtxBox.get_Vec(vname,dim);
  flxVec vt(vp,dim);
  if (istrm->get_vec(vt,dim)==false) {
    std::ostringstream ssV;
    ssV << "There were not enough numbers in the input stream '" << strS << "'.";
    throw FlxException_NeglectInInteractive("FlxObjISread_vec::task",ssV.str());
  }
}

FlxObjISread_vec::~FlxObjISread_vec()
{
  delete VecConstStr;
  if (dimF) delete dimF;
  if (strV) delete strV;
}

FlxObjBase* FlxObjReadISread_vec::read()
{
  FlxString* VecConstStr = new FlxString(false,false);
  FlxFunction* dimF = NULL;
  FlxString* strm = NULL;
  try {
    reader->getChar('(');
    if (reader->whatIsNextChar()!=')') {
      dimF = new FlxFunction(funReader,false);
    }
    reader->getChar(')');
    reader->getChar('=');
    strm = new FlxString(false,false);
    read_optionalPara(false);
    return new FlxObjISread_vec(get_doLog(),VecConstStr,dimF,strm);
  } catch (FlxException &e) {
    FLXMSG("FlxObjReadISread_vec::read",1);
    delete VecConstStr;
    if (dimF) delete dimF;
    if (strm) delete strm;
    throw;
  }
}

void FlxObjSleep::task()
{
  tuint t = time2wait->cast2tuint(false);
  GlobalVar.slogcout(3) << "Sleep: going to sleep for " << t << " seconds." << std::endl;;
  //flx_put_to_sleep(t);
  throw FlxException_NotImplemented("FlxObjSleep::task");
}

FlxObjBase* FlxObjReadSleep::read()
{
  reader->getChar('(');
  FlxFunction* time2wait = new FlxFunction(funReader,false);
  try {
    reader->getChar(')');
    read_optionalPara(false);
    return new FlxObjSleep(get_doLog(),time2wait);
  } catch (FlxException &e) {
    FLXMSG("FlxObjReadSleep::read",1);
    delete time2wait;
    throw;
  }
}


FunIvStream_size::~FunIvStream_size()
{
  if (strV != NULL) delete strV;
}

void FunIvStream_size::set_istrm()
{
  strS = strV->eval_word(true);
    delete strV;
    strV = NULL;
  istrm = dynamic_cast<FlxIstream_vector*>(&(data->IstreamBox.get(strS)));
  if (istrm==NULL) {
    std::ostringstream ssV;
    ssV << "Input-stream '" << strS << "' ins not a vector-input stream!";
    throw FlxException_NeglectInInteractive("FunIvStream_size::set_istrm", ssV.str() );
  }
}

const tdouble FunIvStream_size::calc()
{
  if (istrm==NULL) {
    set_istrm();
  }
  return tdouble(istrm->get_total_size());
}

const std::string FunIvStream_size::write()
{
  if (istrm==NULL) {
    set_istrm();
  }
  return "ivstream_size(" + strS + ")";
}

const bool FunIvStream_size::search_circref(FlxFunction* fcr)
{
  if (strV != NULL) {
    return strV->search_circref(fcr);
  } else {
    return false;
  }
}

FunBase* FunReadFunIvStream_size::read(bool errSerious)
{
  FlxString* strV = new FlxString(false,errSerious);
  return new FunIvStream_size( strV );
}


FunBase* FunReadFunISread::read(bool errSerious)
{
  FlxString* strV = new FlxString(false,errSerious);
  return new FunISread( strV );
}

FunISread::~FunISread()
{
  if (strV) delete strV;
}

const bool FunISread::search_circref(FlxFunction* fcr)
{
  if (strV) {
    return strV->search_circref(fcr);
  } else {
    return false;
  }
}

void FunISread::set_istrm()
{
  strS = strV->eval_word(true);
  istrm = &(data->IstreamBox.get(strS));
  delete strV;
  strV = NULL;
}

const std::string FunISread::write()
{
  if (istrm==NULL) {
    set_istrm();
  }
  return "isread(" + strS + ")";
}

const tdouble FunISread::calc()
{
  if (istrm==NULL) {
    set_istrm();
  }
  tdouble cv;
  istrm->get_value(cv);
  return cv;
}

const tdouble FunLinesInFile::calc()
{
  const std::string fstr = file->eval(false);
  std::ifstream inFile(fstr.c_str());
  return std::count(std::istreambuf_iterator<char>(inFile), 
             std::istreambuf_iterator<char>(), '\n');
}

const bool FunLinesInFile::search_circref(FlxFunction* fcr)
{
  return file->search_circref(fcr);
}

const std::string FunLinesInFile::write()
{
  return "lines_in_file(" + file->write() + ")";
}

FunBase* FunReadFunLinesInFile::read(bool errSerious)
{
  FlxString* file = new FlxString(false,errSerious);
  return new FunLinesInFile( file );
}

const tdouble FunRNDvecID::calc()
{
  tuint N = 0;
  const tdouble* wV = data->ConstMtxBox.get_Vec(wVec->eval(),N);
  flxVec v(wV,N);
  return tdouble(RndCreator->gen_smp_index2(v));
}

const std::string FunRNDvecID::write()
{
  return "rnd_vec_id(" + wVec->write() + ")";
}

const bool FunRNDvecID::search_circref(FlxFunction* fcr)
{
  return wVec->search_circref(fcr);
}

FunBase* FunReadFunRndVecID::read(bool errSerious)
{
  FlxMtxConstFun* wVec = new FlxMtxConstFun(true);
  return new FunRNDvecID(wVec);
}

FunObjExec::~FunObjExec()
{
  if (fun) delete fun; 
  delete InternListLoop;
}

const tdouble FunObjExec::calc()
{
  try {
  InternListLoop->loop_block_exec_1();
  } catch (FlxReturnBreakContinue_baseE &e) {
    InternListLoop->loop_block_exec_2();
    throw;
  }
  const tdouble res = (fun)?(fun->calc()):(*fres);
  InternListLoop->loop_block_exec_2();
  return res;
}

const std::string FunObjExec::write()
{
  std::ostringstream ssV;
  ssV << "objexec(";
  if (fun) {
    ssV << fun->write();
  } else {
    ssV << ':' << data->ConstantBox.get(fres);
  }
  ssV << ",{...}";
  ssV << ")";
  return ssV.str();
}

const bool FunObjExec::search_circref(FlxFunction* fcr)
{
  if (fun) {
    return fun->search_circref(fcr);
  } else {
    return false;
  }
}

const bool FunObjExec::optimize(FunBasePtr& optf, const Fun_OptimizeInfo& foi)
{
  if (fun) child_optimize(fun,foi);
  return false;
}

const bool FunObjExec::dependOn_Const(const tdouble* const thenumber)
{
  throw FlxException_NotImplemented("FunObjExec::dependOn_Const");
}

FunBase* FunReadObjExec::read(bool errSerious)
{
  FunBase* fun = NULL;
  tdouble *fres = NULL;
  if (reader->whatIsNextChar()==':') {
    reader->getChar(':',errSerious);
    fres = ConstantBox->get(reader->getWord(true,errSerious),true);
  } else {
    fun = FunctionList->read(errSerious);
  }
  try {
    reader->getChar(',',errSerious);
    FlxCodeBlock* loopblock = FlxObjReadCodeBlock::read_block(true,errSerious);
    if (fres) loopblock->add_internal_const(fres);
    return new FunObjExec(fun,fres,loopblock);
  } catch (FlxException &e) {
    if (fun) delete fun;
    throw;
  }
}

FunCatchError::FunCatchError(std::vector< FunBase* >* ParaListV)
: FunBaseFun_multPara(ParaListV)
{

}

const tdouble FunCatchError::calc()
{
  try {
    return ParaListP[0]->calc();
  } catch (FlxException_NeglectInInteractive &e) {
    return ParaListP[1]->calc();
  } catch (FlxException &e) {
    if (ParaListP[2]->calc()>ZERO) {
      return ParaListP[1]->calc();
    }
    throw;
  }
}


FunBase* FunReadFunRoot::read(bool errSerious)
{
  FunBase* funR = NULL;
  FunBase *startF=NULL,*endF=NULL,*dy=NULL,*dx=NULL;
  ostreamp streamp = NULL;
  try {
    // get the constant
      tdouble* d1 = read_const_var(errSerious,true);
    reader->getChar(',',errSerious);
    reader->getChar('[',errSerious);
    startF = FunctionList->read(errSerious);
    reader->getChar(',',errSerious);
    endF = FunctionList->read(errSerious);
    reader->getChar(']',errSerious);
    reader->getChar(',',errSerious);
    funR = FunctionList->read(errSerious);
    reader->getChar(',',errSerious);
    const std::string method = reader->getWord(true,errSerious);
    tuint ID;
    // which method do we have to call?
      if (method == "bisec") {
        ID = 0;
      } else if (method == "rgfsi") {
        ID = 1;
      } else {
        std::ostringstream ssV;
        ssV << "Unknown method '" << method << "' for root-serach.";
        throw FlxException("FunReadFunRoot::read_1", ssV.str() );
      }
    while (reader->whatIsNextChar()==',') {
      reader->getChar(',',errSerious);
      const std::string keyw = reader->getWord(true,false);
      reader->getChar('=',errSerious);
      if (keyw=="dy") {
        dy = FunctionList->read(errSerious);
      } else if (keyw=="dx") {
        dx = FunctionList->read(errSerious);
      } else if (keyw=="stream") {
        std::string stream = reader->getWord(true,false);
        streamp = (data->OstreamBox.get(stream));
      } else {
        std::ostringstream ssV;
        ssV << "Unknown parameter-name '" << keyw << "' for root-serach.";
        throw FlxException("FunReadFunRoot::read_2", ssV.str() );
      }
    }
    if (dx==NULL) {
      dx = new FunNumber(1e-6);
    }
    if (dy==NULL) {
      dy = new FunNumber(1e-8);
    }
    return new FunRoot(ID,funR,d1,startF,endF,dy,dx,streamp);  
  } catch (FlxException &e) {
    FLXMSG("FunReadFunRoot::read_3",1);
    if (funR) delete funR;
    if (startF) delete startF;
    if (endF) delete endF;
    if (dx) delete dx;
    if (dy) delete dy;
    throw;
  }
}


// ---------------------------------------------------------------------------------------------------

FlxObjRndSmp::~FlxObjRndSmp()
{
  if (rbrvsets) delete rbrvsets;
  if (RndBox) delete RndBox;
}

void FlxObjRndSmp::task() 
{
  if (rbrvsets) {
    const std::vector<std::string> set_str_vec = parse_strseq_as_vec(rbrvsets->eval(true));
    RndBox = new RBRV_constructor(set_str_vec,data->rbrv_box);
    delete rbrvsets;
    rbrvsets = NULL;
  }
  RndBox->gen_smp();
  //GlobalVar.slog(4) << "Random variable: samples generated." << std::endl;
}

FlxObjBase* FlxObjReadRndSmp::read()
{
  FlxString* rbrvsets = new FlxString(false,false);
  try {
    read_optionalPara(false);
    return new FlxObjRndSmp(get_doLog(),rbrvsets);
  } catch (FlxException &e) {
    delete rbrvsets;
    throw;
  }
}

void FlxObjRndTrackRecord::task()
{
  tuint nmb = numb->cast2int();
  if ( numb == 0) return;
  sout() << rv_normal();
  for (tuint i = 1; i < nmb; ++i) {
    sout() << std::endl << rv_normal();
  }
  sout() << std::endl;
}

void FlxObjRndTrackReplay::task()
{
  data->RndCreator.replay_start( &data->IstreamBox.get(isname->eval_word(true)) );
}

void FlxObjRndSeed::task()
{
  rv_initialize(false, true, seedv->cast2tuint(false), icv->cast2tuintW0(false));
}

FlxObjBase* FlxObjReadRndTrack::read()
{
  std::string keyS2 = reader->getWord(true,false);
  if (keyS2 == "record") {
    FlxFunction* numb = new FlxFunction(funReader,false);
    try {
      read_optionalPara(false);
      return new FlxObjRndTrackRecord(get_doLog(),numb,get_stream(), get_verbose());
    } catch (FlxException &e) {
      FLXMSG("FlxObjReadRndTrack::read_3",1);
      delete numb;
      throw;
    }
  } else if ( keyS2 == "replay" ) {
    FlxString* isname = new FlxString(false,false);
    read_optionalPara(false);
    return new FlxObjRndTrackReplay(get_doLog(),isname);
  } else if ( keyS2 == "stop" ) {
    read_optionalPara(false);
    return new FlxObjRndTrackStop(get_doLog());
  } else {
    std::ostringstream ssV;
    ssV << "Unknown keyword '" << keyS2 << "'.";
    throw FlxException_NeglectInInteractive("FlxObjReadRndTrack::read_4", ssV.str(), reader->getCurrentPos() );
  }
}

FlxObjBase* FlxObjReadRndSeed::read()
{
  FlxFunction* seedv = NULL; FlxFunction* icv = NULL;
  try {
    seedv = new FlxFunction(funReader,false);
    if (reader->whatIsNextChar() == ';') {
      icv = new FlxFunction(new FunNumber(ZERO));
    } else {
      icv = new FlxFunction(funReader,false);
    }
    read_optionalPara(false);
    return new FlxObjRndSeed(get_doLog(),seedv, icv);
  } catch (FlxException &e) {
    FLXMSG("FlxObjReadRnd::read_5",1);
    if (seedv) delete seedv;
    if (icv) delete icv;
    throw;
  }
}








