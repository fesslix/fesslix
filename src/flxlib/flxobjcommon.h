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

#include "flxobjects.h"


class FlxCreateObjReaders_Common : public FlxCreateObjReaders {
  public:
    void createObjReaders (FlxObjectReadBox* objReadBox );
    void createFunReaders (FlxData* dataBox );
};




// -------------------------- OBJECTS ----------------------------------------------------------------------



/**
* @brief object class: define a variable expression
*
* var varName = FlxFunction
*/
class FLXLIB_EXPORT FlxObjVar : public FlxObjBase {
  private:
    const std::string cname;
    FlxFunction *fun;
    const bool put_on_IL;
    
    void task();
  public:
    FlxObjVar ( const bool dolog, const std::string& cnameV, FlxFunction* funV, const bool put_on_IL ) : FlxObjBase(dolog), cname(cnameV), fun(funV), put_on_IL(put_on_IL) {}
    ~FlxObjVar();
};

/**
* @brief object class: define a function
*
* fun funName = FlxFunction
*/
class FLXLIB_EXPORT FlxObjFun : public FlxObjBase {
  private:
    const std::string cname;
    FunReadFunBase *fun;
    void task();
  public:
    FlxObjFun ( const bool dolog, const std::string& cnameV, FunReadFunBase* funV ) : FlxObjBase(dolog), cname(cnameV), fun(funV) {}
    ~FlxObjFun();
};

class FLXLIB_EXPORT FlxObjEnd : public FlxObjBase {
  private:
    void task() { throw FlxEndE(); }
  public:
    FlxObjEnd(const bool dolog) : FlxObjBase(dolog) {};
};

class FLXLIB_EXPORT FlxObjExit : public FlxObjBase {
  private:
    void task() { throw FlxExitE(); }
  public:
    FlxObjExit(const bool dolog) : FlxObjBase(dolog) {};
};

class FLXLIB_EXPORT FlxObjReturn : public FlxObjBase {
  private:
    void task() { throw FlxReturnE(); }
  public:
    FlxObjReturn(const bool dolog) : FlxObjBase(dolog) {};
};

class FLXLIB_EXPORT FlxObjContinue : public FlxObjBase {
  private:
    void task() { throw FlxContinueE(); }
  public:
    FlxObjContinue(const bool dolog) : FlxObjBase(dolog) {};
};

class FLXLIB_EXPORT FlxObjBreak : public FlxObjBase {
  private:
    void task() { throw FlxBreakE(); }
  public:
    FlxObjBreak(const bool dolog) : FlxObjBase(dolog) {};
};

class FLXLIB_EXPORT FlxObjThrow : public FlxObjBase {
  private:
    void task() { throw FlxException_NeglectInInteractive("FlxObjThrow::task","User initiated error call."); }
  public:
    FlxObjThrow(const bool dolog) : FlxObjBase(dolog) {};
};

/**
* @brief object class: calculate an arithmetic expression and output the result
*
* calc FlxFunction
*/
class FLXLIB_EXPORT FlxObjCalc : public FlxObjOutputBase {
  private:
    FlxFunction *fun;
    tdouble *ansptr;
    void task();
  public:
    FlxObjCalc ( const bool dolog, FlxFunction* funV, const std::string& ostreamV, const bool checkTOL );
    ~FlxObjCalc() { delete fun; }
};


/**
* @brief object class: output a string
*
* warranty;
*/
class FLXLIB_EXPORT FlxObjWarranty : public FlxObjOutputBase {
  private:
    void task();
  public:
    FlxObjWarranty ( bool dolog, std::string ostreamV ) : FlxObjOutputBase(dolog,ostreamV) {}
    ~FlxObjWarranty() {}
};


/**
* @brief object class: runs an external command
*
* run_external "COMMAND";
*/
class FLXLIB_EXPORT FlxObjRunExternal : public FlxObjOutputBase {
  private:
    FlxString* ext_cmd;
    const bool bthrow;
    void task();
  public:
    FlxObjRunExternal ( bool dolog, FlxString* ext_cmd, std::string ostreamV, const bool bthrow ) 
      : FlxObjOutputBase(dolog,ostreamV), ext_cmd(ext_cmd), bthrow(bthrow) {}
    ~FlxObjRunExternal();
};


/**
* @brief object class: executes Fesslix 'a second time'
*
* run_files ID FILENAME;
*/
class FLXLIB_EXPORT FlxObjRunExternal_Files : public FlxObjOutputBase {
  private:
    const std::string jid;
    FlxString* filename;
    FlxString* destination;
    void task();
  public:
    FlxObjRunExternal_Files ( bool dolog, const std::string& jid, FlxString* filename, FlxString* destination, std::string ostreamV ) 
      : FlxObjOutputBase(dolog,ostreamV), jid(jid), filename(filename), destination(destination) {}
    ~FlxObjRunExternal_Files();
};

/**
* @brief object class: file-filter for constants/variables
*
* file_filter_cv (filename)
*/
class FLXLIB_EXPORT FlxObjFileFilterCV : public FlxObjOutputBase {
  private:
    const std::string s_init;
    const std::string s_end;
    FlxString* filename;
    const bool tprec;
    void task();
  public:
    FlxObjFileFilterCV ( bool dolog, FlxString* filename, std::string ostreamV, const std::string s_init, const std::string s_end, const bool tprec ) 
      : FlxObjOutputBase(dolog,ostreamV), s_init(s_init), s_end(s_end), filename(filename), tprec(tprec) { }
    ~FlxObjFileFilterCV() { if (filename) delete filename; }
    
    void parse_str(const std::string& str, std::ostream& sot);
};

/**
* @brief object class: file-filter for SOFiSTiK
*/
class FLXLIB_EXPORT FlxObjFileFilterSOFiSTiK : public FlxObjBase {
  private:
    FlxString* filename;
    const std::string syst_stream;
    const std::string mat_stream;
    tdouble& cvar;
    tdouble& cvar2;
    const std::string mat_string;
    FlxObjBase* block;
    FlxObjFileFilterCV* fcv;
    FlxMtxConstFun* midF;
    FlxFunction* mid_startF;
    void task();
  public:
    FlxObjFileFilterSOFiSTiK ( bool dolog, FlxString* filename, const std::string& syst_stream, const std::string& mat_stream, tdouble& cvar, tdouble& cvar2, const std::string& mat_string, FlxObjBase* block, FlxMtxConstFun* midF, FlxFunction* mid_startF );
    ~FlxObjFileFilterSOFiSTiK() { delete filename; delete block; delete fcv; delete midF, delete mid_startF; }
};

/**
* @brief object class: creates a new output-stream (write into a file)
*
* filestream streamname filename
*/
class FLXLIB_EXPORT FlxObjFileStream : public FlxObjBase {
  private:
    FlxString* streamname;
    FlxString* filename;
    const bool trunc;
    void task();
  public:
    /**
    * @param streamnameV MUST BE LOWERCASE
    */
    FlxObjFileStream ( const bool dolog, FlxString* streamname, FlxString* filename, const bool truncV = true ) 
      : FlxObjBase(dolog), streamname(streamname), filename(filename), trunc(truncV) {}
    ~FlxObjFileStream() { delete streamname, delete filename; }
};

/**
* @brief object class: creates a new string-stream (write into a string-buffer)
*
* filestream streamname filename
*/
class FLXLIB_EXPORT FlxObjStringStream : public FlxObjBase {
  private:
    FlxString* streamname;
    void task();
  public:
    /**
    * @param streamname MUST BE LOWERCASE
    */
    FlxObjStringStream ( const bool dolog, FlxString* streamname ) 
      : FlxObjBase(dolog), streamname(streamname) {}
    ~FlxObjStringStream() { delete streamname; }
};

/**
* @brief object class: creates a new output-stream (write into a file)
*
* ostream_close streamname
*/
class FLXLIB_EXPORT FlxObjOStream_close : public FlxObjBase {
  private:
    FlxString* streamname;
    void task();
  public:
    /**
    * @param streamnameV MUST BE LOWERCASE
    */
    FlxObjOStream_close ( const bool dolog, FlxString* streamname ) 
      : FlxObjBase(dolog), streamname(streamname) {}
    ~FlxObjOStream_close() { delete streamname; }
};

/**
* @brief object class: input file stream (read numbers from a file)
*/
class FLXLIB_EXPORT FlxObjInputFileStream : public FlxObjBase {
  protected:
    FlxString* streamname;
    FlxString* filename;
    FlxFunction* bs;
    const bool erreof;
    FlxFunction* CnumbF;
    FlxString* CvecF;
    const bool binary;
    const bool binaryfloat;
    
    const tuint get_Columns(std::vector<tuint> &Cvec, const bool errSerious=true);
    virtual void task();
  public:
    /**
    * @param streamname MUST BE LOWERCASE
    */
    FlxObjInputFileStream ( const bool dolog, FlxString* streamname, FlxString* filename, FlxFunction* bs, FlxFunction* CnumbF, FlxString* CvecF, const bool erreofV = false, const bool binary = false, const bool binaryfloat = false ) 
      : FlxObjBase(dolog), streamname(streamname), filename(filename), bs(bs), erreof(erreofV), CnumbF(CnumbF), CvecF(CvecF), binary(binary), binaryfloat(binaryfloat) {}
    virtual ~FlxObjInputFileStream();
};

class FLXLIB_EXPORT FlxObjInputFileStreamCombine : public FlxObjInputFileStream {
  private:
    std::vector<FlxString*> filename_vec;
    std::vector<FlxFunction*> weightfun_vec;

    void task();
  public:
    FlxObjInputFileStreamCombine ( const bool dolog, FlxString* streamname, std::vector<FlxString*>& filename_vec, std::vector<FlxFunction*>& weightfun_vec, FlxFunction* bs, FlxFunction* CnumbF, FlxString* CvecF, const bool erreofV = false );
    ~FlxObjInputFileStreamCombine();
};

/**
* @brief object class: input vector stream (read numbers from a vector)
*/
class FLXLIB_EXPORT FlxObjInputVectorStream : public FlxObjBase {
  private:
    FlxString* streamname;
    FlxString* inputStreamName;
    FlxFunction* Nreserve;
    const bool erreof;
    void task();
  public:
    /**
    * @param streamname MUST BE LOWERCASE
    */
    FlxObjInputVectorStream ( const bool dolog, FlxString* streamname, FlxString* inputStreamName, FlxFunction* Nreserve, const bool erreofV = false ) 
      : FlxObjBase(dolog), streamname(streamname), inputStreamName(inputStreamName), Nreserve(Nreserve), erreof(erreofV) {}
    ~FlxObjInputVectorStream();
};

/**
* @brief object class: append a number to an input vector stream
*/
class FLXLIB_EXPORT FlxObjivstream_append : public FlxObjBase {
  private:
    FlxIstream_vector* istream;
    FlxFunction* fun;
    FlxMtxConstFun* mcf;
    FlxString* sn;
    void task();
  public:
    FlxObjivstream_append ( bool dolog,FlxString* sn, FlxFunction* fun, FlxMtxConstFun* mcf )
      : FlxObjBase(dolog), istream(NULL), fun(fun), mcf(mcf), sn(sn) {}
    ~FlxObjivstream_append() { if (fun) delete fun; if (mcf) delete mcf; delete sn; }
};

/**
* @brief object class: clear an input vector stream
*/
class FLXLIB_EXPORT FlxObjivstream_clear : public FlxObjBase {
  private:
    FlxString* sn;
    const bool is_reset;
    void task();
  public:
    FlxObjivstream_clear ( bool dolog, FlxString* sn, const bool is_reset ) : FlxObjBase(dolog), sn(sn), is_reset(is_reset) {}
    ~FlxObjivstream_clear() { delete sn; }
};

/**
* @brief object class: writes an input stream to an output stream
*/
class FLXLIB_EXPORT FlxObjistream_write : public FlxObjOutputBase {
  private:
    FlxString* isname;
    void task();
  public:
    FlxObjistream_write ( bool dolog, FlxString* isname, std::string ostreamV )
      : FlxObjOutputBase(dolog,ostreamV), isname(isname) {}
    ~FlxObjistream_write() { delete isname; }
};

/**
* @brief object class: reads input from a file
*
* read "PATHNAME";
*/
class FLXLIB_EXPORT FlxObjReadFile : public FlxObjBase, public FlxEvaluateCmdBase {
  private:
    FlxString* filename;
    void task();
  public:
    FlxObjReadFile ( bool dolog, FlxString* filenameV ) : FlxObjBase(dolog), filename(filenameV) {}
    ~FlxObjReadFile() { delete filename; }
};

/**
* @brief object class: creates a new Distributor-Stream
*
* filestream streamname filename
*/
class FLXLIB_EXPORT FlxObjDistributorStream : public FlxObjBase {
  private:
    FlxString* streamname;
    FlxString* stream1;
    FlxString* stream2;
    void task();
  public:
    /**
    * @note all streamnames MUST BE LOWERCASE
    */
    FlxObjDistributorStream ( const bool dolog, FlxString* streamnameV, FlxString* stream1V, FlxString* stream2V ) 
      : FlxObjBase(dolog), streamname(streamnameV), stream1(stream1V), stream2(stream2V) {}
    ~FlxObjDistributorStream() {delete streamname; delete stream1; delete stream2; };
};

/**
* @brief object class: set default values of optional parameters
*
* default DefSetName = value;
*/
class FLXLIB_EXPORT FlxObjDefault : public FlxObjBase {
  private:
    void* value;
    void task();
    FlxOptionalParaBase* para;
  public:
    FlxObjDefault( bool dolog, void* valueV, FlxOptionalParaBase* paraV) : FlxObjBase(dolog), value(valueV), para(paraV) {};
    ~FlxObjDefault();
};

/**
* @brief object class: control structure - if
*
* if (FlxFunction) { ... } [ else { ... }; ]
*/
class FLXLIB_EXPORT FlxObjIf : public FlxObjBase {
  private:
    FlxFunction* funIf;
    FlxObjBase* InternListThen;
    FlxObjBase* InternListElse;
    void task() {
      if ( funIf->calc() > ZERO ) { if (InternListThen != NULL) InternListThen->exec(); }
      else { if (InternListElse != NULL) InternListElse->exec(); }
    }
  public:
    FlxObjIf ( bool dolog, FlxFunction* funIfV, FlxObjBase* InternListThenV, FlxObjBase* InternListElseV ) : FlxObjBase(dolog), funIf(funIfV), InternListThen(InternListThenV), InternListElse(InternListElseV) {}
    ~FlxObjIf();
};

/**
* @brief object class: control structure - while-loop
*
* while ( FlxFunction ) { ... };
*/
class FLXLIB_EXPORT FlxObjWhile : public FlxObjBase {
  private:
    unsigned int maxCycles;
    FlxFunction* funWhile;
    FlxObjBase* InternListLoop;
    void task();
  public:
    FlxObjWhile ( bool dolog, FlxFunction* funWhileV, FlxObjBase* InternListLoopV, unsigned int MaxCycles=0) 
      : FlxObjBase(dolog), maxCycles(MaxCycles), funWhile(funWhileV), InternListLoop(InternListLoopV) {}
    ~FlxObjWhile();
};

/**
* @brief object class: control structure - for-loop
*
* for ( const-name, flxfunction, flxfunction ) { ... };
*/
class FLXLIB_EXPORT FlxObjFor : public FlxObjBase {
  private:
    unsigned int maxCycles;
    tdouble* theConst;
    FlxFunction* funCond;
    FlxFunction* funConst;
    FlxObjBase* InternListLoop;
    FlxObjConst* ConstDef;
    void task();
  public:
    FlxObjFor ( bool dolog, tdouble* theConstV, FlxObjConst* ConstDefV, FlxFunction* funCondV, FlxFunction* funConstV, FlxObjBase* InternListLoopV, unsigned int MaxCycles=0 ) 
      : FlxObjBase(dolog), maxCycles(MaxCycles), theConst(theConstV), funCond(funCondV), funConst(funConstV), InternListLoop(InternListLoopV), ConstDef(ConstDefV) {}
    ~FlxObjFor();
};

/**
* @brief object class: control structure - for-loop
*
* sfor ( const-name; to; start0 ) { ... };
*/
class FLXLIB_EXPORT FlxObjSFor : public FlxObjBase {
  private:
    tdouble* theConst;
    FlxFunction* funTo;
    bool start0;
    FlxObjBase* InternListLoop;
    void task();
  public:
    FlxObjSFor ( bool dolog, tdouble* theConst, FlxFunction* funTo, bool start0, FlxObjBase* InternListLoop ) 
      : FlxObjBase(dolog), theConst(theConst), funTo(funTo), start0(start0), InternListLoop(InternListLoop) {}
    ~FlxObjSFor() { delete funTo; delete InternListLoop;}
};

/**
* @brief object class: control structure - for-each-loop
*
* for_each VAR in ( ITERSTR; "SEP" ) { ... };
*/
class FLXLIB_EXPORT FlxObjForEach : public FlxObjBase {
  private:
    std::string& loop_var;
    FlxString* iterstr;
    const std::string sep_char;
    FlxObjBase* InternListLoop;
    const bool do_trim;
    void task();
  public:
    FlxObjForEach ( bool dolog, std::string& loop_var, FlxString* iterstr, const std::string& sep_char, FlxObjBase* InternListLoop, const bool do_trim ) 
      : FlxObjBase(dolog), loop_var(loop_var), iterstr(iterstr), sep_char(sep_char), InternListLoop(InternListLoop), do_trim(do_trim) {}
    ~FlxObjForEach() { delete iterstr; delete InternListLoop; }
};

/**
* @brief object class: catches minor errors
*
* catch_error { ... };
*/
class FLXLIB_EXPORT FlxObjCatch_Error : public FlxObjBase {
  private:
    tdouble* theConst;
    FlxFunction* funTo;
    bool start0;
    FlxObjBase* ill;
    const bool errSerious;
    FlxObjBase* catchblock;
    
    void task();
  public:
    FlxObjCatch_Error ( bool dolog, FlxObjBase* ill, FlxObjBase* catchblock, const bool errSerious ) 
      : FlxObjBase(dolog), ill(ill), errSerious(errSerious), catchblock(catchblock) {}
    ~FlxObjCatch_Error() { delete ill; delete catchblock; }
};

/**
* @brief object class: timer
*
* time { CODE_TO_STOP; };
*/
class FLXLIB_EXPORT FlxObjTime : public FlxObjOutputBase {
  private:
    FlxObjBase* InternListTime;
    tdouble rt;                // time needed to read the object
    const bool spt;        // store physical time (and not CPU time)
    void task();
  public:
    FlxObjTime ( bool dolog, FlxObjBase* InternListTime, std::string ostreamV, tdouble rt, const bool spt ) : FlxObjOutputBase(dolog, ostreamV), InternListTime(InternListTime), rt(rt), spt(spt) {}
    ~FlxObjTime();
};

/**
* @brief object class: define a timer
*
* timer define TIMER_NAME
*/
class FLXLIB_EXPORT FlxObjTimerDefine : public FlxObjBase {
  private:
    const std::string tname;
    void task() { data->TimerBox.insert(tname, new FlxTimer() ); 
      GlobalVar.slog(4) << "timer: timer '" << tname << "' defined." << std::endl; }
  public:
    FlxObjTimerDefine ( const bool dolog, const std::string& tnameV ) : FlxObjBase(dolog), tname(tnameV) {}
};

/**
* @brief object class: start a timer
*
* timer start TIMER_NAME
*/
class FLXLIB_EXPORT FlxObjTimerStart : public FlxObjBase {
  private:
    std::string tname;
    void task() { GlobalVar.slog(4) << "timer: timer '" << tname << "' started. (t=" << GlobalVar.Double2String(data->TimerBox.get(tname)->get_time()) << ")" << std::endl;
      data->TimerBox.get(tname)->start(); }
  public:
    FlxObjTimerStart ( bool dolog, std::string tnameV ) : FlxObjBase(dolog), tname(tnameV) {}
};

/**
* @brief object class: stop a timer
*
* timer stop TIMER_NAME
*/
class FLXLIB_EXPORT FlxObjTimerStop : public FlxObjBase {
  private:
    const std::string tname;
    void task() { data->TimerBox.get(tname)->stop();
      GlobalVar.slog(4) << "timer: timer '" << tname << "' stopped. (t=" << GlobalVar.Double2String(data->TimerBox.get(tname)->get_time()) << ")" << std::endl; }
  public:
    FlxObjTimerStop ( const bool dolog, const std::string& tnameV ) : FlxObjBase(dolog), tname(tnameV) {}
};

/**
* @brief object class: delete a timer
*
* timer delete TIMER_NAME
*/
class FLXLIB_EXPORT FlxObjTimerDelete : public FlxObjBase {
  private:
    const std::string tname;
    void task() { data->TimerBox.deleteEl(tname); 
      GlobalVar.slog(4) << "timer: timer '" << tname << "' deleted." << std::endl; }
  public:
    FlxObjTimerDelete ( const bool dolog, const std::string& tnameV ) : FlxObjBase(dolog), tname(tnameV) {}
};

/**
* @brief object class: print time of a timer
*
* timer print TIMER_NAME
*/
class FLXLIB_EXPORT FlxObjTimerPrint : public FlxObjOutputBase {
  private:
    const std::string tname;
    void task();
  public:
    FlxObjTimerPrint ( const bool dolog, const std::string& tnameV, const std::string& ostreamV ) : FlxObjOutputBase(dolog,ostreamV), tname(tnameV)  {}
};

/**
* @brief object class: plot some functions (for GNUplot)
*
* funplot f1,f2,...;
*/
class FLXLIB_EXPORT FlxObjFunPlot : public FlxObjOutputBase {
  protected:
    std::vector<int> B;
    std::vector<FlxFunction*> V;
    std::vector<FlxMtxConstFun*> M;
    std::vector<FlxString*> S;
    const std::string sep_str;
    const std::string init_str;
    const std::string end_str;
    const bool binary;
    const bool binaryfloat;
    virtual void task();
    void check_first(std::ostream& sout_, bool& not_first);
  public:
    FlxObjFunPlot(const bool dolog,const std::vector<int>& B, const std::vector<FlxFunction*>& V, const std::vector<FlxMtxConstFun*>& M, std::vector<FlxString*>& S, const std::string& ostreamV, const bool check_TOL, const int prec, const int fixW, const std::string& boost_str, const std::string& sep_str, const std::string& init_str, const std::string& end_str, const bool binary, const bool binaryfloat): FlxObjOutputBase(dolog,ostreamV,false,check_TOL,prec,fixW,boost_str), B(B), V(V), M(M), S(S), sep_str(sep_str), init_str(init_str), end_str(end_str), binary(binary), binaryfloat(binaryfloat) {}
    ~FlxObjFunPlot();
};

/**
* @brief object class: plot some functions (for GNUplot)
*
* funplot_header "f1","f2",...;
*/
class FLXLIB_EXPORT FlxObjFunPlot_header : public FlxObjOutputBase {
  protected:
    std::vector<std::string> hdr;
    bool executed;
    const bool only_once;
    virtual void task();
  public:
    FlxObjFunPlot_header(const bool dolog,const std::vector<std::string>& hdr, const std::string ostreamV, const int prec, const int fixW, const bool only_once): FlxObjOutputBase(dolog,ostreamV,false,false,prec,fixW), hdr(hdr), executed(false), only_once(only_once) {}
    static void write_entry(std::string hestr, std::ostream& os, const int prec, const int fixW, const bool start);
};


/**
* @brief object class: Counting Intervals
*
* IntervalCount (Number) { ... };
*/
class FLXLIB_EXPORT FlxObjIntervalCount : public FlxObjBase {
  protected:
    /**
    * @brief number of intervals
    */
    tuint Np;
    /**
    * @brief current number of intervals
    */
    tuint countI;
    /**
    * @brief the function used to calculate number of intervals
    */
    FlxFunction* funNp;
    /**
    * @brief the MCI_CODE
    */
    FlxObjBase* InternList;
    void task();
  public:
    FlxObjIntervalCount ( bool dolog, FlxFunction* funNp, FlxObjBase* InternList ) : FlxObjBase(dolog), countI(0),funNp(funNp),InternList(InternList) {}
    virtual ~FlxObjIntervalCount();
};

/**
* @brief object class: Counting Intervals
*
* Filter (i;MtxFun) { ... };
*/
class FLXLIB_EXPORT FlxObjFilter : public FlxObjBase {
  protected:
    /**
    * @brief pointer to the constant (in ConstBox)
    */
    tdouble* cv;
    /**
    * @brief the matrix to run through
    */
    FlxMtxConstFun* seqMtx;
    /**
    * @brief the FILTER_CODE
    */
    FlxObjBase* block;
    void task();
  public:
    FlxObjFilter ( bool dolog, tdouble* cv, FlxMtxConstFun* seqMtx, FlxObjBase* block ) : FlxObjBase(dolog), cv(cv),seqMtx(seqMtx),block(block) {}
    virtual ~FlxObjFilter();
};

/**
* @brief object class: set configuration-options for floating point conversion
*/
class FLXLIB_EXPORT FlxObjFloatingPointConversion : public FlxObjBase {
  private:
    int id;
    FlxFunction *fun;
    void task();
  public:
    FlxObjFloatingPointConversion ( bool dolog, FlxFunction* funV, int idV ) : FlxObjBase(dolog), id(idV), fun(funV) {}
    ~FlxObjFloatingPointConversion() { delete fun; }
};

/**
* @brief object class: set some default values at runtime
*/
class FLXLIB_EXPORT FlxObjSetVariousDefault : public FlxObjBase {
  private:
    const tuint ID;
    FlxFunction *fun;
    void task();
  public:
    /**
    * 1 set loglevel at runtime
    * 2 set default polynomial degree of FlxFunDeg
    */
    FlxObjSetVariousDefault ( const bool dolog, const int ID, FlxFunction* fun ) : FlxObjBase(dolog), ID(ID), fun(fun) {}
    ~FlxObjSetVariousDefault() { delete fun; }
};


// -------------------------- READERS ----------------------------------------------------------------------

/**
* @brief object read class: for FlxObjVar
*/
class FlxObjReadVar : public FlxObjReadFCVbase {
  private:
    bool put_on_IL;
  public:
    FlxObjReadVar(const bool put_on_IL = true) : put_on_IL(put_on_IL) { if (put_on_IL) data->IgnoreBox.set_iL_recur("var"); }
    FlxObjBase* read ();
};

/**
* @brief object read class: for FlxObjFun
*/
class FlxObjReadFun : public FlxObjReadFCVbase {
  protected:
    const std::string get_name();
  public:
    FlxObjReadFun() { data->IgnoreBox.set_iL_recur("fun"); }
    FlxObjBase* read ();
};


class Interpolate_help : public FlxDataBase {
  private:
    FlxMtxConstFun* mtxfun;
    flxVec* xvec;
    flxVec* yvec;
  public:
    Interpolate_help(FlxMtxConstFun* mtxfun) : mtxfun(mtxfun), xvec(NULL), yvec(NULL) {}
    ~Interpolate_help();
    /**
    * @brief read the content of matrix mtxfun and assigns it to vectors xvec and yvec
    */
    void initialize();
    flxVec& get_xvec();
    flxVec& get_yvec();
};

/**
* @brief arithmetic class: cosinus
*/
class FunInterpolate : public FunBaseFun_onePara {
  private:
    const std::string fname;
    Interpolate_help* data;        // memory not managed by this class!
  public:
    FunInterpolate (std::vector<FunBase*> *ParaListV, const std::string fname, Interpolate_help* data) : FunBaseFun_onePara(ParaListV), fname(fname), data(data) {};
    const tdouble calc();
    const std::string write_v() { return fname;}
};

class FunReadFunInterpolate : public FunReadFunBase {
  private:
    const std::string fname;
    Interpolate_help data;
  public:
    FunReadFunInterpolate ( const std::string fname, FlxMtxConstFun* mtxfun) : fname(fname), data(mtxfun) {};
    ~FunReadFunInterpolate() {}
    virtual void initialize();
    FunBase* read ( bool errSerious );
};

/**
* @brief object read class: for FlxObjFun
*/
class FlxObjReadInterpolate : public FlxObjReadFun {
  public:
    FlxObjReadInterpolate() { data->IgnoreBox.set_iL_recur("interpolate"); }
    FlxObjBase* read ();
};

/**
* @brief object read class: for FlxObjEnd
*/
class FlxObjReadEnd : public FlxObjReadBase {
  public:
    FlxObjBase* read ();
};

/**
* @brief object read class: for FlxObjExit
*/
class FlxObjReadExit : public FlxObjReadBase {
  public:
    FlxObjBase* read ();
};

/**
* @brief object read class: for FlxObjReturn
*/
class FlxObjReadReturn : public FlxObjReadBase {
  public:
    FlxObjBase* read ();
};

/**
* @brief object read class: for FlxObjContinue
*/
class FlxObjReadContinue : public FlxObjReadBase {
  public:
    FlxObjBase* read ();
};

/**
* @brief object read class: for FlxObjBreak
*/
class FlxObjReadBreak : public FlxObjReadBase {
  public:
    FlxObjBase* read ();
};

/**
* @brief object read class: for FlxObjThrow
*/
class FlxObjReadThrow : public FlxObjReadBase {
  public:
    FlxObjBase* read ();
};

/**
* @brief object read class: for FlxObjCalc
*/
class FlxObjReadCalc : public FlxObjReadOutputBase {
  public:
    FlxObjBase* read ();
};


/**
* @brief object read class: for FlxObjEcho
*/
class FlxObjReadEcho : public FlxObjReadOutputBase {
  public:
    FlxObjReadEcho();
    FlxObjBase* read ();
};

/**
* @brief object read class: for FlxObjWarranty
*/
class FlxObjReadWarranty : public FlxObjReadOutputBase {
  public:
    FlxObjBase* read ();
};

/**
* @brief object read class: for FlxObjFileFilterCV
*/
class FlxObjReadFileFilterCV : public FlxObjReadOutputBase {
  public:
    FlxObjReadFileFilterCV();
    FlxObjBase* read ();
};

/**
* @brief object read class: for FlxObjFileFilterSOFiSTiK
*/
class FlxObjReadFileFilterSOFiSTiK : public FlxObjReadBase {
  public:
    FlxObjBase* read ();
};

/**
* @brief object read class: for FlxObjRunExternal
*/
class FlxObjReadRunExternal : public FlxObjReadOutputBase {
  public:
    FlxObjReadRunExternal();
    FlxObjBase* read ();
};

/**
* @brief object read class: for FlxObjRunExternal_Files
*/
class FlxObjReadRunExternal_Files : public FlxObjReadOutputBase {
  public:
    FlxObjBase* read ();
};

/**
* @brief object read class: for FlxObjFileStream
*/
class FlxObjReadFileStream : public FlxObjReadBase {
  public:
    FlxObjReadFileStream();
    FlxObjBase* read ();
};

/**
* @brief object read class: for FlxObjFileStream
*/
class FlxObjReadStringStream : public FlxObjReadBase {
  public:
    FlxObjBase* read ();
};

/**
* @brief object read class: for FlxObjFileStream
*/
class FlxObjReadOStream_close : public FlxObjReadBase {
  public:
    FlxObjBase* read ();
};

/**
* @brief object read class: for FlxObjInputFileStream
*/
class FlxObjReadInputFileStream : public FlxObjReadBase {
  public:
    FlxObjReadInputFileStream();
    FlxObjBase* read ();
};

/**
* @brief object read class: for FlxObjInputFileStream
*/
class FlxObjReadInputFileStreamCombine : public FlxObjReadInputFileStream {
  public:
    FlxObjBase* read ();
};

/**
* @brief object read class: for FlxObjInputVectorStream
*/
class FlxObjReadInputVectorStream : public FlxObjReadBase {
  public:
    FlxObjReadInputVectorStream();
    FlxObjBase* read ();
};

/**
* @brief object read class: for FlxObjivstream_append
*/
class FlxObjReadivstream_append : public FlxObjReadBase {
  public:
    FlxObjBase* read ();
};

/**
* @brief object read class: for FlxObjivstream_clear
*/
class FlxObjReadivstream_clear : public FlxObjReadBase {
  private:
    const bool is_reset;
  public:
    FlxObjReadivstream_clear(const bool is_reset=false) : is_reset(is_reset) {}
    FlxObjBase* read ();
};

/**
* @brief object read class: for FlxObjistream_write
*/
class FlxObjReadistream_write : public FlxObjReadOutputBase {
  public:
    FlxObjBase* read ();
};

/**
* @brief object read class: for FlxObjFileStream
*/
class FlxObjReadReadFile : public FlxObjReadBase {
  public:
    FlxObjBase* read ();
};

/**
* @brief object read class: for FlxObjFileStream
*/
class FlxObjReadDistributorStream : public FlxObjReadBase {
  public:
    FlxObjBase* read ();
};

/**
* @brief object read class: for FlxObjDefault
*/
class FlxObjReadDefault : public FlxObjReadBase {
  private:
    FlxObjBase* read_special(std::string& key);
  public:
    FlxObjBase* read();
};

/**
* @brief object read class: for FlxObjIf
*/
class FlxObjReadIf : public FlxObjReadBase {
  public:
    FlxObjBase* read ();
};

/**
* @brief object read class: for FlxObjIf
*/
class FlxObjReadIf_no_read : public FlxObjReadBase {
  public:
    FlxObjBase* read ();
};

/**
* @brief object read class: for FlxObjWhile
*/
class FlxObjReadWhile : public FlxObjReadLoops {
  public:
    FlxObjBase* read ();
};

/**
* @brief object read class: for FlxObjFor
*/
class FlxObjReadFor : public FlxObjReadLoops {
  public:
    FlxObjBase* read ();
};

/**
* @brief object read class: for FlxObjSFor
*/
class FlxObjReadSFor : public FlxObjReadBase {
  public:
    FlxObjBase* read ();
};

/**
* @brief object read class: for FlxObjForEach
*/
class FlxObjReadForEach : public FlxObjReadBase {
  public:
    FlxObjReadForEach();
    FlxObjBase* read ();
};

/**
* @brief object read class: for FlxObjCatch_Error
*/
class FlxObjReadCatch_Error : public FlxObjReadBase {
  public:
    FlxObjReadCatch_Error();
    FlxObjBase* read ();
};

/**
* @brief object read class: for FlxObjSFor
*/
class FlxObjReadTime : public FlxObjReadOutputBase {
  public:
    FlxObjReadTime();
    FlxObjBase* read ();
};

/**
* @brief object read class: for Timer-Objects
*/
class FlxObjReadTimer : public FlxObjReadOutputBase {
  public:
    FlxObjBase* read();
};

/**
* @brief object read class: for FunPlot-Objects
*/
class FlxObjReadFunPlot : public FlxObjReadOutputBase {
  public:
    FlxObjReadFunPlot();
    FlxObjBase* read();
};

/**
* @brief object read class: for FunPlot_header-Objects
*/
class FlxObjReadFunPlot_header : public FlxObjReadOutputBase {
  public:
    FlxObjReadFunPlot_header();
    FlxObjBase* read();
};


/**
* @brief object read class: for FlxObjIntervalCount
*/
class FlxObjReadIntervalCount : public FlxObjReadBase {
  public:
    FlxObjBase* read ();
};

/**
* @brief object read class: for FlxObjFilter
*/
class FlxObjReadFilter : public FlxObjReadBase {
  public:
    FlxObjBase* read ();
};


// ----- isread_vec

/**
* @brief object class: defines a RBRV_vec
*/
class FLXLIB_EXPORT FlxObjISread_vec : public FlxObjBase {
  private:
    FlxString* VecConstStr;                // -> eval -> = vecName
    FlxFunction* dimF;
    FlxString* strV;
    FlxIstream* istrm;
    std::string strS;
    
    void set_istrm();
    void task();
  public:
    FlxObjISread_vec ( const bool dolog, FlxString* VecConstStr, FlxFunction* dimF, FlxString* strV ) : FlxObjBase(dolog), VecConstStr(VecConstStr), dimF(dimF), strV(strV), istrm(NULL) {}
    ~FlxObjISread_vec();
};

/**
* @brief object read class: for FlxObjRBRV_vec
*/
class FlxObjReadISread_vec : public FlxObjReadBase {
  public:
    FlxObjBase* read ();
};



// ----- sleep

/**
* @brief wait for a specified number of seconds
*/
class FLXLIB_EXPORT FlxObjSleep : public FlxObjBase {
  private:
    FlxFunction* time2wait;
    
    void task();
  public:
    FlxObjSleep ( const bool dolog, FlxFunction* time2wait ) : FlxObjBase(dolog), time2wait(time2wait) {}
    ~FlxObjSleep() { delete time2wait; }
};

/**
* @brief object read class: for FlxObjRBRV_vec
*/
class FlxObjReadSleep : public FlxObjReadBase {
  public:
    FlxObjBase* read ();
};



// -------------------------- FUNCTIONS ----------------------------------------------------------------------


/**
* @brief arithmetic class: isread (reads a value from an input stream)
*/
class FLXLIB_EXPORT FunIvStream_size : public FunBase, public FlxDataBase {
  private:
    FlxString* strV;
    FlxIstream_vector* istrm;
    std::string strS;
    void set_istrm();
  public:
    FunIvStream_size (FlxString* strV) : strV(strV), istrm(NULL) {};
    ~FunIvStream_size();
    const tdouble calc();
    const bool search_circref(FlxFunction* fcr);
    const std::string write();
    const tuint precedence() { return 0; }
    const bool dependOn_Const(const tdouble* const thenumber) { return false; }
    const bool optimize(FunBasePtr& optf, const Fun_OptimizeInfo &foi) { return false; }
};

class FunReadFunIvStream_size : public FunReadFunBase {
  public:
    FunBase* read ( bool errSerious );
};


/**
* @brief arithmetic class: isread (reads a value from an input stream)
*/
class FLXLIB_EXPORT FunISread : public FunBase, public FlxDataBase {
  private:
    FlxString* strV;
    FlxIstream* istrm;
    std::string strS;
    void set_istrm();
  public:
    FunISread (FlxString* strV) : strV(strV), istrm(NULL) {};
    ~FunISread();
    const tdouble calc();
    const bool search_circref(FlxFunction* fcr);
    const std::string write();
    const tuint precedence() { return 0; }
    const bool dependOn_Const(const tdouble* const thenumber) { return false; }
    const bool optimize(FunBasePtr& optf, const Fun_OptimizeInfo &foi) { return false; }
};

class FunReadFunISread : public FunReadFunBase {
  public:
    FunBase* read ( bool errSerious );
};


class FLXLIB_EXPORT FunLinesInFile : public FunBase, public FlxDataBase {
  private:
    FlxString* file;
    void set_istrm();
  public:
    FunLinesInFile (FlxString* file) : file(file) {};
    ~FunLinesInFile() { delete file; } 
    const tdouble calc();
    const bool search_circref(FlxFunction* fcr);
    const std::string write();
    const tuint precedence() { return 0; }
    const bool dependOn_Const(const tdouble* const thenumber) { return false; }
    const bool optimize(FunBasePtr& optf, const Fun_OptimizeInfo &foi) { return false; }
};

class FunReadFunLinesInFile : public FunReadFunBase {
  public:
    FunBase* read ( bool errSerious );
};

class FLXLIB_EXPORT FunRNDvecID : public FunBase, public FlxDataBase {
  private:
    FlxMtxConstFun* wVec;
  public:
    FunRNDvecID (FlxMtxConstFun* wVec) : wVec(wVec) {};
    ~FunRNDvecID() { if (wVec) delete wVec; } 
    const tdouble calc();
    const bool search_circref(FlxFunction* fcr);
    const std::string write();
    const tuint precedence() { return 0; }
    const bool dependOn_Const(const tdouble* const thenumber) { return false; }
    const bool optimize(FunBasePtr& optf, const Fun_OptimizeInfo &foi) { return false; }
};

class FunReadFunRndVecID : public FunReadFunBase {
  public:
    FunBase* read ( bool errSerious );
};

/**
* @brief arithmetic class: isread (reads a value from an input stream)
*/
class FLXLIB_EXPORT FunObjExec : public FunBase, public FlxDataBase {
  private:
    FunBase* fun;
    tdouble* fres;
    FlxCodeBlock* InternListLoop;
  public:
    FunObjExec (FunBase* fun, tdouble* fres, FlxCodeBlock* InternListLoop) : fun(fun), fres(fres), InternListLoop(InternListLoop) {};
    ~FunObjExec();
    const tdouble calc();
    const bool search_circref(FlxFunction* fcr);
    const std::string write();
    const tuint precedence() { return 0; }
    const bool dependOn_Const(const tdouble* const thenumber);
    const bool optimize(FunBasePtr& optf, const Fun_OptimizeInfo &foi);
};

class FunReadObjExec : public FunReadFunBase {
  public:FlxString* strV;
    FlxIstream* istrm;
    std::string strS;
    void set_istrm();
    FunBase* read ( bool errSerious );
};


class FunCatchError : public FunBaseFun_multPara {
  public:
    FunCatchError (std::vector<FunBase*> *ParaListV);
    const tdouble calc();
    const std::string write_v() { return "catch_error";}
};

class FunReadFunCatchError : public FunReadFunBase {
  public:
    FunBase* read ( bool errSerious ) { return new FunCatchError( read_parameters(3,errSerious) ); }
};


/**
* @brief the 'calc'-part of this object is defined in flxfunction_fun
*        the reader is defined here to allow implementation of the logging-feature
*/
class FunReadFunRoot : public FunReadFunBase, public FlxDataBase {
  public:
    FunBase* read ( bool errSerious );
};


// -------------------------------------------------------------------------------------------------------------


/**
* @brief object class: create random samples
*
* rnd smp;
*/
class FlxObjRndSmp : public FlxObjBase {
  protected:
    FlxString* rbrvsets;
    RBRV_constructor* RndBox;
    virtual void task();
  public:
    FlxObjRndSmp(const bool dolog,FlxString* rbrvsets): FlxObjBase(dolog),rbrvsets(rbrvsets), RndBox(NULL) {}
    ~FlxObjRndSmp();
};

class FlxObjReadRndSmp : public FlxObjReadOutputBase {
  public:
    FlxObjBase* read();
};

// -------------------------------------------------------------------------------------------------------------

/**
* @brief object class: create random samples
*
* rnd track record NUMB;
*/
class FlxObjRndTrackRecord : public FlxObjOutputBase {
  protected:
    FlxFunction* numb;
    virtual void task();
  public:
    FlxObjRndTrackRecord(bool dolog, FlxFunction* numb, std::string ostreamV, bool verbose): FlxObjOutputBase(dolog,ostreamV,verbose), numb(numb) {}
    ~FlxObjRndTrackRecord() { delete numb; }
};

/**
* @brief object class: stop reading of semi random samples
*
* rnd track replay "FILE";
*/
class FlxObjRndTrackReplay : public FlxObjBase {
  protected:
    FlxString* isname;
    virtual void task();
  public:
    FlxObjRndTrackReplay(bool dolog, FlxString* isname) : FlxObjBase(dolog), isname(isname) {};
    ~FlxObjRndTrackReplay() { delete isname; }
};

/**
* @brief object class: stop reading of semi random samples
*
* rnd track stop;
*/
class FlxObjRndTrackStop : public FlxObjBase {
  protected:
    virtual void task() {data->RndCreator.replay_stop(true);};
  public:
    FlxObjRndTrackStop(bool dolog) : FlxObjBase(dolog) {};
};

/**
* @brief object read class: for managed random variables
*/
class FlxObjReadRndTrack : public FlxObjReadOutputBase {
  public:
    FlxObjBase* read();
};

/**
* @brief object class: set a seed value
*
* rnd seed FlxFunction;
*/
class FlxObjRndSeed : public FlxObjBase {
  protected:
    FlxFunction* seedv;
    FlxFunction* icv;
    virtual void task();
  public:
    FlxObjRndSeed( bool dolog, FlxFunction* seedv, FlxFunction* icv) : FlxObjBase(dolog), seedv(seedv), icv(icv) {}
    ~FlxObjRndSeed() {delete seedv; delete icv;}
};

class FlxObjReadRndSeed : public FlxObjReadOutputBase {
  public:
    FlxObjBase* read();
};



