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

#define FLXLIB_CPP

#define fesslix_flxglobal_CPP

#include "flxmath.h"
#include "flxdefault.h"
#include "flxio.h"

#include <string>
#include <sstream>
#include <fstream>
#include <limits>
#include <iostream>
#include <iomanip>
#include <thread>


size_t bVec_NumbTrue(const bVec& bv)
{
  size_t N=0;
  for (size_t i=0;i<bv.size();++i) {
    if (bv[i]) ++N;
  }
  return N;
}

void FlxError(bool errSerious, std::string errnumber, std::string Titel, std::string Msg)
{
  if (errSerious) {
    throw FlxException(errnumber,Titel,Msg);
  } else {
    throw FlxException_NeglectInInteractive(errnumber,Titel,Msg);
  }
}

const int strGlobalVar::D2S_get_fixW(const int prec, const int fixW) const
{
  if (fixW!=0) return fixW;
  const tuint precTU = ((prec<0)?D2S_Prec:tuint(prec));
  return precTU + 7;
}

const std::string strGlobalVar::Double2String_maxPrec(tdouble dv)
{
  std::ostringstream os;
  os << std::setprecision( std::numeric_limits<double>::digits10 + 2 ) << dv;
  return os.str();
}

const std::string strGlobalVar::Double2String(tdouble dv, const bool checkTOL, const int prec, int fixW)
{
  if (checkTOL) {
    if (fabs(dv) <= GlobalVar.TOL()) dv = ZERO;
  }
  std::ostringstream ssV;
  // perform some settings (output-type)
    const tuint precTU = ((prec<0)?D2S_Prec:tuint(prec));
    ssV.precision(precTU);
    if (D2S_Type == 0) {
      // do nothing
    } else if (D2S_Type == 1) {
      ssV << std::scientific;
    } else if (D2S_Type == 2) {
      if ( fabs(dv)>=D2S_BValU || fabs(dv)<D2S_BValL ) {
        ssV << std::scientific;
      } else {
        if (dv>ONE) ssV.precision(precTU+1);
      }
    } else if (D2S_Type == 3) {
      ssV << std::fixed;
    }
  // create the number
    ssV << dv;
    std::string strV = ssV.str();
  // edit the number
    // delete 0s after the e
    if (D2S_Del0) {
      size_t iI = strV.rfind("e");
      if (iI != std::string::npos) {
        iI+=2;
        tuint n = 0;
        while ( strV[iI+n] == '0' ) {
          ++n;
        }
        strV.erase(iI, n);
      }
    }
    // delete + after e
    if (D2S_DelP) {
      size_t iI = strV.rfind("+");
      if (iI != std::string::npos) {
        strV.erase(iI,1);
      }
    }
    // delete 0s after .XYZ
    if (D2S_Del0 && strV.rfind(".") != std::string::npos) {
      size_t iI = strV.rfind("e");
      tuint n;
      if (iI != std::string::npos) {
        n = iI-1;
      } else {
        n = strV.length()-1;
      }
      tuint m = 0;
      while ( strV[n-m] == '0' ) {
        ++m;
      }
      strV.erase(n-m+1,m);
      iI = strV.rfind("e");
      if (iI != std::string::npos) {
        n = iI-1;
      } else {
        n = strV.length()-1;
      }
      if (strV[n] == '.') {
        strV.erase(n,1);
      }
    }
    // delete e if it is at the end of the string
    if (D2S_Del0) {
      size_t iI = strV.rfind("e");
      if (iI != std::string::npos) {
        if (iI == strV.length()-1) {
          strV.erase(iI,1);
        }
      }
    }
  // apply fixed width setting
    if (fixW==0) fixW = precTU + 7;
    if (fixW>0 && strV.length()<tuint(fixW)) {
      const tuint tufW = tuint(fixW);
      std::size_t it = strV.find(".");
      if (it==std::string::npos) it = strV.length();
      if (it < 2) {
        it = 2-it;
        if (strV.length()>tufW) {
          it = 0;
        } else if (strV.length()+it>tufW) {
          it = tufW-strV.length();
        }
        strV.insert(0,it,' ');
      }
      // align the e's
        it = strV.find(".");
        if (it!=std::string::npos) {
          std::size_t it2 = strV.find("e");
          if (it2!=std::string::npos) {
            if (it2>=it+1) {
              it = it2-(it+1);        // number of digits after the dot
              if (it<precTU) {
                it = precTU - it;
                strV.insert(it2,it,'0');
              }
            }
          }
        }
      // append after the number
        if (tufW>strV.length()) {
          it = tufW - strV.length();
          strV.insert(strV.length(),it,' ');
        }
    }
  return strV;
}

const std::string strGlobalVar::Double2String_sci(tdouble dv, const int prec, int fixW)
{
  std::ostringstream ssV;
  // perform some settings (output-type)
    const tuint precTU = ((prec<0)?D2S_Prec:tuint(prec));
    ssV.precision(precTU);
    ssV << std::scientific;
  // create the number
    ssV << dv;
    std::string strV = ssV.str();
  // edit the number
    // delete 0s after the e
    if (D2S_Del0) {
      size_t iI = strV.rfind("e");
      if (iI != std::string::npos) {
        iI+=2;
        tuint n = 0;
        while ( strV[iI+n] == '0' ) {
          ++n;
        }
        strV.erase(iI, n);
      }
    }
    // delete + after e
    if (D2S_DelP) {
      size_t iI = strV.rfind("+");
      if (iI != std::string::npos) {
        strV.erase(iI,1);
      }
    }
    // delete e if it is at the end of the string
    if (D2S_Del0) {
      size_t iI = strV.rfind("e");
      if (iI != std::string::npos) {
        if (iI == strV.length()-1) {
          strV.erase(iI,1);
        }
      }
    }
  // apply fixed width setting
    if (fixW==0) fixW = precTU + 7;
    if (fixW>0 && strV.length()<tuint(fixW)) {
      const tuint tufW = tuint(fixW);
      std::size_t it = strV.find(".");
      if (it==std::string::npos) it = strV.length();
      if (it < 2) {
        it = 2-it;
        if (strV.length()>tufW) {
          it = 0;
        } else if (strV.length()+it>tufW) {
          it = tufW-strV.length();
        }
        strV.insert(0,it,' ');
      }
      // align the e's
        it = strV.find(".");
        if (it!=std::string::npos) {
          std::size_t it2 = strV.find("e");
          if (it2!=std::string::npos) {
            if (it2>=it+1) {
              it = it2-(it+1);        // number of digits after the dot
              if (it<precTU) {
                it = precTU - it;
                strV.insert(it2,it,'0');
              }
            }
          }
        }
      // append after the number
        if (tufW>strV.length()) {
          it = tufW - strV.length();
          strV.insert(strV.length(),it,' ');
        }
    }
  return strV;
}

void strGlobalVar::Double2String_setPrec(const tuint prec)
{
  D2S_Prec = prec;
}

void strGlobalVar::Double2String_setType(const tuint type)
{
  D2S_Type = type;
}

void strGlobalVar::Double2String_setBValU(const tdouble bv)
{
  D2S_BValU = fabs(bv);
}

void strGlobalVar::Double2String_setBValL(const tdouble bv)
{
  D2S_BValL = fabs(bv);
}

void strGlobalVar::Double2String_setDel0(const bool d0)
{
  D2S_Del0 = d0;
}

void strGlobalVar::Double2String_setDelP(const bool dP)
{
  D2S_DelP = dP;
}

void strGlobalVar::Double2String_log()
{
  slog(4) << "  floating point conversion:" << std::endl;
  slog(4) << "        Precision:              " << D2S_Prec << std::endl;
  slog(4) << "        Type:                   " << D2S_Type;
  if (D2S_Type==2) {
    slog(4) << " [" << D2S_BValL << "; " << D2S_BValU << "]";
  }
  slog(4) << std::endl;
  slog(4) << "        Del0:                   ";
  if (D2S_Del0) {
    slog(4) << "yes";
  } else {
    slog(4) << "no";
  }
  slog(4) << std::endl;
  slog(4) << "        DelP:                   ";
  if (D2S_DelP) { 
    slog(4) << "yes";
  } else {
    slog(4) << "no";
  }
  slog(4) << std::endl;
}

const boost::basic_format<char> strGlobalVar::D2S_totalPrec(const tdouble dv)
{
  return (format("%19.12e") % dv);
}

std::ostream& strGlobalVar::slog(const int logLevel_)
{
  if (slogP == NULL || logLevel_ > logLevel) return *sdummy;
  else return *slogP; 
}

std::ostream& strGlobalVar::slogcout(const int logLevel_)
{
  if (slogcoutP == NULL || logLevel_ > logLevel) return *sdummy;
  else return *slogcoutP; 
}

void strGlobalVar::set_TOL(const tdouble tolV)
{
  TOLval = fabs(tolV);
  if (TOLval > 1e-6) {
    slog(2) << std::endl << "ALERT: TOL (" << TOLval << ") is very large. This leads to wrong results!" << std::endl << std::endl;
  }
}

strGlobalVar::strGlobalVar() 
: true_stdcout(&std::cout),true_cerr(&std::cerr),sdummy(new flxDummyOstream),slogP(NULL),stdcout(true_stdcout),stdcerr(true_cerr),slogcoutP(NULL),
  // NOTE this are not the default values - the default values are set in the settings of the configuration file
  D2S_Prec(DEFAULT_FPC_PREC), D2S_Type(DEFAULT_FPC_TYPE), D2S_BValU(DEFAULT_FPC_BVALU), D2S_BValL(DEFAULT_FPC_BVALL), D2S_Del0(DEFAULT_FPC_DEL0), D2S_DelP(DEFAULT_FPC_DELP),
  TOLval(DEFAULT_TOL),prelog_active(true),prelog_stream(""),
  logLevel_strong(DEFAULT_LOG_LEVEL),leak_check(false),logLev_counter(0),
  logLevel(DEFAULT_LOG_LEVEL),
  MT19937_init_calls(DEFAULT_MT19937_INIT_CALLS), MT19937_init_RAND(DEFAULT_MT19937_INIT_RAND), MT19937_init_seed(DEFAULT_MT19937_INIT_SEED), MT19937_init_seedvalue(DEFAULT_MT19937_INIT_SEEDVALUE),
  LegendrePolyH_init_numb(DEFAULT_LEGENDRE_NUMB), allowInteractiveMode(true), prgBar(DEFAULT_FLX_PRGBAR), configDir(""),configDir_local(""), pwd(""),
  sqrtEps(std::sqrt(std::numeric_limits<tdouble>::epsilon())), cmd1st(false), max_parallel_threads(std::thread::hardware_concurrency())
{
  set_slogcout(true_stdcout,true_stdcout);
  if (max_parallel_threads<=1) {
    max_parallel_threads = 2;
  }
}

strGlobalVar::~strGlobalVar()
{
  using namespace std;
  if (sdummy) delete sdummy;
  // delete slogcoutP
    if (slogcoutP) {
      flxStreamAlloc* fsa = dynamic_cast<flxStreamAlloc*> (slogcoutP);
      if (fsa) delete fsa;
    }
  // delete log-stream
    if (slogP) {
      ofstream* thestream = dynamic_cast<ofstream*> (slogP);
      if ( thestream ) {
        thestream->close();
        delete thestream;
      }
    }
}

void strGlobalVar::set_slogcout(ostreamp slogV, ostreamp stdcoutV)
{
  const bool b_same_slog = (slogP==slogV);
  // update cout
    const bool b_same_cout = (stdcout==stdcoutV);
    if (b_same_slog && b_same_cout) return;
    if (!b_same_cout) {
      const ostreamp old_cout = stdcout;
      stdcout = stdcoutV;
      if (slogP==old_cout) slogP = stdcout;
      if (slogcoutP==old_cout) slogcoutP = stdcout;
    }
  // update slog
    if (!b_same_slog) {
      slogP = slogV;
    }
  // try to delete slogcoutP
    flxStreamAlloc* fsa = dynamic_cast<flxStreamAlloc*> (slogcoutP);
    if (fsa) delete fsa;
    if (slogP==NULL) slogcoutP = NULL;
  // reallocate slogcoutP
    if (slogP==stdcout) {
      slogcoutP = slogP;
    } else {
      slogcoutP = new flxStreamAlloc(stdcout, slogP);
    }
}

void strGlobalVar::logLevel_strong_set(const int logL)
{
  logLevel_strong = logL; 
  logLevel_strong_reset();
}

void strGlobalVar::logLevel_log_deactivate(const bool deact)
{
  if (deact) {                // deactivates logging
    ++logLev_counter;
    if (logLev_counter>1) return;        // log is already deactivated
    if (GlobalVar.logLevel>2) GlobalVar.logLevel=2;
  } else {                // activates logging
    if (logLev_counter==0) {
      throw FlxException_Crude("strGlobalVar::logLevel_log_deactivate");
    }
    --logLev_counter;
    if (logLev_counter==0) {
      logLevel_strong_reset();
    }
  }
}

void strGlobalVar::prelog_activated(const bool b)
{
  if (prelog_active==b) return;
  prelog_active=b;
  prelog_stream = "";
}

void strGlobalVar::prelog_append(char c)
{
  prelog_stream+=c;
}

void strGlobalVar::prelog_append(std::string str)
{
  prelog_stream+=str;
}

void strGlobalVar::prelog_prepare_write()
{
  bool b1=true;
  size_t pos; char c;
  do {
    pos = prelog_stream.size();
    if (prelog_stream == "" || pos <= 0 || prelog_stream.empty() ) return;
    --pos;
    c = prelog_stream[pos];
    if (c == 10 || c == 9 || c == 13 || c == ' ' ) {
      prelog_stream.erase(pos,1);
    } else {
      b1 = false;
    }
  } while (b1);
  b1 = true;
  do {
    if (prelog_stream == "" || prelog_stream.size() <= 0 || prelog_stream.empty() ) return;
    c = prelog_stream[0];
    if (c == 10 || c == 9 || c == 13 || c == ' ' ) {
      prelog_stream.erase(0,1);
    } else {
      b1 = false;
    }
  } while (b1);
}

void strGlobalVar::prelog_clear()
{
  prelog_stream = "";
}

std::string strGlobalVar::prelog_force_write()
{
  prelog_prepare_write();
  std::string strtmp = prelog_stream;
//   prelog_stream = "";
  return strtmp;
}

void strGlobalVar::prelog_write()
{
  if (prelog_active) {
    if ( prelog_isNOTempty() ) {
      prelog_prepare_write();
      slog(4) << prelog_stream << std::endl;
      prelog_stream = "";
    }
  } else {
    prelog_stream = "";
  }
}

void FlxProgress::start(const tuint N_v) {
  if (!active) return;
  N = N_v;
  last_t = clock();
  op << "0%"; op.flush(); op << "\b\b";
  last_p = 0;
  is_running = true;
}

void FlxProgress::tick_t(const tuint idx) {
  const tuint p = int((idx*100)/N);
  if (p==last_p) return;
  op << p << "%"; op.flush();
  op << "\b\b";
  if (p>=10) {
    op << '\b';
    if (p>=100) {
      op << '\b';
    }
  }
  last_p = p;
  last_t = clock();
}

void FlxProgress::stop() {
  if (!active) return;
  op << "    "; op.flush(); op << "\b\b\b\b";
  is_running = false;
}


const std::string bool2string(const bool b)
{
  if (b) return "yes";
  else return "no";
}

void fesslix_logInfo(std::ostream& lout)
{
  lout << " FlxLib: " << std::endl;
  lout << "   version: " << FLX_VERSION << std::endl;
  lout << "   compiled with the options ..." << std::endl;
  // FLXDEBUG
    lout << "     FLX_DEBUG                     ";
    #if FLX_DEBUG
      lout << "ON";
    #else 
      lout << "OFF";
    #endif
    lout << std::endl;
  // FLXDEBUG_COUT
    lout << "     FLX_DEBUG_COUT                ";
    #if FLX_DEBUG_COUT
      lout << "ON";
    #else 
      lout << "OFF";
    #endif
    lout << std::endl;
  // FLX_BOOST_FS
    lout << "     FLX_BOOST_FS                  ";
    #if FLX_BOOST_FS
      lout << "ON";
    #else 
      lout << "OFF";
    #endif  
    lout << std::endl;
  // FLX_BOOST_RX
    lout << "     FLX_BOOST_RX                  ";
    #if FLX_BOOST_RX
      lout << "ON";
    #else 
      lout << "OFF";
    #endif  
    lout << std::endl;
  // FLX_KAHAN_MCI
    lout << "     FLX_KAHAN_MCI                 ";
    #if FLX_KAHAN_MCI
      lout << "ON";
    #else 
      lout << "OFF";
    #endif  
    lout << std::endl;
  // FLX_KAHAN_KDE
    lout << "     FLX_KAHAN_KDE                 ";
    #if FLX_KAHAN_KDE
      lout << "ON";
    #else 
      lout << "OFF";
    #endif  
    lout << std::endl;
  // FLX_KAHAN_DOT
    lout << "     FLX_KAHAN_DOT                 ";
    #if FLX_KAHAN_DOT
      lout << "ON";
    #else 
      lout << "OFF";
    #endif  
    lout << std::endl;
  // FLX_KAHAN_2NORM
    lout << "     FLX_KAHAN_2NORM               ";
    #if FLX_KAHAN_2NORM
      lout << "ON";
    #else 
      lout << "OFF";
    #endif  
    lout << std::endl;
  // FLX_KAHAN_MTX_FULL
    lout << "     FLX_KAHAN_MTX_FULL            ";
    #if FLX_KAHAN_MTX_FULL
      lout << "ON";
    #else 
      lout << "OFF";
    #endif  
    lout << std::endl;
  // FLX_KAHAN_MTX_LTRI
    lout << "     FLX_KAHAN_MTX_LTRI            ";
    #if FLX_KAHAN_MTX_LTRI
      lout << "ON";
    #else 
      lout << "OFF";
    #endif  
    lout << std::endl;
  // FLX_KAHAN_MTX_SLTRI
    lout << "     FLX_KAHAN_MTX_SLTRI           ";
    #if FLX_KAHAN_MTX_SLTRI
      lout << "ON";
    #else 
      lout << "OFF";
    #endif  
    lout << std::endl;
  // FLX_KAHAN_MTX_SSILU
    lout << "     FLX_KAHAN_MTX_SSILU           ";
    #if FLX_KAHAN_MTX_SSILU
      lout << "ON";
    #else 
      lout << "OFF";
    #endif  
    lout << std::endl;
  // FLX_KAHAN_MTX_SYM
    lout << "     FLX_KAHAN_MTX_SYM             ";
    #if FLX_KAHAN_MTX_SYM
      lout << "ON";
    #else 
      lout << "OFF";
    #endif  
    lout << std::endl;
  // FLX_KAHAN_MTX_SSYM
    lout << "     FLX_KAHAN_MTX_SSYM            ";
    #if FLX_KAHAN_MTX_SSYM
      lout << "ON";
    #else 
      lout << "OFF";
    #endif  
    lout << std::endl;
  // FLX_KAHAN_MTX_SSFEM
    lout << "     FLX_KAHAN_MTX_SSFEM           ";
    #if FLX_KAHAN_MTX_SSFEM
      lout << "ON";
    #else 
      lout << "OFF";
    #endif  
    lout << std::endl;
  // FLX_KAHAN_MTX_PFEM
    lout << "     FLX_KAHAN_MTX_PFEM            ";
    #if FLX_KAHAN_MTX_PFEM
      lout << "ON";
    #else 
      lout << "OFF";
    #endif
    lout << std::endl;
  // FLX_KAHAN_CG_STAGE
    lout << "     FLX_KAHAN_CG_STAGE            ";
    lout << FLX_KAHAN_CG_STAGE << std::endl;    
}

