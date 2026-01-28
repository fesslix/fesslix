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

#define FLXLIB_CPP

#include "flxfunction_fun_read.h"
#include "flxfunction_fun.h"

using namespace std;

void flxfunction_fun_insert(FlxFunctionBox& FunBox)
{
  FunBox.insert("sin", new FunReadFunSin() );
  FunBox.insert("arcsin", new FunReadFunArcsin() );
  FunBox.insert("cos", new FunReadFunCos() );
  FunBox.insert("arccos", new FunReadFunArccos() );
  FunBox.insert("tan", new FunReadFunTan() );
  FunBox.insert("arctan", new FunReadFunArctan() );
  FunBox.insert("ln", new FunReadFunLn() );
  FunBox.insert("log", new FunReadFunLog() );
  FunBox.insert("loga", new FunReadFunLoga() );
  FunBox.insert("exp", new FunReadFunExp() );
  FunBox.insert("cot", new FunReadFunCot() );
  FunBox.insert("arccot", new FunReadFunArccot() );
  FunBox.insert("sinh", new FunReadFunSinh() );
  FunBox.insert("cosh", new FunReadFunCosh() );
  FunBox.insert("tanh", new FunReadFunTanh() );
  FunBox.insert("arctanh", new FunReadFunArctanh() );
  FunBox.insert("sqrt", new FunReadFunSqrt() );
  FunBox.insert("abs", new FunReadFunAbs() );
  FunBox.insert("frac", new FunReadFunFrac() );
  FunBox.insert("round", new FunReadFunRound() );
  FunBox.insert("mod", new FunReadFunMod() );
  FunBox.insert("sig", new FunReadFunSig() );
  FunBox.insert("factorial", new FunReadFunFactorial() );
  FunBox.insert("factorialln", new FunReadFunFactorialLn() );
  FunBox.insert("pdfn", new FunReadFunPDFn() );
  FunBox.insert("pdfn_ln", new FunReadFunPDFn_ln() );
  FunBox.insert("cdfn", new FunReadFunCDFn() );
  FunBox.insert("cdfn_inv", new FunReadFunInvCDFn() );
  FunBox.insert("relidx", new FunReadFunRelIdx() );
  FunBox.insert("faillsf", new FunReadFunFailLSF() );
  FunBox.insert("hvi", new FunReadFunHvi() );
  FunBox.insert("time0", new FunReadFunTime0() );
  FunBox.insert("integ", new FunReadFunInteg() );
  FunBox.insert("sum", new FunReadFunSum() );
  FunBox.insert("optimize1d", new FunReadFunOptimize1D() );
  FunBox.insert("deg2gauss", new FunReadFundeg2gauss() );
  FunBox.insert("binomialcoeff", new FunReadFunBinomialCoeff() );
  FunBox.insert("binomialcoeff_ln", new FunReadFunBinomialCoeff_ln() );
  FunBox.insert("pdfn2", new FunReadFunPDFn2() );
  FunBox.insert("pdfn2_ln", new FunReadFunPDFn2_ln() );
  FunBox.insert("pmf_beta_binomial_ln", new FunReadFunPMF_beta_binomial_ln() );
  FunBox.insert("evalr", new FunReadFunEvalR() );
  FunBox.insert("evalc", new FunReadFunEvalC() );
  FunBox.insert("evalw", new FunReadFunEvalW() );
  FunBox.insert("autocorr_exp", new FunReadFunAutoCorrExp() );
  FunBox.insert("autocorr_exp2", new FunReadFunAutoCorrExp2() );
  FunBox.insert("tgamma", new FunReadFunGamma() );
  FunBox.insert("lngamma", new FunReadFunLnGamma() );
  FunBox.insert("igamma", new FunReadFunIGamma() );
  FunBox.insert("igammal", new FunReadFunIGammaL() );
  FunBox.insert("rgamma", new FunReadFunRGamma() );
  FunBox.insert("rgammal", new FunReadFunRGammaL() );
  FunBox.insert("rgamma_inv", new FunReadFunRGamma_inv() );
  FunBox.insert("rgammal_inv", new FunReadFunRGammaL_inv() );
  FunBox.insert("ibeta", new FunReadFunIBeta() );
  FunBox.insert("ibeta_inv", new FunReadFunIBeta_inv() );
  FunBox.insert("betafun", new FunReadFunBeta() );
  FunBox.insert("lnbetafun", new FunReadFunLnBeta() );
  FunBox.insert("erf", new FunReadFunErf() );
  FunBox.insert("erf_inv", new FunReadFunErf_inv() );
  FunBox.insert("error", new FunReadFunError() );
  FunBox.insert("isnan", new FunReadFunIsNaN() );
}


FunBase* FunReadFunPDFn2::read( bool errSerious )
{
  std::vector<FunBase*>* tv = read_parameters(-1,errSerious);
  if (tv->size()==3) {
    return new FunPDFn2( tv );
  } else if (tv->size()==7) {
    return new FunPDFn2_general( tv );
  } else {
    ostringstream ssV;
    ssV << "'pdfn2' expects either 3 or 7 parameters - and not " << tv->size() << ".";
    throw FlxException("FunReadFunPDFn2::read", ssV.str());
  }
}

FunBase* FunReadFunPDFn2_ln::read( bool errSerious )
{
  std::vector<FunBase*>* tv = read_parameters(-1,errSerious);
  if (tv->size()==3) {
    return new FunPDFn2_ln( tv );
  } else if (tv->size()==7) {
    return new FunPDFn2_ln_general( tv );
  } else {
    ostringstream ssV;
    ssV << "'pdfn2_ln' expects either 3 or 7 parameters - and not " << tv->size() << ".";
    throw FlxException("FunReadFunPDFn2_ln::read", ssV.str());
  }
}

FunBase* FunReadFunPMF_beta_binomial_ln::read( bool errSerious )
{
  std::vector<FunBase*>* tv = read_parameters(4,errSerious);
  return new FunPMF_beta_binomial_ln( tv );
}

FunBase* FunReadFunPDFn::read( bool errSerious )
{
  std::vector<FunBase*>* tv = read_parameters(-1,errSerious);
  if (tv->size()==1) {
    return new FunPDFn( tv );
  } else if (tv->size()==3) {
    return new FunPDFn_general( tv );
  } else {
    ostringstream ssV;
    ssV << "'pdfn' expects either 1 or 3 parameters - and not " << tv->size() << ".";
    throw FlxException("FunReadFunPDFn::read", ssV.str());
  }
}

FunBase* FunReadFunPDFn_ln::read( bool errSerious )
{
  std::vector<FunBase*>* tv = read_parameters(-1,errSerious);
  if (tv->size()==1) {
    return new FunPDFn_ln( tv );
  } else if (tv->size()==3) {
    return new FunPDFn_ln_general( tv );
  } else {
    ostringstream ssV;
    ssV << "'pdfn_ln' expects either 1 or 3 parameters - and not " << tv->size() << ".";
    throw FlxException("FunReadFunPDFn_ln::read", ssV.str());
  }
}

FunBase* FunReadFunCDFn::read( bool errSerious )
{
  std::vector<FunBase*>* tv = read_parameters(-1,errSerious);
  if (tv->size()==1) {
    return new FunCDFn( tv );
  } else if (tv->size()==3) {
    return new FunCDFn_general( tv );
  } else {
    ostringstream ssV;
    ssV << "'cdfn' expects either 1 or 3 parameters - and not " << tv->size() << ".";
    throw FlxException("FunReadFunCDFn::read", ssV.str());
  }
}

const tdouble FunLoga::calc()
{
  const tdouble x = ParaListP[0]->calc();
  const tdouble base = ParaListP[1]->calc();
  if (base<=ZERO || base==ONE) {
    std::ostringstream ssV;
    ssV << "The base of a logarithm must be a positive number not equal to one (and not " << GlobalVar.Double2String(base) << ").";
    throw FlxException("FunLoga::calc", ssV.str() );
  }
  return std::log10(x)/std::log10(base);
}

const tdouble FunRound::calc()
{
  const tuint Narg = ParaList->size();
  if (Narg==1) {
    return round_flx( ParaListP[0]->calc() );
  } else if (Narg==2) {
    return round_flx( ParaListP[0]->calc(), tuint_from(ParaListP[1]->calc(),"Number of parameters") );
  } else {
    std::ostringstream ssV;
    ssV << "Invalid number of arguments: " << ParaList->size() << ".";
    throw FlxException("FunRound::calc", ssV.str() );
  }
}

FunInteg::~FunInteg()
{
  delete funI;
  delete startF;
  delete endF;
  delete gaussF;
  if (intF) delete intF;
}

const bool FunInteg::search_circref(FlxFunction* fcr)
{
  return funI->search_circref(fcr) || startF->search_circref(fcr) || endF->search_circref(fcr) || gaussF->search_circref(fcr) || intF->search_circref(fcr);
}

const bool FunInteg::dependOn_Const(const tdouble*const thenumber)
{
  if (startF->dependOn_Const(thenumber)) return true;
  if (endF->dependOn_Const(thenumber)) return true;
  if (gaussF->dependOn_Const(thenumber)) return true;
  if (intF->dependOn_Const(thenumber)) return true;
  if (theconst==thenumber) { 
    return false;
  } else {
    return funI->dependOn_Const(thenumber);
  }
}

const bool FunInteg::optimize(FunBasePtr& optf, const Fun_OptimizeInfo &foi) {
  child_optimize(funI, foi);
  child_optimize(startF, foi);
  child_optimize(endF, foi);
  child_optimize(gaussF, foi);
  if (intF) child_optimize(intF, foi);
  return false;
}

const std::string FunInteg::write()
{
  string str1 = "integ(";
  str1+=ConstantBox->get(theconst);
  str1+="=[";
  str1+=startF->write();
  str1+=",";
  str1+=endF->write();
  str1+="],";
  str1+=funI->write();
  str1+=",gp=";
  str1+=gaussF->write();
  if (intF!=NULL) {
    str1+=",int=";
    str1+=intF->write();
  }
  str1+=")";
  return str1;
}

const tdouble FunInteg::calc()
{
  // get number of integration points
    const tuint N = tnlong_from(gaussF->calc(),"number of Gauss-points",true,false,gaussF);
  // determine start and end
    const tdouble start = startF->calc();
    const tdouble end = endF->calc();
    if (end<=start) {
      if (fabs(start-end)<=GlobalVar.TOL()) return ZERO;
      std::ostringstream ssV;
      ssV << "Error with integration boundaries.";
      throw FlxException("FunInteg::calc_01", ssV.str() );
    }
  if (intF==NULL) {
    return FlxFun_GaussIntegration(funI,theconst,start,end,N,*GI);
  } else {
    // get the number of intervals
    const tuint I = tnlong_from(intF->calc(),"the number of intervals",true,false,intF);
    if (is_log) {
      if (start<=ZERO) {
        std::ostringstream ssV;
        ssV << "Log-scaled interval spacing is only allowed for positive intervals.";
        throw FlxException("FunInteg::calc_02", ssV.str() );
      }
      const tdouble dx = (log(end) - log(start))/I;
      tdouble res = ZERO;
      tdouble t = log(start);
      for (tuint i=0;i<I;++i) {
        res+=FlxFun_GaussIntegration(funI,theconst,exp(t),exp(t+dx),N,*GI);
        t+=dx;
      }
      return res;
    } else {
      const tdouble dx = (end - start)/I;
      tdouble res = ZERO;
      tdouble t = start;
      for (tuint i=0;i<I;++i) {
        res+=FlxFun_GaussIntegration(funI,theconst,t,t+dx,N,*GI);
        t+=dx;
      }
      return res;
    }
  }
}


FunBase* FunReadFunInteg::read(bool errSerious)
{
  FunBase* funI = NULL;
  FunBase *startF=NULL,*endF=NULL,*gaussF=NULL,*intF=NULL;
  bool is_log = false;
  try {
    // get the constant
      tdouble* d1 = read_const_var(errSerious,true);
    // integration bounds
      reader->getChar('=',errSerious);
      reader->getChar('[',errSerious);
      startF = FunctionList->read(errSerious);
      reader->getChar(',',errSerious);
      endF = FunctionList->read(errSerious);
      reader->getChar(']',errSerious);
      reader->getChar(',',errSerious);
    // function to integrate
      funI = FunctionList->read(errSerious);
    reader->getChar(',',errSerious);
    reader->getWord("gp",errSerious);
    reader->getChar('=',errSerious);
    gaussF = FunctionList->read(errSerious);
    intF = NULL;
    if (reader->whatIsNextChar()==',') {
      reader->getChar(',',errSerious);
      reader->getWord("int",errSerious);
      reader->getChar('=',errSerious);
      intF = FunctionList->read(errSerious);
      if (reader->whatIsNextChar()==',') {
        reader->getChar(',',errSerious);
        reader->getWord("log",errSerious);
        is_log = true;
      }
    }
    return new FunInteg(funI,d1,startF,endF,gaussF,intF,is_log);  
  } catch (FlxException &e) {
    FLXMSG("FunReadFunInteg::read",1);
    if (funI) delete funI;
    if (startF!=NULL) delete startF;
    if (endF!=NULL) delete endF;
    if (gaussF!=NULL) delete gaussF;
    if (intF!=NULL) delete intF;
    throw;
  }
}

FunSum::~FunSum()
{
  delete funS;
  delete startF;
  delete cond;
  delete step;
}

const tdouble FunSum::calc()
{
  const tdouble STORE_VAL = *theconst;
  tdouble res = 0.0;
  *theconst = startF->calc();
  while (cond->calc()>0.0) {
    res += funS->calc();
    *theconst = step->calc();
  }
  *theconst = STORE_VAL;
  return res;
}

const bool FunSum::search_circref(FlxFunction* fcr)
{
  return funS->search_circref(fcr) || startF->search_circref(fcr) || cond->search_circref(fcr) || step->search_circref(fcr);
}

const bool FunSum::dependOn_Const(const tdouble*const thenumber)
{
  if (startF->dependOn_Const(thenumber)) return true;
  if (cond->dependOn_Const(thenumber)) return true;
  if (step->dependOn_Const(thenumber)) return true;
  if (theconst==thenumber) { 
    return false;
  } else {
    return funS->dependOn_Const(thenumber);
  }
}

const bool FunSum::optimize(FunBasePtr& optf, const Fun_OptimizeInfo &foi) {
  child_optimize(funS, foi);
  child_optimize(startF, foi);
  child_optimize(cond, foi);
  child_optimize(step, foi);
  return false;
}

const std::string FunSum::write()
{
  string str1 = "sum(";
  str1+=ConstantBox->get(theconst);
  str1+="=";
  str1+=startF->write();
  str1+=",";
  str1+=cond->write();
  str1+=",";
  str1+=step->write();
  str1+=",";
  str1+=funS->write();
  str1+=")";
  return str1;
}

FunBase* FunReadFunSum::read(bool errSerious)
{
  FunBase* funS = NULL;
  FunBase *startF=NULL,*cond=NULL,*step=NULL;
  try {
    // get the constant
      tdouble* d1 = read_const_var(errSerious,true);
    // start-value
      reader->getChar('=',errSerious);
      startF = FunctionList->read(errSerious);
      reader->getChar(',',errSerious);
    // condition
      cond = FunctionList->read(errSerious);
      reader->getChar(',',errSerious);
    // step-expression
      step = FunctionList->read(errSerious);
      reader->getChar(',',errSerious);
    // function to integrate
      funS = FunctionList->read(errSerious);
    return new FunSum(funS,d1,startF,cond,step);  
  } catch (FlxException &e) {
    FLXMSG("FunReadFunSum::read",1);
    if (startF) delete startF;
    if (cond) delete cond;
    if (step) delete step;
    if (funS) delete funS;
    throw;
  }
}

FunRoot::~FunRoot()
{
  delete funR;
  delete startF;
  delete endF;
  delete dx;
  delete dy;
}

const bool FunRoot::search_circref(FlxFunction* fcr)
{
  return funR->search_circref(fcr) || startF->search_circref(fcr) || endF->search_circref(fcr) || dx->search_circref(fcr) || dy->search_circref(fcr);
}

const bool FunRoot::dependOn_Const(const tdouble*const thenumber)
{
  if (startF->dependOn_Const(thenumber)) return true;
  if (endF->dependOn_Const(thenumber)) return true;
  if (dx->dependOn_Const(thenumber)) return true;
  if (dy->dependOn_Const(thenumber)) return true;
  if (theconst==thenumber) { 
    return false;
  } else {
    return funR->dependOn_Const(thenumber);
  }
}

const bool FunRoot::optimize(FunBasePtr& optf, const Fun_OptimizeInfo &foi) {
  child_optimize(funR, foi);
  child_optimize(startF, foi);
  child_optimize(endF, foi);
  child_optimize(dx, foi);
  child_optimize(dy, foi);
  return false;
}

const std::string FunRoot::write()
{
  string str1 = "root(";
  str1+=ConstantBox->get(theconst);
  str1+="=[";
  str1+=startF->write();
  str1+=",";
  str1+=endF->write();
  str1+="],";
  str1+=funR->write();
  str1+=",";
  switch (ID) {
    case 0:
      str1+="bisec";
      break;
    case 1:
      str1+="rgfsi";
      break;
    default:
      throw FlxException_Crude("FunRoot::write");
  }
  str1+=",dy=";
  str1+=dy->write();
  str1+=",dx=";
  str1+=dx->write();
  str1+=")";
  return str1;
}

const tdouble FunRoot::calc()
{
  switch (ID) {
    case 0:
      return FlxFun_RootSearch_Bisec(funR,theconst,startF->calc(),endF->calc(),dx->calc(),dy->calc(),streamp);
    case 1:
      return FlxFun_RootSearch_RegulaFalsi(funR,theconst,startF->calc(),endF->calc(),dx->calc(),dy->calc(),streamp);
    default:
      throw FlxException_Crude("FunRoot::calc");
  }
}

FunOptimize1D::~FunOptimize1D()
{
  delete funR;
  delete startF;
  delete endF;
  delete eps1;
  delete eps2;
  delete niter;
  delete nex;
}

const bool FunOptimize1D::search_circref(FlxFunction* fcr)
{
  return funR->search_circref(fcr) || startF->search_circref(fcr) || endF->search_circref(fcr) || eps1->search_circref(fcr) || eps2->search_circref(fcr) || niter->search_circref(fcr) || nex->search_circref(fcr);
}

const bool FunOptimize1D::dependOn_Const(const tdouble*const thenumber)
{
  if (startF->dependOn_Const(thenumber)) return true;
  if (endF->dependOn_Const(thenumber)) return true;
  if (eps1->dependOn_Const(thenumber)) return true;
  if (eps2->dependOn_Const(thenumber)) return true;
  if (niter->dependOn_Const(thenumber)) return true;
  if (nex->dependOn_Const(thenumber)) return true;
  if (theconst==thenumber) { 
    return false;
  } else {
    return funR->dependOn_Const(thenumber);
  }
}

const bool FunOptimize1D::optimize(FunBasePtr& optf, const Fun_OptimizeInfo &foi) {
  child_optimize(funR, foi);
  child_optimize(startF, foi);
  child_optimize(endF, foi);
  child_optimize(eps1, foi);
  child_optimize(eps2, foi);
  child_optimize(niter, foi);
  child_optimize(nex, foi);
  return false;
}

const std::string FunOptimize1D::write()
{
  string str1 = "optimize1D(";
  str1+=ConstantBox->get(theconst);
  str1+=",[";
  str1+=startF->write();
  str1+=",";
  str1+=endF->write();
  str1+="],";
  str1+=funR->write();
  if (eps1) str1+=",eps1="+eps1->write();
  if (eps2) str1+=",eps2="+eps2->write();
  if (niter) str1+=",niter="+niter->write();
  if (nex) str1+=",nex="+nex->write();
  str1+=",method=";
  str1+= (use_brent?"brent":"golden");
  str1+=")";
  return str1;
}

const tdouble FunOptimize1D::calc()
{
  const tdouble storedF = *theconst;
  flx_fun_void_data ds;
    ds.theconst = theconst;
    ds.funR = funR;
  const tuint NiterV = tuint_from(niter->calc(),"niter",true,false,niter);  
  const tuint NexV = tuint_from(nex->calc(),"nex",true,false,nex);
  const tdouble eps1V = eps1->calc();
  const tdouble eps2V = eps2->calc();
  tdouble res = ZERO;
  flx_optim(startF->calc(),endF->calc(),res,&flx_fun_void,&ds,use_brent,false,NiterV,NexV,eps1V,eps2V);
  *theconst = storedF;
  return res;
}

FunBase* FunReadFunOptimize1D::read(bool errSerious)
{
  FunBase *startF=NULL,*endF=NULL,*eps1=NULL,*eps2=NULL,*niter=NULL,*nex=NULL;
  // get the constant
    tdouble* d1 = read_const_var(errSerious,true);
  try {
      reader->getChar(',',errSerious);
        reader->getChar('[',errSerious);
          startF = FunctionList->read(errSerious);
          reader->getChar(',',errSerious);
          endF = FunctionList->read(errSerious);
        reader->getChar(']',errSerious);
      reader->getChar(',',errSerious);
        FunBase* funR = FunctionList->read(errSerious);
      bool use_brent = true;
      while (reader->whatIsNextChar()==',') {
        reader->getChar(',',false);
        const std::string pid = reader->getWord(true,false);
        reader->getChar('=',false);
        if (pid=="eps1") {
          if (eps1) throw FlxException("FunReadFunOptimize1D::read_01","Parameter '"+pid+"' was already specified");
          eps1 = FunctionList->read(errSerious);
        } else if (pid=="eps2") {
          if (eps2) throw FlxException("FunReadFunOptimize1D::read_02","Parameter '"+pid+"' was already specified");
          eps2 = FunctionList->read(errSerious);
        } else if (pid=="niter") {
          if (niter) throw FlxException("FunReadFunOptimize1D::read_03","Parameter '"+pid+"' was already specified");
          niter = FunctionList->read(errSerious);
        } else if (pid=="nex") {
          if (nex) throw FlxException("FunReadFunOptimize1D::read_04","Parameter '"+pid+"' was already specified");
          nex = FunctionList->read(errSerious);
        } else if (pid=="meth") {
          const std::string methstr = reader->getWord(true,false);
          if (methstr=="brent") use_brent = true;
          else if (methstr=="golden") use_brent = false;
          else {
            throw FlxException("FunReadFunOptimize1D::read_05","Unknown method identifier name '"+methstr+"'.");
          }
        } else {
          throw FlxException("FunReadFunOptimize1D::read_06","Unknown parameter name '"+pid+"'.");
        }
      }
      if (eps1==NULL) eps1 = new FunNumber(1e-6);
      if (eps2==NULL) eps2 = new FunNumber(1e-6);
      if (niter==NULL) niter = new FunNumber(1e3);
      if (nex==NULL) nex = new FunNumber(1e2);
      return new FunOptimize1D(funR,d1,startF,endF,eps1,eps2,niter,nex,use_brent);  
  } catch (FlxException &e) {
    FLXMSG("FunReadFunOptimize1D::read",1);
    if (startF) delete startF;
    if (endF) delete endF;
    if (eps1) delete eps1;
    if (eps2) delete eps2;
    if (niter) delete niter;
    if (nex) delete nex;
    throw;
  }
}

FunBase* FunReadFunEvalR::read(bool errSerious)
{
  FunBase* funI = FunctionList->read(errSerious);
  const tdouble cv = funI->calc();
  delete funI;
  return new FunNumber(cv);
}

FunBase* FunReadFunEvalC::read(bool errSerious)
{
  FunBase* funI = FunctionList->read(errSerious);
  Fun_OptimizeInfo foi(true);
  FunBase::child_optimize(funI,foi);
  return funI;
}

FunBase* FunReadFunAutoCorrExp::read(bool errSerious)
{
  tdouble* d1 = read_const_var(errSerious);
  reader->getChar(',');
  tdouble* d2 = read_const_var(errSerious);
  return new FunAutoCorrExp(d1,d2);
}

FunBase* FunReadFunAutoCorrExp2::read(bool errSerious)
{
  tdouble* d1 = read_const_var(errSerious);
  reader->getChar(',');
  tdouble* d2 = read_const_var(errSerious);
  return new FunAutoCorrExp2(d1,d2);
}

const tdouble FunPDFn_general::calc()
{
  const tdouble x = ParaListP[0]->calc();
  const tdouble m = ParaListP[1]->calc();
  const tdouble s = ONE/(ParaListP[2]->calc());
  return  rv_phi( (x-m)*s )*s; 
}

const tdouble FunPDFn_ln_general::calc()
{
  const tdouble x = ParaListP[0]->calc();
  const tdouble m = ParaListP[1]->calc();
  const tdouble s = ONE/(ParaListP[2]->calc());
  return -pow2((x-m)*s)/2-log(sqrt(2*PI)/s);
}

const tdouble FunCDFn_general::calc()
{
  const tdouble x = ParaListP[0]->calc();
  const tdouble m = ParaListP[1]->calc();
  const tdouble s = ONE/(ParaListP[2]->calc());
  return  rv_Phi( (x-m)*s ); 
}

const tdouble FunDeg2Gauss::calc()
{
  return GI->pDegree2GPs(tuint_from(child_1->calc(),"Degree of a polynomial",true,true,child_1));
}

const tdouble FunBinomialCoeff::calc()
{
  const int i1 = int(round_flx(ParaListP[0]->calc()));
  const int i2 = int(round_flx(ParaListP[1]->calc()));
  return BinomCoeff(i1,i2);
}

const tdouble FunBinomialCoeff_ln::calc()
{
  const int i1 = int(round_flx(ParaListP[0]->calc()));
  const int i2 = int(round_flx(ParaListP[1]->calc()));
  return BinomCoeff_ln(i1,i2);
}

const tdouble FunPDFn2::calc()
{
  const tdouble u1 = ParaListP[0]->calc();
  const tdouble u2 = ParaListP[1]->calc();
  const tdouble r = ParaListP[2]->calc();
  if (r<=-ONE || r>=ONE) {
    std::ostringstream ssV;
    ssV << "The specified correlation (" << GlobalVar.Double2String(r) << ") must be within the interval ]-1;1[";
    throw FlxException("FunPDFn2::calc_1",ssV.str());
  }
  const tdouble t = ONE-pow2(r);
  const tdouble res = exp(-(pow2(u1)+pow2(u2)-2*r*u1*u2)/(2*t))/(2*PI*sqrt(t));
  #if FLX_DEBUG
    if (std::isnan(res)) throw FlxException_Crude("FunPDFn2::calc_2");
  #endif
  return res;
}

const tdouble FunPDFn2_general::calc()
{
  const tdouble x1 = ParaListP[0]->calc();
  const tdouble x2 = ParaListP[1]->calc();
  const tdouble r = ParaListP[2]->calc();
  const tdouble m1 = ParaListP[3]->calc();
  const tdouble m2 = ParaListP[4]->calc();
  const tdouble s1 = ParaListP[5]->calc();
  const tdouble s2 = ParaListP[6]->calc();
  const tdouble t = ONE-pow2(r);
  const tdouble u1 = (x1-m1)/s1;
  const tdouble u2 = (x2-m2)/s2;
  return exp(-(pow2(u1)+pow2(u2)-2*r*u1*u2)/(2*t))/(2*PI*s1*s2*sqrt(t));
}

const tdouble FunPDFn2_ln::calc()
{
  const tdouble u1 = ParaListP[0]->calc();
  const tdouble u2 = ParaListP[1]->calc();
  const tdouble r = ParaListP[2]->calc();
  if (r<=-ONE || r>=ONE) {
    std::ostringstream ssV;
    ssV << "The specified correlation (" << GlobalVar.Double2String(r) << ") must be within the interval ]-1;1[";
    throw FlxException("FunPDFn2_ln::calc_1",ssV.str());
  }
  const tdouble t = ONE-pow2(r);
  const tdouble res = -(pow2(u1)+pow2(u2)-2*r*u1*u2)/(2*t) - log(2*PI) - log(t)/2;
  #if FLX_DEBUG
    if (std::isnan(res)) throw FlxException_Crude("FunPDFn2_ln::calc_2");
  #endif
  return res;
}

const tdouble FunPDFn2_ln_general::calc()
{
  const tdouble x1 = ParaListP[0]->calc();
  const tdouble x2 = ParaListP[1]->calc();
  const tdouble r = ParaListP[2]->calc();
  const tdouble m1 = ParaListP[3]->calc();
  const tdouble m2 = ParaListP[4]->calc();
  const tdouble s1 = ParaListP[5]->calc();
  const tdouble s2 = ParaListP[6]->calc();
  const tdouble t = ONE-pow2(r);
  const tdouble u1 = (x1-m1)/s1;
  const tdouble u2 = (x2-m2)/s2;
  return -(pow2(u1)+pow2(u2)-2*r*u1*u2)/(2*t) - log(2*PI*s1*s2) - log(t)/2;
}

const tdouble FunPMF_beta_binomial_ln::calc()
{
  const tuint M = tuint_from(ParaListP[0]->calc(),"number of hits",true,true,ParaListP[0]);
  const tuint N = tuint_from(ParaListP[1]->calc(),"number of trials",true,false,ParaListP[1]);
  const tdouble alpha = ParaListP[2]->calc();
  const tdouble beta = ParaListP[3]->calc();
  if (alpha<=ZERO) {
    std::ostringstream ssV;
    ssV << "Parameter 'alpha' must be positive.";
    throw FlxException("FunPMF_beta_binomial_ln::calc_01", ssV.str() );
  }
  if (beta<=ZERO) {
    std::ostringstream ssV;
    ssV << "Parameter 'beta' must be positive.";
    throw FlxException("FunPMF_beta_binomial_ln::calc_02", ssV.str() );
  }
  return rv_pmf_beta_binomial_ln(M,N,alpha,beta);
}

const tdouble FunAutoCorrExp::calc()
{
  return exp(-ONE*(*d)/(*l));
}

const std::string FunAutoCorrExp::write()
{
  string str1 = "autocov_exp(";
  str1+=ConstantBox->get(d);
  str1+=",";
  str1+=ConstantBox->get(l);
  str1+=")";
  return str1;
}

const bool FunAutoCorrExp::optimize(FunBasePtr& optf, const Fun_OptimizeInfo &foi) {
  if (foi.computeConstants) {
    calc_me(optf);
    return true;
  } else {
    return false;
  }
}

const tdouble FunAutoCorrExp2::calc()
{
  return exp(-ONE*pow2((*d)/(*l)));
}

const tdouble FunGamma::calc()
{
  const tdouble c1 = child_1->calc();
  try {
    return flxgamma( c1 );
  } catch (...) {
    FLXMSG("FunGamma::calc",1);
    std::ostringstream ssV;
    ssV << "Error evaluating gamma-function at '" << c1 << "'.";
    throw FlxException("FunGamma::calc", ssV.str() );
  }
}

const tdouble FunLnGamma::calc()
{
  const tdouble c1 = child_1->calc();
  try {
    return GammaLn( c1 );
  } catch (...) {
    FLXMSG("FunLnGamma::calc",1);
    std::ostringstream ssV;
    ssV << "Error evaluating lngamma-function at '" << c1 << "'.";
    throw FlxException("FunLnGamma::calc", ssV.str() );
  }
}

const tdouble FunIGamma::calc()
{
  const tdouble c1 = ParaListP[0]->calc();
  const tdouble c2 = ParaListP[1]->calc();
  try {
    return flxgamma( c1, c2 );
  } catch (...) {
    FLXMSG("tdouble FunIGamma::calc_1",1);
    std::ostringstream ssV;
    ssV << "Error evaluating incomplete upper gamma-function at '(" << c1 << "," << c2 << ")'.";
    throw FlxException("FunIGamma::calc_2", ssV.str() );
  }
}

const tdouble FunIGammaL::calc()
{
  const tdouble c1 = ParaListP[0]->calc();
  const tdouble c2 = ParaListP[1]->calc();
  try {
    return flxgamma_l( c1, c2 );
  } catch (...) {
    FLXMSG("FunIGammaL::calc_1",1);
    std::ostringstream ssV;
    ssV << "Error evaluating incomplete lower gamma-function at '(" << c1 << "," << c2 << ")'.";
    throw FlxException("FunIGammaL::calc_2", ssV.str() );
  }
}

const tdouble FunRGamma::calc()
{
  const tdouble c1 = ParaListP[0]->calc();
  const tdouble c2 = ParaListP[1]->calc();
  try {
    return flxgamma_ru( c1, c2 );
  } catch (...) {
    FLXMSG("FunRGamma::calc_1",1);
    std::ostringstream ssV;
    ssV << "Error evaluating regularized incomplete upper gamma-function at '(" << c1 << "," << c2 << ")'.";
    throw FlxException("FunIGamma::calc_2", ssV.str() );
  }
}

const tdouble FunRGammaL::calc()
{
  const tdouble c1 = ParaListP[0]->calc();
  const tdouble c2 = ParaListP[1]->calc();
  try {
    return flxgamma_rl( c1, c2 );
  } catch (...) {
    FLXMSG("FunRGammaL::calc_1",1);
    std::ostringstream ssV;
    ssV << "Error evaluating regularized incomplete lower gamma-function at '(" << c1 << "," << c2 << ")'.";
    throw FlxException("FunIGammaL::calc_2", ssV.str() );
  }
}

const tdouble FunRGamma_inv::calc()
{
  const tdouble c1 = ParaListP[0]->calc();
  const tdouble c2 = ParaListP[1]->calc();
  try {
    return flxgamma_ru_inv( c1, c2 );
  } catch (...) {
    FLXMSG("FunRGamma::calc_1",1);
    std::ostringstream ssV;
    ssV << "Error evaluating inverse regularized incomplete upper gamma-function at '(" << c1 << "," << c2 << ")'.";
    throw FlxException("FunRGamma_inv::calc_2", ssV.str() );
  }
}

const tdouble FunRGammaL_inv::calc()
{
  const tdouble c1 = ParaListP[0]->calc();
  const tdouble c2 = ParaListP[1]->calc();
  try {
    return flxgamma_rl_inv( c1, c2 );
  } catch (...) {
    FLXMSG("FunRGammaL::calc_1",1);
    std::ostringstream ssV;
    ssV << "Error evaluating inverse regularized incomplete lower gamma-function at '(" << c1 << "," << c2 << ")'.";
    throw FlxException("FunRGammaL_inv::calc_2", ssV.str() );
  }
}

const tdouble FunBeta::calc()
{
  const tdouble alpha = ParaListP[0]->calc();
  const tdouble beta = ParaListP[1]->calc();
  return BetaFun( alpha,beta);
}

const tdouble FunLnBeta::calc()
{
  const tdouble alpha = ParaListP[0]->calc();
  const tdouble beta = ParaListP[1]->calc();
  return BetaFunLn( alpha,beta);
}

const tdouble FunIBeta::calc()
{
  const tdouble alpha = ParaListP[0]->calc();
  const tdouble beta = ParaListP[1]->calc();
  const tdouble x = ParaListP[2]->calc();
  try {
    return iBeta_reg( alpha,beta,x );
  } catch (...) {
    FLXMSG("FunIBeta::calc_1",1);
    std::ostringstream ssV;
    ssV << "Error evaluating regularized incomplete beta-function at '(" << alpha << "," << beta << "," << x << ")'.";
    throw FlxException("FunIBeta::calc_2", ssV.str() );
  }
}

const tdouble FunIBeta_inv::calc()
{
  const tdouble alpha = ParaListP[0]->calc();
  const tdouble beta = ParaListP[1]->calc();
  const tdouble p = ParaListP[2]->calc();
  try {
    return iBeta_reg_inv( alpha,beta,p );
  } catch (FlxException &e) {
    throw;
  } catch (...) {
    FLXMSG("FunIBeta_inv::calc_1",1);
    std::ostringstream ssV;
    ssV << "Error evaluating inverse of the regularized incomplete beta-function at '(" << alpha << "," << beta << "," << p << ")'.";
    throw FlxException("FunIBeta::calc_2", ssV.str() );
  }
}

const tdouble FunErf_inv::calc()
{
  const tdouble dc1 = child_1->calc();
  try {
    return flxerf_inv( dc1 );
  } catch (...) {
    FLXMSG("FunErf_inv::calc_1",1);
    std::ostringstream ssV;
    ssV << "Error evaluating inverse error function for p=" << GlobalVar.Double2String(dc1) << ".";
    throw FlxException("FunErf_inv::calc_2", ssV.str() );
  }
}

const tdouble FunError::calc()
{
  const tdouble c1 = child_1->calc();
  if (c1>ZERO) {
    std::ostringstream ssV;
    ssV << "error-function caused an error: '" << GlobalVar.Double2String(c1) << "'.";
    throw FlxException_NeglectInInteractive("FunError::calc", ssV.str() );
  }
  return c1;
}

const bool FunError::optimize(FunBasePtr& optf, const Fun_OptimizeInfo &foi)
{
  child_optimize(child_1, foi);
  return false;
}
