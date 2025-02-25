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

#ifndef fesslix_flxfunction_fun_calc_H
#define fesslix_flxfunction_fun_calc_H

#include "flxfunction_ext.h"



/**
* @brief arithmetic class: FunEvalW
*/
class FunEvalW : public FunBaseFun_onePara {
  public:
    FunEvalW (std::vector<FunBase*> *ParaListV) : FunBaseFun_onePara(ParaListV) {};
    const tdouble calc() { return child_1->calc(); }
    const std::string write_v() { return "evalw";}
    virtual const bool evalw() { return true; }
};

/**
* @brief arithmetic class: sinus
*/
class FunSin : public FunBaseFun_onePara {
  public:
    FunSin (std::vector<FunBase*> *ParaListV) : FunBaseFun_onePara(ParaListV) {};
    const tdouble calc() { return std::sin( child_1->calc() ); }
    const std::string write_v() { return "sin";}
};

/**
* @brief arithmetic class: arcsin
*/
class FunArcsin : public FunBaseFun_onePara {
  public:
    FunArcsin (std::vector<FunBase*> *ParaListV) : FunBaseFun_onePara(ParaListV) {};
    const tdouble calc() { return std::asin( child_1->calc() ); }
    const std::string write_v() { return "arcsin";}
};

/**
* @brief arithmetic class: cosinus
*/
class FunCos : public FunBaseFun_onePara {
  public:
    FunCos (std::vector<FunBase*> *ParaListV) : FunBaseFun_onePara(ParaListV) {};
    const tdouble calc() { return std::cos( child_1->calc() ); }
    const std::string write_v() { return "cos";}
};

/**
* @brief arithmetic class: arccos
*/
class FunArccos : public FunBaseFun_onePara {
  public:
    FunArccos (std::vector<FunBase*> *ParaListV) : FunBaseFun_onePara(ParaListV) {};
    const tdouble calc() { return std::acos( child_1->calc() ); }
    const std::string write_v() { return "arccos";}
};

/**
* @brief arithmetic class: tangens
*/
class FunTan : public FunBaseFun_onePara {
  public:
    FunTan (std::vector<FunBase*> *ParaListV) : FunBaseFun_onePara(ParaListV) {};
    const tdouble calc() { return std::tan( child_1->calc() ); }
    const std::string write_v() { return "tan";}
};

/**
* @brief arithmetic class: arctan
*/
class FunArctan : public FunBaseFun_onePara {
  public:
    FunArctan (std::vector<FunBase*> *ParaListV) : FunBaseFun_onePara(ParaListV) {};
    const tdouble calc() { return std::atan( child_1->calc() ); }
    const std::string write_v() { return "arctan";}
};

/**
* @brief arithmetic class: ln
*/
class FunLn : public FunBaseFun_onePara {
  public:
    FunLn (FunBase *child_1) : FunBaseFun_onePara(child_1) {};
    FunLn (std::vector<FunBase*> *ParaListV) : FunBaseFun_onePara(ParaListV) {};
    const tdouble calc() { return std::log( child_1->calc() ); }
    const std::string write_v() { return "ln";}
};

/**
* @brief arithmetic class: log
*/
class FunLog : public FunBaseFun_onePara {
  public:
    FunLog (std::vector<FunBase*> *ParaListV) : FunBaseFun_onePara(ParaListV) {};
    const tdouble calc() { return std::log10( child_1->calc() ); }
    const std::string write_v() { return "log";}
};

/**
* @brief arithmetic class: log to arbitrary base
*/
class FunLoga : public FunBaseFun_multPara {
  public:
    FunLoga (std::vector<FunBase*> *ParaListV) : FunBaseFun_multPara(ParaListV) {};
    const tdouble calc();
    const std::string write_v() { return "loga";}
};

/**
* @brief arithmetic class: exp = e^x
*/
class FunExp : public FunBaseFun_onePara {
  public:
    FunExp (FunBase *child_1) : FunBaseFun_onePara(child_1) {};
    FunExp (std::vector<FunBase*> *ParaListV) : FunBaseFun_onePara(ParaListV) {};
    const tdouble calc() { return std::exp( child_1->calc() ); }
    const std::string write_v() { return "exp";}
};

/**
* @brief arithmetic class: cot
*/
class FunCot : public FunBaseFun_onePara {
  public:
    FunCot (std::vector<FunBase*> *ParaListV) : FunBaseFun_onePara(ParaListV) {};
    const tdouble calc() { return ONE / std::tan( child_1->calc() ); }
    const std::string write_v() { return "cot";}
};

/**
* @brief arithmetic class: arccot
*/
class FunArccot : public FunBaseFun_onePara {
  public:
    FunArccot (std::vector<FunBase*> *ParaListV) : FunBaseFun_onePara(ParaListV) {};
    const tdouble calc() { return arccot( child_1->calc() ); }
    const std::string write_v() { return "arccot";}
};

/**
* @brief arithmetic class: sinh
*/
class FunSinh : public FunBaseFun_onePara {
  public:
    FunSinh (std::vector<FunBase*> *ParaListV) : FunBaseFun_onePara(ParaListV) {};
    const tdouble calc() { return std::sinh( child_1->calc() ); }
    const std::string write_v() { return "sinh";}
};

/**
* @brief arithmetic class: cosh
*/
class FunCosh : public FunBaseFun_onePara {
  public:
    FunCosh (std::vector<FunBase*> *ParaListV) : FunBaseFun_onePara(ParaListV) {};
    const tdouble calc() { return std::cosh( child_1->calc() ); }
    const std::string write_v() { return "cosh";}
};

/**
* @brief arithmetic class: tanh
*/
class FunTanh : public FunBaseFun_onePara {
  public:
    FunTanh (std::vector<FunBase*> *ParaListV) : FunBaseFun_onePara(ParaListV) {};
    const tdouble calc() { return std::tanh( child_1->calc() ); }
    const std::string write_v() { return "tanh";}
};

/**
* @brief arithmetic class: arctanh
*/
class FunArctanh : public FunBaseFun_onePara {
  public:
    FunArctanh (std::vector<FunBase*> *ParaListV) : FunBaseFun_onePara(ParaListV) {};
    const tdouble calc() { return arctanh( child_1->calc() ); }
    const std::string write_v() { return "arctanh";}
};

/**
* @brief arithmetic class: sqrt
*/
class FunSqrt : public FunBaseFun_onePara {
  public:
    FunSqrt (std::vector<FunBase*> *ParaListV) : FunBaseFun_onePara(ParaListV) {};
    const tdouble calc() { return std::sqrt( child_1->calc() ); }
    const std::string write_v() { return "sqrt";}
};

/**
* @brief arithmetic class: abs
*/
class FunAbs : public FunBaseFun_onePara {
  public:
    FunAbs (std::vector<FunBase*> *ParaListV) : FunBaseFun_onePara(ParaListV) {};
    const tdouble calc() { return std::fabs( child_1->calc() ); }
    const std::string write_v() { return "abs";}
};

/**
* @brief arithmetic class: frac
*/
class FunFrac : public FunBaseFun_onePara {
  public:
    FunFrac (std::vector<FunBase*> *ParaListV) : FunBaseFun_onePara(ParaListV) {};
    const tdouble calc() { return frac( child_1->calc() ); }
    const std::string write_v() { return "frac";}
};

/**
* @brief arithmetic class: round
*/
class FunRound : public FunBaseFun_multPara {
  public:
    FunRound (std::vector<FunBase*> *ParaListV) : FunBaseFun_multPara(ParaListV) {};
    const tdouble calc();
    const std::string write_v() { return "round";}
};

/**
* @brief arithmetic class: mod
*/
class FunMod : public FunBaseFun_multPara {
  public:
    FunMod (std::vector<FunBase*> *ParaListV) : FunBaseFun_multPara(ParaListV) {};
    const tdouble calc() { return int(round_flx(ParaListP[0]->calc()))%int(round_flx(ParaListP[1]->calc())); }
    const std::string write_v() { return "mod";}
};

/**
* @brief arithmetic class: sig
*/
class FunSig : public FunBaseFun_onePara {
  public:
    FunSig (std::vector<FunBase*> *ParaListV) : FunBaseFun_onePara(ParaListV) {};
    const tdouble calc() { 
      const tdouble res = child_1->calc();
      if (std::fabs(res)<=GlobalVar.TOL()) return ZERO;
      else if ( res > ZERO ) return ONE;
      else return -ONE;
    }
    const std::string write_v() { return "sig";}
};

/**
* @brief arithmetic class: factorial
*/
// TODO: tread overflow
class FunFactorial : public FunBaseFun_onePara {
  public:
    FunFactorial (std::vector<FunBase*> *ParaListV) : FunBaseFun_onePara(ParaListV) {};
    const tdouble calc() { 
      const int n=int(round_flx(child_1->calc()));
      return Factorial(n);
    }
    const std::string write_v() { return "factorial";}
};

/**
* @brief arithmetic class: factorialLn
*/
// TODO: tread overflow
class FunFactorialLn : public FunBaseFun_onePara {
  public:
    FunFactorialLn (std::vector<FunBase*> *ParaListV) : FunBaseFun_onePara(ParaListV) {};
    const tdouble calc() { 
      const int n=int(round_flx(child_1->calc()));
      return FactorialLn(n);
    }
    const std::string write_v() { return "factorialln";}
};

/**
* @brief arithmetic class: pdfn
*/
class FunPDFn_general : public FunBaseFun_multPara {
  public:
    FunPDFn_general (std::vector<FunBase*> *ParaListV) : FunBaseFun_multPara(ParaListV) {};
    const tdouble calc();
    const std::string write_v() { return "pdfn";}
};

/**
* @brief arithmetic class: pdfn
*/
class FunPDFn : public FunBaseFun_onePara {
  public:
    FunPDFn (std::vector<FunBase*> *ParaListV) : FunBaseFun_onePara(ParaListV) {};
    const tdouble calc() { return  rv_phi( child_1->calc() ); }
    const std::string write_v() { return "pdfn";}
};

/**
* @brief arithmetic class: pdfn
*/
class FunPDFn_ln_general : public FunBaseFun_multPara {
  public:
    FunPDFn_ln_general (std::vector<FunBase*> *ParaListV) : FunBaseFun_multPara(ParaListV) {};
    const tdouble calc();
    const std::string write_v() { return "pdfn_ln";}
};

/**
* @brief arithmetic class: pdfn
*/
class FunPDFn_ln : public FunBaseFun_onePara {
  public:
    FunPDFn_ln (std::vector<FunBase*> *ParaListV) : FunBaseFun_onePara(ParaListV) {};
    const tdouble calc() { return  -pow2(child_1->calc())/2-log(2*PI)/2; }
    const std::string write_v() { return "pdfn_ln";}
};

/**
* @brief arithmetic class: cdfn
*/
class FunCDFn_general : public FunBaseFun_multPara {
  public:
    FunCDFn_general (std::vector<FunBase*> *ParaListV) : FunBaseFun_multPara(ParaListV) {};
    const tdouble calc();
    const std::string write_v() { return "cdfn";}
};

/**
* @brief arithmetic class: cdfn
*/
class FunCDFn : public FunBaseFun_onePara {
  public:
    FunCDFn (std::vector<FunBase*> *ParaListV) : FunBaseFun_onePara(ParaListV) {};
    const tdouble calc() { return rv_Phi( child_1->calc() ); }
    const std::string write_v() { return "cdfn";}
};

/**
* @brief arithmetic class: invcdfn
*/
class FunInvCDFn : public FunBaseFun_onePara {
  public:
    FunInvCDFn (std::vector<FunBase*> *ParaListV) : FunBaseFun_onePara(ParaListV) {};
    const tdouble calc() { return rv_InvPhi( child_1->calc() ); }
    const std::string write_v() { return "cdfn_inv";}
};

/**
* @brief arithmetic class: relidx
*/
class FunRelIdx : public FunBaseFun_onePara {
  public:
    FunRelIdx (std::vector<FunBase*> *ParaListV) : FunBaseFun_onePara(ParaListV) {};
    const tdouble calc() { return -rv_InvPhi( child_1->calc() ); }
    const std::string write_v() { return "relidx";}
};

/**
* @brief arithmetic class: faillsf
*/
class FunFailLSF : public FunBaseFun_onePara {
  public:
    FunFailLSF (std::vector<FunBase*> *ParaListV) : FunBaseFun_onePara(ParaListV) {};
    const tdouble calc() { return (child_1->calc() <= ZERO); }
    const std::string write_v() { return "faillsf";}
};

/**
* @brief arithmetic class: hvi
*/
class FunHvi : public FunBaseFun_onePara {
  public:
    FunHvi (std::vector<FunBase*> *ParaListV) : FunBaseFun_onePara(ParaListV) {};
    const tdouble calc() { return ( !(child_1->calc() < ZERO));  }
    const std::string write_v() { return "hvi";}
};

/**
* @brief arithmetic class: time0
*/
class FunTime0 : public FunBaseFun_multPara {
  public:
    FunTime0 (std::vector<FunBase*> *ParaListV) : FunBaseFun_multPara(ParaListV) {};
    const tdouble calc() { return tdouble(std::time(0)); }
    const std::string write() { return "time0()";}
    const std::string write_v() { return "time0";}
    const bool optimize(FunBasePtr& optf, const Fun_OptimizeInfo &foi) { return false; }
};


/**
* @brief arithmetic class: numerical integration of a function
*/
class FLXLIB_EXPORT FunInteg : public FunBase, public FlxBoxBase {
  private:
    /**
    * @brief the function to integrate
    */
    FunBase* funI;
    /**
    * @brief the integration variable
    */
    tdouble* theconst;
    /**
    * @brief start of the integration
    */
    FunBase* startF;
    /**
    * @brief end of the integration
    */
    FunBase* endF;
    /**
    * @brief the number of Gauss-points to use
    */
    FunBase* gaussF;
    /**
    * @brief the number of intervals to use
    */
    FunBase* intF;
    /**
    * @brief integrate in log-scaled intervals
    */
    const bool is_log;
  public:
    FunInteg(FunBase* funI, tdouble* theconst, FunBase* startF, FunBase* endF, FunBase* gaussF, FunBase* intF, const bool is_log=false) : funI(funI),theconst(theconst),startF(startF),endF(endF),gaussF(gaussF),intF(intF), is_log(is_log) {}
    ~FunInteg();
    const tdouble calc();
    const bool search_circref(FlxFunction* fcr);
    const std::string write();
    const tuint precedence() { return 0; }
    const bool dependOn_Const(const tdouble* const thenumber);
    const bool optimize(FunBasePtr& optf, const Fun_OptimizeInfo &foi);
};

/**
* @brief arithmetic class: numerical integration of a function
*/
class FunSum : public FunBase, public FlxBoxBase {
  private:
    /**
    * @brief the function to sum
    */
    FunBase* funS;
    /**
    * @brief the integration variable
    */
    tdouble* theconst;
    /**
    * @brief start of the integration
    */
    FunBase* startF;
    /**
    * @brief condition
    */
    FunBase* cond;
    /**
    * @brief step-expression
    */
    FunBase* step;
  public:
    FunSum(FunBase* funS, tdouble* theconst, FunBase* startF, FunBase* cond, FunBase* step) : funS(funS),theconst(theconst),startF(startF),cond(cond),step(step) {}
    ~FunSum();
    const tdouble calc();
    const bool search_circref(FlxFunction* fcr);
    const std::string write();
    const tuint precedence() { return 0; }
    const bool dependOn_Const(const tdouble* const thenumber);
    const bool optimize(FunBasePtr& optf, const Fun_OptimizeInfo &foi);
};

/**
* @brief arithmetic class: root search (a single one) of a function
*/
class FLXLIB_EXPORT FunRoot : public FunBase, public FlxBoxBase {
  protected:
    /**
    * @brief the ID of the method
    * 0: bisection
    * 1: regula falsi
    */
    const int ID;
    /**
    * @brief the function to investigate
    */
    FunBase* funR;
    /**
    * @brief the root-search variable
    */
    tdouble* theconst;
    /**
    * @brief start of the serach-interval
    */
    FunBase* startF;
    /**
    * @brief end of the search-interval
    */
    FunBase* endF;
    /**
    * @brief stop criteria
    */
    FunBase* dx;
    FunBase* dy;
    /**
    * @brief stream for logging
    */
    ostreamp streamp;
  public:
    FunRoot(const int ID, FunBase* funR, tdouble* theconst, FunBase* startF, FunBase* endF, FunBase* dx, FunBase* dy, ostreamp streamp) : ID(ID), funR(funR),theconst(theconst),startF(startF),endF(endF),dx(dx),dy(dy),streamp(streamp) {}
    ~FunRoot();
    const tdouble calc();
    const bool search_circref(FlxFunction* fcr);
    const tuint precedence() { return 0; }
    const bool dependOn_Const(const tdouble* const thenumber);
    const bool optimize(FunBasePtr& optf, const Fun_OptimizeInfo &foi);
    const std::string write();
};

/**
* @brief arithmetic class: root search (a single one) of a function
*/
class FunOptimize1D : public FunBase, public FlxBoxBase {
  protected:
    /**
    * @brief the function to investigate
    */
    FunBase* funR;
    /**
    * @brief the root-search variable
    */
    tdouble* theconst;
    /**
    * @brief start of the serach-interval
    */
    FunBase* startF;
    /**
    * @brief end of the search-interval
    */
    FunBase* endF;
    /**
    * @brief stop criteria
    */
    FunBase* eps1;
    FunBase* eps2;
    FunBase* niter;
    FunBase* nex;
    const bool use_brent;
  public:
    FunOptimize1D(FunBase* funR, tdouble* theconst, FunBase* startF, FunBase* endF, FunBase* eps1, FunBase* eps2, FunBase* niter, FunBase* nex, const bool use_brent) : funR(funR),theconst(theconst),startF(startF),endF(endF),eps1(eps1),eps2(eps2),niter(niter),nex(nex), use_brent(use_brent) {}
    ~FunOptimize1D();
    const tdouble calc();
    const bool search_circref(FlxFunction* fcr);
    const tuint precedence() { return 0; }
    const bool dependOn_Const(const tdouble* const thenumber);
    const bool optimize(FunBasePtr& optf, const Fun_OptimizeInfo &foi);
    const std::string write();
};

/**
* @brief arithmetic class: sinus
*/
class FunDeg2Gauss : public FunBaseFun_onePara {
  public:
    FunDeg2Gauss (std::vector<FunBase*> *ParaListV) : FunBaseFun_onePara(ParaListV) {};
    const tdouble calc();
    const std::string write_v() { return "deg2gauss";}
};

/**
* @brief arithmetic class: binomial coefficient
*/
class FunBinomialCoeff : public FunBaseFun_multPara {
  public:
    FunBinomialCoeff (std::vector<FunBase*> *ParaListV) : FunBaseFun_multPara(ParaListV) {};
    const tdouble calc();
    const std::string write_v() { return "binomialcoeff";}
};

/**
* @brief arithmetic class: log of the binomial coefficient
*/
class FunBinomialCoeff_ln : public FunBaseFun_multPara {
  public:
    FunBinomialCoeff_ln (std::vector<FunBase*> *ParaListV) : FunBaseFun_multPara(ParaListV) {};
    const tdouble calc();
    const std::string write_v() { return "binomialcoeff_ln";}
};

/**
* @brief arithmetic class: pdfn2
*/
class FLXLIB_EXPORT FunPDFn2 : public FunBaseFun_multPara {
  public:
    FunPDFn2 (std::vector<FunBase*> *ParaListV) : FunBaseFun_multPara(ParaListV) {};
    const tdouble calc();
    const std::string write_v() { return "pdfn2";}
};
class FunPDFn2_general : public FunBaseFun_multPara {
  public:
    FunPDFn2_general (std::vector<FunBase*> *ParaListV) : FunBaseFun_multPara(ParaListV) {};
    const tdouble calc();
    const std::string write_v() { return "pdfn2";}
};

/**
* @brief arithmetic class: pdfn2_ln
*/
class FLXLIB_EXPORT FunPDFn2_ln : public FunBaseFun_multPara {
  public:
    FunPDFn2_ln (std::vector<FunBase*> *ParaListV) : FunBaseFun_multPara(ParaListV) {};
    const tdouble calc();
    const std::string write_v() { return "pdfn2_ln";}
};
class FunPDFn2_ln_general : public FunBaseFun_multPara {
  public:
    FunPDFn2_ln_general (std::vector<FunBase*> *ParaListV) : FunBaseFun_multPara(ParaListV) {};
    const tdouble calc();
    const std::string write_v() { return "pdfn2_ln";}
};

/**
* @brief arithmetic class: pmf_beta_binomial_ln
*/
class FLXLIB_EXPORT FunPMF_beta_binomial_ln : public FunBaseFun_multPara {
  public:
    FunPMF_beta_binomial_ln (std::vector<FunBase*> *ParaListV) : FunBaseFun_multPara(ParaListV) {};
    const tdouble calc();
    const std::string write_v() { return "pmf_beta_binomial_ln";}
};

/**
* @brief arithmetic class: numerical integration of a function
*/
class FunAutoCorrExp : public FunBase, public FlxBoxBase {
  protected:
    const tdouble* d;
    const tdouble* l;
  public:
    FunAutoCorrExp(const tdouble* d, const tdouble* l) : d(d), l(l) {}
    const tdouble calc();
    const bool search_circref(FlxFunction* fcr) { return false; }
    const std::string write();
    const tuint precedence() { return 0; }
    const bool dependOn_Const(const tdouble* const thenumber) { return thenumber==d || thenumber==l;} 
    const bool optimize(FunBasePtr& optf, const Fun_OptimizeInfo &foi);
};

/**
* @brief arithmetic class: numerical integration of a function
*/
class FunAutoCorrExp2 : public FunAutoCorrExp {
  public:
    FunAutoCorrExp2(const tdouble* d, const tdouble* l) : FunAutoCorrExp(d,l) {}
    const tdouble calc();
};

/**
* @brief arithmetic class: Gamma-Function
*/
class FunGamma : public FunBaseFun_onePara {
  public:
    FunGamma (std::vector<FunBase*> *ParaListV) : FunBaseFun_onePara(ParaListV) {};
    const tdouble calc();
    const std::string write_v() { return "tgamma";}
};

/**
* @brief arithmetic class: LnGamma-Function
*/
class FunLnGamma : public FunBaseFun_onePara {
  public:
    FunLnGamma (std::vector<FunBase*> *ParaListV) : FunBaseFun_onePara(ParaListV) {};
    const tdouble calc();
    const std::string write_v() { return "lngamma";}
};
/**
* @brief arithmetic class: incomplete Gamma-Function
*/
class FunIGamma : public FunBaseFun_multPara {
  public:
    FunIGamma (std::vector<FunBase*> *ParaListV) : FunBaseFun_multPara(ParaListV) {};
    const tdouble calc();
    const std::string write_v() { return "igamma";}
};

/**
* @brief arithmetic class: incomplete Gamma-Function
*/
class FunIGammaL : public FunBaseFun_multPara {
  public:
    FunIGammaL (std::vector<FunBase*> *ParaListV) : FunBaseFun_multPara(ParaListV) {};
    const tdouble calc();
    const std::string write_v() { return "igammal";}
};

/**
* @brief arithmetic class: incomplete upper Gamma-Function
*/
class FunRGamma : public FunBaseFun_multPara {
  public:
    FunRGamma (std::vector<FunBase*> *ParaListV) : FunBaseFun_multPara(ParaListV) {};
    const tdouble calc();
    const std::string write_v() { return "rgamma";}
};

/**
* @brief arithmetic class: incomplete lower Gamma-Function
*/
class FunRGammaL : public FunBaseFun_multPara {
  public:
    FunRGammaL (std::vector<FunBase*> *ParaListV) : FunBaseFun_multPara(ParaListV) {};
    const tdouble calc();
    const std::string write_v() { return "rgammal";}
};

/**
* @brief arithmetic class: inverse incomplete upper Gamma-Function
*/
class FunRGamma_inv : public FunBaseFun_multPara {
  public:
    FunRGamma_inv (std::vector<FunBase*> *ParaListV) : FunBaseFun_multPara(ParaListV) {};
    const tdouble calc();
    const std::string write_v() { return "rgamma_inv";}
};

/**
* @brief arithmetic class: inverse incomplete lower Gamma-Function
*/
class FunRGammaL_inv : public FunBaseFun_multPara {
  public:
    FunRGammaL_inv (std::vector<FunBase*> *ParaListV) : FunBaseFun_multPara(ParaListV) {};
    const tdouble calc();
    const std::string write_v() { return "rgammal_inv";}
};

/**
* @brief arithmetic class: regularized incomplete Beta-Function
*/
class FunBeta : public FunBaseFun_multPara {
  public:
    FunBeta (std::vector<FunBase*> *ParaListV) : FunBaseFun_multPara(ParaListV) {};
    const tdouble calc();
    const std::string write_v() { return "betafun";}
};

/**
* @brief arithmetic class: regularized incomplete Beta-Function (log-transform)
*/
class FunLnBeta : public FunBaseFun_multPara {
  public:
    FunLnBeta (std::vector<FunBase*> *ParaListV) : FunBaseFun_multPara(ParaListV) {};
    const tdouble calc();
    const std::string write_v() { return "lnbetafun";}
};

/**
* @brief arithmetic class: regularized incomplete Beta-Function
*/
class FunIBeta : public FunBaseFun_multPara {
  public:
    FunIBeta (std::vector<FunBase*> *ParaListV) : FunBaseFun_multPara(ParaListV) {};
    const tdouble calc();
    const std::string write_v() { return "ibeta";}
};

/**
* @brief arithmetic class: regularized incomplete Beta-Function
*/
class FunIBeta_inv : public FunBaseFun_multPara {
  public:
    FunIBeta_inv (std::vector<FunBase*> *ParaListV) : FunBaseFun_multPara(ParaListV) {};
    const tdouble calc();
    const std::string write_v() { return "ibeta_inv";}
};

/**
* @brief arithmetic class: erf
*/
class FunErf : public FunBaseFun_onePara {
  public:
    FunErf (std::vector<FunBase*> *ParaListV) : FunBaseFun_onePara(ParaListV) {};
    const tdouble calc() { return flxerf( child_1->calc() ); }
    const std::string write_v() { return "erf";}
};

/**
* @brief arithmetic class: erf_inv
*/
class FunErf_inv : public FunBaseFun_onePara {
  public:
    FunErf_inv (std::vector<FunBase*> *ParaListV) : FunBaseFun_onePara(ParaListV) {};
    const tdouble calc();
    const std::string write_v() { return "erf_inv";}
};

/**
* @brief arithmetic class: error-Function
*         throws an error if the expression is true
*/
class FunError : public FunBaseFun_onePara {
  public:
    FunError (std::vector<FunBase*> *ParaListV) : FunBaseFun_onePara(ParaListV) {};
    const tdouble calc();
    const std::string write_v() { return "error";}
    virtual const bool optimize(FunBasePtr& optf, const Fun_OptimizeInfo &foi);
};

/**
* @brief arithmetic class: isNaN
*/
class FunIsNaN : public FunBaseFun_onePara {
  public:
    FunIsNaN (std::vector<FunBase*> *ParaListV) : FunBaseFun_onePara(ParaListV) {};
    const tdouble calc() { const tdouble t=child_1->calc(); return ((t!=t)?ONE:ZERO); }
    const std::string write_v() { return "isnan";}
};


#endif // fesslix_flxfunction_fun_calc_H
