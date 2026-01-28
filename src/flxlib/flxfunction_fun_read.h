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

#ifndef fesslix_flxfunction_fun_read_H
#define fesslix_flxfunction_fun_read_H

#include "flxfunction_fun_calc.h"


class FunReadFunSin : public FunReadFunBase {
  public:
    FunBase* read ( bool errSerious ) { return new FunSin( read_parameters(1,errSerious) ); }
};

class FunReadFunArcsin : public FunReadFunBase {
  public:
    FunBase* read ( bool errSerious ) { return new FunArcsin( read_parameters(1,errSerious) ); }
};

class FunReadFunCos : public FunReadFunBase {
  public:
    FunBase* read ( bool errSerious ) { return new FunCos( read_parameters(1,errSerious) ); }
};

class FunReadFunArccos : public FunReadFunBase {
  public:
    FunBase* read ( bool errSerious ) { return new FunArccos( read_parameters(1,errSerious) ); }
};

class FunReadFunTan : public FunReadFunBase {
  public:
    FunBase* read ( bool errSerious ) { return new FunTan( read_parameters(1,errSerious) ); }
};

class FunReadFunArctan : public FunReadFunBase {
  public:
    FunBase* read ( bool errSerious ) { return new FunArctan( read_parameters(1,errSerious) ); }
};

class FunReadFunLn : public FunReadFunBase {
  public:
    FunBase* read ( bool errSerious ) { return new FunLn( read_parameters(1,errSerious) ); }
};

class FunReadFunLog : public FunReadFunBase {
  public:
    FunBase* read ( bool errSerious ) { return new FunLog( read_parameters(1,errSerious) ); }
};

class FunReadFunLoga : public FunReadFunBase {
  public:
    FunBase* read ( bool errSerious ) { return new FunLoga( read_parameters(2,errSerious) ); }
};

class FunReadFunExp : public FunReadFunBase {
  public:
    FunBase* read ( bool errSerious ) { return new FunExp( read_parameters(1,errSerious) ); }
};

class FunReadFunCot : public FunReadFunBase {
  public:
    FunBase* read ( bool errSerious ) { return new FunCot( read_parameters(1,errSerious) ); }
};

class FunReadFunArccot : public FunReadFunBase {
  public:
    FunBase* read ( bool errSerious ) { return new FunArccot( read_parameters(1,errSerious) ); }
};

class FunReadFunSinh : public FunReadFunBase {
  public:
    FunBase* read ( bool errSerious ) { return new FunSinh( read_parameters(1,errSerious) ); }
};

class FunReadFunCosh : public FunReadFunBase {
  public:
    FunBase* read ( bool errSerious ) { return new FunCosh( read_parameters(1,errSerious) ); }
};

class FunReadFunTanh : public FunReadFunBase {
  public:
    FunBase* read ( bool errSerious ) { return new FunTanh( read_parameters(1,errSerious) ); }
};

class FunReadFunArctanh : public FunReadFunBase {
  public:
    FunBase* read ( bool errSerious ) { return new FunArctanh( read_parameters(1,errSerious) ); }
};

class FunReadFunSqrt : public FunReadFunBase {
  public:
    FunBase* read ( bool errSerious ) { return new FunSqrt( read_parameters(1,errSerious) ); }
};

class FunReadFunPDFn2 : public FunReadFunBase {
  public:
    FunBase* read ( bool errSerious );
};

class FunReadFunPDFn2_ln : public FunReadFunBase {
  public:
    FunBase* read ( bool errSerious );
};

class FunReadFunPMF_beta_binomial_ln : public FunReadFunBase {
  public:
    FunBase* read ( bool errSerious );
};

class FunReadFunAbs : public FunReadFunBase {
  public:
    FunBase* read ( bool errSerious ) { return new FunAbs( read_parameters(1,errSerious) ); }
};

class FunReadFunFrac : public FunReadFunBase {
  public:
    FunBase* read ( bool errSerious ) { return new FunFrac( read_parameters(1,errSerious) ); }
};

class FunReadFunRound : public FunReadFunBase {
  public:
    FunBase* read ( bool errSerious ) { return new FunRound( read_parameters(-1,errSerious) ); }
};

class FunReadFunMod : public FunReadFunBase {
  public:
    FunBase* read ( bool errSerious ) { return new FunMod( read_parameters(2,errSerious) ); }
};

class FunReadFunSig : public FunReadFunBase {
  public:
    FunBase* read ( bool errSerious ) { return new FunSig( read_parameters(1,errSerious) ); }
};

class FunReadFunFactorial : public FunReadFunBase {
  public:
    FunBase* read ( bool errSerious ) { return new FunFactorial( read_parameters(1,errSerious) ); }
};

class FunReadFunFactorialLn : public FunReadFunBase {
  public:
    FunBase* read ( bool errSerious ) { return new FunFactorialLn( read_parameters(1,errSerious) ); }
};

class FunReadFunPDFn : public FunReadFunBase {
  public:
    FunBase* read ( bool errSerious );
};

class FunReadFunPDFn_ln : public FunReadFunBase {
  public:
    FunBase* read ( bool errSerious );
};

class FunReadFunCDFn : public FunReadFunBase {
  public:
    FunBase* read ( bool errSerious );
};

class FunReadFunInvCDFn : public FunReadFunBase {
  public:
    FunBase* read ( bool errSerious ) { return new FunInvCDFn( read_parameters(1,errSerious) ); }
};

class FunReadFunRelIdx : public FunReadFunBase {
  public:
    FunBase* read ( bool errSerious ) { return new FunRelIdx( read_parameters(1,errSerious) ); }
};

class FunReadFunFailLSF : public FunReadFunBase {
  public:
    FunBase* read ( bool errSerious ) { return new FunFailLSF( read_parameters(1,errSerious) ); }
};

class FunReadFunHvi : public FunReadFunBase {
  public:
    FunBase* read ( bool errSerious ) { return new FunHvi( read_parameters(1,errSerious) ); }
};

class FunReadFunTime0 : public FunReadFunBase {
  public:
    FunBase* read ( bool errSerious ) { return new FunTime0( read_parameters(0,errSerious) ); }
};

class FunReadFunInteg : public FunReadFunBase {
  public:
    FunBase* read ( bool errSerious );
};

class FunReadFunSum : public FunReadFunBase {
  public:
    FunBase* read ( bool errSerious );
};

class FunReadFunOptimize1D : public FunReadFunBase {
  public:
    FunBase* read ( bool errSerious );
};

class FunReadFundeg2gauss : public FunReadFunBase {
  public:
    FunBase* read ( bool errSerious ) { return new FunDeg2Gauss( read_parameters(1,errSerious) ); }
};

class FunReadFunBinomialCoeff : public FunReadFunBase {
  public:
    FunBase* read ( bool errSerious ) { return new FunBinomialCoeff( read_parameters(2,errSerious) ); }
};

class FunReadFunBinomialCoeff_ln : public FunReadFunBase {
  public:
    FunBase* read ( bool errSerious ) { return new FunBinomialCoeff_ln( read_parameters(2,errSerious) ); }
};

class FunReadFunEvalR : public FunReadFunBase {
  public:
    FunBase* read ( bool errSerious );
};

class FunReadFunEvalC : public FunReadFunBase {
  public:
    FunBase* read ( bool errSerious );
};

class FunReadFunEvalW : public FunReadFunBase {
  public:
    FunBase* read ( bool errSerious ) { return new FunEvalW( read_parameters(1,errSerious) ); }
};

class FunReadFunAutoCorrExp : public FunReadFunBase {
  public:
    FunBase* read ( bool errSerious );
};

class FunReadFunAutoCorrExp2 : public FunReadFunBase {
  public:
    FunBase* read ( bool errSerious );
};

class FunReadFunGamma : public FunReadFunBase {
  public:
    FunBase* read ( bool errSerious ) { return new FunGamma( read_parameters(1,errSerious) ); }
};

class FunReadFunLnGamma : public FunReadFunBase {
  public:
    FunBase* read ( bool errSerious ) { return new FunLnGamma( read_parameters(1,errSerious) ); }
};

class FunReadFunIGamma : public FunReadFunBase {
  public:
    FunBase* read ( bool errSerious ) { return new FunIGamma( read_parameters(2,errSerious) ); }
};

class FunReadFunIGammaL : public FunReadFunBase {
  public:
    FunBase* read ( bool errSerious ) { return new FunIGammaL( read_parameters(2,errSerious) ); }
};

class FunReadFunRGamma : public FunReadFunBase {
  public:
    FunBase* read ( bool errSerious ) { return new FunRGamma( read_parameters(2,errSerious) ); }
};

class FunReadFunRGammaL : public FunReadFunBase {
  public:
    FunBase* read ( bool errSerious ) { return new FunRGammaL( read_parameters(2,errSerious) ); }
};

class FunReadFunRGamma_inv : public FunReadFunBase {
  public:
    FunBase* read ( bool errSerious ) { return new FunRGamma_inv( read_parameters(2,errSerious) ); }
};

class FunReadFunRGammaL_inv : public FunReadFunBase {
  public:
    FunBase* read ( bool errSerious ) { return new FunRGammaL_inv( read_parameters(2,errSerious) ); }
};

class FunReadFunBeta : public FunReadFunBase {
  public:
    FunBase* read ( bool errSerious ) { return new FunBeta( read_parameters(2,errSerious) ); }
};

class FunReadFunLnBeta : public FunReadFunBase {
  public:
    FunBase* read ( bool errSerious ) { return new FunLnBeta( read_parameters(2,errSerious) ); }
};

class FunReadFunIBeta : public FunReadFunBase {
  public:
    FunBase* read ( bool errSerious ) { return new FunIBeta( read_parameters(3,errSerious) ); }
};

class FunReadFunIBeta_inv : public FunReadFunBase {
  public:
    FunBase* read ( bool errSerious ) { return new FunIBeta_inv( read_parameters(3,errSerious) ); }
};

class FunReadFunErf : public FunReadFunBase {
  public:
    FunBase* read ( bool errSerious ) { return new FunErf( read_parameters(1,errSerious) ); }
};

class FunReadFunErf_inv : public FunReadFunBase {
  public:
    FunBase* read ( bool errSerious ) { return new FunErf_inv( read_parameters(1,errSerious) ); }
};

class FunReadFunRV_y2x : public FunReadFunBase {
  public:
    FunBase* read ( bool errSerious );
};

class FunReadFunRV_x2y : public FunReadFunBase {
  public:
    FunBase* read ( bool errSerious );
};

class FunReadFunError : public FunReadFunBase {
  public:
    FunBase* read ( bool errSerious ) { return new FunError( read_parameters(1,errSerious) ); }
};

class FunReadFunIsNaN : public FunReadFunBase {
  public:
    FunBase* read ( bool errSerious ) { return new FunIsNaN( read_parameters(1,errSerious) ); }
};


#endif // fesslix_flxfunction_fun_read_H
