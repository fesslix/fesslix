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

#include "flxfunction_ext.h"
#include "flxfunction_fun.h"
#include "flxfunction_ope.h"
#include "flxfunction_ope_calc.h"



FlxFunPoint::FoR FlxFunPoint::get_FoR(char strlFoR, bool errSerious)
{
  switch (strlFoR) {
    case 'c':
      return cartesian;
    case 'p':
      return cylindrical;
    case 's':
      return spherical;
    default:
      std::ostringstream ssV;
      ssV << "Unkown frame of reference '" << strlFoR << "'.";
      FlxError(errSerious,"FunReadPara::set_NumbOfPara", ssV.str() );
      return cartesian;        // DUMMY-RETURN
  }
}

char FlxFunPoint::get_FoR(FlxFunPoint::FoR elFoR)
{
  switch(elFoR) {
    case cartesian:
      return 'c';
      break;
    case cylindrical:
      return 'p';
      break;
    case spherical:
      return 's';
      break;
    default:
      return ' ';
  }
}

FlxFunPoint::FlxFunPoint(const FlxFunPoint& ref)
: FrORef(ref.FrORef), d1(new FlxFunction(*(ref.d1))), d2(new FlxFunction(*(ref.d2))), d3(new FlxFunction(*(ref.d3)))
{
  
}

FlxFunPoint& FlxFunPoint::operator=(const FlxFunPoint& ref)
{
  if (this != & ref) {
    FrORef = ref.FrORef;
    delete d1; delete d2; delete d3;
    d1 = new FlxFunction(*(ref.d1));
    d2 = new FlxFunction(*(ref.d2));
    d3 = new FlxFunction(*(ref.d3));
  }
  return *this;
}

FlxFunPoint::FlxFunPoint(ReadStream* reader, FlxFunctionReader* funReader, bool errSerious)
:FrORef(cartesian),d1(NULL),d2(NULL),d3(NULL)
{
  char c = reader->getChar();
  // Cartesian frame of reference
    if (c=='[' || c=='c' || c=='C') {
      FrORef = cartesian;
    } 
  // Cylindrical frame of reference
    else if (c=='p' || c=='P') {
      FrORef = cylindrical;
    }
  // Spherical frame of reference
    else if (c=='s' || c=='S') {
      FrORef = spherical;
    }
  // Syntax-Error
    else {
      std::ostringstream ssV;
      ssV << "Expected input of a point and not '" << c << "'.";
      FlxError(errSerious,"FlxFunPoint::FlxFunPoint_1", ssV.str() );
    }
  if (c!='[') reader->getChar('[',errSerious);
  // input coordinate-information
    d1 = new FlxFunction(funReader,errSerious);
    if (reader->whatIsNextChar()==']') {
      d2 = new FlxFunction(new FunNumber(0.0) );
    } else {
      reader->getChar(',');
      d2 = new FlxFunction(funReader,errSerious);
    }
    if (reader->whatIsNextChar()==']') {
      d3 = new FlxFunction(new FunNumber(0.0) );
    } else {
      reader->getChar(',');
      d3 = new FlxFunction(funReader,errSerious);
    }
  // Finish the input of the point
    reader->getChar(']',errSerious);
}

FlxFunPoint::FlxFunPoint(const FlxFunPoint::FoR FrORef, FlxFunction* d1, FlxFunction* d2, FlxFunction* d3)
: FrORef(FrORef), d1(d1), d2(d2), d3(d3) 
{
   
}

FlxFunPoint::~FlxFunPoint()
{
  delete d1; 
  delete d2; 
  delete d3;
}

const flxPoint FlxFunPoint::calc() const
{
  tdouble rr;
  tdouble pp;
  tdouble tt;
  switch (FrORef) {
    case cartesian:
      return flxPoint(d1->calc(), d2->calc(), d3->calc());
    case cylindrical:
      rr = d1->calc();
      pp = d2->calc();
      return flxPoint(rr*cos(pp), rr*sin(pp), d3->calc());
    case spherical:
      rr = d1->calc();
      pp = d2->calc();
      tt = d3->calc();
      return flxPoint(rr*cos(pp)*sin(tt), rr*sin(pp)*sin(tt), rr*cos(tt));      
    default:
      #if FLX_DEBUG
        throw FlxException_Crude("FlxFunPoint::calc");
      #else
        return flxPoint(ZERO,ZERO,ZERO);
      #endif
  }
}

const bool FlxFunPoint::is_direction_zero(const char dir)
{
  switch (dir) {
    case 'x':
      return d1->is_zero();
    case 'y':
      return d2->is_zero();
    case 'z':
      return d3->is_zero();
    default:
      throw FlxException_Crude("FlxFunPoint::is_direction_zero");
  }
}

std::ostream& operator<<(std::ostream& os, const FlxFunPoint& val)
{
  return os << val.get_FoR(val.FrORef) << "[" << val.d1->write() << "," << val.d2->write() << "," << val.d3->write() << "]";
}


const tdouble FlxFun_GaussIntegration(FunBase* funI, tdouble* theconst, const tdouble& start, const tdouble& end, const tuint& GPN, GaussIntegration& GI)
{
  const tdouble theconst_orig = *theconst;        // remember the value of theconst 
  // get number of integration points
    GI.check_GA(GPN);
    const tdouble* GP = GI.get_GP(GPN);
    const tdouble* GW = GI.get_GW(GPN);
  // determine start and end
    const tdouble c1 = (end-start)/2;
    const tdouble c2 = (end+start)/2;
  // perform the integration
    tdouble sum = ZERO;
    for (tuint p=0;p<GPN;++p) {
      *theconst = c1*GP[p] + c2;
      sum += funI->calc() * GW[p];
    }
    sum*=c1;
    #if FLX_DEBUG
      if ( std::isnan(sum) ) throw FlxException_Crude("FlxFun_GaussIntegration");
    #endif
  *theconst = theconst_orig;
  return sum;
}

const tuint FlxFun_EstimateGaussPoint(FunBase* funI, tdouble* theconst, const tdouble& start, const tdouble& end, const tuint& GPmax, const tuint& GPtest, const tdouble& err, GaussIntegration& GI)
{
  #if FLX_DEBUG
    if (GPmax+2 > GPtest) {
      std::ostringstream ssV;
      ssV << "ERROR";
      throw FlxException("FlxFun_EstimateGaussPoint", ssV.str() );
    }
  #endif
  // check if the function is constant
    const tdouble theconst_orig = *theconst;
    *theconst = start; 
      tdouble t = funI->calc();
    *theconst = end; 
      tdouble t2 = funI->calc();
    if (t==t2) {
      *theconst = (start+end)/2; 
        t2 = funI->calc();
      if (t==t2) {
        const tdouble c1 = (end-start)/2;
        const tdouble c2 = (end+start)/2;
        *theconst = c2+c1/sqrt(3.); 
          t2 = funI->calc();
        if (t==t2) {
          *theconst = c2-c1/sqrt(3.); 
            t2 = funI->calc();
          if (t==t2) {
            *theconst = theconst_orig;
            return 0;
          }
        }
      }
    }
    *theconst = theconst_orig;
  
  // compute the reference value
    const tdouble refRes = FlxFun_GaussIntegration(funI,theconst,start,end,GPtest,GI);
  
  // do the actual iteration
    for (tuint i=1;i<=GPmax;++i) {
      t = FlxFun_GaussIntegration(funI,theconst,start,end,i,GI);
      if (fabs((t-refRes)/refRes)<=err) {
        t = FlxFun_GaussIntegration(funI,theconst,start,end,i+1,GI);
        if (fabs((t-refRes)/refRes)<=err) {
          return i;
        }
      }
    }
  return GPmax+1;  
}




const tdouble FlxFun_RootSearch_Bisec(FunBase* funR, tdouble* theconst, tdouble start, tdouble end, const tdouble dx, const tdouble dy, ostreamp streamp)
{
  const tdouble storedF = *theconst;
  flx_fun_void_data ds;
    ds.theconst = theconst;
    ds.funR = funR;
  const tdouble res = flx_RootSearch_Bisec(&flx_fun_void,&ds,start,end,dx,dy,streamp);
  *theconst = storedF;
  return res;
}

const tdouble FlxFun_RootSearch_RegulaFalsi(FunBase* funR, tdouble* theconst, tdouble start, tdouble end, const tdouble dx, const tdouble dy, ostreamp streamp)
{
  const tdouble storedF = *theconst;
  flx_fun_void_data ds;
    ds.theconst = theconst;
    ds.funR = funR;
  const tdouble res = flx_RootSearch_RegulaFalsi(&flx_fun_void,&ds,start,end,dx,dy,streamp);
  *theconst = storedF;
  return res;
}

