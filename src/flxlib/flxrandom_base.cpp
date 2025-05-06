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

#include "flxrandom_base.h"

#include "flxfunction_fun_calc.h"
#include "flxfunction_ope_calc.h"

#if FLX_DEBUG
  int FlxRndCreator::Cinst = 0;
#endif
FlxRndCreator* FlxBoxBaseR::RndCreator = NULL;
GaussIntegration* FlxBoxBaseR::GI = NULL;


FlxRndCreator* RndCreator_ptr = nullptr;


FlxRndCreator& get_RndCreator()
{
  if (RndCreator_ptr) {
    return *RndCreator_ptr;
  } else {
    throw FlxException("get_RndCreator","Please start the engine.");
  }
}

void set_RndCreator_ptr(FlxRndCreator* RndCreator_ptr_)
{
  RndCreator_ptr = RndCreator_ptr_;
}



void FlxBoxBaseR::set_Boxes(FlxRndCreator* RndCreatorV)
{
  RndCreator = RndCreatorV;
}


void FlxRndCreator::replay_start(FlxIstream* rndReaderV)
{
  replay_stop();
  rndReader = rndReaderV;
  GlobalVar.slog(3) << "rnd track: started replay of semi random values." << std::endl;
}

FlxRndCreator::FlxRndCreator(rng_type* rngp) : rngp(rngp), rndReader(NULL)
{
  if (rngp==NULL) {
    #if FLX_DEBUG
      ++Cinst;
      if (Cinst > 1) {
        std::ostringstream ssV;
        ssV << "More than one instance of 'FlxRndBox' created ...";
        throw FlxException("FlxRndBox::FlxRndBox", ssV.str() );
      }
    #endif
    FlxBoxBaseR::set_Boxes(this);
  }
}

const tdouble FlxRndCreator::gen_smp()
{
  if (rndReader == NULL) {
    if (rngp) {
      return rv_normal(*rngp);
    } else {
      return rv_normal();
    }
  } else {
    tdouble t;
    if (!rndReader->get_value(t)) {
      replay_stop();
      GlobalVar.alert.alert("FlxRndCreator::gen_smp","Replay of semi random values stopped.");
      return gen_smp();
    }
    return t;
  }
}

const tdouble FlxRndCreator::gen_smp_uniform()
{
  if (rndReader == NULL) {
    if (rngp) {
      return rv_uniform(*rngp);
    } else {
      return rv_uniform();
    }
  } else {
    return rv_Phi(gen_smp());
  }
}

const tuint FlxRndCreator::gen_smp_index(const flxVec& v)
{
  const tdouble p = gen_smp_uniform();
  tuint s = 0;        // start of block
  tuint l = v.get_N();        // length of block
  const tdouble* pV = v.get_tmp_vptr_const();
  while (l>1) {
    const tuint ln = (l+1)/2;                // estimate the separation index
    if ( p <= pV[s+ln-1] ) {                // we have to search in the first part
      l = ln;
    } else {                                // we have to search in the second part
      s = s + ln;
      l -= ln;
    }
  }
  return s;  
}

const tuint FlxRndCreator::gen_smp_index2(const flxVec& v)
{
  return gen_smp_index2_help(gen_smp_uniform(),v);
}

const tuint FlxRndCreator::gen_smp_index2_help(const tdouble p, const flxVec& v)
{
  const tdouble pt = p*v.get_sum();
  tuint N = v.get_N();
  // pdouble sum;  -> too ineficient!!!
  tdouble sum = ZERO;
  for (tuint i=0;i<N;++i) {
    sum += v[i];
    if (pt<=sum) return i;
  }
  throw FlxException_Crude("FlxRndCreator::gen_smp_index2_help");
}


void FlxRndCreator::gen_smp(flxVec& y)
{
  if (rndReader == NULL) {
    if (rngp) {
      rv_normal(y,*rngp);
    } else {
      rv_normal(y);
    }
  } else {
    tuint L;
    if (!rndReader->get_vec(y,L)) {
      replay_stop();
      GlobalVar.alert.alert("FlxRndCreator::gen_smp","Replay of semi random values stopped.");
      for (tuint i = L; i<y.get_N(); ++i) {
        y[i] = gen_smp();
      }
    }
  }
}

void FlxRndCreator::shuffle(tuint* vp, const tuint N) {
  // initialize vector
  for (tuint i=0;i<N;++i) {
    vp[i] = i;
  }
  // random sampling
  for (tuint i=0;i<N;++i) {
    const tuint ip = gen_smp_index(N);
    if (ip!=i) {
      const tuint k = vp[N-i-1];
      vp[N-i-1] = vp[ip];
      vp[ip] = k;
    }
  }
}

void FlxRndCreator::latin_hypercube(tdouble *dp, const tuint N, const tuint M)
{
  std::valarray<tuint> ivec(N);
  const tdouble dx = ONE/N;
  // handle each dimension separately
  for (tuint i=0;i<M;++i) {
    shuffle(&ivec[0],N);
    for (tuint j=0;j<N;++j) {
      tdouble u = gen_smp_uniform();
      u += j;
      u*=dx;
      dp[ivec[j]*M+i] = u;
    }
  }
}

void FlxRndCreator::replay_stop(const bool doAlert)
{
  if (rndReader==NULL) {
    if (doAlert) GlobalVar.alert.alert("FlxRndCreator::replay_stop","Replay of semi random values already stopped.");
  }
  if (rndReader!=NULL) {
    rndReader = NULL;
    GlobalVar.slog(3) << "rnd track: stopped replay of semi random values." << std::endl;
  }
}
