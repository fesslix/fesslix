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
#include "flxrbrv_rvs_read.h"




FlxRndSamplingSpace_base::FlxRndSamplingSpace_base(RBRV_constructor& RndBox)
: RndBox(RndBox), DIM(RndBox.get_NRV()), y_rem(DIM), z_rem(DIM)
{

}

void FlxRndSamplingSpace_base::gen_smp(flxVec& y)
{
  RndBox.propose_y(y);
}

void FlxRndSamplingSpace_base::calc_foverh(tdouble& foverh, const flxVec& z)
{
  // get h
    tdouble h;
    get_h(h,z);
  // get f
    tdouble f = ONE;
    for (size_t i=0;i<DIM;++i) {
      f*=rv_phi(z[i]);
    }
  foverh = f/h;
}

void FlxRndSamplingSpace_base::transform_space(tdouble& foverh)
{
  // generate uncorrelated standard normal samples
    gen_smp(y_rem);
  // get vector in z-space (sampling space)
    y2z(y_rem,z_rem);
  // calc foverh
    calc_foverh(foverh,z_rem);
  // set z as new y-vector
    RndBox.set_smp(z_rem);
}

FlxRndSamplingSpace_uni::~FlxRndSamplingSpace_uni()
{
  delete rv;
}

void FlxRndSamplingSpace_uni::y2z(flxVec& y, flxVec& z)
{
  rv->eval_para();
  for (tuint i=0;i<z.get_N();++i) {
    z[i] = rv->transform_y2x(y[i]);
  }
}

void FlxRndSamplingSpace_uni::get_h(tdouble& h, const flxVec& z) const
{
  h = ONE;
  for (tuint i = 0; i<z.get_N();++i) {
    h*=rv->calc_pdf_x(z[i]);
  }
}

void FlxRndSamplingSpace_uni::print_info(std::ostream& sout, const bool verbose)
{
  sout << "uniform " ;
  tuint c = 0;
  rv->eval_para();
  rv->print(sout,"",c,false);
}

FlxRndSamplingSpace_Generator_base* FlxRndSamplingSpace_Generator_base::createSS(FlxRndSamplingSpace_base::sstype sst, bool errSerious)
{
  switch(sst) {
    case FlxRndSamplingSpace_base::ssuni:
      return new FlxRndSamplingSpace_Generator_Uni(errSerious);
      break;
    case FlxRndSamplingSpace_base::ssnormal:
      return new FlxRndSamplingSpace_Generator_Normal(errSerious);
      break;
    case FlxRndSamplingSpace_base::sstail:
      return new FlxRndSamplingSpace_Generator_TailStdN(errSerious);
      break;
    default:
      std::ostringstream ssV;
      ssV << "ERROR.";
      throw FlxException("FlxRndSamplingSpace_Generator_base::createSS", ssV.str() );
  }
}

std::string FlxRndSamplingSpace_Generator_base::get_rvt(FlxRndSamplingSpace_base::sstype sst)
{
  switch(sst) {
    case FlxRndSamplingSpace_base::ssuni:
      return "uni";
      break;
    case FlxRndSamplingSpace_base::ssnormal:
      return "normal";
      break;
    case FlxRndSamplingSpace_base::sstail:
      return "tailstdn";
      break;
    default:
      return "";
  }
}

FlxRndSamplingSpace_base::sstype FlxRndSamplingSpace_Generator_base::get_sst(std::string strsst, bool errSerious)
{
  std::transform(strsst.begin(), strsst.end(), strsst.begin(), (int(*)(int)) std::tolower);
  if (strsst == "uni" ) {
    return FlxRndSamplingSpace_base::ssuni;
  } else if (strsst == "normal" ) {
    return FlxRndSamplingSpace_base::ssnormal;
  } else if (strsst == "tailstdn" ) {
    return  FlxRndSamplingSpace_base::sstail;
  } else {
    std::ostringstream ssV;
    ssV << "Unkown type of sampling space '" << strsst << "'.";
    FlxError(errSerious,"FlxRndSamplingSpace_Generator_base::get_sst", ssV.str() ); return FlxRndSamplingSpace_base::ssuni; // DUMMY-RETURN
  }
}

FlxRndSamplingSpace_Generator_Uni::FlxRndSamplingSpace_Generator_Uni(bool errSerious)
{
  rvG = RBRV_entry_read_base::read_entry(false,true);
}

FlxRndSamplingSpace_Generator_Uni::~FlxRndSamplingSpace_Generator_Uni()
{
  delete rvG;
}

FlxRndSamplingSpace_base* FlxRndSamplingSpace_Generator_Uni::generate_SS(RBRV_constructor& RndBox)
{
  tuint riID = 0;
  RBRV_entry* rvtmp0 = rvG->generate_entry("dummy",riID);
  RBRV_entry_RV_base* rvtmp = dynamic_cast<RBRV_entry_RV_base*> (rvtmp0);
  if (rvtmp==NULL) {
    throw FlxException("FlxRndSamplingSpace_Generator_Uni::generate_SS","The specified random variable cannot be sampled from directly.");
  }
  return new FlxRndSamplingSpace_uni( rvtmp, RndBox );
}

FlxRndSamplingSpace_Generator_TailStdN::FlxRndSamplingSpace_Generator_TailStdN(bool errSerious)
{
  reader->getWord("betadp",false);
  betaDP_F = new FlxFunction(funReader,errSerious);
}

FlxRndSamplingSpace_base* FlxRndSamplingSpace_Generator_TailStdN::generate_SS(RBRV_constructor& RndBox)
{
  return new FlxRndSamplingSpace_TailStdN(betaDP_F->calc(),RndBox);
}

FlxRndSamplingSpace_base* FlxRndSamplingSpace_Generator_Normal::generate_SS(RBRV_constructor& RndBox)
{
  // prepare mu
    std::string str = muMF->eval();
    tuint N = 0;
    const flxVec mu_V(ConstMtxBox->get_Vec(str,N),N);
  // prepare sigma
    str = sigmaMF->eval();
    N = 0;
    const flxVec sigma_V(ConstMtxBox->get_Vec(str,N),N);
  return new FlxRndSamplingSpace_normal(mu_V,sigma_V,(betaTrunc_F?betaTrunc_F->calc():-ONE),NinitScale->cast2tulong(),RndBox);  
}

FlxRndSamplingSpace_Generator_Normal::FlxRndSamplingSpace_Generator_Normal(bool errSerious)
{
  reader->getWord("mu",errSerious);
  muMF = new FlxMtxConstFun(true);
  sigmaMF = NULL;
  betaTrunc_F = NULL;
  NinitScale = NULL;
  try {
    reader->getWord("sd",errSerious);
    sigmaMF = new FlxMtxConstFun(true);
    if (reader->whatIsNextString(9,true)=="betatrunc") {
      reader->getWord("betatrunc",errSerious);
      betaTrunc_F = new FlxFunction(funReader,errSerious);
      if (reader->whatIsNextString(5,true)=="ninit") {
        reader->getWord("ninit",errSerious);
        NinitScale = new FlxFunction(funReader,errSerious);
      }
    }
  } catch (FlxException &e) {
    FLXMSG("FlxRndSamplingSpace_Generator_Normal::FlxRndSamplingSpace_Generator_Normal",1);
    delete muMF;
    if (sigmaMF) delete sigmaMF;
    if (betaTrunc_F) delete betaTrunc_F;
    if (NinitScale) delete NinitScale;
    throw;
  }
  if (NinitScale==NULL) {
    NinitScale = new FlxFunction(new FunNumber(1e6));
  }
}

void FlxRndSamplingSpace_normal::print_info(std::ostream& sout, const bool verbose)
{
  sout << "Normal";
  if (verbose) {
    sout << " - mean=" << mu << "; sigma=" << sigma;
    if (betaTrunc>ZERO) {
      sout << "; betaTrunc=" << GlobalVar.Double2String(betaTrunc) << " (p=" << GlobalVar.Double2String(ONE/h_scale) << " with ninit=" << GlobalVar.Double2String(nInit) << ")";
    }
  }
}

FlxRndSamplingSpace_normal::FlxRndSamplingSpace_normal(const flxVec& mu, const flxVec& sigma, const tdouble betaTrunc, const tulong nInit, RBRV_constructor& RndBox)
: FlxRndSamplingSpace_base(RndBox), mu(mu),sigma(sigma),betaTrunc(betaTrunc),nInit(nInit)
{
  if (mu.get_N() != sigma.get_N()) {
    std::ostringstream ssV;
    ssV << "Vector sizes do not match.";
    throw FlxException("FlxRndSamplingSpace_normal::FlxRndSamplingSpace_normal_1", ssV.str() );
  }
  if (DIM != mu.get_N() ) {
    std::ostringstream ssV;
    ssV << "Size of vectors 'mu' and 'sigma' (" << DIM << ") out of range (" << mu.get_N() << ").";
    throw FlxException("FlxRndSamplingSpace_normal::FlxRndSamplingSpace_normal_2", ssV.str() );
  }
  if (betaTrunc<=ZERO) {
    h_scale = ONE;
  } else {
    std::size_t DIM = RndBox.get_NRV();
    // generate uncorrelated standard normal samples
      const tdouble* yp = y_rem.get_tmp_vptr();
    tulong Nout = 0;
    for (tulong i=0;i<nInit;++i) {
      gen_smp(y_rem);
      tdouble l = ZERO;
      for (std::size_t i=0;i<DIM;++i) {
        l += pow2(yp[i]*sigma[i]+mu[i]);
      }
      l = sqrt(l);
      if (l>betaTrunc) ++Nout;      
    }
    h_scale = tdouble(nInit)/tdouble(Nout);
  }
}

void FlxRndSamplingSpace_normal::y2z(flxVec& y, flxVec& z)
{
  do {
    tdouble l = ZERO;
    for (size_t i=0;i<z.get_N();++i) {
      z[i] = y[i]*sigma[i]+mu[i];
      l += pow2(z[i]);
    }
    l = sqrt(l);
    if (l<betaTrunc) gen_smp(y);
    else break;
  } while (true);
}

void FlxRndSamplingSpace_normal::get_h(tdouble& h, const flxVec& z) const
{
  h = h_scale;
  for (tuint i = 0; i<z.get_N();++i) {
    h*=rv_phi((z[i]-mu[i])/sigma[i])/sigma[i];
  }
}

void FlxRndSamplingSpace_TailStdN::calc_foverh(tdouble& foverh, const flxVec& z)
{
  // calc f/h
    foverh=cF;
}

FlxRndSamplingSpace_TailStdN::FlxRndSamplingSpace_TailStdN(const tdouble betaDP, RBRV_constructor& RndBox)
: FlxRndSamplingSpace_base(RndBox), betaDP2(pow2(betaDP)),cF(ONE-rv_cdf_ChiSquare(DIM,pow2(betaDP)))
{
  if (betaDP < ZERO) {
    std::ostringstream ssV;
    ssV << "'betaDP' (" << betaDP << ") has to be larger than zero.";
    throw FlxException("FlxRndSamplingSpace_TailStdN::FlxRndSamplingSpace_TailStdN", ssV.str() );
  }
}

void FlxRndSamplingSpace_TailStdN::print_info(std::ostream& sout, const bool verbose)
{
  sout << "Tails of the Normal - Normal defined outside the interval [-" << sqrt(betaDP2) << "; " << sqrt(betaDP2) << "]";
  if (verbose) {
    sout << std::endl << "   to get the same accuracy with Monte Carlo, you would need " << GlobalVar.Double2String(ONE/cF) << " times more samples.";
  }
}

void FlxRndSamplingSpace_TailStdN::y2z(flxVec& y, flxVec& z)
{
  const tdouble ly2 = y.get_Norm2_NOroot();
  const tdouble lz2 = rv_InvCDF_ChiSquare(DIM,ONE-cF*rv_cdf_ChiSquare(DIM,ly2));
  z = y;
  z *= sqrt(lz2/ly2);
}

FlxRndKernel_base* FlxRndKernel_base::read(const bool errSerious)
{
  return FlxRndKernel_base::read(reader->getWord(true,errSerious),errSerious);
}

FlxRndKernel_base* FlxRndKernel_base::read(const std::string& kernelStr, const bool errSerious)
{
  if (kernelStr == "gauss") {
    return new FlxRndKernel_Gauss();
  } else if (kernelStr == "uniform") {
    return new FlxRndKernel_Uniform();
  } else {
    std::ostringstream ssV;
    ssV << "Unknown keyword '" << kernelStr << "'.";
    FlxError(errSerious,"FlxRndKernel_base::read", ssV.str() ); return NULL; // DUMMY-RETURN
  }
}

void FlxRndKernel_base::set_h(const tdouble h_value)
{
  if (h_value <= ZERO) {
    std::ostringstream ssV;
    ssV << "Bandwidth has to be positive.";
    throw FlxException("FlxRndKernel_base::set_h", ssV.str() );
  } else {
    h = h_value;
    hinv = ONE/h_value;
  }    
}

const tdouble FlxRndKernel_Uniform::cdf(const tdouble x, const tdouble xs) const
{
  const tdouble xt = (x-xs)*hinv;
  if (xt<=-ONE) return ZERO;
  if (xt>=ONE) return ONE;
  return (xt+ONE)/2;
}


