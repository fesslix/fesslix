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

#include "flxmtxfun.h"


//======================= Preliminary Definitions ============================================

class FLXLIB_EXPORT FlxBoxBase2 : public FlxBoxBase, public FlxMtxBoxBase {
  public:
    virtual ~FlxBoxBase2() {};
};




//======================= Sampling Spaces ============================================

class RBRV_constructor;
class FLXLIB_EXPORT FlxRndSamplingSpace_base {
  protected:
    RBRV_constructor& RndBox;
    const tuint DIM;
    flxVec y_rem;
    flxVec z_rem;
    
    /**
    * @brief generates samples in uncorrelated standard normal space
    */
    void gen_smp(flxVec& y);
    /**
    * @brief transform from uncorrelated standard normal space to uncorrelated sampling space
    */
    virtual void y2z(flxVec& y, flxVec& z) = 0;
    /**
    * @brief calculates the h-coefficient for Importance Sampling
    */
    virtual void get_h(tdouble &h, const flxVec& z) const = 0;
    /**
    * @brief calculates the f/h-factor for Importance Sampling
    */
    virtual void calc_foverh(tdouble& foverh, const flxVec& z);
    
  public:
    FlxRndSamplingSpace_base(RBRV_constructor& RndBox);
    virtual ~FlxRndSamplingSpace_base() {}
    /**
    * @brief generates y-realizations in the sampling sapce
    */
    void transform_space(tdouble& foverh);
    virtual void print_info(std::ostream &sout, const bool verbose) = 0;
    virtual const tdouble rely_confidence() const { return -ONE; }
    // Types of SamplingSpaces
      enum sstype { 
        ssuni, ssnormal,sstail
      };
};

class FLXLIB_EXPORT FlxRndSamplingSpace_Generator_base : public FlxReaderBase2, protected FlxBoxBase2 {
  public:
    virtual ~FlxRndSamplingSpace_Generator_base() {};
    /**
    * @brief create an instance of a sampling space
    */
    virtual FlxRndSamplingSpace_base* generate_SS(RBRV_constructor& RndBox) = 0;
    
    // Types of sampling spaces
      static FlxRndSamplingSpace_base::sstype get_sst(std::string strsst, bool errSerious=true);
      static std::string get_rvt(FlxRndSamplingSpace_base::sstype sst);
    static FlxRndSamplingSpace_Generator_base* createSS(FlxRndSamplingSpace_base::sstype sst, bool errSerious=true);
};

class FLXLIB_EXPORT FlxRndSamplingSpace_normal : public FlxRndSamplingSpace_base {
  private:
    const flxVec mu;
    const flxVec sigma;
    const tdouble betaTrunc;
    tdouble h_scale;
    const tulong nInit;
  protected:
    void y2z(flxVec& y, flxVec& z);
    void get_h(tdouble &h, const flxVec& z) const;
  public:
    FlxRndSamplingSpace_normal(const flxVec& mu, const flxVec& sigma, const tdouble betaTrunc, const tulong nInit, RBRV_constructor& RndBox);
    void print_info(std::ostream &sout, const bool verbose);
};

class FLXLIB_EXPORT FlxRndSamplingSpace_Generator_Normal : public FlxRndSamplingSpace_Generator_base {
  private:
    FlxMtxConstFun* muMF;
    FlxMtxConstFun* sigmaMF;
    FlxFunction* betaTrunc_F;
    FlxFunction* NinitScale;
  public:
    FlxRndSamplingSpace_Generator_Normal(bool errSerious);
    ~FlxRndSamplingSpace_Generator_Normal() { delete muMF; delete sigmaMF; if (betaTrunc_F) delete betaTrunc_F; delete NinitScale; }
    FlxRndSamplingSpace_base* generate_SS(RBRV_constructor& RndBox);
};

class FLXLIB_EXPORT FlxRndSamplingSpace_TailStdN : public FlxRndSamplingSpace_base {
  private:
    const tdouble betaDP2;
    tdouble cF;
  protected:
    void y2z(flxVec& y, flxVec& z);
    void get_h(tdouble &h, const flxVec& z) const {};  // not explicitely needed
    void calc_foverh(tdouble& foverh, const flxVec& z);
  public:
    FlxRndSamplingSpace_TailStdN(const tdouble betaDP, RBRV_constructor& RndBox);
    void print_info(std::ostream &sout, const bool verbose); 
    const tdouble rely_confidence() const { return (ONE/cF); }
};

class FLXLIB_EXPORT FlxRndSamplingSpace_Generator_TailStdN : public FlxRndSamplingSpace_Generator_base {
  private:
    FlxFunction* betaDP_F;
  public:
    FlxRndSamplingSpace_Generator_TailStdN(bool errSerious);
    ~FlxRndSamplingSpace_Generator_TailStdN() { delete betaDP_F; }
    FlxRndSamplingSpace_base* generate_SS(RBRV_constructor& RndBox);
};

class RBRV_entry_RV_base;
class FLXLIB_EXPORT FlxRndSamplingSpace_uni : public FlxRndSamplingSpace_base {
  private:
    RBRV_entry_RV_base* const rv;
  public:
    FlxRndSamplingSpace_uni(RBRV_entry_RV_base* rv, RBRV_constructor& RndBox) : FlxRndSamplingSpace_base(RndBox), rv(rv) {}
    ~FlxRndSamplingSpace_uni();
    void y2z(flxVec& y, flxVec& z);
    void get_h(tdouble &h, const flxVec& z) const;
    void print_info(std::ostream &sout, const bool verbose);
};

class RBRV_entry_read_base;
class FLXLIB_EXPORT FlxRndSamplingSpace_Generator_Uni : public FlxRndSamplingSpace_Generator_base {
  private:
    RBRV_entry_read_base* rvG;
  public:
    FlxRndSamplingSpace_Generator_Uni(bool errSerious);
    ~FlxRndSamplingSpace_Generator_Uni();
    FlxRndSamplingSpace_base* generate_SS(RBRV_constructor& RndBox);
};


//======================= KERNELS ============================================

class FLXLIB_EXPORT FlxRndKernel_base : public FlxReaderBase2 {
  protected:
    tdouble h;
    tdouble hinv;        // bandwidth of the kernel
  public:
    FlxRndKernel_base() : h(ONE), hinv(ONE) {}
    virtual ~FlxRndKernel_base() {};
    
    /**
    * @brief set the bandwidth of the kernel
    */
    void set_h(const tdouble h_value);
    void set_h_fast(const tdouble h_value) { h = h_value; hinv = ONE/h_value; }
    const tdouble get_h() const { return h;}
    /**
    * @brief kernel-PDF centered at xs evaluated at x
    */
    virtual const tdouble pdf(const tdouble x, const tdouble xs) const = 0;
    /**
    * @brief kernel-PDF centered at xs evaluated at x
    */
    virtual const tdouble cdf(const tdouble x, const tdouble xs) const = 0;
    /**
    * @brief transfrom from standard normal space to kernel-space (centered around zero)
    */
    virtual const tdouble transform_y2x(const tdouble y) const = 0;
    virtual void print_info(std::ostream& sout) const = 0;
    
    static FlxRndKernel_base* read(const bool errSerious=true);
    static FlxRndKernel_base* read(const std::string& kernelStr, const bool errSerious=true);
};

class FLXLIB_EXPORT FlxRndKernel_Gauss : public FlxRndKernel_base {
  public:
    const tdouble pdf(const tdouble x, const tdouble xs) const { return hinv*rv_phi((x-xs)*hinv); }
    const tdouble cdf(const tdouble x, const tdouble xs) const { return rv_Phi((x-xs)*hinv); }
    const tdouble transform_y2x(const tdouble y) const { return y*h; }
    void print_info(std::ostream& sout) const { sout << "Gaussian kernel; h=" << GlobalVar.Double2String(h); }
};

class FLXLIB_EXPORT FlxRndKernel_Uniform : public FlxRndKernel_base {
  public:
    const tdouble pdf(const tdouble x, const tdouble xs) const { return hinv/2; }
    const tdouble cdf(const tdouble x, const tdouble xs) const;
    const tdouble transform_y2x(const tdouble y) const { return ((y>0)?(ONE-rv_Phi(-y)*2):(rv_Phi(y)*2-ONE))*h; }
    void print_info(std::ostream& sout) const { sout << "uniform kernel; h=" << GlobalVar.Double2String(h); }
};

