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

#ifndef fesslix_flxrandom_base_H
#define fesslix_flxrandom_base_H

#include "flxfunction_data.h"
#include "flxMtx.h"


/**
* @brief A class for managing the generation of random number
*/
class FLXLIB_EXPORT FlxRndCreator {
  private:
    #if FLX_DEBUG
      /**
      * @brief number of instances of this class - only one class may exist
      */
      static int Cinst;
    #endif
    /**
    * @brief random sequence input: reader
    */
    rng_type* rngp;
    FlxIstream* rndReader;
  public:
    FlxRndCreator(rng_type* rngp=NULL);
    ~FlxRndCreator() {}
    /**
    * @brief generate new random samples - standard normal distributed (uncorrelated)
    */
    const tdouble gen_smp();
    /**
    * @brief generate new random samples - uniform distributed on [0;1]
    */
    const tdouble gen_smp_uniform();
    /**
    * @brief generate new random samples - binary RV with Probability p to be true
    */
    const bool gen_smp_binary(const tdouble p) { return (gen_smp_uniform()<=p); }
    /**
    * @brief decide on an index of a CDF-vector
    * @param v a CDF-vector (probability that the index is smaller or equal than the current index)
    */
    const tuint gen_smp_index(const flxVec& v);
    /**
    * @brief decide on an index of a vector of probabilities
    * @param v a unnormalized vector of probability masses
    */
    const tuint gen_smp_index2(const flxVec& v);
      const tuint gen_smp_index2_help(const tdouble p, const flxVec& v);
    /**
    * @brief decide on an index of a vector if each compoent is equally likely to occur
    * @param N size of the vector -> index 0...N-1
    */
    const tuint gen_smp_index(const tuint N) { const tuint id(gen_smp_uniform()*N); return (id<N?id:N-1); }
    void gen_smp(flxVec& y);
    /**
    * @brief shuffle indices 1,2,...,N randomly around
    * @param N length of vp
    * @param vp pointer to result vector of length N
    */
    void shuffle(tuint* vp, const tuint N);
    /**
    * @brief generates N latin hypercube samples of dimension M, the samples are uniform between [0,1]
    *
    * @param dp double-pointer of size N*M
    * @param N number of samples
    * @param M dimension of parameter vector
    * @returns dp-matrix with M columns and N rows
    */
    void latin_hypercube(tdouble *dp, const tuint N, const tuint M);
    /**
    * @brief random sequence input: stop
    */
    void replay_stop(const bool doAlert=false);
    /**
    * @brief random sequence input: start
    */
    void replay_start(FlxIstream* rndReaderV);
};


class FLXLIB_EXPORT FlxBoxBaseR {
  protected:
    static FlxRndCreator* RndCreator;
    static GaussIntegration* GI;
  public:
    ~FlxBoxBaseR() {};
    static void set_Boxes( FlxRndCreator* RndCreatorV);
    static void set_GI(GaussIntegration* GIv) { GI = GIv;}
};




#endif // fesslix_flxrandom_base_H

