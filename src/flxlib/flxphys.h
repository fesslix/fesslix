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

#ifndef fesslix_phys_H
#define fesslix_phys_H

#include "flxobjects.h"



class FLXLIB_EXPORT FlxCreateObjReaders_FlxPhys : public FlxCreateObjReaders {
  public:
    void createObjReaders (FlxObjectReadBox* objReadBox );
    void createFunReaders (FlxData* dataBox );
};


// -------------------------- saturation vapor pressure ----------------------------------------------------------------------

FLXLIB_EXPORT const tdouble phys_temp2svp(const tdouble temp);

// -------------------------- phys_dewpoint ----------------------------------------------------------------------

FLXLIB_EXPORT const tdouble flxPhys_dewpoint(const tdouble temp, const tdouble phi);

class FunPhys_dewpoint : public FunBaseFun_multPara {
  public:
    FunPhys_dewpoint(std::vector<FunBase*> *ParaListV) : FunBaseFun_multPara(ParaListV) {}
    const tdouble calc();
    const std::string write_v() { return "phys_dewpoint"; }
};

class FunReadPhys_dewpoint : public FunReadFunBase {
  public:
    FunReadPhys_dewpoint() :  FunReadFunBase( true ) {}
    FunBase* read ( bool errSerious );
};

// -------------------------- FunPhys_tauphi2temp ----------------------------------------------------------------------

FLXLIB_EXPORT const tdouble flxPhys_tauphi2temp(const tdouble temp_dewp, const tdouble phi);

class FunPhys_tauphi2temp : public FunBaseFun_multPara {
  public:
    FunPhys_tauphi2temp(std::vector<FunBase*> *ParaListV) : FunBaseFun_multPara(ParaListV) {}
    const tdouble calc();
    const std::string write_v() { return "phys_tauphi2temp"; }
};

class FunReadPhys_tauphi2temp : public FunReadFunBase {
  public:
    FunReadPhys_tauphi2temp() :  FunReadFunBase( true ) {}
    FunBase* read ( bool errSerious );
};

// -------------------------- FunPhys_tautemp2phi ----------------------------------------------------------------------

FLXLIB_EXPORT const tdouble flxPhys_tautemp2phi(const tdouble temp_dewp, const tdouble temp);

class FunPhys_tautemp2phi : public FunBaseFun_multPara {
  public:
    FunPhys_tautemp2phi(std::vector<FunBase*> *ParaListV) : FunBaseFun_multPara(ParaListV) {}
    const tdouble calc();
    const std::string write_v() { return "phys_tautemp2phi"; }
};

class FunReadPhys_tautemp2phi : public FunReadFunBase {
  public:
    FunReadPhys_tautemp2phi() :  FunReadFunBase( true ) {}
    FunBase* read ( bool errSerious );
};

// -------------------------- absolut humidity ----------------------------------------------------------------------

FLXLIB_EXPORT const tdouble flxPhys_abs_humidity(const tdouble temp,  const tdouble phi);


#endif // fesslix_phys_H

