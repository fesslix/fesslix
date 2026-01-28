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

#include "flxphys.h"


const tdouble* magnus_k1 = NULL;  // = 6.112 hPa
const tdouble* magnus_k2 = NULL;  // = 17.62
const tdouble* magnus_k3 = NULL;  // = 243.12 °C

const tdouble t_absolut_ZERO = -273.15;

/**
* @brief Berechnet den Sättigungsdampfdruck abhängig von der Temperatur (temp)
*/
inline tdouble flx_fun_magnus_base(const tdouble temp) {
  const tdouble es =  (*magnus_k1) * exp( ((*magnus_k2) * temp) / ((*magnus_k3) + temp) );
  return es;  // in [hPa]
}

/**
* @brief Berechnet den Taupunkt abhängig von der Temperatur (temp) und relativer Luftfeuchtigkeit (phi)
*/
inline tdouble flx_fun_magnus(const tdouble temp, const tdouble phi) {
  if (phi<=ZERO) return t_absolut_ZERO;
  const tdouble lphi = log(phi);
  const tdouble help = (*magnus_k2)/(*magnus_k3+temp);
  return (*magnus_k3)*(help*temp+lphi)/(help*(*magnus_k3)-lphi);
}


// ********************************************
// OBSOLETE IMPLEMENTATION
// /**
// * @brief Datenkontainer für flx_fun_magnus_rh
// */
// struct flx_fun_magnus_data {
//   const tdouble d1;
//   const tdouble d2;
//   flx_fun_magnus_data(const tdouble d1, const tdouble d2) : d1(d1), d2(d2) {}
// };
//
// /**
// * @brief Hilfsfunktion um von Taupunkt und Luftfeuchtigkeit auf Temperatur zu schließen
// */
// inline tdouble flx_fun_magnus_rh(const tdouble temp, void *dp) {
//   flx_fun_magnus_data *p = (flx_fun_magnus_data *)dp;
//   return flx_fun_magnus(temp,p->d2)-p->d1;
// }
//
// /**
// * @brief Hilfsfunktion um von Taupunkt und Temperatur auf Luftfeuchtigkeit zu schließen
// */
// inline tdouble flx_fun_magnus_tp(const tdouble phi, void *dp) {
//   flx_fun_magnus_data *p = (flx_fun_magnus_data *)dp;
//   return flx_fun_magnus(p->d2,phi)-p->d1;
// }
// ********************************************

// ------------------------------------------------------------------------------------------------

void FlxCreateObjReaders_FlxPhys::createObjReaders(FlxObjectReadBox* objReadBox) {

}

void FlxCreateObjReaders_FlxPhys::createFunReaders(FlxData* dataBox)
{
  dataBox->FunBox.insert("phys_dewpoint", new FunReadPhys_dewpoint() );
  dataBox->FunBox.insert("phys_tauphi2temp", new FunReadPhys_tauphi2temp() );
  dataBox->FunBox.insert("phys_tauptemp2phi", new FunReadPhys_tautemp2phi() );

  // constants (https://de.wikipedia.org/wiki/S%C3%A4ttigungsdampfdruck#Magnus-Formel)
    magnus_k1 = dataBox->ConstantBox.insert("phys_magnus_k1",6.112);	// [hPa]
    magnus_k2 = dataBox->ConstantBox.insert("phys_magnus_k2",17.62);
    magnus_k3 = dataBox->ConstantBox.insert("phys_magnus_k3",243.12);	// [°C]
}

// ------------------------------------------------------------------------------------------------

const tdouble phys_temp2svp(const tdouble temp)
{
  return flx_fun_magnus_base(temp);
}

// ------------------------------------------------------------------------------------------------

const tdouble flxPhys_dewpoint(const tdouble temp, const tdouble phi)
{
  return flx_fun_magnus(temp,phi);
}

const tdouble FunPhys_dewpoint::calc()
{
  const tdouble temp = (*ParaList)[0]->calc();
  tdouble phi = (*ParaList)[1]->calc();
  return flx_fun_magnus(temp,phi);
}

FunBase* FunReadPhys_dewpoint::read(bool errSerious)
{
  return new FunPhys_dewpoint( read_parameters(2,errSerious) );
}

// ------------------------------------------------------------------------------------------------

const tdouble flxPhys_tauphi2temp(const tdouble temp_dewp, const tdouble phi)
{
  // ********************************************
  // OBSOLETE IMPLEMENTATION
  // flx_fun_magnus_data dp(temp_dewp,phi);
  // return flx_RootSearch_RegulaFalsi(&flx_fun_magnus_rh,&dp,ONE,30*ONE);
  // ********************************************
  const tdouble vp = flx_fun_magnus_base(temp_dewp);
  const tdouble gamma = log(vp/((*magnus_k1)*phi));
  return (*magnus_k3)*gamma/((*magnus_k2)-gamma);
}

const tdouble FunPhys_tauphi2temp::calc()
{
  return flxPhys_tauphi2temp((*ParaList)[0]->calc(),(*ParaList)[1]->calc());
}

FunBase* FunReadPhys_tauphi2temp::read(bool errSerious)
{
  return new FunPhys_tauphi2temp( read_parameters(2,errSerious) );
}

// ------------------------------------------------------------------------------------------------

const tdouble flxPhys_tautemp2phi(const tdouble temp_dewp, const tdouble temp)
{
  // ********************************************
  // OBSOLETE IMPLEMENTATION
  // flx_fun_magnus_data dp(temp_dewp,temp);
  // return flx_RootSearch_RegulaFalsi(&flx_fun_magnus_tp,&dp,ONE,30*ONE);
  // ********************************************
  const tdouble vp = flx_fun_magnus_base(temp_dewp);
  const tdouble svp = flx_fun_magnus_base(temp);
  return vp/svp;
}

const tdouble FunPhys_tautemp2phi::calc()
{
  return flxPhys_tautemp2phi((*ParaList)[0]->calc(),(*ParaList)[1]->calc());
}

FunBase* FunReadPhys_tautemp2phi::read(bool errSerious)
{
  return new FunPhys_tautemp2phi( read_parameters(2,errSerious) );
}

// ------------------------------------------------------------------------------------------------

/**
 * @brief Berechnet die absolute Feuchte [g/m³]
 *        aus Temperatur (temp, °C) und relativer Feuchte (phi, 0..1)
 */
const tdouble flxPhys_abs_humidity(const tdouble temp,  const tdouble phi)
{
    if (phi <= ZERO) return ZERO;

    // Sättigungsdampfdruck nach Magnus [hPa]
    const tdouble es =  flx_fun_magnus_base(temp);

    // Tatsächlicher Dampfdruck [hPa]
    const tdouble e = phi * es;

    // Absolute Feuchte [g/m³]
    return 216.7 * e / (temp + 273.15);
}


