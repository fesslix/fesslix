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

#include "flxglobal.h"

#include <ostream>

const ostreamp* FlxAlert::cerr = NULL;


void FlxAlert::alert(const std::string& alertID, const std::string& alertStr) const
{
  *(*cerr) << std::endl << "ALERT (" << alertID << ")" << std::endl << "  " << alertStr << std::endl << std::endl;
  
  if (GlobalVar.check_logNOTcout()) {
    GlobalVar.slog(2) << std::endl << "ALERT (" << alertID << ")" << std::endl << "  " << alertStr << std::endl << std::endl;
  }
}

const std::string FlxException::what()
{
  std::ostringstream ssV;
  ssV << "ERROR - an error occurred (" << errNumber << ")" << std::endl << "\t" << titel << std::endl << "\t" << msg << std::endl;
  if (GlobalVar.prelog_isNOTempty()) {
    ssV << "Last parsed input:" << std::endl << GlobalVar.prelog_force_write() << std::endl;
  }
  return ssV.str();
}



