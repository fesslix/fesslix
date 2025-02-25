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

#include "flxMemoryManager.h"
#include "flxexception.h"

#if FLX_DEBUG
  int FlxMemoryManager::Cinst = 0;
#endif

const tuint MEM_DBL_N = 1000;


FlxMemoryManager::FlxMemoryManager()
: cntnr_id(0)
{
  #if FLX_DEBUG
    ++Cinst;
    if (Cinst > 1) {
      std::ostringstream ssV;
      ssV << "More than one instance of 'FlxMemoryManager' created ...";
      throw FlxException("FlxMemoryManager::FlxMemoryManager", ssV.str() );
    }
  #endif
  // define the first container
    cntnr.push_back(new tdouble[MEM_DBL_N]);
    cntnr_entries.push_back(0);
}

FlxMemoryManager::~FlxMemoryManager()
{
  // cntnr
    const std::size_t N = cntnr.size();
    for (std::size_t i=0;i<N;++i) {
      delete [] cntnr[i];
    }
}

tdouble* FlxMemoryManager::new_double()
{
  const tuint N = cntnr.size();
  if ( cntnr_id>=N ) {                                // allocate new memory
    cntnr.push_back(new tdouble[MEM_DBL_N]);
    cntnr_entries.push_back(0);
  } else  if (cntnr_entries[cntnr_id] >= MEM_DBL_N) {        // current container is full
    ++cntnr_id;
    return new_double();
  }
  tuint &Nentries = cntnr_entries[cntnr_id];
  tdouble* res = &(cntnr[cntnr_id][Nentries]);
  ++Nentries;
  return res;
}

tdouble* FlxMemoryManager::new_vector(const tuint N)
{
  const tuint Nc = cntnr.size();
  tuint cid;        // current id of the container
  for (cid = cntnr_id;cid<Nc;++cid) {        // determine where enough storage is available
    if (cntnr_entries[cid]+N <= MEM_DBL_N) break;
  }
  if (cid>=Nc) {        // allocate new container
    const tuint cS = ((N<MEM_DBL_N)?MEM_DBL_N:N);
    cntnr.push_back(new tdouble[cS]);
    cntnr_entries.push_back(0);
    cid = Nc;
  } 
  tuint &Nentries = cntnr_entries[cid];
  tdouble* res = &(cntnr[cid][Nentries]);
  Nentries += N;
  return res;
}





