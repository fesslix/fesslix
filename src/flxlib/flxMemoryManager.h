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

#ifndef fesslix_flxMemoryManager_H
#define fesslix_flxMemoryManager_H

#include "flxglobaldef.h"
#include <vector>


/**
* @brief manages memory that is not deleted anymore ... allocated memory is stored in blocks
*/
class FlxMemoryManager {
  private:
    #if FLX_DEBUG
      /**
      * @brief number of instances of this calss - only one class may exist
      */
      static int Cinst;
    #endif
    tuint cntnr_id;                // the current container id (all previous containers are already full)
    /**
    * @brief a container that contains tdouble values
    * a single entry contains (at least) MEM_DBL_N tdouble values
    */
    std::vector<tdouble*> cntnr;
    std::vector<tuint> cntnr_entries;                // the number of entries in the container
      
  public:
    FlxMemoryManager();
    ~FlxMemoryManager ();


    /**
    * @brief 'allocates' memory for a double entry
    */
    tdouble* new_double();
    /**
    * @brief 'allocates' memory for a double vector of size N
    */
    tdouble* new_vector(const tuint N);
};






#endif // fesslix_flxMemoryManager_H

