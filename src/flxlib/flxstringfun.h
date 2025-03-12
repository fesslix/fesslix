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

#include "flxfunction_ext.h"


/**
* @brief The base class for all other arithmetic classes.
*/
class StringFunBase {
  public:
    virtual void eval(std::ostream& os) = 0;
    virtual ~StringFunBase () {};
    virtual const std::string write() = 0;
};

class FunReadFlxStringFunBase : public FlxReaderBase, public FlxBoxBase {
  protected:
  public:
    virtual StringFunBase* read ( bool errSerious ) = 0;
    virtual ~FunReadFlxStringFunBase() {};
};



/**
* @brief A class for storing pre-defined flxString-functions
*/
class FLXLIB_EXPORT FlxStringFunBox {
  private:
    #if FLX_DEBUG
      /**
      * @brief number of instances of this calss - only one class may exist
      */
      static int Cinst;
    #endif
    std::map<std::string, FunReadFlxStringFunBase*> box;
  public:
    FlxStringFunBox ();
    ~FlxStringFunBox ();
    void insert ( const std::string& name, FunReadFlxStringFunBase* fR);
    FunReadFlxStringFunBase* get ( const std::string& name) ;
    StringFunBase* read(ReadStream* reader, const bool errSerious);
};




