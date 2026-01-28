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

#ifndef fesslix_flxfunction_ope_read_H
#define fesslix_flxfunction_ope_read_H

#include "flxfunction_ope_calc.h"


/**
* @brief arithmetic read class for 'words' (variables, funnctions, constants)
*/
class FunReadWord : public FunReadBase {
  private:
    FunBase* read ( const bool errSerious );
  public:
    /**
    * @brief Create a FunReadWord class.
    *
    * An instance of each arithmetic function read class has to be created !!!
    *
    */
    FunReadWord () : FunReadBase(1010) {};
};

/**
* @brief considers brackets in arithmetic functions
*/
class FunReadBracket : public FunReadBase {
  private:
    FunBase* read ( bool errSerious );
  public:
    FunReadBracket () : FunReadBase(1100) {}
};

/**
* @brief arithmetic read class of FunAdd
*/
class FunReadAdd : public FunReadBase {
  private:
    FunBase* read ( bool errSerious );
  public:
    FunReadAdd () : FunReadBase(40) {}
};

/**
* @brief arithmetic read class of FunMult
*/
class FunReadMult : public FunReadBase {
  private:
    FunBase* read ( bool errSerious );
  public:
    FunReadMult () : FunReadBase(50) {}
};

/**
* @brief arithmetic read class of FunPower
*/
class FunReadPower : public FunReadBase {
  private:
    FunBase* read ( bool errSerious );
  public:
    FunReadPower () : FunReadBase(60) {}
};

/**
* @brief arithmetic read class of FunLessThan
*/
class FunReadLessThan : public FunReadBase {
  private:
    FunBase* read ( bool errSerious );
  public:
    FunReadLessThan () : FunReadBase(31) {}
};

/**
* @brief arithmetic read class of FunAnd
*/
class FunReadAnd : public FunReadBase {
  private:
    FunBase* read ( bool errSerious );
  public:
    FunReadAnd () : FunReadBase(20) {}
};

/**
* @brief arithmetic read class of FunAnd
*/
class FunReadOr : public FunReadBase {
  private:
    FunBase* read ( bool errSerious );
  public:
    FunReadOr () : FunReadBase(19) {}
};

/**
* @brief arithmetic read class of FunNot
*/
class FunReadNot : public FunReadBase {
  private:
    FunBase* read ( bool errSerious );
  public:
    FunReadNot () : FunReadBase(900) {}
};

/**
* @brief arithmetic read class of FunTernary
*/
class FunReadTernary : public FunReadBase {
  private:
    FunBase* read ( bool errSerious );
  public:
    FunReadTernary () : FunReadBase(10) {}
};

/**
* @brief arithmetic read class of FunEqual
*/
class FunReadEqual : public FunReadBase {
  private:
    FunBase* read ( bool errSerious );
  public:
    FunReadEqual () : FunReadBase(30) {}
};

#endif // fesslix_flxfunction_ope_read_H
