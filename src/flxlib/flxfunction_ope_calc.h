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

#ifndef fesslix_flxfunction_ope_calc_H
#define fesslix_flxfunction_ope_calc_H

#include "flxfunction_ext.h"


/**
* @brief arithmetic class: stores a constant from the ConstantBox
*/
class FLXLIB_EXPORT FunConst : public FunBase, public FlxBoxBase {
  private:
    const tdouble* thenumber;
  public:
    FunConst (const tdouble* TheNumber) : thenumber(TheNumber) {};
    const tdouble calc() { return *thenumber; }
    const bool search_circref(FlxFunction* fcr) { return false; }
    const tuint precedence() { return 0; }
    const std::string write();
    const bool dependOn_Const(const tdouble* const theconst) { return thenumber==theconst; }
    const bool optimize(FunBasePtr& optf, const Fun_OptimizeInfo &foi);
};

/**
* @brief arithmetic class: stores a arithmetic expression from the VarBox
*/
class FLXLIB_EXPORT FunVar : public FunBase, public FlxBoxBase {
  private:
    FlxFunction* thefun;
  public:
    FunVar (FlxFunction* TheFun) : thefun(TheFun) {};
    const tdouble calc() { return thefun->calc(); }
    const bool search_circref(FlxFunction* fcr);
    const tuint precedence() { return 0; }
    const std::string write();
    const bool dependOn_Const(const tdouble* const thenumber) { return thefun->dependOn_Const(thenumber); }
    const bool optimize(FunBasePtr& optf, const Fun_OptimizeInfo &foi);
};

/**
* @brief arithmetic class: add two values
*/
class FLXLIB_EXPORT FunAdd : public FunBaseOperat2 {
  public:
    FunAdd (FunBase *child1, FunBase *child2) : FunBaseOperat2(child1, child2) {};
    const tdouble calc() { 
      const tdouble d1 = child_1->calc(); const tdouble d2 = child_2->calc(); return d1 + d2; }
    const std::string write();
    const tuint precedence() { return 12; }
    const bool optimize(FunBasePtr& optf, const Fun_OptimizeInfo &foi);
};

/**
* @brief arithmetic class: subtract two values
*/
class FLXLIB_EXPORT FunSub : public FunBaseOperat2 {
  public:
    FunSub (FunBase *child1, FunBase *child2) : FunBaseOperat2(child1, child2) {};
    const tdouble calc() { const tdouble d1 = child_1->calc(); const tdouble d2 = child_2->calc(); return d1 - d2; }
    const std::string write();
    const tuint precedence() { return 12; }
    const bool optimize(FunBasePtr& optf, const Fun_OptimizeInfo &foi);
};

/**
* @brief arithmetic class: Negation
*/
class FLXLIB_EXPORT FunNot : public FunBaseOperat1 {
  public:
    FunNot (FunBase *child1) : FunBaseOperat1(child1) {};
    const tdouble calc() { 
      if ( fabs(child_1->calc())<=GlobalVar.TOL() ) { return ONE; }
      else { return ZERO; }
    }
    const std::string write();
    const tuint precedence() { return 0; }
};

/**
* @brief arithmetic class: multiply two values
*/
class FLXLIB_EXPORT FunMult : public FunBaseOperat2 {
  public:
    FunMult (FunBase *child1, FunBase *child2) : FunBaseOperat2(child1, child2) {};
    const tdouble calc() { const tdouble d1 = child_1->calc(); const tdouble d2 = child_2->calc(); return d1*d2; }
    const std::string write();
    const tuint precedence() { return 11; }
    const bool optimize(FunBasePtr& optf, const Fun_OptimizeInfo &foi);
};

/**
* @brief arithmetic class: divide two values
*/
class FLXLIB_EXPORT FunMult_Div : public FunBaseOperat2 {
  public:
    FunMult_Div (FunBase *child1, FunBase *child2) : FunBaseOperat2(child1, child2) {};
    const tdouble calc() {  const tdouble d1 = child_1->calc(); const tdouble d2 = child_2->calc(); return d1 / d2; }
    const std::string write();
    const tuint precedence() { return 11; }
    const bool optimize(FunBasePtr& optf, const Fun_OptimizeInfo &foi);
};

/**
* @brief arithmetic class: raise a value to a given power
*/
class FLXLIB_EXPORT FunPower : public FunBaseOperat2 {
  public:
    FunPower (FunBase *child1, FunBase *child2) : FunBaseOperat2(child1, child2) {};
    const tdouble calc() { const tdouble d1 = child_1->calc(); const tdouble d2 = child_2->calc(); return std::pow (d1, d2); }
    const std::string write();
    const tuint precedence() { return 10; }
    const bool optimize(FunBasePtr& optf, const Fun_OptimizeInfo &foi);
};

/**
* @brief arithmetic class: compares two values ((not) equal to)
*/
class FLXLIB_EXPORT FunEqual : public FunBaseOperat2 {
  private:
    const bool isEqual;
  public:
    FunEqual (FunBase *child1, FunBase *child2, const bool IsEqual) : FunBaseOperat2(child1, child2), isEqual(IsEqual) {};
    const tdouble calc() { 
      const tdouble d1 = child_1->calc(); const tdouble d2 = child_2->calc();
      if ( d1 == d2 ) { return isEqual; }
      else { return !isEqual; }
    }
    const std::string write();
    const tuint precedence() { return 14; }
};

/**
* @brief arithmetic class: compares two values (<, <=, >=, >)
*/
class FLXLIB_EXPORT FunLessThan : public FunBaseOperat2 {
  private:
    bool isEqual;
    bool isLess;
  public:
    FunLessThan (FunBase *child1, FunBase *child2, bool IsEqual, bool IsLess) : FunBaseOperat2(child1, child2), isEqual(IsEqual), isLess(IsLess) {};
    const tdouble calc();
    const std::string write();
    const tuint precedence() { return 13; }
};

/**
* @brief arithmetic class: Logical AND
*/
class FLXLIB_EXPORT FunAnd : public FunBaseOperat2 {
  public:
    FunAnd (FunBase *child1, FunBase *child2) : FunBaseOperat2(child1, child2) {};
    const tdouble calc() { 
      if (!(fabs(child_1->calc()) > GlobalVar.TOL())) return false;
      return (fabs(child_2->calc()) > GlobalVar.TOL());
    }
    const std::string write();
    const tuint precedence() { return 15; }
};

/**
* @brief arithmetic class: Logical AND
*/
class FLXLIB_EXPORT FunOr : public FunBaseOperat2 {
  public:
    FunOr (FunBase *child1, FunBase *child2) : FunBaseOperat2(child1, child2) {};
    const tdouble calc() { 
      if (fabs(child_1->calc()) > GlobalVar.TOL()) return true;
      return (fabs(child_2->calc()) > GlobalVar.TOL());
    }
    const std::string write();
    const tuint precedence() { return 16; }
};

/**
* @brief arithmetic class: ternary operator
*/
class FLXLIB_EXPORT FunTernary : public FunBaseOperat3 {
  public:
    FunTernary (FunBase *child1, FunBase *child2, FunBase *child3) : FunBaseOperat3(child1, child2, child3) {};
    const tdouble calc() { 
      if ( fabs(child_1->calc()) > GlobalVar.TOL() ) {
        return child_2->calc();
      } else {
        return child_3->calc();
      }
    }
    const std::string write();
    const tuint precedence() { return 17; }
};


#endif // fesslix_flxfunction_ope_calc_H
