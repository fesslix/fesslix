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

#define FLXLIB_CPP

#include "flxfunction_ope_read.h"
#include "flxfunction_ope.h"

using namespace std;

void flxfunction_ope_insert(FunReadSTART* FunList)
{
  FunList->insert( new FunReadAdd() );                // 40
  FunList->insert( new FunReadMult() );                // 50
  FunList->insert( new FunReadPower() );        // 60
  FunList->insert( new FunReadEqual() );        // 30
  FunList->insert( new FunReadLessThan() );        // 31
  FunList->insert( new FunReadAnd() );                // 20
  FunList->insert( new FunReadOr() );                // 19
  FunList->insert( new FunReadTernary() );        // 10
  FunList->insert( new FunReadNot() );                // 900
  FunList->insert( new FunReadWord() );                // 1010
  FunList->insert( new FunReadBracket() );        // 1100
}

const bool FunVar::search_circref(FlxFunction* fcr) {
  if ( fcr == thefun ) {
    return true;
  }
  return thefun->search_circref(fcr);
}

const string FunConst::write() {
  return ConstantBox->get(thenumber);
}

const bool FunConst::optimize(FunBasePtr& optf, const Fun_OptimizeInfo& foi)
{
  if (foi.computeConstants) {
    calc_me(optf);
    return true;
  } else {
    return false;
  }
}

const std::string FunVar::write() {
  return VarBox->get(thefun);
}

const bool FunVar::optimize(FunBasePtr& optf, const Fun_OptimizeInfo& foi)
{
  if (foi.computeConstants) {
    calc_me(optf);
    return true;
  } else {
    return false;
  }
}

const string FunNot::write() {
  string str1 = "!";
  if (child_1->precedence() > 0 ) {
    str1 += "(" + child_1->write() + ")";
  } else {
    str1 += child_1->write();
  }
  return str1;
}

const bool FunPower::optimize(FunBasePtr& optf, const Fun_OptimizeInfo& foi)
{
  if (FunBaseOperat2::optimize(optf,foi)) return true;
  else if (is_number(child_1) ) {
    const tdouble t = child_1->calc();
    if (t==ZERO) {
      optf = new FunNumber(ZERO);
      return true;
    } else if (t==ONE) {
      optf = new FunNumber(ONE);
      return true;
    }
  } else if (is_number(child_2) ) {
    const tdouble t = child_2->calc();
    if (t==ZERO) {
      optf = new FunNumber(ONE);
      return true;
    } else if (t==ONE) {
      optf = child_1;
      child_1 = new FunDummy();
      return true;
    }
  }
  return false;
}

const std::string FunPower::write() {
  string str1 = "";
  if (child_1->precedence() > precedence()) {
    str1 += "(" + child_1->write() + ")";
  } else {
    str1 += child_1->write();
  }
  str1+="^";
  if (child_2->precedence() > precedence()) {
    str1 += "(" + child_2->write() + ")";
  } else {
    str1 += child_2->write();
  }
  return str1;
}

const std::string FunTernary::write() {
  string str1 = "";
  if (child_1->precedence() >= precedence()) {
    str1 += "(" + child_1->write() + ")";
  } else {
    str1 += child_1->write();
  }
  str1+="?";
  if (child_2->precedence() >= precedence()) {
    str1 += "(" + child_2->write() + ")";
  } else {
    str1 += child_2->write();
  }
  str1+=":";
  if (child_3->precedence() >= precedence()) {
    str1 += "(" + child_3->write() + ")";
  } else {
    str1 += child_3->write();
  }
  return str1;
}

const std::string FunAnd::write() {
  string str1 = "";
  // child_1
    if (child_1->precedence() > precedence()) { 
      str1 += "(" + child_1->write() + ")";
    } else {
      str1 += child_1->write();
    }
    str1 += "&&";
  // child_2
    if (child_2->precedence() > precedence()) { 
      str1 += "(" + child_2->write() + ")";
    } else {
      str1 += child_2->write();
    }
  return str1;
}

const std::string FunOr::write() {
  string str1 = "";
  // child_1
    if (child_1->precedence() > precedence()) { 
      str1 += "(" + child_1->write() + ")";
    } else {
      str1 += child_1->write();
    }
    str1 += "||";
  // child_2
    if (child_2->precedence() > precedence()) { 
      str1 += "(" + child_2->write() + ")";
    } else {
      str1 += child_2->write();
    }
  return str1;
}

const std::string FunEqual::write() {
  string str1 = "";
  if (child_1->precedence() > precedence() ) {
    str1 += "(" + child_1->write() + ")";
  } else {
    str1 += child_1->write();
  }
  if (isEqual) {
    str1+="==";
  } else {
    str1+="!=";
  }
  if (child_2->precedence() >= precedence() ) {
    str1 += "(" + child_2->write() + ")";
  } else {
    str1 += child_2->write();
  }
  return str1;
}

const bool FunMult::optimize(FunBasePtr& optf, const Fun_OptimizeInfo& foi)
{
  if (FunBaseOperat2::optimize(optf,foi)) return true;
  else if (is_number(child_1)) {
    const tdouble t = child_1->calc();
    if (t==ZERO) {
      optf = new FunNumber(ZERO);
      return true;
    } else if (t==ONE) {
      optf = child_2;
      child_2 = new FunDummy();
      return true;
    }
  } else if (is_number(child_2)) {
    const tdouble t = child_2->calc();
    if (t==ZERO) {
      optf = new FunNumber(ZERO);
      return true;
    } else if (t==ONE) {
      optf = child_1;
      child_1 = new FunDummy();
      return true;
    }
  }
  return false;
}

const std::string FunMult::write() {
  string str1 = "";
  if (child_1->precedence() > precedence() ) {
    str1 += "(" + child_1->write() + ")";
  } else {
    str1 += child_1->write();
  }
  str1+="*";
  if (child_2->precedence() > precedence() ) {
    str1 += "(" + child_2->write() + ")";
  } else {
    str1 += child_2->write();
  }
  return str1;  
}

const bool FunMult_Div::optimize(FunBasePtr& optf, const Fun_OptimizeInfo& foi)
{
  if (FunBaseOperat2::optimize(optf,foi)) return true;
  else if (is_number(child_1)) {
    const tdouble t = child_1->calc();
    if (t==ZERO) {
      optf = new FunNumber(ZERO);
      return true;
    }
  }
  return false;
}

const std::string FunMult_Div::write() {
  string str1 = "";
  if (child_1->precedence() > precedence() ) {
    str1 += "(" + child_1->write() + ")";
  } else {
    str1 += child_1->write();
  }
  str1+="/";
  if (child_2->precedence() >= precedence() ) {
    str1 += "(" + child_2->write() + ")";
  } else {
    str1 += child_2->write();
  }
  return str1;  
}

const string FunAdd::write() {
  string str1 = "";
  if (child_1->precedence() > precedence() ) {
    str1 += "(" + child_1->write() + ")";
  } else {
    str1 += child_1->write();
  }
  str1+="+";
  if (child_2->precedence() > precedence() ) {
    str1 += "(" + child_2->write() + ")";
  } else {
    str1 += child_2->write();
  }
  return str1;  
}

const bool FunAdd::optimize(FunBasePtr& optf, const Fun_OptimizeInfo &foi)
{
  if (FunBaseOperat2::optimize(optf,foi)) return true;
  else if (is_number(child_1) ) {
    const tdouble t = child_1->calc();
    if (t==ZERO) {
      optf = child_2;
      child_2 = new FunDummy();
      return true;
    }
  } else if (is_number(child_2) ) {
    const tdouble t = child_2->calc();
    if (t==ZERO) {
      optf = child_1;
      child_1 = new FunDummy();
      return true;
    }
  }
  return false;
}

const string FunSub::write() {
  string str1 = "";
  if (child_1->precedence() > precedence() ) {
    str1 += "(" + child_1->write() + ")";
  } else {
    str1 += child_1->write();
  }
  str1+="-";
  if (child_2->precedence() >= precedence() ) {
    str1 += "(" + child_2->write() + ")";
  } else {
    str1 += child_2->write();
  }
  return str1;  
}

const bool FunSub::optimize(FunBasePtr& optf, const Fun_OptimizeInfo& foi)
{
  if (FunBaseOperat2::optimize(optf,foi)) return true;
  else if (is_number(child_2) ) {
    const tdouble t = child_2->calc();
    if (t==ZERO) {
      optf = child_1;
      child_1 = new FunDummy();
      return true;
    }
  }
  return false;
}

const tdouble FunLessThan::calc() {
  const tdouble r1 = child_1->calc();
  const tdouble r2 = child_2->calc();
  if ( isEqual) {
    if (isLess) {
      if (r1 <= r2) return ONE;
    } else {
      if (r1 >= r2) return ONE;
    }
  } else {
    if (isLess) {
      if (r1 < r2) return ONE;
    } else {
      if (r1 > r2) return ONE;
    }
  }
  return ZERO;
}
    
const std::string FunLessThan::write() {
  string str1 = "";
  if (child_1->precedence() > precedence() ) {
    str1 += "(" + child_1->write() + ")";
  } else {
    str1 += child_1->write();
  }
  if (isLess) {
    str1+="<";
  } else {
    str1+=">";
  }
  if (isEqual) {
    str1+="=";
  }
  if (child_2->precedence() >= precedence() ) {
    str1 += "(" + child_2->write() + ")";
  } else {
    str1 += child_2->write();
  }
  return str1;
}

FunBase* FunReadWord::read (const bool errSerious) {
  if ( reader->getNextType() == ReadStream::STRING ) {
    const string strV = reader->getWord(true,errSerious);
    
    // predefined functions
      FunReadFunBase* fr = funBox->get(strV);
      if ( fr ) {
        FunBase* fb = NULL;
        try {
        reader->getChar('(',errSerious);
        fb = fr->read(errSerious);
        reader->getChar(')',errSerious);
        } catch (FlxException& e) {
          FLXMSG("FunReadWord::read",1);
          if (fb) delete fb;
          throw;
        }
        return fb;
      }

    // user-defined constants
      tdouble* d1 = ConstantBox->get(strV);
      if ( d1 ) {
        return new FunConst( d1 );
      }

    // user-defined variables
      FlxFunction* f1 = VarBox->get(strV);
      if ( f1 ) {
        return new FunVar(f1);
      }

    // Word is unknown ...
      ostringstream ssV_2;
      ssV_2 << "'" << strV << "' was not defined yet.";
      FlxError(errSerious,"FunReadWord::read_1", ssV_2.str(), reader->getCurrentPos());
      return NULL;
  }
  else return Next->read(errSerious);
}

FunBase* FunReadBracket::read (bool errSerious) {
  if ( reader->whatIsNextChar() == '(' ) {
    reader->getChar();
    FunBase *theNode = startLink->read(errSerious);
    if ( reader->getChar() == ')' ) {
      return theNode;
    } else {
      ostringstream ssV_2;
      ssV_2 << "Right parenthesis '(' expected.";
      FlxError(errSerious,"FunReadBracket::read_1", ssV_2.str(), reader->getCurrentPos());
      return NULL;
    }
  } else {
    return Next->read(errSerious);
  }
}

FunBase* FunReadAdd::read (bool errSerious) {
  FunBase *theNode = Next->read(errSerious);
  try {
    while ( reader->whatIsNextChar() == '+' || reader->whatIsNextChar() == '-' ) {
      if ( reader->getChar() == '+' ) {
        theNode = new FunAdd(theNode, Next->read(errSerious));
      } else {
        theNode = new FunSub(theNode, Next->read(errSerious));
      }
    }
  } catch (FlxException& e) {
    FLXMSG("FunReadAdd::read",1);
    if ( theNode) { delete theNode; theNode=NULL; }
    throw;
  }
  return theNode;
}

FunBase* FunReadEqual::read (bool errSerious) {
  FunBase *theNode = Next->read(errSerious);
  bool IsEqual;
  try {
    while ( reader->whatIsNextString(2) == "==" || reader->whatIsNextString(2) == "!=" ) {
      if ( reader->getChar() == '=' ) {
        IsEqual = true;
      }
      else {
        IsEqual = false;
      }
      reader->getChar('=');
      theNode = new FunEqual(theNode, Next->read(errSerious), IsEqual);
    }
  } catch (FlxException& e) {
    FLXMSG("FunReadEqual::read",1);
    if ( theNode) { delete theNode; theNode=NULL; }
    throw;
  }
  return theNode;
}

FunBase* FunReadLessThan::read(bool errSerious) {
  FunBase *theNode = Next->read(errSerious);
  bool IsEqual;
  bool IsLess;
  try {
    while ( reader->whatIsNextString(2) == "<=" || reader->whatIsNextString(2) == ">=" || reader->whatIsNextChar() == '>' || reader->whatIsNextChar() == '<' ) {
      if ( reader->getChar() == '<' ) {
        IsLess = true;
      } else {
        IsLess = false;
      }
      if ( reader->whatIsNextChar() == '=' ) {
        reader->getChar();
        IsEqual = true;
      } else {
        IsEqual = false;
      }

      theNode = new FunLessThan(theNode, Next->read(errSerious), IsEqual, IsLess);
    }
  } catch (FlxException& e) {
    FLXMSG("FunReadLessThan::read",1);
    if ( theNode!=NULL) { delete theNode; theNode=NULL; }
    throw;
  }
  return theNode;  
}

FunBase* FunReadAnd::read(bool errSerious) {
  FunBase *theNode = Next->read(errSerious);
  try {
    while ( reader->whatIsNextString(2) == "&&") {
      reader->getChar(); reader->getChar();
      theNode = new FunAnd(theNode, Next->read(errSerious));
    }
  } catch (FlxException& e) {
    FLXMSG("FunReadAnd::read",1);
    if ( theNode!=NULL) { delete theNode; theNode=NULL; }
    throw;
  }
  return theNode; 
}

FunBase* FunReadOr::read(bool errSerious) {
  FunBase *theNode = Next->read(errSerious);
  try {
    while ( reader->whatIsNextString(2) == "||") {
      reader->getChar(); reader->getChar();
      theNode = new FunOr(theNode, Next->read(errSerious));
    }
  } catch (FlxException& e) {
    FLXMSG("FunReadAnd::read",1);
    if ( theNode!=NULL) { delete theNode; theNode=NULL; }
    throw;
  }
  return theNode; 
}

FunBase* FunReadNot::read(bool errSerious) {
  if ( reader->whatIsNextChar() == '!' ) {
    reader->getChar();
    FunBase* rtn = NULL;
    try {
      rtn = new FunNot(Next->read(errSerious));
    } catch (FlxException& e) {
      FLXMSG("FunReadNot::read",1);
      if ( rtn!=NULL) { delete rtn; rtn=NULL; }
      throw;
    }
    return rtn;
  } else {
    return Next->read(errSerious);
  }
}

FunBase* FunReadTernary::read(bool errSerious) {
  FunBase *theNode = Next->read(errSerious);
  FunBase *theNode2 = NULL;
  try {
    while ( reader->whatIsNextChar() == '?' ) {
      reader->getChar();
      theNode2 = Next->read(errSerious);
      reader->getChar(':');
      theNode = new FunTernary(theNode, theNode2, Next->read(errSerious));
    }
  } catch (FlxException& e) {
    FLXMSG("FunReadTernary::read",1);
    if ( theNode!=NULL) { delete theNode; theNode=NULL; }
    if ( theNode2!=NULL) { delete theNode2; theNode2=NULL; }
    throw;
  }
  return theNode;
}


FunBase* FunReadMult::read (bool errSerious) {
  FunBase *theNode = Next->read(errSerious);
  try {
    while ( reader->whatIsNextChar() == '*' || reader->whatIsNextChar() == '/' ) {
      if ( reader->getChar() == '/' ) {
        theNode = new FunMult_Div(theNode, Next->read(errSerious));
      } else { 
        theNode = new FunMult(theNode, Next->read(errSerious));
      }
    }
  } catch (FlxException& e) {
    FLXMSG("FunReadMult::read",1);
    if ( theNode!=NULL) { delete theNode; theNode=NULL; }
    throw;
  }
  return theNode;
}

FunBase* FunReadPower::read (bool errSerious) {
  FunBase *theNode = Next->read(errSerious);
  try {
    while ( reader->whatIsNextChar() == '^' ) {
      reader->getChar();
      theNode = new FunPower(theNode, Next->read(errSerious));
    }
  } catch (FlxException& e) {
    FLXMSG("FunReadPower::read",1);
    if ( theNode!=NULL) { delete theNode; theNode=NULL; }
    throw;
  }
  return theNode;
}

