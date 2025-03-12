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

#include "flxstring.h"


#if FLX_DEBUG
  int FlxStringFunBox::Cinst = 0;
#endif


FlxStringFunBox::FlxStringFunBox()
{
  #if FLX_DEBUG
    ++Cinst;
    if (Cinst > 1) {
      std::ostringstream ssV;
      ssV << "More than one instance of 'FlxStringFunBox' created ...";
      throw FlxException("FlxStringFunBox::FlxStringFunBox", ssV.str() );
    }
  #endif
}

FlxStringFunBox::~FlxStringFunBox()
{
  for (std::map<std::string, FunReadFlxStringFunBase*>::iterator pos = box.begin(); pos != box.end(); ++pos) delete pos->second;
}

FunReadFlxStringFunBase* FlxStringFunBox::get(const std::string& name)
{
  std::map<std::string, FunReadFlxStringFunBase*>::iterator pos;
  pos = box.find(name);
  if ( pos != box.end() ) {
    return pos->second;
  } else {
    return NULL;
  }
}

void FlxStringFunBox::insert(const std::string& name, FunReadFlxStringFunBase* fR)
{
  std::pair<std::string, FunReadFlxStringFunBase*>Element(name, fR);
  if ( ! box.insert(Element).second ) {
    std::map<std::string, FunReadFlxStringFunBase*>::iterator pos;
    pos = box.find(name);
    delete pos->second;
    pos->second = fR;
  }
}

StringFunBase* FlxStringFunBox::read(ReadStream* reader, const bool errSerious)
{
  const std::string fn = reader->getWord(true,true);
  FunReadFlxStringFunBase* resR = get(fn);
  if (resR==NULL) {
    std::ostringstream ssV;
    ssV << "FlxString-function '" << fn << "' does not exist.";
    throw FlxException("FlxStringFunBox::read", ssV.str() );
  }
  reader->getChar('(');
  StringFunBase *res = resR->read(errSerious);
  reader->getChar(')');
  return res;
}





