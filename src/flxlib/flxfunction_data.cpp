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

#include "flxfunction.h"
#include "flxdefault.h"

tdouble *flxPoint::GX = NULL;
tdouble *flxPoint::GY = NULL;
tdouble *flxPoint::GZ = NULL;
tdouble *flxPoint::GX2 = NULL;
tdouble *flxPoint::GY2 = NULL;
tdouble *flxPoint::GZ2 = NULL;
tdouble *flxPoint::DELTAX = NULL;
tdouble *flxPoint::DELTAY = NULL;
tdouble *flxPoint::DELTAZ = NULL;
tdouble *flxPoint::DELTAP = NULL;

const tuint CONST_CNTNR_SIZE = 1000;
tuint GaussIntegration::GaussPointArraySize = DEFAULT_GAUSS_NUMB;
tuint GaussIntegration::GaussPointMaxArraySize = DEFAULT_GAUSS_MAXNUMB;

#if FLX_DEBUG
  int FlxConstantBox::Cinst = 0;
  int GaussIntegration::Cinst = 0;
#endif

using namespace std;


const tdouble flxPoint::get_phi(const flxPoint& rhs) const
{
  tdouble c = (*this * rhs)/(this->length()*rhs.length());
  return std::acos(c);
}

flxPoint& flxPoint::normalize()
{
  tdouble il = 1/length();
  x*=il; y*=il; z*=il;
  return *this;
}

const tdouble& flxPoint::operator[](const int index) const
{
  switch (index) {
    case 0:
      return x;
      break;
    case 1:
      return y;
      break;
    case 2:
      return z;
      break;
    default:
      std::ostringstream ssV;
      ssV << "Index out of range (" << index << ").";
      throw FlxException("flxPoint::operator[]", ssV.str() );
  }
}

flxPoint flxPoint::operator%(const flxPoint& rhs) const
{
  return flxPoint(y*rhs.z-z*rhs.y,z*rhs.x-x*rhs.z,x*rhs.y-y*rhs.x);
}

void flxPoint::get_d(const flxPoint& p, tdouble& lx, tdouble& ly, tdouble& lz, tdouble& l) const
{
  lx = p.x-x;
  ly = p.y-y;
  lz = p.z-z;
  l = sqrt(pow2(lx)+pow2(ly)+pow2(lz));
}

void flxPoint::set_global() const
{
  *(GX+0) = x;
  *(GX+1) = y;
  *(GX+2) = z;
}

void flxPoint::set_global(const flxPoint& p) const
{
  *(GX+0) = x;
  *(GX+1) = y;
  *(GX+2) = z;
  *(GX+3) = p.x;
  *(GX+4) = p.y;
  *(GX+5) = p.z;
  *(GX+6) = fabs(*(GX+0) - *(GX+3));
  *(GX+7) = fabs(*(GX+1) - *(GX+4));
  *(GX+8) = fabs(*(GX+2) - *(GX+5));
  *(GX+9) = sqrt(pow2(*(GX+6)) + pow2(*(GX+7)) + pow2(*(GX+8)));
}

void flxPoint::set_global(const flxPoint& p, const tdouble& factor) const
{
  *(GX+0) = x + factor*(p.x-x);
  *(GX+1) = y + factor*(p.y-y);
  *(GX+2) = z + factor*(p.z-z);
}

void flxPoint::set_global_dist(const flxPoint& p, const tdouble& factor) const
{
  *(GX2-3) = x + factor*(p.x-x);
  *(GX2-2) = y + factor*(p.y-y);
  *(GX2-1) = z + factor*(p.z-z);
  *(GX2+3) = fabs(*(GX2+0) - *(GX2-3));
  *(GX2+4) = fabs(*(GX2+1) - *(GX2-2));
  *(GX2+5) = fabs(*(GX2+2) - *(GX2-1));
  *(GX2+6) = sqrt(pow2(*(GX2+3)) + pow2(*(GX2+4)) + pow2(*(GX2+5)));
}

void flxPoint::set_global2() const
{
  *(GX2+0) = x;
  *(GX2+1) = y;
  *(GX2+2) = z;
}

void flxPoint::set_global_dist() const
{
  *(GX2-3) = x;
  *(GX2-2) = y;
  *(GX2-1) = z;
  *(GX2+3) = fabs(*(GX2+0) - *(GX2-3));
  *(GX2+4) = fabs(*(GX2+1) - *(GX2-2));
  *(GX2+5) = fabs(*(GX2+2) - *(GX2-1));
  *(GX2+6) = sqrt(pow2(*(GX2+3)) + pow2(*(GX2+4)) + pow2(*(GX2+5)));
}

void flxPoint::set_global2_dist() const
{
  *(GX2+0) = x;
  *(GX2+1) = y;
  *(GX2+2) = z;
  *(GX2+3) = fabs(*(GX2+0) - *(GX2-3));
  *(GX2+4) = fabs(*(GX2+1) - *(GX2-2));
  *(GX2+5) = fabs(*(GX2+2) - *(GX2-1));
  *(GX2+6) = sqrt(pow2(*(GX2+3)) + pow2(*(GX2+4)) + pow2(*(GX2+5)));
}

void flxPoint::set_Const(FlxConstantBox& ConstantBox)
{
  GX = ConstantBox.get("gx",true);
  GY = ConstantBox.get("gy",true);
  GZ = ConstantBox.get("gz",true);
  GX2 = ConstantBox.get("gx2",true);
  GY2 = ConstantBox.get("gy2",true);
  GZ2 = ConstantBox.get("gz2",true);
  DELTAX = ConstantBox.get("deltax",true);
  DELTAY = ConstantBox.get("deltay",true);
  DELTAZ = ConstantBox.get("deltaz",true);
  DELTAP = ConstantBox.get("deltap",true);
}

std::ostream& operator<<(std::ostream& os, const flxPoint& val)
{
  return os << "(" << GlobalVar.Double2String(val.x) << "," << GlobalVar.Double2String(val.y) << "," << GlobalVar.Double2String(val.z) << ")";
}

FlxConstantBox::FlxConstantBox() 
{
  #if FLX_DEBUG
    ++Cinst;
    if (Cinst > 1) {
      std::ostringstream ssV;
      ssV << "More than one instance of 'FlxConstantBox' created ...";
      throw FlxException("FlxConstantBox::FlxConstantBox", ssV.str() );
    }
  #endif
  // define some geometric constants
    fastPtr = GlobalVar.MemMngr.new_vector(13);
    lz = fastPtr+0; ly = fastPtr+1; lx = fastPtr+2; 
    gx = fastPtr+3; gy = fastPtr+4; gz = fastPtr+5;
    gx2 = fastPtr+6; gy2 = fastPtr+7; gz2 = fastPtr+8;
    deltax = fastPtr+9; deltay = fastPtr+10; deltaz = fastPtr+11; deltap = fastPtr+12;
    pair<std::string, tdouble*>Element0("lz", lz); box.insert(Element0);
    pair<std::string, tdouble*>Element1("ly", ly); box.insert(Element1);
    pair<std::string, tdouble*>Element2("lx", lx); box.insert(Element2);
    pair<std::string, tdouble*>Element3("gx", gx); box.insert(Element3);
    pair<std::string, tdouble*>Element4("gy", gy); box.insert(Element4);
    pair<std::string, tdouble*>Element5("gz", gz); box.insert(Element5);
    pair<std::string, tdouble*>Element6("gx2", gx2); box.insert(Element6);
    pair<std::string, tdouble*>Element7("gy2", gy2); box.insert(Element7);
    pair<std::string, tdouble*>Element8("gz2", gz2); box.insert(Element8);
    pair<std::string, tdouble*>Element9("deltax", deltax); box.insert(Element9);
    pair<std::string, tdouble*>Element10("deltay", deltay); box.insert(Element10);
    pair<std::string, tdouble*>Element11("deltaz", deltaz); box.insert(Element11);
    pair<std::string, tdouble*>Element12("deltap", deltap); box.insert(Element12);
  FlxBoxBase::set_constBox(this);
  // insert more constants
    insert("pi", PI );
    insert("gamma", GAMMA );
    insert("e", std::exp(ONE) );
    insert("true", ONE);
    insert("false", ZERO);
    insert("ans", ZERO);
    insert("leak_check", GlobalVar.is_leak_check());
  #ifdef __WINDOWS__
    insert("is_win", ONE);
  #else
    insert("is_win", ZERO);
  #endif
}

tdouble* FlxConstantBox::insert(const std::string& name, const tdouble& value) {
  tdouble *d1 = get(name,false);
  if (d1) {        // constant does already exist
    *d1 = value;
    return d1;
  } 
  d1 = GlobalVar.MemMngr.new_double();
  *d1 = value;
  pair<std::string, tdouble*>Element(name, d1);
  if ( ! box.insert(Element).second ) {
    throw FlxException_Crude("FlxConstantBox::insert");
  }
  return d1;
}

tdouble* FlxConstantBox::declareC(const string& name, const tdouble default_val) {
  tdouble* res = get(name);
  if (res) return res;
  return insert(name, default_val);
}

tdouble* FlxConstantBox::get(const std::string& name, const bool doDeclare) {
  map<std::string, tdouble*>::iterator pos;
  pos = box.find(name);
  if ( pos != box.end() ) {
    return pos->second;
  } else {
    if (doDeclare) {
      return declareC(name);
    } else {
      return NULL;
    }
  }
}

tdouble& FlxConstantBox::getRef(const string& name)
{
  tdouble* p = get(name,false);
  if (p==NULL) {
    std::ostringstream ssV;
    ssV << "A constant with name '" << name << "' does not exist.";
    throw FlxException("FlxConstantBox::getRef", ssV.str() );
  }
  return *p;
}

const string FlxConstantBox::get(const tdouble* dv)
{
  for (map<string, tdouble*>::iterator pos = box.begin(); pos != box.end(); ++pos) {
    if (pos->second == dv) {
      return pos->first;
    }
  }
  std::ostringstream ssV;
  ssV << "Pointer not part of the list.";
  throw FlxException("FlxConstantBox::get", ssV.str() );
}

// void FlxConstantBox::update_globalConst(const flxPoint& p, const string& add_str)
// {
//   if (add_str=="2") {
//     *gx2 = p.get_x();
//     *gy2 = p.get_y();
//     *gz2 = p.get_z();
//   } else {
//     insert("gx"+add_str, p.get_x() );
//     insert("gy"+add_str, p.get_y() );
//     insert("gz"+add_str, p.get_z() );
//   }
// }

// void FlxConstantBox::update_globalConst(const flxPoint& p1, const flxPoint& p2)
// {
//   p1.set_global();
//   update_globalConst(p2,"2");
//   *deltap = p1.get_d(p2);
//   *deltax = fabs(p2.get_x()-p1.get_x());
//   *deltay = fabs(p2.get_y()-p1.get_y());
//   *deltaz = fabs(p2.get_z()-p1.get_z());
// }

void FlxConstantBox::update_globalConst(const flxPoint& p1, const flxPoint& p2, const tdouble& coord_lx_start, const tdouble& coord_lx_end)
{
  const flxPoint p1t = p1+(p2-p1)*((coord_lx_start+1)/2.0);
  const flxPoint p2t = p1+(p2-p1)*((coord_lx_end+1)/2.0);
  p1t.set_global(p2t);
}

const bool FlxConstantBox::dependOn_GlobalSpatialVar(FunBase* fun) const
{
  if (fun->dependOn_Const(gx)) return true;
  if (fun->dependOn_Const(gy)) return true;
  if (fun->dependOn_Const(gz)) return true;
  if (fun->dependOn_Const(gx2)) return true;
  if (fun->dependOn_Const(gy2)) return true;
  if (fun->dependOn_Const(gz2)) return true;
  if (fun->dependOn_Const(deltap)) return true;
  if (fun->dependOn_Const(deltax)) return true;
  if (fun->dependOn_Const(deltay)) return true;
  if (fun->dependOn_Const(deltaz)) return true;
  if (fun->dependOn_Const(lx)) return true;
  if (fun->dependOn_Const(ly)) return true;
  if (fun->dependOn_Const(lz)) return true;
  return false;
}

const bool FlxConstantBox::is_SpatialVar(const tdouble* const thenumber) const
{
  if (thenumber>=fastPtr && thenumber <= fastPtr+12) {
    return true;
  } else {
    return false;
  }
}

void GaussIntegration::check_GA(const tuint i)
{
  if (i <= numbGP) return;
  ReadGP(i);
  if (i > numbGP) {
    std::ostringstream ssV;
    ssV << "Not enough Gauss points available.";
    throw FlxException("GaussIntegration::check_GA", ssV.str() );
  }
}

void GaussIntegration::check_GA_polynomialDegree(const tuint i)
{
  check_GA(pDegree2GPs(i));
}

void GaussIntegration::ReadGP(tuint numbR, const std::string gaussFile)
{
  if (gaussFile!="") {
    open_GaussFile(gaussFile);    
  }
  if (GaussPointMaxArraySize <= numbGP) return;
  if (gaussRS==NULL) {
    if (numbR==0) {
      GlobalVar.alert.alert("GaussIntegration::ReadGP_0","No Gauss-file loaded.");
      return;
    } else {
      throw FlxException("GaussIntegration::ReadGP_1","No Gauss-file loaded.");
    }
  }
  
  // the number of Gauss points to read
  if (numbR == 0) numbR = GaussPointArraySize;
  
  if (!fileRead) {        // check if file has already been read
    const tuint fnumbU = gaussRS->get_UInt<tuint>();
    fileRead=true;
    if (fnumbU < GaussPointMaxArraySize) {
      GlobalVar.alert.alert("GaussIntegration::ReadGP_3", "Gauss-file contains less Gauss points than specified in 'gauss.maxnumb'.");
      GaussPointMaxArraySize = fnumbU;
    }
  }
  if (numbR > GaussPointMaxArraySize) {
    std::ostringstream ssV;
    ssV << "It is not possilbe to load more than " << GaussPointMaxArraySize << " Gauss points (" << numbR << ").";
    throw FlxException("GaussIntegration::ReadGP_4", ssV.str() );
  }
  
  tdouble* tmpVec1 = new tdouble[tuint((GaussPointMaxArraySize+1)/2)];
  tdouble* tmpVec2 = new tdouble[tuint((GaussPointMaxArraySize+1)/2)];
  tuint tmpV;
  for (tuint i = filepos; i <= numbR; ++i) {   // ... von 1 bis ...
    tmpV = gaussRS->get_UInt<tuint>();
    if (tmpV != i) {
      GlobalVar.alert.alert("GaussIntegration::ReadGP_5", "Gauss-file: error while reading file.");
    }
    tmpV = (i+1)/2;
    for (tuint j = 0; j < tmpV; ++j) {
      tmpVec1[j] = gaussRS->get_Double();
    }
    for (tuint j = 0; j < tmpV; ++j) {
      tmpVec2[j] = gaussRS->get_Double();
    }

    if (i > numbGP) {
      GP[numbGP] = createGPWvector(tmpVec1,i);
      GW[numbGP] =  createGPWvector(tmpVec2,i,true);
      ++numbGP;
    }
  }
  filepos = numbGP + 1;
  
  delete[] tmpVec1;
  delete[] tmpVec2;
  GlobalVar.slog(4) << "Gauss-points: up to " << numbGP << " Gauss points inserted." << std::endl;
}

const tdouble GaussIntegration::get_Point(const tdouble* vec, const tuint& index, const tuint& GA, const bool weight)
{
  const tuint h = tuint((GA+1)/2)-GA%2;
  if (index < h) { return vec[tuint((GA+1)/2)-index-1]*(weight?1:-1); } 
  else if (index==h && GA%2==1) { return vec[0]; }
  else if (index < GA) { return vec[index-h]; } 
  else { 
    std::ostringstream ssV; ssV << "Index '" << index << "' out of range '" << GA << "'.";
    throw FlxException("GaussIntegration::get_Point", ssV.str() );
  }
}

tdouble* GaussIntegration::createGPWvector(tdouble* symVec, const tuint& GA, const bool weight)
{
  tdouble* tmpV = new tdouble[GA];
  for (tuint i = 0; i < GA; ++i) {
    tmpV[i] = get_Point(symVec,i,GA,weight);
  }
  return tmpV;
}

void GaussIntegration::open_GaussFile(string gaussFile)
{
  if (gaussRS) throw FlxException_Crude("GaussIntegration::open_GaussFile_1");
  if ( gaussFile == "{no}" ) return;
  if ( gaussFile == "{default}") {
    gaussFile = GlobalVar.get_exe_dir();
    if (!gaussFile.empty()) {
      const char last_char = gaussFile.back();
      if (last_char!='/' || last_char!='\\') {
        gaussFile.append("/");
      }
    }
    gaussFile.append("gausspoints.dat");
    try {
      gaussRS = new ReadStream(gaussFile.c_str());
    } catch (FlxException &e) {
      FLXMSG("GaussIntegration::open_GaussFile_2",1);
    }
    return;
  }
  gaussRS = new ReadStream(gaussFile.c_str());
}

GaussIntegration::GaussIntegration() : filepos(1), fileRead(false), gaussRS(NULL)
{
  tdouble* tmpVec = new tdouble[tuint((GaussPointMaxArraySize+1)/2)];
  #if FLX_DEBUG
    ++Cinst;
    if (Cinst > 1) {
      std::ostringstream ssV;
      ssV << "More than one instance of 'GaussIntegration' created ...";
      throw FlxException("GaussIntegration::GaussIntegration", ssV.str() );
    }
  #endif
  GP = new tdouble*[GaussPointMaxArraySize];
  GW = new tdouble*[GaussPointMaxArraySize];
  // 1 Gauss Points
    tmpVec[0] = 0.0;
    GP[0] = createGPWvector(tmpVec,1);
    tmpVec[0] = 2.0;
    GW[0] = createGPWvector(tmpVec,1,true);
  // 2 Gauss Points
    tmpVec[0] = 1.0/sqrt(3.0);
    GP[1] = createGPWvector(tmpVec,2);
    tmpVec[0] = 1.0;
    GW[1] = createGPWvector(tmpVec,2,true);
  // 3 Gauss Points    
    tmpVec[0] = 0.0;
    tmpVec[1] = sqrt(3.0/5.0);
    GP[2] = createGPWvector(tmpVec,3);
    tmpVec[0] = 8.0/9.0;
    tmpVec[1] = 5.0/9.0;
    GW[2] = createGPWvector(tmpVec,3,true);
  // 4 Gauss Points
    tmpVec[0] = sqrt((3.0-2.0*sqrt(6.0/5.0))/7.0);
    tmpVec[1] = sqrt((3.0+2.0*sqrt(6.0/5.0))/7.0);
    GP[3] = createGPWvector(tmpVec,4);
    tmpVec[0] = (18.0+sqrt(30.0))/36.0;
    tmpVec[1] = (18.0-sqrt(30.0))/36.0;
    GW[3] = createGPWvector(tmpVec,4,true);
  // 5 Gauss Points
    tmpVec[0] = 0.0;
    tmpVec[1] = sqrt(5.0-2.0*sqrt(10.0/7.0))/3;
    tmpVec[2] = sqrt(5.0+2.0*sqrt(10.0/7.0))/3;
    GP[4] = createGPWvector(tmpVec,5);
    tmpVec[0] = 128.0/225.0;
    tmpVec[1] = (322.0+13.0*sqrt(70.0))/900.0;
    tmpVec[2] = (322.0-13.0*sqrt(70.0))/900.0;
    GW[4] = createGPWvector(tmpVec,5,true);
  numbGP = 5;    
  for (tuint i=numbGP;i<GaussPointMaxArraySize;++i) {
    GP[i]=NULL;
    GW[i]=NULL;
  }
  delete[] tmpVec;
}

const tdouble* GaussIntegration::get_GP(const tuint i)
{
  #if FLX_DEBUG
    if (i==0 || i>numbGP) {
      std::ostringstream ssV;
      ssV << "Number of Gauss points '" << i << "' is not valid.";
      throw FlxException("GaussIntegration::get_GP_1", ssV.str() );
    }
  #endif
  return GP[i-1];  
}

const tdouble* GaussIntegration::get_GW(const tuint i)
{
  #if FLX_DEBUG
    if (i==0 || i>numbGP) {
      std::ostringstream ssV;
      ssV << "Number of Gauss points '" << i << "' is not valid.";
      throw FlxException("GaussIntegration::get_GW_1", ssV.str() );
    }
  #endif
  return GW[i-1];  
}

GaussIntegration::~GaussIntegration()
{
  for (tuint i=0;i<numbGP;i++) {
    if (GP[i]!=NULL) {
      delete [] GP[i];
      delete [] GW[i];
    }
  }
  delete [] GP;
  delete [] GW;
  if (gaussRS) {
    delete gaussRS;
    gaussRS = NULL;
  }
}







