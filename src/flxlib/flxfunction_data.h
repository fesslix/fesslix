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

#ifndef fesslix_flxfunction_data_H
#define fesslix_flxfunction_data_H

#include "flxio.h"
#include "flxmath.h"
#include "flxmath_rnd.h"

#include <map>
#include <cctype>
#include <ctime>


class FlxConstantBox;

/**
* @brief A class for storing a point in space
*/
class FLXLIB_EXPORT flxPoint {
  private:
    static tdouble *GX;
    static tdouble *GY;
    static tdouble *GZ;
    static tdouble *GX2;
    static tdouble *GY2;
    static tdouble *GZ2;
    static tdouble *DELTAX;
    static tdouble *DELTAY;
    static tdouble *DELTAZ;
    static tdouble *DELTAP;
    tdouble x;
    tdouble y;
    tdouble z;
  public:
    flxPoint(const tdouble& x, const tdouble& y, const tdouble& z) :x(x),y(y),z(z) {};
    flxPoint() : x(ZERO),y(ZERO),z(ZERO) {};
    const tdouble& get_x() const { return x;}
    const tdouble& get_y() const { return y;}
    const tdouble& get_z() const { return z;}
    flxPoint& operator+=(const flxPoint& rhs) { x+=rhs.x; y+=rhs.y; z+=rhs.z; return *this; }
    flxPoint& operator*=(const tdouble& rhs) { x*=rhs;y*=rhs;z*=rhs; return *this; }
    const flxPoint operator+(const flxPoint& rhs) const {return flxPoint(x+rhs.x,y+rhs.y,z+rhs.z);}
    const flxPoint operator-(const flxPoint& rhs) const {return flxPoint(x-rhs.x,y-rhs.y,z-rhs.z);}
    const flxPoint operator*(const tdouble& rhs) const {return flxPoint(rhs*x,rhs*y,rhs*z);}
    const tdouble operator*(const flxPoint& rhs) const { return x*rhs.x+y*rhs.y+z*rhs.z;}
    const bool operator==(const flxPoint& rhs) const {return fabs(x-rhs.x)<=GlobalVar.TOL() && fabs(y-rhs.y)<=GlobalVar.TOL() && fabs(z-rhs.z)<=GlobalVar.TOL();}
    const bool operator!=(const flxPoint& rhs) const {return !operator==(rhs);}
    const tdouble& operator[](const int index) const;
    tdouble& id_fast(const int i) { return *(&x+i); }        // 0,1,2
    /**
    * @brief cross product
    */
    flxPoint operator%(const flxPoint& rhs) const;
    const tdouble length() const {return sqrt(pow2(x)+pow2(y)+pow2(z));}
    flxPoint& normalize();
    const tdouble get_phi(const flxPoint& rhs) const;
    FLXLIB_EXPORT friend std::ostream& operator<<(std::ostream& os, const flxPoint& val);
    /**
    * @brief get distance to another point
    */
    const tdouble get_d(const flxPoint& p) const { return sqrt(pow2(x-p.x)+pow2(y-p.y)+pow2(z-p.z));};
    void get_d(const flxPoint& p, tdouble& lx, tdouble &ly, tdouble &lz, tdouble &l) const;
    /**
    * @brief updates global constants GX,GY,GZ
    */
    void set_global() const;
    /**
    * @brief updates global constants GX,GY,GZ,GX2,GY2,GZ2,DELTAX,DELTAY,DELTAZ,DELTAP
    */
    void set_global(const flxPoint &p) const;
    /**
    * @brief updates global constants GX,GY,GZ - point on line !!!
    * @param factor [0;1]
    */
    void set_global(const flxPoint &p, const tdouble &factor) const;
    /**
    * @brief updates global constants GX,GY,GZ,DELTAX,DELTAY,DELTAZ,DELTAP
    */
    void set_global_dist() const;
    /**
    * @brief updates global constants GX,GY,GZ - point on line !!! and DELTAX,DELTAY,DELTAZ,DELTAP
    */
    void set_global_dist(const flxPoint &p, const tdouble &factor) const;
    /**
    * @brief updates global constants GX2,GY2,GZ2
    */
    void set_global2() const;    
    /**
    * @brief updates global constants GX2,GY2,GZ2,DELTAX,DELTAY,DELTAZ,DELTAP
    */
    void set_global2_dist() const;
    
    static void set_Const(FlxConstantBox& ConstantBox);
};

FLXLIB_EXPORT std::ostream& operator<<(std::ostream& os, const flxPoint& val);


struct p2d {
  tdouble xi;
  tdouble nu;
};

class FunBase;
/**
* @brief A class for storing constant values (the result of the evaluated function is stored).
*/
class FLXLIB_EXPORT FlxConstantBox {
  private:
    #if FLX_DEBUG
      /**
      * @brief number of instances of this calss - only one class may exist
      */
      static int Cinst;
    #endif
    std::map<std::string, tdouble*> box;
    tdouble *fastPtr;
      tdouble *lz, *ly, *lx;
      tdouble *gx, *gy, *gz;
      tdouble *gx2, *gy2, *gz2;
      tdouble *deltax, *deltay, *deltaz, *deltap;
  public:
    FlxConstantBox();
    /**
    * @brief insert: the reference of a variable is not changed by 'insert'
    */
    tdouble* insert ( const std::string& name, const tdouble& value);
    tdouble* declareC ( const std::string& name, const tdouble default_val=ZERO );
    tdouble* get ( const std::string& name, const bool doDeclare = false );
    tdouble& getRef ( const std::string& name );
    const std::string get( const tdouble* dv);
    /**
    * @brief update the global constants: gx, gy, gz
    */
//     void update_globalConst(const flxPoint& p, const std::string& add_str);
    /**
    * @brief update the global constants: gx, gy, gz, gx2, gy2, gz2, deltap, deltax, deltay, deltaz
    */
//     void update_globalConst(const flxPoint& p1, const flxPoint& p2);
    void update_globalConst(const flxPoint& p1, const flxPoint& p2, const tdouble& coord_lx_start, const tdouble& coord_lx_end);
    
    /**
    * @brief Determines if the expression 'fun' depends on a global spatial const-variable
    */
    const bool dependOn_GlobalSpatialVar(FunBase* fun) const;
    /**
    * @brief checks if the variable is a spatial random variable
    */
    const bool is_SpatialVar(const tdouble* const thenumber) const;
};

/**
* @brief provides the Gauss points and weights for an integration on the domain [-1;1]
*/
class FLXLIB_EXPORT GaussIntegration {
  private:
    #if FLX_DEBUG
      /**
      * @brief number of instances of this calss - only one class may exist
      */
      static int Cinst;
    #endif
    /**
    * @brief At which Gauss point is the file pointer located?
    */
    tuint filepos;
     /**
    * @brief Has ReadGP() been called at least once?
    */
    bool fileRead;
    /**
    * @brief Vector containing the Gauss points (takes advantage of symmetry)
    */
    tdouble** GP;
    /**
    * @brief Vector containing the Gauss weights (takes advantage of symmetry)
    */
    tdouble** GW;
    /**
    * @brief Number of Gauss points in GP and GW
    */
    tuint numbGP;
    /**
    * @brief the input stream to read the Gauss points
    */
    ReadStream* gaussRS;
    /**
    * @brief returns the Gauss point/weight from a list of points/weights
    *  resolves the symmetry problem
    * @param vec return of either get_GP or get_GW
    * @param index index to extract [0;GA-1]
    * @param GA number of Gauss points/weights (length of vec)
    * @param weight if a weight should be returned, this parameter has to be true
    * @todo maybe faster solution possible - by neglecting symmetry and store everything
    */
    const tdouble get_Point(const tdouble* vec, const tuint& index, const tuint& GA, const bool weight=false);
    /**
    * @brief resolves the symmetry from the input
    * @param vec a vector with either Gauss points or weights
    * @param GA number of Gauss points/weights (length of symVec)
    * @param weight if a weight should be returned, this parameter has to be true
    */
    tdouble* createGPWvector(tdouble* symVec, const tuint& GA, const bool weight=false);
    
    void open_GaussFile(std::string gaussFile);
    
  public:
    /**
    * @brief contains the default number of entries in GP and GW
    */
    static tuint GaussPointArraySize;
    /**
    * @brief contains the maximum number of entries in GP and GW
    */
    static tuint GaussPointMaxArraySize;
    
    GaussIntegration();
    ~GaussIntegration();
    /**
    * @brief reads 'numbR' Gauss points from a file
    */
    void ReadGP(tuint numbR = 0, const std::string gaussFile="");
    /**
    * @brief check if enough Gauss points are loaded (returns an error if it was not possible to load enough points)
    * @param i Required number of Gauss points
    * @throw FlxException - in case not enough Gauss points could be loaded
    */
    void check_GA(const tuint i);
    /**
    * @brief check if enough Gauss points are loaded (returns an error if it was not possible to load enough points)
    * @param i Degree of the polynomial to integrate
    * @throw FlxException - in case not enough Gauss points could be loaded
    */
    void check_GA_polynomialDegree(const tuint i);
    /**
    * @brief returns the corresponding vector for i Gauss points (takes advantage of symmetry)
    *  please make sure that i is within the valid range
    * @see check_GA
    * @return a vector with the Gauss points (takes advantage of symmetry)
    */
    const tdouble* get_GP(const tuint i);
    /**
    * @brief returns the corresponding vector for i Gauss weights (takes advantage of symmetry)
    *  please make sure that i is within the valid range
    * @see get_GP
    */
    const tdouble* get_GW(const tuint i);
    /**
    * @brief transforms a polynomial degree into a number of Gauss points
    */
    const tuint pDegree2GPs(const tuint i) const {
      tuint n = i+1;
      if (n%2 == 1) {
        n = (n+1)/2;
      } else {
        n = n / 2;
      }
      return n;
    }
    /**
    * @brief transforms a number of Gauss points into a polynomial degree
    * @param i the number of Gauss points
    * @param upper if true, the upper polynomial degree is returned; if false, the lower one 
    */
    const tuint GPs2pDegree(const tuint i, bool upper=true) const {
      if (upper) return 2*i-1;
      else return 2*i-2;
    }
};



#endif // fesslix_flxfunction_data_H

