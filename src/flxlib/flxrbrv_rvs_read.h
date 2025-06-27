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

#include "flxrbrv_rvs.h"
#include "flxstring.h"


class FlxObjRBRV_set_creator_box;

class FLXLIB_EXPORT RBRV_entry_read_base : public FlxReaderBase2 {
  protected:
    FlxString* nameF;
    // correlation with a single RBRV from the set
      FlxString* corrName;
      FlxFunction* corrVal;
      bool corrFixed;        // evaluated only once!
      bool eval_once;        // evaluate parameters only once (not supported by every type!)
      
      static FlxObjRBRV_set_creator_box* rbrv_set_creator;
      
      /**
      * @brief reads the value of optional parameter eval_once
      */
      void read_eval_once();
      /**
      * @brief reads the value of an optional parameters
      * @returns NULL if pName is not found next
      */
      FlxFunction* read_opt_para(const std::string& pName);
  public:
    RBRV_entry_read_base(const bool readName, const bool readBrakets, const bool hasParas=true);
    virtual ~RBRV_entry_read_base();
    
    virtual RBRV_entry* generate_entry(const std::string& family, tuint& running_iID) = 0;
    RBRV_entry_RV_base* generate_entry_rv(const bool errSerious=true);
    void activate_eval_once() { eval_once = true; }
    
    void read_corr(const bool errSerious);
    /**
    * @note generates correlations for Rosenblatt random variables
    * @param ensureEmpty no correlation definition must exist - otherwise an error is thrown (needed for Nataf transformation)
    */
    void generate_corr(std::vector<RBRV_entry*>& entries, const tuint IDinSet, const bool ensureEmpty);
    
    /**
    * @note if readName=false: dummy will be assigned as name
    */
    static RBRV_entry_read_base* read_entry(const bool readName=true, const bool readBrakets=true);                // reads a single entry
    /**
    * @note reads an entry and generates the corresponding random variable at the same time (used for stand-alone random variables - or if the type of a distribution is required as input to a function)
    */
    static RBRV_entry_RV_base* read_gen_entry(bool errSerious);
    static void read_parents(std::vector<FlxString*>& set_parents, const bool errSerious); // reads the parents
    static void read(std::vector<RBRV_entry_read_base*>& set_entries, std::vector<FlxString*>& set_parents, const bool errSerious); // reads an entire set of entries
    /**
    * @brief ensure that no set of random variables already uses the proposed name
    */
    static void generate_set_base_check_name(RBRV_set_box& box, const std::string& name);
    /**
    * @brief generates the vector with the parent sets and checks if the proposed name is valid and unique
    */
    static void generate_set_base(RBRV_set_box& box, const std::string& name, std::vector<FlxString*> set_parents, RBRV_set_baseDPtr& parents);
    static void generate_set_base(RBRV_set_box& box, std::vector<std::string> set_parents, RBRV_set_baseDPtr& parents);
    
    /**
    * @brief this map manages that that are currently being defined ... map defined in 'flxobjrbrv.cpp'
    */
    static void set_rbrv_set_creator(FlxObjRBRV_set_creator_box* rbrv_set_creatorV) { rbrv_set_creator=rbrv_set_creatorV; }
};


// -------------------------------------------------------------------------------------------


class FLXLIB_EXPORT FlxObjRBRV_set_creator {
  private:
    const std::string set_name;
    const bool is_Nataf;
    const bool is_Nataf_evalOnce;
    RBRV_set_baseDPtr parents; 
    const tuint Nparents;
    std::vector<RBRV_entry*> set_entries;
    const bool allow_x2y;
    tuint rID;
    std::map<std::pair<std::string,std::string>,tdouble> cormap;

    RBRV_entry_RV_base* get_rv(const std::string rvn, const bool throwErr);
    const tuint get_rvID(const std::string rvn);
    
    RBRV_set_Nataf* register_set_Nataf(RBRV_set_box& box, const bool doreg);
  public:
    FlxObjRBRV_set_creator ( const std::string& set_name, RBRV_set_baseDPtr parents, const tuint Nparents, const bool allow_x2y );
    FlxObjRBRV_set_creator ( RBRV_set_box& box, const std::string& set_name, RBRV_set_baseDPtr parents, const tuint Nparents, const bool allow_x2y, std::vector<RBRV_entry_read_base*>& set_entriesV );
    FlxObjRBRV_set_creator ( const std::string& set_name, const bool eval_once );
    ~FlxObjRBRV_set_creator();
    
    const bool get_Nataf_evalOnce() const;
    /**
    * @brief adds an entry to the set
    * @note the memory of 'entry' is not managed!!!
    */
    void add_entry(RBRV_set_box& box, RBRV_entry_read_base* entry);
    /**
    * @brief adds entry to the set
    * @note the memory of 'ep' and 'csVal' is managed!!! (also in case of an error)
    */
    void add_entry(RBRV_set_box& box, RBRV_entry* ep, FlxFunction* csVal=nullptr, const std::string csNam="", const bool csFix=false);
    void add_corr(const std::string& rv1, const std::string& rv2, const tdouble rho, const bool corr_approx, const bool rhogauss, const bool dolog);
    RBRV_set_base* register_set(RBRV_set_box& box, const bool doreg);
    RBRV_set* register_set_rbrv(RBRV_set_box& box, const bool doreg);
};


class FlxObjRBRV_set_creator_box {
  private:
    std::map<std::string,FlxObjRBRV_set_creator*> cmap;
    
  public:
    FlxObjRBRV_set_creator_box ( );
    ~FlxObjRBRV_set_creator_box( );
    
    /**
    * @brief create a new 'set of random variables' that still has to be specified
    * @param name: It must be ensured beforehand, that the used name is unique
    */
    void create_new(const std::string& name, FlxObjRBRV_set_creator* ele);
    /**
    * @brief returns the creator with name 'set_name'
    */
    FlxObjRBRV_set_creator* get_creator( const std::string& set_name, const bool throwErr=true);
    /**
    * @brief registers the completely specified set 'name' as rbrv_set
    */   
    void register_set(const std::string& name, RBRV_set_box& box);
};


// -------------------------------------------------------------------------------------------


class FLXLIB_EXPORT RBRV_entry_read_fun : public RBRV_entry_read_base {
  protected:
    FlxFunction* fun;
  public:
    RBRV_entry_read_fun(const bool readName, const bool readBrakets=true);
    virtual ~RBRV_entry_read_fun();
    
    virtual RBRV_entry* generate_entry(const std::string& family, tuint& running_iID);
};

class FLXLIB_EXPORT RBRV_entry_read_normal : public RBRV_entry_read_base {
  protected:
    int pid;
    FlxFunction* p1;
    FlxFunction* p2;
    FlxFunction* p3;
    FlxFunction* p4;
  public:
    RBRV_entry_read_normal(const bool readName, const bool readBrakets=true);
    virtual ~RBRV_entry_read_normal();
    
    virtual RBRV_entry* generate_entry(const std::string& family, tuint& running_iID);
};

class FLXLIB_EXPORT RBRV_entry_read_stdn : public RBRV_entry_read_base {
  protected:
  public:
    RBRV_entry_read_stdn(const bool readName, const bool readBrakets=true);
    virtual ~RBRV_entry_read_stdn();
    
    virtual RBRV_entry* generate_entry(const std::string& family, tuint& running_iID);
};

class FLXLIB_EXPORT RBRV_entry_read_logn : public RBRV_entry_read_base {
  protected:
    int pid;
    FlxFunction* p1;
    FlxFunction* p2;
    FlxFunction* p3;
    FlxFunction* p4;
    FlxFunction* epsilon;
  public:
    RBRV_entry_read_logn(const bool readName, const bool readBrakets=true);
    virtual ~RBRV_entry_read_logn();
    
    virtual RBRV_entry* generate_entry(const std::string& family, tuint& running_iID);
};

class FLXLIB_EXPORT RBRV_entry_read_uniform : public RBRV_entry_read_base {
  protected:
    FlxFunction* a;
    FlxFunction* b;
  public:
    RBRV_entry_read_uniform(const bool readName, const bool readBrakets=true);
    virtual ~RBRV_entry_read_uniform();
    
    virtual RBRV_entry* generate_entry(const std::string& family, tuint& running_iID);
};

class FLXLIB_EXPORT RBRV_entry_read_Gumbel : public RBRV_entry_read_base {
  protected:
    int methID;
    FlxFunction* p1;
    FlxFunction* p2;
    FlxFunction* p3;
    FlxFunction* p4;
  public:
    RBRV_entry_read_Gumbel(const bool readName, const bool readBrakets=true);
    virtual ~RBRV_entry_read_Gumbel();
    
    virtual RBRV_entry* generate_entry(const std::string& family, tuint& running_iID);
};

class FLXLIB_EXPORT RBRV_entry_read_normal_trunc : public RBRV_entry_read_base {
  protected:
    FlxFunction* m;
    FlxFunction* s;
    FlxFunction* a;
    FlxFunction* b;
  public:
    RBRV_entry_read_normal_trunc(const bool readName, const bool readBrakets=true);
    virtual ~RBRV_entry_read_normal_trunc();
    
    virtual RBRV_entry* generate_entry(const std::string& family, tuint& running_iID);
};

class FLXLIB_EXPORT RBRV_entry_read_beta : public RBRV_entry_read_base {
  protected:
    bool is_mean;
    FlxFunction* p1;
    FlxFunction* p2;
    FlxFunction* a;
    FlxFunction* b;
  public:
    RBRV_entry_read_beta(const bool readName, const bool readBrakets=true);
    virtual ~RBRV_entry_read_beta();
    
    virtual RBRV_entry* generate_entry(const std::string& family, tuint& running_iID);
};

class FLXLIB_EXPORT RBRV_entry_read_exponential : public RBRV_entry_read_base {
  protected:
    FlxFunction* lambda;
    FlxFunction* epsilon;
  public:
    RBRV_entry_read_exponential(const bool readName, const bool readBrakets=true);
    virtual ~RBRV_entry_read_exponential();
    
    virtual RBRV_entry* generate_entry(const std::string& family, tuint& running_iID);
};

class FLXLIB_EXPORT RBRV_entry_read_gamma : public RBRV_entry_read_base {
  protected:
    bool is_mean;
    FlxFunction* p1;
    FlxFunction* p2;
    FlxFunction* epsilon;
  public:
    RBRV_entry_read_gamma(const bool readName, const bool readBrakets=true);
    virtual ~RBRV_entry_read_gamma();
    
    virtual RBRV_entry* generate_entry(const std::string& family, tuint& running_iID);
};

class FLXLIB_EXPORT RBRV_entry_read_Poisson : public RBRV_entry_read_base {
  protected:
    FlxFunction* mean;
  public:
    RBRV_entry_read_Poisson(const bool readName, const bool readBrakets=true);
    virtual ~RBRV_entry_read_Poisson();
    
    virtual RBRV_entry* generate_entry(const std::string& family, tuint& running_iID);
};

class FLXLIB_EXPORT RBRV_entry_read_Binomial : public RBRV_entry_read_base {
  protected:
    FlxFunction* p;
    FlxFunction* N;
  public:
    RBRV_entry_read_Binomial(const bool readName, const bool readBrakets=true);
    virtual ~RBRV_entry_read_Binomial();
    
    virtual RBRV_entry* generate_entry(const std::string& family, tuint& running_iID);
};

class FLXLIB_EXPORT RBRV_entry_read_Cauchy : public RBRV_entry_read_base {
  protected:
    FlxFunction* loc;
    FlxFunction* scale;
  public:
    RBRV_entry_read_Cauchy(const bool readName, const bool readBrakets=true);
    virtual ~RBRV_entry_read_Cauchy();
    
    virtual RBRV_entry* generate_entry(const std::string& family, tuint& running_iID);
};

class FLXLIB_EXPORT RBRV_entry_read_Weibull : public RBRV_entry_read_base {
  protected:
    bool is_mean;
    FlxFunction* p1;
    FlxFunction* p2;
    FlxFunction* epsilon;
  public:
    RBRV_entry_read_Weibull(const bool readName, const bool readBrakets=true);
    virtual ~RBRV_entry_read_Weibull();
    
    virtual RBRV_entry* generate_entry(const std::string& family, tuint& running_iID);
};

class FLXLIB_EXPORT RBRV_entry_read_ChiSquared : public RBRV_entry_read_base {
  protected:
    const bool isSquared;
    FlxFunction* p1;
  public:
    RBRV_entry_read_ChiSquared(const bool isSquared, const bool readName, const bool readBrakets=true);
    virtual ~RBRV_entry_read_ChiSquared();
    
    virtual RBRV_entry* generate_entry(const std::string& family, tuint& running_iID);
};

class FLXLIB_EXPORT RBRV_entry_read_StudentsT : public RBRV_entry_read_base {
  protected:
    FlxFunction* p1;
  public:
    RBRV_entry_read_StudentsT(const bool readName, const bool readBrakets=true);
    virtual ~RBRV_entry_read_StudentsT();
    
    virtual RBRV_entry* generate_entry(const std::string& family, tuint& running_iID);
};

class FLXLIB_EXPORT RBRV_entry_read_StudentsT_generalized : public RBRV_entry_read_base {
  protected:
    FlxFunction* nu;
    FlxFunction* loc;
    FlxFunction* scale;
  public:
    RBRV_entry_read_StudentsT_generalized(const bool readName, const bool readBrakets=true);
    virtual ~RBRV_entry_read_StudentsT_generalized();

    virtual RBRV_entry* generate_entry(const std::string& family, tuint& running_iID);
};

class FLXLIB_EXPORT RBRV_entry_read_Laplace : public RBRV_entry_read_base {
  protected:
    FlxFunction* loc;
    FlxFunction* scale;
  public:
    RBRV_entry_read_Laplace(const bool readName, const bool readBrakets=true);
    virtual ~RBRV_entry_read_Laplace();

    virtual RBRV_entry* generate_entry(const std::string& family, tuint& running_iID);
};

class FLXLIB_EXPORT RBRV_entry_read_UserTransform : public RBRV_entry_read_base {
  protected:
    bool is_z2x;
    FlxFunction* t1;                   // z2x or y2z
    FlxFunction* t2;                 // x2z or z2y
    FlxFunction* dh;                 // dx2z/dx or dz2y/dz
    FlxFunction* checkXf;        // returns true if x is valid
    RBRV_entry_read_base* rv_z;
  public:
    RBRV_entry_read_UserTransform(const bool readName, const bool readBrakets=true);
    virtual ~RBRV_entry_read_UserTransform();
    
    virtual RBRV_entry* generate_entry(const std::string& family, tuint& running_iID);
};

class FLXLIB_EXPORT RBRV_entry_read_Truncated : public RBRV_entry_read_base {
  protected:
    FlxFunction* a;                   // z2x or y2z
    FlxFunction* b;                 // x2z or z2y
    RBRV_entry_read_base* rv_z;
  public:
    RBRV_entry_read_Truncated(const bool readName, const bool readBrakets=true);
    virtual ~RBRV_entry_read_Truncated();

    virtual RBRV_entry* generate_entry(const std::string& family, tuint& running_iID);
};

class FLXLIB_EXPORT RBRV_entry_read_maxminTransform : public RBRV_entry_read_base {
  protected:
    bool is_max;
    FlxFunction* n;
    RBRV_entry_read_base* rv_z;
  public:
    RBRV_entry_read_maxminTransform(const bool readName, const bool readBrakets=true);
    virtual ~RBRV_entry_read_maxminTransform();
    
    virtual RBRV_entry* generate_entry(const std::string& family, tuint& running_iID);
};



