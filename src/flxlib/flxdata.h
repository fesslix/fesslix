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

#include "flxmtxfun.h"
#include "flxrbrv.h"
#include "config.h"


class FlxObjBase;
class FlxCodeBlock;

class FLXLIB_EXPORT FlxReadManager : public FlxReaderBase2 {
  private:
    #if FLX_DEBUG
      /**
      * @brief number of instances of this calss - only one class may exist
      */
      static int Cinst;
    #endif
    std::stack<ReadStream*> s;
  public:
    FlxReadManager();
    /**
    * @brief Set a new reader on top of the stack
    */
    void push(ReadStream* readerV);
    /**
    * @brief Removes the reader from the stack -> memory is NOT deallocated
    */
    void pop();
    
    FlxFunction* parse_function(const std::string& funStr);
    FlxFunction* parse_function(py::object pyobj, std::string descr="");
    FlxMtxFun_base* parse_FlxMtxFun(const tuint N, py::object pyobj, std::string descr="");
    FlxCodeBlock* parse_code(const std::string& codeStr);
    
    static void set_funReader( FlxFunctionReader* funReaderV) {funReader = funReaderV;};
};


FLXLIB_EXPORT void set_ReadManager(FlxReadManager* readManager_ptr_);
FLXLIB_EXPORT FlxReadManager* get_ReadManager();

/**
* @brief A class for storing output-streams (ostream&)
*/
class FLXLIB_EXPORT FlxOstreamBox {
  private:
    #if FLX_DEBUG
      /**
      * @brief number of instances of this class - only one class may exist
      */
      static int Cinst;
    #endif
    std::map<std::string, ostreamp*> box;
    /**
    * @brief deletes an existing stream
    * @return true, if stream could be deleted; false: stream not allowed for delete
    */
    const bool delete_stream(ostreamp& strm);
  public:
    FlxOstreamBox();
    ~FlxOstreamBox();
    /**
    * @brief Insert a ostream.
    * @param name An unique name (MUSTE BE LOWERCASE)
    * @param value The ostream to store (a pointer to that ostream is stored)
    * @return true: new object inserted; false: stream overwritten
    */
    const bool insert ( const std::string& name, const ostreamp value);
    /**
    * @brief close the stream - point a dummy-stream to it
    * @param name (MUSTE BE LOWERCASE)
    * @param err throw error if stream does not exist
    */
    void close( const std::string& name, const bool err=true );
    /**
    * @brief Get a ostream by name.
    * @param name The name of the ostream to return. (MUSTE BE LOWERCASE)
    * @return ostreamp&
    */
    const ostreamp& get ( const std::string& name);
    
};


/**
* @brief A class for storing input-streams (ostream&)
*/
class FLXLIB_EXPORT FlxIstreamBox {
  private:
    #if FLX_DEBUG
      /**
      * @brief number of instances of this class - only one class may exist
      */
      static int Cinst;
    #endif
    std::map<std::string, FlxIstream*> box;
  public:
    FlxIstreamBox();
    ~FlxIstreamBox();
    /**
    * @brief Insert a istream.
    * @param name An unique name
    * @param value The istream to store (a pointer to that ostream is stored)
    */
    void insert ( const std::string& name, FlxIstream* value, bool errSerious=true);
    /**
    * @brief Get a FlxIstream by name.
    * @param name The name of the istream to return.
    * @return FlxIstream&
    */
    FlxIstream& get ( const std::string& name );
    /**
    * @brief returns NULL if an vector-input-stream with that name does not exist
    */
    FlxIstream_vector* get_isVector ( const std::string& name );
};

class FLXLIB_EXPORT FlxSubBox {
  private:
    #if FLX_DEBUG
      /**
      * @brief number of instances of this calss - only one class may exist
      */
      static int Cinst;
    #endif
    std::map<std::string, FlxObjBase*> box;
  public:
    FlxSubBox();
    ~FlxSubBox();
    void insert (const std::string& name, FlxObjBase* value);
    FlxObjBase* get(const std::string& name);
};

class FLXLIB_EXPORT FlxIgnoreBox {
  private:
    #if FLX_DEBUG
      /**
      * @brief number of instances of this calss - only one class may exist
      */
      static int Cinst;
    #endif
    std::vector<std::string> ignoreList_recur;
    int iL_recur_active;
  public:
    FlxIgnoreBox ();
    /**
    * @param objName MUSTE BE LOWERCASE
    */
    void set_iL_recur(const std::string& objName ) { ignoreList_recur.push_back(objName); }
    void activate_iL_recur() { iL_recur_active++; }
    void deactivate_iL_recur() { iL_recur_active--; }
    const bool isActive_iL_recur() { return (iL_recur_active > 0)?true:false; }
    /**
    * @param objName MUSTE BE LOWERCASE
    */
    const bool isOnIgnoreList_recur( const std::string& objName );
};


class FLXLIB_EXPORT FlxTimer {
  private:
    bool is_running;
    clock_t startV, delta;
  
  public:
    FlxTimer() : is_running(false), delta(0) {}
    void start() { if (!is_running) { is_running = true; startV = clock(); } }
    void stop() { if (is_running) { is_running = false; delta += clock()-startV; } }
    clock_t get_ticks () { return delta; }
    const tdouble get_time() { return tdouble(delta)/CLOCKS_PER_SEC; }
    bool IsRunning () { return is_running; }
    void reset() { delta = 0; };
};

class FLXLIB_EXPORT FlxTimerBox {
  private:
    #if FLX_DEBUG
      /**
      * @brief number of instances of this calss - only one class may exist
      */
      static int Cinst;
    #endif
    std::map<std::string, FlxTimer*> box;
  public:
    FlxTimerBox();
    ~FlxTimerBox();
    void insert (const std::string& name, FlxTimer* value);
    FlxTimer* get( const std::string& name);
    void deleteEl ( const std::string& name );  
};


/**
* @brief A class for storing string expressions
*/
class FLXLIB_EXPORT flxStrConstBox {
  private:
    #if FLX_DEBUG
      /**
      * @brief number of instances of this class - only one class may exist
      */
      static int Cinst;
    #endif
    std::map<std::string, std::string> box;
  public:
    void insert ( const std::string& name, const std::string& value);
    const std::string& get ( const std::string& name );
    std::string& get_ref ( const std::string& name );
};


/**
* @brief A class for storing data
*/
class FLXLIB_EXPORT FlxData {
  private:
    #if FLX_DEBUG
      /**
      * @brief number of instances of this calss - only one class may exist
      */
      static int Cinst;
    #endif
  public:
    FlxOstreamBox OstreamBox;
    FlxIstreamBox IstreamBox;
    FlxConstantBox ConstantBox;
    FlxVarBox VarBox;
    FlxFunctionBox FunBox;
    FlxConstMtxBox ConstMtxBox;
    FlxSubBox SubBox;
    FlxIgnoreBox IgnoreBox;
    FlxTimerBox TimerBox;
    GaussIntegration GaussInt;
    FlxRndCreator RndCreator;
    FlxReadManager ReadManager;
    FlxStringFunBox StrFunBox;
    RBRV_set_box rbrv_box;
    flxStrConstBox strConstBox;
    
    FlxData();
    ~FlxData();
};

class FLXLIB_EXPORT FlxDataBase {
  protected:
    static FlxData* data;
  public:
    ~FlxDataBase() {};
     /**
    * @brief Set the DataBox 
    * @param dataV a pointer to the DataBox to use
    */
    static void set_data ( FlxData *dataV );;
    static FlxData& get_data() { return *data; }
};


// ****************** flxOBJECTS *************************************************

/**
* @brief The base class for all other 'object classes'
*/
class FLXLIB_EXPORT FlxObjBase : public FlxDataBase {
  private:
    /**
    * @brief The next 'object class'
    *
    * If the command is NOT inside of a loop or an if-block, 'Next' will be NULL.
    */
    FlxObjBase* Next;
  protected:
    /**
    * @brief true if this object should not log
    */
    const bool NOTdolog;
  private:
    void NOTdolog_start() { GlobalVar.logLevel_log_deactivate(true); }
    void NOTdolog_stop() { GlobalVar.logLevel_log_deactivate(false); }
  protected:
    /**
    * @brief This function has to be overriden by all 'object classes' (that is the hart)
    */
    virtual void task() = 0;
    //virtual void task() {};
  public:
    FlxObjBase (const bool dolog) : Next(NULL),NOTdolog(!dolog) {}
    virtual ~FlxObjBase();
    /**
    * @brief Execute this and all following 'object classes'
    */
    virtual void exec();
    /**
    * @brief Attach an 'object class' to the end of the list of 'object classes'
    */
    void attach_obj ( FlxObjBase* NextV );
    
};

/**
* @brief The base class for all other 'object classes'
*/
class FLXLIB_EXPORT FlxCodeBlock : public FlxObjBase {
  private:
    /**
    * contains the const-variables that are marked as 'internal' within the code block
    */
    std::vector<tdouble*> cvec;    
    std::valarray<tdouble> dvec;
    bool catch_return;
    bool catch_continue;
    
    /**
    * @brief This function has to be overriden by all 'object classes' (that is the hart)
    */
    virtual void task() {}
    
  public:
    FlxCodeBlock (const bool dolog) : FlxObjBase(dolog), catch_return(true), catch_continue(false) {}
    virtual ~FlxCodeBlock() {}
    /**
    * @brief Execute this and all following 'object classes'
    */
    virtual void exec();
    
    void loop_block_exec_1();
    void loop_block_exec_2();
    /**
    * @brief marks a const-variable as internal
    */
    void add_internal_const(tdouble* icp);
    void deactivate_return_catch() { catch_return = false; }
    void activate_continue_catch() { catch_continue = true; }
};


/**
* @brief object class: dummy class - does nothing
*/
class FLXLIB_EXPORT FlxObjDummy : public FlxObjBase {
  private:
    void task() {};
  public:
    FlxObjDummy ( ) : FlxObjBase(true) {}
};


/**
* @brief object class: dummy class - does nothing
*/
template <typename T> class FlxVoidBox {
  private:
      #if FLX_DEBUG
        /**
        * @brief number of instances of this calss - only one class may exist
        */
        static int Cinst;
      #endif
      std::map<std::string, T*> box;
  public:
      FlxVoidBox()
      {
        #if FLX_DEBUG
          ++Cinst;
          if (Cinst > 1) {
            std::ostringstream ssV;
            ssV << "More than one instance of 'FlxVoidBox'<" << typeid(T).name() << "> created ...";
            throw FlxException("FlxVoidBox::FlxVoidBox", ssV.str() );
          }
        #endif
      }
      ~FlxVoidBox() {
        for (auto &pos : box) delete pos.second;
      }
      void insert (const std::string& name, T* value) {
        std::pair<std::string, T*>Element(name, value);
        if ( ! box.insert(Element).second ) {;
          delete value;
          std::ostringstream ssV;
          ssV << "Element<" << typeid(T).name() << "> '" << name << "' is already defined.";
          throw FlxException("FlxVoidBox::insert_1", ssV.str() );
        }
      }
      T& get( const std::string& name) {
        auto pos = box.find(name);
        if ( pos != box.end() ) {
          return *(pos->second);
        } else {
          std::ostringstream ssV;
          ssV << "Element<" << typeid(T).name() << "> '" << name << "' does not exist.";
          throw FlxException("FlxVoidBox::get_1", ssV.str() );
        }
      }
      void deleteEl ( const std::string& name)
      {
        T* tt = get(name);
        delete tt;
        box.erase(name);
      }
};
#if FLX_DEBUG
template <typename T>
int FlxVoidBox<T>::Cinst = 0;
#endif




//=================== flxMtxFun ==========================


class FLXLIB_EXPORT FlxMtxFun_MtxConst : public FlxMtxFun_base, public FlxDataBase {
  protected:
    FlxMtxConstFun mtxConstFun;
  public:
    FlxMtxFun_MtxConst(const tuint N, const char* mtxName_strV, FlxObjBase* block=NULL);
    FlxMtxFun_MtxConst(const tuint N, FlxMtxConstFun& mtxConstFun_);
    virtual void eval();
};




