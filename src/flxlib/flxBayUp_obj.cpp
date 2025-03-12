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

#define fesslix_flxBayUp_obj_CPP

#include "flxBayUp_obj.h"

#include "flxBayDA.h"
#include "flxobjmtx.h"

#include <boost/random.hpp>

void FlxCreateObjReaders_BU::createObjReaders(FlxObjectReadBox* objReadBox) {
  objReadBox->insert("bayup_new", new FlxObjReadBayUp_new());
  objReadBox->insert("bayup_likelihood", new FlxObjReadBayUp_likelihood());
  objReadBox->insert("bayup_uncertobsv", new FlxObjReadBayUp_uncertobsv());
  objReadBox->insert("bayup_glbllikelihood", new FlxObjReadBayUp_glbllikelihood(flxBayUp::UNDEFINED));
  objReadBox->insert("bayup_abcmetric", new FlxObjReadBayUp_glbllikelihood(flxBayUp::ABC));
  objReadBox->insert("bayup_ralsf", new FlxObjReadBayUp_glbllikelihood(flxBayUp::RA));
  objReadBox->insert("bayup_update", new FlxObjReadBayUp_update());
  objReadBox->insert("bayup_set", new FlxObjReadBayUp_Set());
  objReadBox->insert("bayup_reset_smpls", new FlxObjReadBayUp_Reset_Smpls());
  objReadBox->insert("bayda_new", new FlxObjReadBayDA_new());
  objReadBox->insert("bayda_sample", new FlxObjReadBayDA_sample());
}

void FlxCreateObjReaders_BU::createFunReaders(FlxData* dataBox)
{
  dataBox->FunBox.insert("bayup_prop", new FunReadFunBayUp_Prop() );
  dataBox->FunBox.insert("bayup_lsf", new FunReadFunBayUp_lsf() );
  dataBox->FunBox.insert("convexp", new FunReadFunConvExp() );
}

//----------------------------------------------------------------------------------------------


FlxObjBayUp_new::FlxObjBayUp_new( const bool dolog, const std::string& ostreamV, FlxString* nameID, FlxString* rbrvsets, FlxFunction* cStart, FlxFunction* scaleconst, const bool cStart_isLog )
: FlxObjOutputBase(dolog,ostreamV), nameID(nameID), rbrvsets(rbrvsets), cStart(cStart), scaleconst(scaleconst), cStart_isLog(cStart_isLog)
{
  
}

FlxObjBayUp_new::~FlxObjBayUp_new()
{
  if (nameID) delete nameID;
  if (rbrvsets) delete rbrvsets;
  if (cStart) delete cStart;
  if (scaleconst) delete scaleconst;
}

void FlxObjBayUp_new::task()
{
  const std::string buID = nameID->eval_word(true);
  tdouble cStartV = ZERO;
  if (cStart_isLog) {
    cStartV = cStart->calc();
  } else {
    cStartV = log(cStart->cast2positive_or0(false));
  }
  flxBayUp* bu = new flxBayUp(buID,scaleconst->cast2positive(false),cStartV,rbrvsets->eval(true));
  GlobalVar.slogcout(4) << "BayUp: new updating object '" << buID << "' created." << std::endl;
  // add 'bu' to the box
    try {
      BayUpBox->insert(buID,bu,false);
    } catch (FlxException &e) {
      FLXMSG("FlxObjBayUp_new::task",1);
      if (bu) delete bu;
      throw;
    }
}

FlxObjReadBayUp_new::FlxObjReadBayUp_new()
{
  AllDefParaBox->insert(new FlxOptionalParaFlxString("nataf","bayup::rbrvsets",true));
    ParaBox.insert("rbrvsets", "bayup::rbrvsets" );
  AllDefParaBox->insert(new FlxOptionalParaFun(ZERO,"bayup::cstart"));
    ParaBox.insert("cstart", "bayup::cstart" );
  AllDefParaBox->insert(new FlxOptionalParaFun(ONE,"bayup::scaleconst"));
    ParaBox.insert("scaleconst", "bayup::scaleconst" );
  AllDefParaBox->insert(new FlxOptionalParaBool(false,"bayup::cstart_log"));
    ParaBox.insert("cstart_log", "bayup::cstart_log" );
}

FlxObjBase* FlxObjReadBayUp_new::read() {
  FlxString* nameID = new FlxString(false,false);
  try {
    read_optionalPara(false);
    return new FlxObjBayUp_new(get_doLog(),get_stream(),nameID, get_optPara_FlxString("rbrvsets"), get_optPara_FlxFunction("cstart"), get_optPara_FlxFunction("scaleconst"), get_optPara_bool("cstart_log") );
  } catch (FlxException &e) {
    FLXMSG("FlxObjReadBayUp_new::read",1);
    delete nameID;
    throw;
  }
}

FlxObjBayUp_likelihood::FlxObjBayUp_likelihood(const bool dolog, const std::string& ostreamV, FlxString* nameID, FlxFunction* lfun, const bool is_log)
: FlxObjOutputBase(dolog,ostreamV), nameID(nameID), lfun(lfun), is_log(is_log)
{

}

FlxObjBayUp_likelihood::~FlxObjBayUp_likelihood()
{
  delete nameID;
  delete lfun;
}

void FlxObjBayUp_likelihood::task()
{
  const std::string name = nameID->eval_word(true);
  flxBayUp& bu = BayUpBox->get(name);
  // create the name of the entry
    std::ostringstream ssV;
    ssV << name << "::" << bu.get_N_localLkl();
    const std::string nameE = ssV.str();
  RBRV_entry* lklEntry = NULL;
  if (is_log) {
    lklEntry = new RBRV_entry_fun_log(nameE,new FlxFunction(*lfun));
  } else {
    lklEntry = new RBRV_entry_fun(nameE,new FlxFunction(*lfun));
  }
  try {
    bu.add_localLkl(lklEntry);
  } catch (FlxException &e) {
    FLXMSG("FlxObjBayUp_likelihood::task",1);
    delete lklEntry;
    throw;
  }
}

FlxObjBayUp_likelihood_data::FlxObjBayUp_likelihood_data(const bool dolog, const std::string& ostreamV, FlxString* nameID, const tuint paraN, FlxString* isname, FlxFunction* lfun, const bool is_log)
: FlxObjOutputBase(dolog,ostreamV), nameID(nameID), paraN(paraN), isname(isname), lfun(lfun), is_log(is_log)
{
  #if FLX_DEBUG
    if (paraN==0) throw FlxException_Crude("FlxObjBayUp_likelihood_data::FlxObjBayUp_likelihood_data");
  #endif
}

FlxObjBayUp_likelihood_data::~FlxObjBayUp_likelihood_data()
{
  delete nameID;
  delete isname;
  delete lfun;
}

void FlxObjBayUp_likelihood_data::task()
{
  const std::string name = nameID->eval_word(true);
  flxBayUp& bu = BayUpBox->get(name);
  // create the name of the entry
    std::ostringstream ssV;
    ssV << name << "::" << bu.get_N_localLkl();
    const std::string nameE = ssV.str();
  // get the vector input stream
    const std::string isname_str = isname->eval_word(true);
    FlxIstream_vector *isv = dynamic_cast<FlxIstream_vector*>(&(data->IstreamBox.get(isname_str)));
    if (isv==NULL) {
      std::ostringstream ssV;
      ssV << "The input stream'" << isname_str << "' is not a vector-input stream.";
      throw FlxException_NeglectInInteractive("FlxObjBayUp_likelihood_data::task_1", ssV.str() );
    }
    isv->reset_stream();
    const tulong N = isv->get_total_size();
    if (N==0) {
      std::ostringstream ssV;
      ssV << "The vector input stream '" << isname_str << "' is empty.";
      throw FlxException_NeglectInInteractive("FlxObjBayUp_likelihood_data::task_2", ssV.str() );
    }
    if (N%paraN!=0) {
      std::ostringstream ssV;
      ssV << "The number of entries in the vector input stream '" << isname_str << "' cannot be divided by the number of parameters of the local likelihood function.";
      throw FlxException_NeglectInInteractive("FlxObjBayUp_likelihood_data::task_3", ssV.str() );
    }
  RBRV_entry* lklEntry = new RBRV_entry_fun_data(nameE,new FlxFunction(*lfun),paraN,isv,is_log);
  try {
    bu.add_localLkl(lklEntry);
  } catch (FlxException &e) {
    FLXMSG("FlxObjBayUp_likelihood_data::task_4",1);
    delete lklEntry;
    throw;
  }
}

FlxObjReadBayUp_likelihood::FlxObjReadBayUp_likelihood(): FlxObjReadOutputBase()
{
  AllDefParaBox->insert(new FlxOptionalParaBool(false,"bayup::log_likeli"));
    ParaBox.insert("log_likeli", "bayup::log_likeli" );
}

FlxObjBase* FlxObjReadBayUp_likelihood::read()
{
  FlxString* nameID = new FlxString(false,false);
  FlxFunction* lfun = NULL;
  FlxString* flxStr = NULL;
  try {
    // number of parameters
      tuint i1=0;
      if (reader->whatIsNextChar()=='(') {
        reader->getChar('(',false);
        FlxFunction* i1f = NULL;
        if (reader->whatIsNextChar()!=')') {
          try {
            i1f = new FlxFunction(funReader,false);
            i1 = i1f->cast2tuintW0(false);
            delete i1f;
          } catch (FlxException &e) {
            FLXMSG("FlxObjReadBayUp_likelihood::read_1",1);
            if (i1f) delete i1f; 
            throw;
          }
          if (i1>0) {                // data!!!!
            reader->getChar(',',false);
            flxStr = new FlxString(false,false);
          }
        }
        reader->getChar(')',false);
      }
    reader->getChar('=',false);
    if (i1==0) {                                // no data
      lfun = new FlxFunction(funReader,false);
      read_optionalPara(false);
      return new FlxObjBayUp_likelihood(get_doLog(),get_stream(),nameID, lfun, get_optPara_bool("log_likeli") );
    } else {                                        // with data
      FunReadPara::set_NumbOfPara(i1);
      try {
        lfun = new FlxFunction(funReader,false);
        read_optionalPara(false);
        FunReadPara::set_NumbOfPara(0);
      } catch (FlxException &e) {
        FLXMSG("FlxObjReadBayUp_likelihood::read_2",1);
        FunReadPara::set_NumbOfPara(0);
        throw;
      }
      read_optionalPara(false);
      return new FlxObjBayUp_likelihood_data(get_doLog(),get_stream(),nameID, i1, flxStr, lfun, get_optPara_bool("log_likeli") );
    }
  } catch (FlxException &e) {
    FLXMSG("FlxObjReadBayUp_likelihood::read_3",1);
    delete nameID;
    if (lfun) delete lfun;
    if (flxStr) delete flxStr;
    throw;
  }
}

FlxObjBayUp_uncertobsv::FlxObjBayUp_uncertobsv(const bool dolog, const std::string& ostreamV, FlxString* nameID_BU, const tuint paraN, FlxString* isname, FlxFunction* lfun, FlxString* set_name, std::vector< RBRV_entry_read_base* >& set_entries, const bool is_log)
: FlxObjOutputBase(dolog,ostreamV), nameID_BU(nameID_BU), paraN(paraN), isname(isname), lfun(lfun), set_name(set_name), set_entries(set_entries), is_log(is_log)
{

}

FlxObjBayUp_uncertobsv::~FlxObjBayUp_uncertobsv()
{
  delete nameID_BU;
  delete isname;
  delete lfun;
  delete set_name;
  for (tuint i=0;i<set_entries.size();++i) delete set_entries[i];
}

void FlxObjBayUp_uncertobsv::task()
{
  const std::string nameBU = nameID_BU->eval_word(true);
  flxBayUp& bu = BayUpBox->get(nameBU);
  // get the vector input stream
    const std::string isname_str = isname->eval_word(true);
    FlxIstream_vector *isv = dynamic_cast<FlxIstream_vector*>(&(data->IstreamBox.get(isname_str)));
    if (isv==NULL) {
      std::ostringstream ssV;
      ssV << "The input stream'" << isname_str << "' is not a vector-input stream.";
      throw FlxException_NeglectInInteractive("FlxObjBayUp_uncertobsv::task_1", ssV.str() );
    }
    isv->reset_stream();
    const tulong N = isv->get_total_size();
    if (N==0) {
      std::ostringstream ssV;
      ssV << "The vector input stream '" << isname_str << "' is empty.";
      throw FlxException_NeglectInInteractive("FlxObjBayUp_likelihood_data::task_2", ssV.str() );
    }
    if (N%paraN!=0) {
      std::ostringstream ssV;
      ssV << "The number of entries in the vector input stream '" << isname_str << "' cannot be divided by the number of parameters of the local likelihood function.";
      throw FlxException_NeglectInInteractive("FlxObjBayUp_likelihood_data::task_3", ssV.str() );
    }
    const tuint Nobsv = N/paraN;
  // Define the RBRV set (the single set)
    const std::string setnameStr = nameBU + "::" + set_name->eval_word(true);
    std::vector<FlxString*> set_parents;                        // just a dummy! -> no parents are allowed
    FlxObjRBRV_set_creator crtr(setnameStr,NULL,0,false,set_entries);
    RBRV_set* ts_single = crtr.register_set_rbrv(data->rbrv_box,false);
    try {
      data->rbrv_box.register_set(ts_single);
    } catch (FlxException& e) {
      FLXMSG("FlxObjBayUp_uncertobsv::task_4",1);
      delete ts_single;
      throw;
    }
    const std::string setnameStr2 = setnameStr+"::lkli";
    flxBayUp_uncertobsv_set* ts = new flxBayUp_uncertobsv_set(setnameStr2,ts_single,new FlxFunction(*lfun),Nobsv,paraN,isv,is_log);
    try {
      data->rbrv_box.register_set(ts);
    } catch (FlxException &e) {
      FLXMSG("FlxObjBayUp_uncertobsv::task_5",1);
      delete ts;
      throw;
    }
    bu.add_localLkl(ts);
}

FlxObjReadBayUp_uncertobsv::FlxObjReadBayUp_uncertobsv(): FlxObjReadOutputBase()
{
  AllDefParaBox->insert(new FlxOptionalParaBool(false,"bayup::log_likeli"));
    ParaBox.insert("log_likeli", "bayup::log_likeli" );
}

FlxObjBase* FlxObjReadBayUp_uncertobsv::read()
{
  FlxString* nameID_BU = new FlxString(false,false);        // ID of the Bayesian updating object
  FlxFunction* lfun = NULL;                                // the Likelihood function
  FlxString* IString = NULL;                                // name of the input stream
  FlxString* set_name = NULL;                                // name of the current RBRV set
  std::vector<FlxString*> set_parents;                        // just a dummy! -> no parents are allowed
  std::vector<RBRV_entry_read_base*> set_entries;        // the RBRVs
  try {
    reader->getChar('(',false);
    set_name = new FlxString(false,false);
    reader->getChar(',',false);
    FlxFunction* i1f = new FlxFunction(funReader,false);
    tuint i1 = 0;
    try {
      i1 = i1f->cast2tuint(false);
      delete i1f;
    } catch (FlxException &e) {
      FLXMSG("FlxObjReadBayUp_uncertobsv::read_1",1);
      delete i1f;
      throw;
    }
    reader->getChar(',',false);
    IString = new FlxString(false,false);   
    reader->getChar(')',false);
    FunReadPara::set_NumbOfPara(i1);
    try {
    // read the RBRV-set
      RBRV_entry_read_base::read(set_entries,set_parents,false);
      if (set_parents.empty()==false) {
        std::ostringstream ssV;
        ssV << "The definition of RBRV-parent-sets is not allowed within this framework.";
        throw FlxException_NeglectInInteractive("FlxObjReadBayUp_uncertobsv::read_2", ssV.str() );
      }
    // read the actual likelihood
      reader->getChar('=',false);
        lfun = new FlxFunction(funReader,false);
        FunReadPara::set_NumbOfPara(0);
    } catch (FlxException &e) {
      FLXMSG("FlxObjReadBayUp_uncertobsv::read_3",1);
      FunReadPara::set_NumbOfPara(0);
      throw;
    }
    read_optionalPara(false);
    return new FlxObjBayUp_uncertobsv(get_doLog(),get_stream(),nameID_BU,i1,IString,lfun,set_name,set_entries,get_optPara_bool("log_likeli"));
  } catch (FlxException &e) {
    FLXMSG("FlxObjReadBayUp_uncertobsv::read_4",1);
    delete nameID_BU;
    if (lfun) delete lfun;
    if (IString) delete IString;
    if (set_name) delete set_name;
    for (tuint i=0;i<set_parents.size();++i) delete set_parents[i];
    for (tuint i=0;i<set_entries.size();++i) delete set_entries[i];
    throw;
  }
}


FlxObjBayUp_glbllikelihood::FlxObjBayUp_glbllikelihood(const bool dolog, const std::string& ostreamV, FlxString* nameID, FlxFunction* lfun, const bool is_log, const flxBayUp::MethCategory methCat)
: FlxObjOutputBase(dolog,ostreamV), nameID(nameID), lfun(lfun), is_log(is_log), methCat(methCat)
{

}

FlxObjBayUp_glbllikelihood::~FlxObjBayUp_glbllikelihood()
{
  delete nameID;
  delete lfun;
}

void FlxObjBayUp_glbllikelihood::task()
{
  const std::string name = nameID->eval_word(true);
  flxBayUp& bu = BayUpBox->get(name);
  bu.set_globalLkl(*lfun,is_log,methCat);
}

FlxObjReadBayUp_glbllikelihood::FlxObjReadBayUp_glbllikelihood(const flxBayUp::MethCategory methCat)
: FlxObjReadOutputBase(), methCat(methCat)
{
  AllDefParaBox->insert(new FlxOptionalParaBool(false,"bayup::log_likeli"));
    ParaBox.insert("log_likeli", "bayup::log_likeli" );
}

FlxObjBase* FlxObjReadBayUp_glbllikelihood::read()
{
  FlxString* nameID = new FlxString(false,false);
  FlxFunction* lfun = NULL;
  try {
    reader->getChar('(');
    lfun = new FlxFunction(funReader,false);
    reader->getChar(')');
    read_optionalPara(false);
    bool log_likeli;
    if (methCat==flxBayUp::RA) {
      log_likeli = false;
    } else {
      log_likeli = get_optPara_bool("log_likeli");
    }
    return new FlxObjBayUp_glbllikelihood(get_doLog(),get_stream(),nameID, lfun, log_likeli, methCat );
  } catch (FlxException &e) {
    FLXMSG("FlxObjReadBayUp_glbllikelihood::read_2",1);
    delete nameID;
    if (lfun) delete lfun;
    throw;
  }
}


FlxObjBayUp_update::FlxObjBayUp_update(const bool dolog, const std::string& ostreamV, FlxString* nameID, FlxFunction* Nc, FlxFunction* Ncl, FlxFunction* Nburn, FlxFunction* Ns_final, FlxFunction* max_runs, const FlxBayUp_Update_List::randomizeTec randomize, flxBayUp_adaptive_ctrl_base* adpt_ctrl, const bool use_cStart, const FlxBayUp_Update_List::MethType meth_id, const bool log_LSF, const Flx_SuS_Control& susControl, SuS_csm_evalStorage* csm_eval)
:FlxObjSuS(dolog,ostreamV,Nc,Ncl,max_runs,randomize,adpt_ctrl,susControl,csm_eval,NULL,NULL),
  nameID(nameID),Ns_final(Ns_final), Nburn(Nburn), use_cStart(use_cStart), meth_id(meth_id), log_LSF(log_LSF)
{
}

FlxObjBayUp_update::~FlxObjBayUp_update()
{
  delete nameID;
  delete Ns_final;
  if (Nburn) delete Nburn;
}

void FlxObjBayUp_update::task()
{
  const std::string name = nameID->eval_word(true);
  flxBayUp& bu = BayUpBox->get(name);
    if (meth_id==FlxBayUp_Update_List::TMCMC) bu.set_TMCMC();
  FlxBayUp_Update_List* list = NULL;
  const tuint NsfV = Ns_final->cast2tuintW0(false);
  // allocate list
    if (adpt_ctrl) {
      const tuint mr = max_runs->cast2tuint(false);
      const tuint NcV = Nc->cast2tuint(false);
      const tuint NclV = Ncl->cast2tuint(false);
      const tuint NburnV = Nburn?(Nburn->cast2tuintW0(false)):0;
      list = new FlxBayUp_Update_List(bu,NcV,NclV,NsfV,NburnV,randomize,adpt_ctrl->copy(),mr,use_cStart,meth_id,log_LSF,susControl.find_multiples);
    } else {
      tuint Nburn = 0;
      if (meth_id==FlxBayUp_Update_List::LS) {
        Nburn = Nc->cast2tuint(false);
      }
      list = new FlxBayUp_Update_List(bu,0,0,NsfV,Nburn,randomize,NULL,0,use_cStart,meth_id,log_LSF,susControl.find_multiples);
    }
  // allocate csm
    FlxBayUP_csm_base* csm = NULL;
    if (
      meth_id==FlxBayUp_Update_List::BUS
      ||meth_id==FlxBayUp_Update_List::UBUS
      ||meth_id==FlxBayUp_Update_List::BUST
      ||meth_id==FlxBayUp_Update_List::ABCSUBSIM
      ||meth_id==FlxBayUp_Update_List::RASUBSIM
      ||meth_id==FlxBayUp_Update_List::TMCMC
    )
    {
      csm = csm_eval->eval(list);
    }
  bu.updater.update(list,csm,use_cStart,susControl);
}

FlxObjReadBayUp_update::FlxObjReadBayUp_update()
{  
  FlxBayUp_Update::define_constants();
  AllDefParaBox->insert(new FlxOptionalParaBool(false,"bayup::use_cstart"));
    ParaBox.insert("use_cstart", "bayup::use_cstart" );
  AllDefParaBox->insert(new FlxOptionalParaBool(false,"bayup::log_lsf"));
    ParaBox.insert("log_lsf", "bayup::log_lsf" );
}

FlxObjBase* FlxObjReadBayUp_update::read()
{
  FlxFunction* funNc=NULL; FlxFunction* funNcl=NULL; FlxFunction* funNburn=NULL; FlxFunction* Ns_final=NULL;
  flxBayUp_adaptive_ctrl_base* adpt_ctrl = NULL;
  SuS_csm_evalStorage* csm_eval = NULL;
  FlxString* nameID = new FlxString(false,false);
  try {
    reader->getChar(':',false);
    // get method-ID
      FlxBayUp_Update_List::MethType meth_id;
      {
        const std::string methStr = reader->getWord(true,false);
        if (methStr=="bus") meth_id = FlxBayUp_Update_List::BUS;
        else if (methStr=="ubus") meth_id = FlxBayUp_Update_List::UBUS;
        else if (methStr=="bust") meth_id = FlxBayUp_Update_List::BUST;
        else if (methStr=="abcsubsim") meth_id = FlxBayUp_Update_List::ABCSUBSIM;
        else if (methStr=="rasubsim") meth_id = FlxBayUp_Update_List::RASUBSIM;
        else if (methStr=="tmcmc") meth_id = FlxBayUp_Update_List::TMCMC;
        else if (methStr=="rs") meth_id = FlxBayUp_Update_List::RS;
        else if (methStr=="mhrs") meth_id = FlxBayUp_Update_List::MHRS;
        else if (methStr=="abcrs") meth_id = FlxBayUp_Update_List::ABCRS;
        else if (methStr=="ramci") meth_id = FlxBayUp_Update_List::RAMCI;
        else if (methStr=="ls") meth_id = FlxBayUp_Update_List::LS;
        else {
          std::ostringstream ssV;
          ssV << "Unknown method-ID (" << methStr << ") for updating method.";
          throw FlxException_NeglectInInteractive("FlxObjReadBayUp_update::read_1", ssV.str() );
        }
      }
    // get number of samples
      reader->getChar('(',false);
      switch (meth_id) {
        case (FlxBayUp_Update_List::BUS):
        case (FlxBayUp_Update_List::UBUS):
        case (FlxBayUp_Update_List::BUST):
        case (FlxBayUp_Update_List::ABCSUBSIM):
        case (FlxBayUp_Update_List::RASUBSIM):
        case (FlxBayUp_Update_List::TMCMC):
        case (FlxBayUp_Update_List::LS):
          funNc = new FlxFunction(funReader,false);
          reader->getChar(';',false);
          if (meth_id!=FlxBayUp_Update_List::LS) {
            if (meth_id==FlxBayUp_Update_List::TMCMC) {
              funNcl = new FlxFunction(new FunNumber(ONE));
              funNburn = new FlxFunction(funReader,false);
            } else {
              funNcl = new FlxFunction(funReader,false);
            }
            reader->getChar(';',false);
          } else {
            funNcl = new FlxFunction(new FunNumber(ONE));
          }
        case (FlxBayUp_Update_List::MHRS):
        case (FlxBayUp_Update_List::RS):
        case (FlxBayUp_Update_List::ABCRS):
        case (FlxBayUp_Update_List::RAMCI):
          Ns_final = new FlxFunction(funReader,false);
          break;
        default:
          throw FlxException_Crude("FlxObjReadBayUp_update::read_2");
      };
      reader->getChar(')',false);
    read_optionalPara(false);
    // get randomize-ID
      FlxBayUp_Update_List::randomizeTec rnd_id;
      switch (meth_id) {
        case (FlxBayUp_Update_List::BUS):
        case (FlxBayUp_Update_List::UBUS):
        case (FlxBayUp_Update_List::BUST):
        case (FlxBayUp_Update_List::ABCSUBSIM):
        case (FlxBayUp_Update_List::RASUBSIM):
        case (FlxBayUp_Update_List::TMCMC):
        {
          rnd_id = get_randomize_id();
          break;
        }
        case (FlxBayUp_Update_List::MHRS):
        {
          rnd_id = FlxBayUp_Update_List::FlxBayUp_Update_List::INIT;
          break;
        }
        case (FlxBayUp_Update_List::RS):
        case (FlxBayUp_Update_List::ABCRS):
        case (FlxBayUp_Update_List::RAMCI):
        case (FlxBayUp_Update_List::LS):
          rnd_id = FlxBayUp_Update_List::NONE;
          break;
        default:
          throw FlxException_Crude("FlxObjReadBayUp_update::read_4");
      };
    // get adaptive information
      switch (meth_id) {
        case (FlxBayUp_Update_List::BUS):
        case (FlxBayUp_Update_List::UBUS):
        case (FlxBayUp_Update_List::BUST):
        case (FlxBayUp_Update_List::ABCSUBSIM):
        case (FlxBayUp_Update_List::RASUBSIM):
        case (FlxBayUp_Update_List::TMCMC):
        {
          adpt_ctrl = get_adpt_ctrl();
          break;
        } 
        default:
          break;
      };
    csm_eval = get_csm_eval();
    
    Flx_SuS_Control susControl = get_susControl();
    
    return new FlxObjBayUp_update(get_doLog(),get_stream(),nameID,funNc,funNcl,funNburn,Ns_final,get_optPara_FlxFunction("max_runs"),
                                  rnd_id,adpt_ctrl,get_optPara_bool("use_cstart"),meth_id,get_optPara_bool("log_lsf"),susControl,csm_eval);
  } catch (FlxException &e) {
    FLXMSG("FlxObjReadBayUp_update::read_7",1);
    delete nameID;
    if (funNc) delete funNc;
    if (funNcl) delete funNcl;
    if (funNburn) delete funNburn;
    if (Ns_final) delete Ns_final;
    if (adpt_ctrl) delete adpt_ctrl;
    if (csm_eval) delete csm_eval;
    throw;
  }
}

const tdouble FunBayUp_Prop::calc_help(const tuint PID)
{
  switch (PID) {
    case 1:
      return buo.updater.get_PrMod();
    case 2:
      return exp(calc_help(5));
    case 3:
      return buo.updater.get_total_LSFc_sim();
    case 4:
      return buo.updater.get_LogPrMod();
    case 5:
      if (buo.updater.get_list().meth_id==FlxBayUp_Update_List::MHRS) {
        return buo.updater.get_list().get_maxL();
      } else {
        return buo.get_cStart();
      }
    case 6:
      return buo.updater.get_LogDataFit();
    case 7:
      return buo.updater.get_LogRelEntr();
    default:
    {
      std::ostringstream ssV;
      ssV << "Unknown property ID (" << PID << ")";
      throw FlxException("FunBayUp_Prop::calc", ssV.str() ); 
    }
  }
}

const tdouble FunBayUp_Prop::calc()
{
  const tuint PID = tuint_from(pidf->calc(),"Property ID",true,true,pidf);
  return calc_help(PID);
}

const std::string FunBayUp_Prop::write()
{
  return "bayup_prop(" + buo.get_name() + "," + pidf->write() + ")";
}

FunBayUp_Prop::~FunBayUp_Prop()
{
  delete pidf;
}

FunBase* FunReadFunBayUp_Prop::read(bool errSerious)
{
  flxBayUp& buo = BayUpBox->get(reader->getWord(true,errSerious));
  reader->getChar(',');
  FunBase* pidf = FunctionList->read(errSerious);
  return new FunBayUp_Prop(buo,pidf);
}

FunBase* FunReadFunBayUp_lsf::read(bool errSerious)
{
  flxBayUp& buo = BayUpBox->get(reader->getWord(true,errSerious));
  buo.freeze();
  return new FunBayUp_lsf(buo);
}

const tdouble FunBayUp_lsf::calc()
{
  switch (buo.get_methCat()) {
    case (flxBayUp::BUS):
    {
      const tdouble L = buo.eval_Likelihood();
      const tdouble p = buo.get_p();
      const tdouble c = buo.get_cStart();
      if (L>c) {
        throw FlxException("FunBayUp_lsf::calc_1","Likelihood larger than scaling constant encountered.");
      }
      return p-rv_InvPhi_noAlert(exp(L-c));
    }
    case (flxBayUp::RA):
      return buo.eval_RAlsf();
    default:
      throw FlxException("FunBayUp_lsf::calc_2","Operation not allowed for this type of updating problem");
  }
}

const std::string FunBayUp_lsf::write()
{
  return "bayup_lsf(" + buo.get_name() + ")";
}

FlxObjBayUp_Set::FlxObjBayUp_Set(const bool dolog, FlxString* setname, const std::vector< FlxString* >& BUname_vec, const std::vector< FlxFunction* >& priorWeight_vec, std::vector< FlxString* > model_res_list_Str, std::vector< FlxFunction** > model_res_map)
: FlxObjBase(dolog), setname(setname), BUname_vec(BUname_vec), priorWeight_vec(priorWeight_vec), Nmodels(BUname_vec.size()), 
  model_res_list_Str(model_res_list_Str), model_res_map(model_res_map), Nvalues(model_res_list_Str.size())
{
  
}

FlxObjBayUp_Set::~FlxObjBayUp_Set()
{
  delete setname;
  for (tuint i=0;i<Nmodels;++i) {
    delete BUname_vec[i];
    delete priorWeight_vec[i];
  }
  for (tuint i=0;i<Nvalues;++i) {
    delete model_res_list_Str[i];
    for (tuint j=0;j<Nmodels;++j) {
      delete model_res_map[i][j];
    }
    delete [] model_res_map[i];
  }
}

void FlxObjBayUp_Set::task()
{
  const std::string name = setname->eval_word(true);
  flxBayUp_mProb_set* bu_set = NULL;
  flxVec v(Nmodels);
  std::vector<std::string> mrl;
  FlxFunction** mrm = NULL;
  flxBayUp_RBRV_set** modelVec = new flxBayUp_RBRV_set*[Nmodels];
  bool no_free = false;
  try {
    for (tuint i=0;i<Nmodels;++i) {
      const std::string sn = BUname_vec[i]->eval_word(true);
      RBRV_set_base* e = data->rbrv_box.get_set(sn,true);
      flxBayUp_RBRV_set* be = dynamic_cast<flxBayUp_RBRV_set*>(e);
      if (be==NULL) {
        std::ostringstream ssV;
        ssV << "The set '" << sn << "' is not (at least not directly) derived from an updating object.";
        throw FlxException("FlxObjBayUp_Set::task_1", ssV.str() );
      }
      modelVec[i] = be;
      v[i] = priorWeight_vec[i]->cast2positive();
    }
    tuint c = 0;
    mrm = new FlxFunction*[Nvalues*Nmodels];
      for (tuint i=0;i<Nvalues;++i) {
        for (tuint j=0;j<Nmodels;++j) {
          mrm[c] = NULL;
          ++c;
        }
      }
    c = 0;
    for (tuint i=0;i<Nvalues;++i) {
      mrl.push_back(model_res_list_Str[i]->eval_word(true));
      for (tuint j=0;j<Nmodels;++j) {
        mrm[c] = new FlxFunction( *(model_res_map[i][j]) );
        ++c;
      }
    }
    no_free = true;
    bu_set = new flxBayUp_mProb_set(false,name,Nmodels,modelVec,v,Nvalues,mrl,mrm);
    data->rbrv_box.register_set(bu_set);
  } catch (FlxException& e) {
    FLXMSG("FlxObjBayUp_Set::task_2",1);
    if (!no_free) delete [] modelVec;
    if (bu_set) delete bu_set;
    if (!no_free) if (mrm) {
      tuint c = 0;
      for (tuint i=0;i<Nvalues;++i) {
        for (tuint j=0;j<Nmodels;++j) {
          if (mrm[c]) delete mrm[c];
          ++c;
        }
      }
      delete [] mrm;
    }
    throw;
  }
}

FlxObjBase* FlxObjReadBayUp_Set::read()
{
  FlxString* setname = new FlxString(false,false);
  std::vector<FlxString*> BUname_vec;
  std::vector<FlxFunction*> priorWeight_vec;
  std::vector<FlxString*> model_res_list_Str;
  std::vector<FlxFunction**> model_res_map;    
  FlxFunction** tfv = NULL;
  try {    
    reader->getChar('{',false);
    while (true) {
      BUname_vec.push_back(new FlxString(false,false));
      if (reader->whatIsNextChar()=='(') {
        reader->getChar('(');
        priorWeight_vec.push_back(new FlxFunction(funReader,false));
        reader->getChar(')');
      } else {
        priorWeight_vec.push_back(new FlxFunction(new FunNumber(1.)));
      }
      if (reader->whatIsNextChar()=='}') {
        break;
      } else {
        reader->getChar(',');
      }
    }
    reader->getChar('}',false);
    reader->getChar('{',false);
    while (true) {
      model_res_list_Str.push_back(new FlxString(false,false));
      reader->getChar('=');
      reader->getChar('(');
      tfv = new FlxFunction*[BUname_vec.size()];
      for (tuint i=0;i<BUname_vec.size();++i) tfv[i] = NULL;
      for (tuint i=0;i<BUname_vec.size();++i) {
        if (i>0) reader->getChar(',');
        tfv[i] = new FlxFunction(funReader,false);
      }
      model_res_map.push_back(tfv);
      tfv = NULL;
      reader->getChar(')');
      if (reader->whatIsNextChar()=='}') {
        break;
      } else {
        reader->getChar(',');
      }
    }
    reader->getChar('}',false);
    read_optionalPara(false);
    return new FlxObjBayUp_Set(get_doLog(),setname,BUname_vec,priorWeight_vec,model_res_list_Str,model_res_map);
  } catch (FlxException &e) {
    FLXMSG("FlxObjReadBayUp_Set::read",1);
    delete setname;
    for (tuint i=0;i<BUname_vec.size();++i) delete BUname_vec[i];
    for (tuint i=0;i<priorWeight_vec.size();++i) delete priorWeight_vec[i];
    for (tuint i=0;i<model_res_list_Str.size();++i) delete model_res_list_Str[i];
    if (tfv) {
      for (tuint i=0;i<BUname_vec.size();++i) if (tfv[i]) delete tfv[i];
      delete [] tfv;
    }
    for (tuint i=0;i<model_res_map.size();++i) {
      for (tuint j=0;j<BUname_vec.size();++j) {
        delete model_res_map[i][j];
      }
      delete [] model_res_map[i];
    }
    throw;
  }
}

void FlxObjBayUp_Reset_Smpls::task()
{
  const std::string name = nameID->eval_word(true);
  flxBayUp& bu = BayUpBox->get(name);
  bu.updater.reset_finalized_smpls();
}

FlxObjBase* FlxObjReadBayUp_Reset_Smpls::read()
{
  FlxString* nameID = new FlxString(false,false);
  read_optionalPara(false);
  return new FlxObjBayUp_Reset_Smpls(get_doLog(),nameID);
}


FunConvExp::FunConvExp(FlxString* GM_name, FlxString* DO_name, FlxString* name_dist_set, FlxString* name_pulse_set, const tuint seed, const tuint Ninteg, const tdouble eps, const tuint Nrgir)
: dist_set(NULL), pulse_set(NULL), GM_name(GM_name), DO_name(DO_name), name_dist_set(name_dist_set), name_pulse_set(name_pulse_set), seed(seed), Ninteg(Ninteg), eps(eps), Nrgir(Nrgir),
  N(0), cv(NULL), y(NULL), eps_Mo(NULL), eps_me(NULL), Odif(NULL), DOp(NULL), GMp(NULL)
{

}

FunConvExp::~FunConvExp()
{
  if (GM_name) delete GM_name;
  if (DO_name) delete DO_name;
  if (name_dist_set) delete name_dist_set;
  if (name_pulse_set) delete name_pulse_set;
  if (cv) delete cv;
  if (y) delete y;
  if (eps_Mo) delete eps_Mo;
  if (eps_me) delete eps_me;
  if (Odif) delete Odif;
}

const tdouble FunConvExp::get_pulse_log()
{
  *eps_me = *Odif;
  *eps_me -= *eps_Mo;
  pulse_set->set_is_valid(false);
  pulse_set->set_x_only_this(eps_me->get_tmp_vptr_const());
  return pulse_set->get_pdf_x_eval_log();
}

const tdouble FunConvExp::compute_cv()
{
  // define the Golden section ;)
    const double resphi = 2*ONE - (ONE+sqrt(5*ONE))/2;
  // assign the starting configuration
    tdouble pa = ZERO;
      dist_set->set_is_valid(false);
      y->set_zero();
      dist_set->set_y_only_this(y->get_tmp_vptr_const());
      dist_set->transform_y2x();
      dist_set->get_x_only_this(eps_Mo->get_tmp_vptr());
      tdouble fa = dist_set->get_pdf_x_eval_log();
      fa += get_pulse_log();
    tdouble pc = ONE;
      dist_set->set_is_valid(false);
      *eps_Mo = *Odif;
      dist_set->set_x_only_this(eps_Mo->get_tmp_vptr_const());
      dist_set->transform_x2y();
      dist_set->get_y_only_this(cv->get_tmp_vptr());
      tdouble fc = dist_set->get_pdf_x_eval_log();
      fc += get_pulse_log();
      tdouble cvl = cv->get_Norm2();
    tdouble pb = pa + resphi*(pc-pa);
      dist_set->set_is_valid(false);
      *y = *cv; *y *= pb;
      dist_set->set_y_only_this(y->get_tmp_vptr_const());
      dist_set->transform_y2x();
      dist_set->get_x_only_this(eps_Mo->get_tmp_vptr());
      tdouble fb = dist_set->get_pdf_x_eval_log();
      fb += get_pulse_log();
      tdouble px, fx;
  while (true) {
    // check if we can already stop
      if ((pc-pa)*cvl<=eps) {
        pb = (pa+pc)/2;
        *cv *= pb;
        // transform cv to original space
          dist_set->set_is_valid(false);
          dist_set->set_y_only_this(cv->get_tmp_vptr_const());
          dist_set->transform_y2x();
          dist_set->get_x_only_this(cv->get_tmp_vptr());
        return pb;
      }
    // get a new point
      if ( (pc-pb) > (pb-pa) ) {
        px = pb + resphi*(pc-pb);
      } else {
        px = pb - resphi*(pb-pa);
      }
      *y = *cv; *y *= px;
      dist_set->set_is_valid(false);
      dist_set->set_y_only_this(y->get_tmp_vptr_const());
      dist_set->transform_y2x();
      dist_set->get_x_only_this(eps_Mo->get_tmp_vptr());
      fx = dist_set->get_pdf_x_eval_log();
      fx += get_pulse_log();
    // choose the next interval
      if (fx>fb) {
        if ( (pc-pb) > (pb-pa) ) {
          pa = pb; fa = fb;
        } else {
          pc = pb; fc = fb;
        }
        pb = px; fb = fx;
      } else {
        if ( (pc-pb) > (pb-pa) ) {
          pc = px; fc = fx;
        } else {
          pa = px; fa = fx;
        }
      }
  }
}

const tdouble FunConvExp::calc()
{
  // retrieve the sets
    if (dist_set==NULL) {
      dist_set = data->rbrv_box.get_set(name_dist_set->eval_word(true),true);
      if (pulse_set) throw FlxException_Crude("FlxObjConvExp::task_1");
      pulse_set = data->rbrv_box.get_set(name_pulse_set->eval_word(true),true);
      if (dist_set->get_NOX_only_this()!=pulse_set->get_NOX_only_this()) {
        std::ostringstream ssV;
        ssV << "The dimensions of the two sets are not identical.";
        throw FlxException("FlxObjConvExp::task_2", ssV.str() );
      }
      if (
        dist_set->get_NRV_only_this()!=dist_set->get_NOX_only_this()
        || pulse_set->get_NRV_only_this()!=pulse_set->get_NOX_only_this()
      ) {
        std::ostringstream ssV;
        ssV << "Such a RBRV-set is not allowed in this context.";
        throw FlxException("FlxObjConvExp::task_3", ssV.str() );
      }
      delete name_dist_set; name_dist_set = NULL;
      delete name_pulse_set; name_pulse_set = NULL;
      N = dist_set->get_NOX_only_this();
      cv = new flxVec(N);
      y = new flxVec(N);
      eps_Mo = new flxVec(N);
      eps_me = new flxVec(N);
      Odif = new flxVec(N);
    }
  // obtain the vectors we have to work with
    GMp = data->ConstMtxBox.get_Vec(GM_name->eval_word(true),N,true);
    DOp = data->ConstMtxBox.get_Vec(DO_name->eval_word(true),N,true);
    *Odif = flxVec(DOp,N);
    *Odif -= flxVec(GMp,N);
  // compute the center of the IS density
    compute_cv();
  // approximate the covariance matrix of the IS density
    // dist_set
      FlxMtxSym* cov_dist = new FlxMtxSym(N);
      try {
      dist_set->add_covMTX(*cov_dist);
      cov_dist->Invert();
    // pulse_set
      FlxMtxSym cov_pulse(N);
      pulse_set->add_covMTX(cov_pulse);
      cov_pulse.Invert();
      (*cov_dist) += cov_pulse;
    cov_dist->Invert();
      } catch (FlxException& e) {
        FLXMSG("FlxObjConvExp::task_4",1);
        delete cov_dist;
        throw;
      }
  // Define the importance sampling density
    RBRV_set_MVN mvn(false,N,static_cast<tuint>(0),std::string("ipstmpinternal"),true,new flxVec(*cv),cov_dist,2);
  // initialize the random number generator
    boost::random::mt19937 rgL;
    rv_initialize(false,true,seed,Nrgir,&rgL,false);
  // perform the sampling
    qdouble ipv(Ninteg,false);
    tdouble* yp = y->get_tmp_vptr();
    for (tuint i=0;i<Ninteg;++i) {
      // generate sample in standard normal space
        for (tuint j=0;j<N;++j) {
          yp[j] = rv_normal(rgL);
        }
      // transform sample to original space
        mvn.set_is_valid(false);
        mvn.set_y_only_this(yp);
        mvn.transform_y2x();
        mvn.get_x_only_this(eps_Mo->get_tmp_vptr());
      // evaluate the integrand
        tdouble t = get_pulse_log();
        dist_set->set_is_valid(false);
        dist_set->set_x_only_this(eps_Mo->get_tmp_vptr_const());
        t += dist_set->get_pdf_x_eval_log();
        t -= mvn.get_pdf_x_eval_log();
        t = exp(t);
      ipv.operator+=(t);
    }
  return ((ipv.cast2double())/Ninteg);
}

const bool FunConvExp::search_circref(FlxFunction* fcr)
{
  return name_dist_set->search_circref(fcr) || name_pulse_set->search_circref(fcr) || GM_name->search_circref(fcr) || DO_name->search_circref(fcr);
}

const bool FunConvExp::dependOn_Const(const tdouble*const thenumber)
{
  return false;
}

const std::string FunConvExp::write()
{
  std::ostringstream ssV;
  ssV << "convexp(";
  if (name_dist_set) ssV << name_dist_set->write();
  else ssV << dist_set->name;
  ssV << ",";
  if (name_pulse_set) ssV << name_pulse_set->write();
  else ssV << pulse_set->name;
  ssV << ","
  << GM_name->write() << ","
  << DO_name->write() << ","
  << "seed=" << GlobalVar.Double2String(seed) << ","
  << "n=" << GlobalVar.Double2String(Ninteg) << ","
  << "eps=" << GlobalVar.Double2String(eps)
  << ")";
  return ssV.str();
}


FunBase* FunReadFunConvExp::read(bool errSerious)
{
  FlxString* GM_name = NULL;
  FlxString* DO_name = NULL;
  FlxString* name_dist_set = NULL;
  FlxString* name_pulse_set = NULL;
  FlxFunction* f = NULL;
  // optional parameters
    tuint seed =        602737839; 
    tuint Ninteg =        1e4; 
    tdouble eps =        0.01;
    tuint rgir =        1000;
  try {
//     reader->getChar('(',false);
    GM_name = new FlxString(false,false);
    reader->getChar(',',false);
    DO_name = new FlxString(false,false);
    reader->getChar(',',false);
    name_dist_set = new FlxString(false,false);
    reader->getChar(',',false);
    name_pulse_set = new FlxString(false,false);
    // seed
    if (reader->whatIsNextChar()==',') {
      reader->getChar(',',errSerious);
      reader->getWord("seed",errSerious);
      reader->getChar('=');
      f = new FlxFunction(FunctionList->read(errSerious));
      seed = f->cast2tuintW0(false);
      delete f; f = NULL;
      // N
      if (reader->whatIsNextChar()==',') {
        reader->getChar(',',errSerious);
        reader->getWord("n",errSerious);
        reader->getChar('=');
        f = new FlxFunction(FunctionList->read(errSerious));
        Ninteg = f->cast2tuint(false);
        delete f; f = NULL;
        // eps
        if (reader->whatIsNextChar()==',') {
          reader->getChar(',',errSerious);
          reader->getWord("eps",errSerious);
          reader->getChar('=');
          f = new FlxFunction(FunctionList->read(errSerious));
          eps = f->cast2positive(false);
          delete f; f = NULL;
          // rgir
          if (reader->whatIsNextChar()==',') {
            reader->getChar(',',errSerious);
            reader->getWord("rgir",errSerious);
            reader->getChar('=');
            f = new FlxFunction(FunctionList->read(errSerious));
            rgir = f->cast2tuintW0(false);
            delete f; f = NULL;
          }
        }
      }
    }
//     reader->getChar(')',false);
  } catch (FlxException& e) {
    FLXMSG("FunReadFunConvExp::read",1);
    if (GM_name) delete GM_name;
    if (DO_name) delete DO_name;
    if (name_dist_set) delete name_dist_set;
    if (name_pulse_set) delete name_pulse_set;
    if (f) delete f;
    throw;
  }
  return new FunConvExp(GM_name,DO_name,name_dist_set,name_pulse_set,seed,Ninteg,eps,rgir);
}

//----------------------------------------------------------------------------------------------

FlxObjBayDA_new::FlxObjBayDA_new( const bool dolog, const std::string& ostreamV, FlxString* nameID, FlxMtxConstFun* mcf_data, FlxFunction* id_transform,
                        FlxFunction* Nchain, FlxFunction* Nburn, FlxFunction* Ntune, FlxFunction* Npost, FlxFunction* N_adapt, FlxFunction* t_plaus, FlxMtxConstFun* dtf, FlxString* pvec_str, FlxString* distid_str )
: FlxObjOutputBase(dolog,ostreamV), nameID(nameID), mcf_data(mcf_data), id_transform(id_transform), Nchain(Nchain), Nburn(Nburn), Ntune(Ntune), Npost(Npost), N_adapt(N_adapt), t_plaus(t_plaus), dtf(dtf), pvec_str(pvec_str), distid_str(distid_str)
{

}

FlxObjBayDA_new::~FlxObjBayDA_new() {
  delete nameID;
  delete mcf_data;
  delete id_transform;
  delete Nchain;
  delete Nburn;
  delete Ntune;
  delete Npost;
  delete N_adapt;
  delete t_plaus;
  delete dtf;
  delete pvec_str;
  delete distid_str;
}
void FlxObjBayDA_new::task()
{
  const std::string name = nameID->eval_word(true);
  const std::string data_name = mcf_data->eval();
  FlxSMtx* mtx_data = data->ConstMtxBox.get(data_name,true);

  const tuint id_transform_ = id_transform->cast2tuintW0(false);
  const tuint Nchain_ = Nchain->cast2tuint(false);
  const tuint Nburn_ = Nburn->cast2tuint(false);
  const tuint Ntune_ = Ntune->cast2tuint(false);
  const tuint Npost_ = Npost->cast2tuint(false);
  const tuint N_adapt_ = N_adapt->cast2tuint(false);
  const tdouble tplaus = t_plaus->cast2positive(false);
  // transform types
    FlxSMtx* dtf_mtx = data->ConstMtxBox.get(dtf->eval(),true);
    std::valarray<int> dtfv(dtf_mtx->get_Ncoeff());
    const tuint Nr = dtf_mtx->get_nrows();
    const tuint Nc = dtf_mtx->get_nrows();
    tuint c = 0;
    for (tuint i=0;i<Nr;++i) {
      for (tuint j=0;j<Nc;++j) {
        dtfv[c++] = dtf_mtx->operator()(i,j);
      }
    }
  // mtxConst to store outcomes of parameter vectors
    const std::string pvec_name = pvec_str->eval_word(true,true);
    const std::string distid_name = distid_str->eval_word(true,true);
  flxBayDA* objDA = new flxBayDA(name,id_transform_,*mtx_data,data->RndCreator,Nchain_,Nburn_,Ntune_,Npost_,N_adapt_,tplaus,dtfv,pvec_name,distid_name);
  try {
    objDA->gen_samples();
    BayUpBox->insert_DA(name,objDA,false);
  } catch (FlxException& e) {
    delete objDA;
    throw;
  }
}

FlxObjReadBayDA_new::FlxObjReadBayDA_new() {
  #if FLX_DEBUG
    if (BayUpBox==NULL) throw FlxException_Crude("FlxObjReadBayDA_new::FlxObjReadBayDA_new");
  #endif
  RBRV_entry_read_bayDA::BayUpBox = BayUpBox;

  AllDefParaBox->insert(new FlxOptionalParaFun(ONE*20,"bayda::nchain"));
    ParaBox.insert("nchain", "bayda::nchain" );
  AllDefParaBox->insert(new FlxOptionalParaFun(ONE*100000,"bayda::nburn"));
    ParaBox.insert("nburn", "bayda::nburn" );
  AllDefParaBox->insert(new FlxOptionalParaFun(ONE*10000,"bayda::ntune"));
    ParaBox.insert("ntune", "bayda::ntune" );
  AllDefParaBox->insert(new FlxOptionalParaFun(ONE*10000,"bayda::npost"));
    ParaBox.insert("npost", "bayda::npost" );
  AllDefParaBox->insert(new FlxOptionalParaFun(ONE*10,"bayda::nadapt"));
    ParaBox.insert("nadapt", "bayda::nadapt" );
  AllDefParaBox->insert(new FlxOptionalParaFun(ZERO,"bayda::id_transform"));
    ParaBox.insert("id_transform", "bayda::id_transform" );
  AllDefParaBox->insert(new FlxOptionalParaFun(ONE*15,"bayda::plausthresh"));
    ParaBox.insert("plausthresh", "bayda::plausthresh" );
  // distribution types to fit
    std::vector<FlxFunction*> vecV;
    vecV.push_back(new FlxFunction(new FunNumber(-ONE)));
    FlxObjBase* blockV = new FlxObjMtxConstNewU(false, new FlxMtxConstFun("internal_baydatypes"), vecV, 1, 1);
    FlxMtxConstFun* mtxvalue = new FlxMtxConstFun( "internal_baydatypes", blockV );
    AllDefParaBox->insert(new FlxOptionalParaMtxFun(mtxvalue,"bayda::types"));
    ParaBox.insert("types", "bayda::types" );
  // name of variable where to store parameter vector
    AllDefParaBox->insert(new FlxOptionalParaFlxString("","bayda::pvec",false));
    ParaBox.insert("pvec", "bayda::pvec" );
  // name of variable where to store distribution id
      AllDefParaBox->insert(new FlxOptionalParaFlxString("","bayda::distid",false));
      ParaBox.insert("distid", "bayda::distid" );

}

FlxObjBase* FlxObjReadBayDA_new::read()
{
  FlxString* nameID = new FlxString(false,false);
  FlxMtxConstFun* mcf_data = NULL;
  try {
    reader->getChar('(',false);
    FlxMtxConstFun* mcf_data = new FlxMtxConstFun(true);
    reader->getChar(')',false);
    read_optionalPara(false);
    return new FlxObjBayDA_new(get_doLog(),get_stream(),nameID,mcf_data,get_optPara_FlxFunction("id_transform"),
                                  get_optPara_FlxFunction("nchain"),get_optPara_FlxFunction("nburn"),get_optPara_FlxFunction("ntune"),get_optPara_FlxFunction("npost"),get_optPara_FlxFunction("nadapt"),get_optPara_FlxFunction("plausthresh"),get_optPara_FlxMtxFun("types"),get_optPara_FlxString("pvec"),get_optPara_FlxString("distid"));
  } catch (FlxException& e) {
    delete nameID;
    if (mcf_data) delete mcf_data;
    throw;
  }
}

//----------------------------------------------------------------------------------------------

FlxObjBayDA_sample::FlxObjBayDA_sample( const bool dolog, const std::string& ostreamV, FlxString* nameID )
: FlxObjOutputBase(dolog,ostreamV), nameID(nameID)
{

}

void FlxObjBayDA_sample::task()
{
  const std::string name = nameID->eval_word(true);
  flxBayDA& objDA = BayUpBox->get_DA(name);
  objDA.sample();
}

FlxObjBase* FlxObjReadBayDA_sample::read()
{
  FlxString* nameID = new FlxString(false,false);
  try {
    read_optionalPara(false);
    return new FlxObjBayDA_sample(get_doLog(),get_stream(),nameID);
  } catch (FlxException& e) {
    delete nameID;
    throw;
  }
}

