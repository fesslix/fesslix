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

#include "flxobjrbrv.h"
#include "flxfunction_ope_calc.h"
#include "flxfunction_fun_calc.h"


FlxObjRBRV_set_creator_box rbrv_set_creator;



void FlxCreateObjReaders_RBRV::createObjReaders(FlxObjectReadBox* objReadBox) {
  objReadBox->insert("rbrv_set_new", new FlxObjReadRBRV_set_new());
  objReadBox->insert("rbrv_set_addrv", new FlxObjReadRBRV_set_addRV());
  objReadBox->insert("rbrv_set_addcorr", new FlxObjReadRBRV_set_addCorr());
  objReadBox->insert("rbrv_set_create", new FlxObjReadRBRV_set_create());
  objReadBox->insert("rbrv_set", new FlxObjReadRBRV_set());
  objReadBox->insert("rbrv_noise", new FlxObjReadRBRV_noise());
  objReadBox->insert("rbrv_proc", new FlxObjReadRBRV_proc());
  objReadBox->insert("rbrv_mvn", new FlxObjReadRBRV_mvn());
  objReadBox->insert("rbrv_mvn_cond_obsv", new FlxObjReadRBRV_mvn_cond_obsv());
  objReadBox->insert("rbrv_psd", new FlxObjReadRBRV_psd());
  objReadBox->insert("rbrv_sphere", new FlxObjReadRBRV_sphere());
  objReadBox->insert("rbrv_dirichlet", new FlxObjReadRBRV_vfset(1));
  objReadBox->insert("rbrv_multinomial", new FlxObjReadRBRV_vfset(2));
  objReadBox->insert("rbrv_vfset", new FlxObjReadRBRV_vfset(0));
  objReadBox->insert("rbrv_print", new FlxObjReadRBRV_print());
  objReadBox->insert("rbrv_info", new FlxObjReadRBRV_info());
  objReadBox->insert("rbrv_vec_get", new FlxObjReadRBRV_vec_get());
  objReadBox->insert("rbrv_vec_set", new FlxObjReadRBRV_vec_set());
}

void FlxCreateObjReaders_RBRV::createFunReaders(FlxData* dataBox)
{
  dataBox->FunBox.insert("rbrv", new FunReadFunRBRV() );
  dataBox->FunBox.insert("rbrv_prob", new FunReadFunRBRV_prob() );
  dataBox->FunBox.insert("rbrv_rp", new FunReadFunRBRV_rp() );
  dataBox->FunBox.insert("pdf", new FunReadFunPDF() );
  dataBox->FunBox.insert("pdf_ln", new FunReadFunPDF_log() );
  dataBox->FunBox.insert("cdf", new FunReadFunCDF() );
  dataBox->FunBox.insert("sf", new FunReadFunSF() );
  dataBox->FunBox.insert("entropy", new FunReadFunEntropy() );
  dataBox->FunBox.insert("cdf_inv", new FunReadFunCDF_inv() );
  dataBox->FunBox.insert("mean", new FunReadFunRBRV_mean() );
  dataBox->FunBox.insert("stddev", new FunReadFunRBRV_sd() );
  dataBox->FunBox.insert("coeffofvar", new FunReadFunRBRV_coeffofvar() );
  dataBox->FunBox.insert("median", new FunReadFunRBRV_median() );
  dataBox->FunBox.insert("mode", new FunReadFunRBRV_mode() );
  dataBox->FunBox.insert("hpd", new FunReadFunHPD() );
  dataBox->FunBox.insert("rnd_y2x", new FunReadFunRBRV_y2x() );
  dataBox->FunBox.insert("rnd_x2y", new FunReadFunRBRV_x2y() );
  dataBox->FunBox.insert("rnd_sample", new FunReadFunRndSample() );
  dataBox->FunBox.insert("expectation_1d", new FunReadFunExpectation_1d() );
  dataBox->FunBox.insert("expectation_mci", new FunReadFunExpectation_mci() );
  
  
  dataBox->ConstantBox.insert("rbrv_corr_eps_x",tdouble(1e-6));
  dataBox->ConstantBox.insert("rbrv_corr_eps_y",tdouble(1e-8));
  dataBox->ConstantBox.insert("rbrv_corr_eps_ubound",tdouble(9.));
  dataBox->ConstantBox.insert("rbrv_corr_eps_intn",tdouble(150.));
}

FlxObjRBRV_set_new::FlxObjRBRV_set_new(const bool dolog, FlxString* set_name, std::vector< FlxString* > set_parents, const bool allow_x2y, const bool is_Nataf, const bool is_Nataf_evalOnce)
: FlxObjBase(dolog), set_name(set_name), set_parents(set_parents), allow_x2y(allow_x2y), is_Nataf(is_Nataf), is_Nataf_evalOnce(is_Nataf_evalOnce)
{

}

FlxObjRBRV_set_new::~FlxObjRBRV_set_new()
{
  delete set_name;
  for (tuint i=0;i<set_parents.size();++i) {
    delete set_parents[i];
  }
}

void FlxObjRBRV_set_new::task()
{
  const std::string name = set_name->eval_word(true);
  const tuint Nparents = set_parents.size();
  if (is_Nataf) {
    rbrv_set_creator.create_new(name,new FlxObjRBRV_set_creator(name,is_Nataf_evalOnce));
  } else {
    RBRV_set_baseDPtr parents = nullptr;
    RBRV_entry_read_base::generate_set_base(data->rbrv_box,name,set_parents,parents);
    rbrv_set_creator.create_new(name,new FlxObjRBRV_set_creator(name,parents,Nparents,allow_x2y));
  }
}

FlxObjReadRBRV_set_new::FlxObjReadRBRV_set_new(): FlxObjReadBase()
{
  AllDefParaBox->insert(new FlxOptionalParaBool(false,"rbrv_set::allow_x2y"));
    ParaBox.insert("allow_x2y", "rbrv_set::allow_x2y" );
  AllDefParaBox->insert(new FlxOptionalParaBool(false,"rbrv_set::is_nataf"));
    ParaBox.insert("is_nataf", "rbrv_set::is_nataf" );
  AllDefParaBox->insert(new FlxOptionalParaBool(true,"rbrv_set::is_nataf_only_once"));
    ParaBox.insert("is_nataf_only_once", "rbrv_set::is_nataf_only_once" );
}

FlxObjBase* FlxObjReadRBRV_set_new::read()
{
  FlxString* set_name = new FlxString(false,false);
  std::vector<FlxString*> set_parents;
  try {
    RBRV_entry_read_base::read_parents(set_parents,false);
    read_optionalPara(false);
    const bool is_Nataf = get_optPara_bool("is_nataf");
    if (is_Nataf) {
      if (set_parents.size()>0) {
        std::ostringstream ssV;
        ssV << "A Nataf set is not allowed to have parents.";
        throw FlxException("FlxObjReadRBRV_set_new::read", ssV.str(), reader->getCurrentPos() );
      }
    }
    return new FlxObjRBRV_set_new(get_doLog(), set_name, set_parents, get_optPara_bool("allow_x2y"), is_Nataf, get_optPara_bool("is_nataf_only_once") );
  } catch (FlxException &e) {
    FLXMSG("FlxObjReadRBRV_set_create::task",1);
    delete set_name;
    for (tuint i=0;i<set_parents.size();++i) {
      delete set_parents[i];
    }
    throw;
  }
}


void FlxObjRBRV_set_addRV::task()
{
  const std::string name = set_name->eval_word(true);
  FlxObjRBRV_set_creator* crtr = rbrv_set_creator.get_creator(name);
  crtr->add_entry(data->rbrv_box, entry);
}


FlxObjBase* FlxObjReadRBRV_set_addRV::read()
{
  FlxString* set_name = new FlxString(false,false);
  RBRV_entry_read_base* entry = NULL;
  try {
    reader->getChar('+');
    reader->getChar('=');
    entry = RBRV_entry_read_base::read_entry();
    entry->read_corr(false);
    read_optionalPara(false);
    return new FlxObjRBRV_set_addRV(get_doLog(),set_name,entry);
  } catch (FlxException& e) {
    delete set_name;
    if (entry) delete entry;
    throw;
  }
}


void FlxObjRBRV_set_addCorr::task()
{
  const std::string name = set_name->eval_word(true);
  const std::string rv1 = name + "::" + name_rv1->eval_word(true,false,true);
  const std::string rv2 = name + "::" + name_rv2->eval_word(true,false,true);
  const tdouble rho = corrVal->calc();
  FlxObjRBRV_set_creator* crtr = rbrv_set_creator.get_creator(name);
  crtr->add_corr(rv1,rv2,rho,corr_approx,rhogauss,!NOTdolog);
}

FlxObjReadRBRV_set_addCorr::FlxObjReadRBRV_set_addCorr(): FlxObjReadBase()
{
  AllDefParaBox->insert(new FlxOptionalParaBool(true,"rbrv_set::corr_approx"));
    ParaBox.insert("corr_approx", "rbrv_set::corr_approx" );
  AllDefParaBox->insert(new FlxOptionalParaBool(false,"rbrv_set::rhogauss"));
    ParaBox.insert("rhogauss", "rbrv_set::rhogauss" );
}

FlxObjBase* FlxObjReadRBRV_set_addCorr::read()
{
  FlxString* set_name = new FlxString(false,false);
  FlxString* name_rv1 = NULL;
  FlxString* name_rv2 = NULL;
  FlxFunction* corrVal = NULL;
  try {
    reader->getChar('(');
    FlxString* name_rv1 = new FlxString(false,false);
    reader->getChar(',');
    FlxString* name_rv2 = new FlxString(false,false);
    reader->getChar(')');
    reader->getChar('=');
    corrVal = new FlxFunction(funReader->read(false));
    read_optionalPara(false);
    return new FlxObjRBRV_set_addCorr(get_doLog(),set_name,name_rv1,name_rv2,corrVal,get_optPara_bool("corr_approx"),get_optPara_bool("rhogauss"));
  } catch (FlxException& e) {
    delete set_name;
    if (name_rv1) delete name_rv1;
    if (name_rv2) delete name_rv2;
    if (corrVal) delete corrVal;
    throw;
  }
}


void FlxObjRBRV_set_create::task()
{
  const std::string name = set_name->eval_word(true);
  rbrv_set_creator.register_set(name,data->rbrv_box);
}

FlxObjBase* FlxObjReadRBRV_set_create::read()
{
  FlxString* set_name = new FlxString(false,false);
  try {
    read_optionalPara(false);
    return new FlxObjRBRV_set_create(get_doLog(),set_name);
  } catch (FlxException& e) {
    delete set_name;
    throw;
  }
}



FlxObjRBRV_set::FlxObjRBRV_set(const bool dolog, FlxString* set_name, std::vector< FlxString* > set_parents, std::vector< RBRV_entry_read_base* > set_entries, const bool allow_x2y)
: FlxObjBase(dolog), set_name(set_name), set_parents(set_parents), set_entries(set_entries), allow_x2y(allow_x2y)
{
  
}

FlxObjRBRV_set::~FlxObjRBRV_set()
{
  delete set_name;
  for (tuint i=0;i<set_parents.size();++i) {
    delete set_parents[i];
  }
  for (tuint i=0;i<set_entries.size();++i) {
    delete set_entries[i];
  }
}

void FlxObjRBRV_set::task()
{
  const std::string name = set_name->eval_word(true);
  const tuint Nparents = set_parents.size();
  RBRV_set_baseDPtr parents = nullptr;
  RBRV_entry_read_base::generate_set_base(data->rbrv_box,name,set_parents,parents);
  FlxObjRBRV_set_creator crtr(data->rbrv_box,name,parents,Nparents,allow_x2y,set_entries);
  crtr.register_set(data->rbrv_box,true);
}

FlxObjReadRBRV_set::FlxObjReadRBRV_set(): FlxObjReadBase()
{
  AllDefParaBox->insert(new FlxOptionalParaBool(false,"rbrv_set::allow_x2y"));
    ParaBox.insert("allow_x2y", "rbrv_set::allow_x2y" );
}

FlxObjBase* FlxObjReadRBRV_set::read() {
  FlxString* set_name = new FlxString(false,false);
  std::vector<FlxString*> set_parents;
  std::vector<RBRV_entry_read_base*> set_entries;
  try {
    RBRV_entry_read_base::read(set_entries,set_parents,false);
    read_optionalPara(false);
    return new FlxObjRBRV_set(get_doLog(), set_name, set_parents, set_entries, get_optPara_bool("allow_x2y") );
  } catch (FlxException &e) {
    FLXMSG("FlxObjRBRV_set::task",1);
    delete set_name;
    for (tuint i=0;i<set_parents.size();++i) {
      delete set_parents[i];
    }
    for (tuint i=0;i<set_entries.size();++i) {
      delete set_entries[i];
    }
    throw;
  }
}

FlxObjRBRV_noise::FlxObjRBRV_noise(const bool dolog, FlxString* set_name, std::vector< FlxString* > set_parents, FlxFunction* Nfun, RBRV_entry_read_base* transf)
: FlxObjBase(dolog), set_name(set_name), set_parents(set_parents), Nfun(Nfun), transf(transf)
{
  
}

FlxObjRBRV_noise::~FlxObjRBRV_noise()
{
  delete set_name;
  for (tuint i=0;i<set_parents.size();++i) {
    delete set_parents[i];
  }
  delete Nfun;
  delete transf;
}

void FlxObjRBRV_noise::task()
{
  const std::string name = set_name->eval_word(true);
  RBRV_set_baseDPtr parents = nullptr;
  RBRV_entry_read_base::generate_set_base(data->rbrv_box,name,set_parents,parents);
  const tuint Nparents = set_parents.size();
  RBRV_set_noise* ts = NULL;
  RBRV_entry* transfe = NULL;
  try {
    const tuint Ndim = Nfun->cast2tuint(false);
    const std::string family = name + "::";
    tuint riID = 0;
    transfe = transf->generate_entry(family,riID);
    ts = new RBRV_set_noise(false,Ndim,name,false,transfe,Nparents,parents);
    parents = NULL;
    transfe = NULL;
    data->rbrv_box.register_set(ts);
    GlobalVar.slog(4) << "rbrv_noise: created new set '" << name << "'." << std::endl;
  } catch (FlxException& e) {
    FLXMSG("FlxObjRBRV_noise::task",1);
    if (parents) delete [] parents;
    if (ts) delete ts;
    if (transfe) delete transfe;
    throw;
  }
}

FlxObjBase* FlxObjReadRBRV_noise::read() {
  FlxString* set_name = new FlxString(false,false);
  FlxFunction* Nfun = NULL;
  std::vector<FlxString*> set_parents;
  try {
    RBRV_entry_read_base::read_parents(set_parents,false);
    reader->getChar('{');
    Nfun = new FlxFunction(funReader->read(false));
    reader->getChar(';');
    RBRV_entry_read_base* transf = RBRV_entry_read_base::read_entry(false);
    reader->getChar('}');
    read_optionalPara(false);
    return new FlxObjRBRV_noise(get_doLog(), set_name, set_parents, Nfun, transf );
  } catch (FlxException& e) {
    FLXMSG("FlxObjReadRBRV_noise::read",1);
    delete set_name;
    if (Nfun) delete Nfun;
    for (tuint i=0;i<set_parents.size();++i) {
      delete set_parents[i];
    }
    throw;
  }  
}

FlxObjRBRV_proc::FlxObjRBRV_proc(const bool dolog, FlxString* set_name, std::vector< FlxString* > set_parents, FlxFunction* Nfun, RBRV_entry_read_base* transf, FlxFunction* rho, FlxFunction* dxf, const tuint M, const tuint evtype, const bool only_once, const bool rhoGauss)
: FlxObjRBRV_noise(dolog,set_name,set_parents,Nfun,transf), rho(rho), dxf(dxf), M(M), evtype(evtype), only_once(only_once), rhoGauss(rhoGauss)
{
  
}

FlxObjRBRV_proc::~FlxObjRBRV_proc()
{
  delete dxf;
  delete rho;
}

void FlxObjRBRV_proc::task()
{
  const std::string name = set_name->eval_word(true);
  RBRV_set_baseDPtr parents = nullptr;
  RBRV_entry_read_base::generate_set_base(data->rbrv_box,name,set_parents,parents);
  const tuint Nparents = set_parents.size();
  RBRV_set_proc* ts = NULL;
  RBRV_entry* transfe = NULL;
  try {
    const tuint Ndim = Nfun->cast2tuint(false);
    const std::string family = name + "::";
    tuint riID = 0;
    transfe = transf->generate_entry(family,riID);
    const tdouble dx = dxf->cast2positive(false);
    ts = new RBRV_set_proc(false,Ndim,M,name,false,transfe,new FlxFunction(*rho),dx,Nparents,parents,evtype,only_once,rhoGauss);
    parents = NULL;
    data->rbrv_box.register_set(ts);
    GlobalVar.slog(4) << "rbrv_proc: created new set '" << name << "'." << std::endl;
  } catch (FlxException& e) {
    FLXMSG("FlxObjRBRV_proc::task",1);
    if (parents) delete [] parents;
    if (ts) {
      delete ts;
    } else {
      if (transfe) delete transfe;
    }
    throw;
  }
}

FlxObjReadRBRV_proc::FlxObjReadRBRV_proc()
{
  AllDefParaBox->insert(new FlxOptionalParaFun(1.,"rbrv_proc::dx"));
    ParaBox.insert("dx", "rbrv_proc::dx" ); 
  AllDefParaBox->insert(new FlxOptionalParaFun(0.,"rbrv_proc::m"));
    ParaBox.insert("m", "rbrv_proc::m" ); 
  AllDefParaBox->insert(new FlxOptionalParaFun(2.,"rbrv_proc::evtype"));
    ParaBox.insert("evtype", "rbrv_proc::evtype" ); 
  AllDefParaBox->insert(new FlxOptionalParaBool(true,"rbrv_proc::only_once"));
    ParaBox.insert("only_once", "rbrv_proc::only_once" );
  AllDefParaBox->insert(new FlxOptionalParaBool(false,"rbrv_proc::rhogauss"));
    ParaBox.insert("rhogauss", "rbrv_proc::rhogauss" );
}

FlxObjBase* FlxObjReadRBRV_proc::read() {
  FlxString* set_name = new FlxString(false,false);
  std::vector<FlxString*> set_parents;
  try {
    RBRV_entry_read_base::read_parents(set_parents,false);
    reader->getChar('{');
    FlxFunction* Nfun = new FlxFunction(funReader->read(false));
    reader->getChar(';');
    RBRV_entry_read_base* transf = RBRV_entry_read_base::read_entry(false);
    reader->getChar(';');
    FlxFunction* rho = new FlxFunction(funReader->read(false));
    reader->getChar('}');
    read_optionalPara(false);
    // get M
      const tuint M = get_optPara_tuint_from_FlxFunction("m",true,false);
    // get evtype
      const tuint evtype = get_optPara_tuint_from_FlxFunction("evtype",false,false);
    return new FlxObjRBRV_proc(get_doLog(), set_name, set_parents, Nfun, transf, rho, get_optPara_FlxFunction("dx"),M,evtype,get_optPara_bool("only_once"),get_optPara_bool("rhogauss"));
  } catch (FlxException& e) {
    FLXMSG("FlxObjReadRBRV_proc::read",1);
    delete set_name;
    for (tuint i=0;i<set_parents.size();++i) {
      delete set_parents[i];
    }
    throw;
  }  
}

FlxObjRBRV_mvn_conv::FlxObjRBRV_mvn_conv(const bool dolog, FlxString* set_name, FlxString* set1_name, FlxString* set2_name, const tuint M, const tuint evtype)
: FlxObjBase(dolog), set_name(set_name), set1_name(set1_name), set2_name(set2_name), M(M), evtype(evtype)
{
  
}

FlxObjRBRV_mvn_conv::~FlxObjRBRV_mvn_conv()
{
  delete set_name;
  delete set1_name;
  delete set2_name;
}

void FlxObjRBRV_mvn_conv::task()
{
  const std::string name = set_name->eval_word(true);
  // obtain the two relevant sets
    const std::string set1n = set1_name->eval_word(true);
    const std::string set2n = set2_name->eval_word(true);
    RBRV_set_base* set1 = data->rbrv_box.get_set(set1n,true);
    RBRV_set_base* set2 = data->rbrv_box.get_set(set2n,true);
    if (set1->get_NOX()!=set1->get_NOX_only_this()
      || set2->get_NOX()!=set2->get_NOX_only_this()
      || set1->get_NRV()!=set1->get_NRV_only_this()
      || set2->get_NRV()!=set2->get_NRV_only_this()
      || set1->get_NOX()!=set1->get_NRV()
      || set2->get_NOX()!=set2->get_NRV()
      || set1->get_NOX()!=set2->get_NRV()
    ) {
      std::ostringstream ssV;
      ssV << "Invalid sets '" << set1n << "' and '" << set2n << "'.";
      throw FlxException("FlxObjRBRV_mvn_conv::task_1", ssV.str() );
    }
    const tuint Ndim = set1->get_NRV();
  // get the covariance matrix
    RBRV_set_base* sb = data->rbrv_box.get_set(name,false);
    RBRV_set_MVN* ts = NULL;
    FlxMtxSym* CovM = NULL;
    if (sb) {
      ts = dynamic_cast<RBRV_set_MVN*>(sb);
    } 
    if (ts) {
      CovM = ts->get_CovM();
      CovM->set_zeroMtx();
    } else {
      CovM = new FlxMtxSym(Ndim);
    }
    try {
    set1->add_covMTX(*CovM);
    set2->add_covMTX(*CovM);
    } catch (FlxException &e) {
      FLXMSG("FlxObjRBRV_mvn_conv::task_2",1);
      if (ts==NULL) delete CovM;
      throw;
    }
  if (ts) {        // update the set
    ts->update_EVP();
  } else {        // generate the set
    try {
      flxVec* mu = new flxVec(Ndim);        // zero mean vector
      ts = new RBRV_set_MVN(false,Ndim,M,name,false,mu,CovM,evtype);
      data->rbrv_box.register_set(ts);
      GlobalVar.slog(4) << "rbrv_mvn: created new set '" << name << "' (as a convolution integral)." << std::endl;
    } catch (FlxException& e) {
      FLXMSG("FlxObjRBRV_mvn_conv::task_3",1);
      if (ts) delete ts;
      throw;
    }
  }
}

FlxObjRBRV_mvn_post::FlxObjRBRV_mvn_post(const bool dolog, FlxString* set_name, FlxString* set1_name, FlxString* set2_name, std::string ov_name, const bool only_obsv, const tuint M, const tuint evtype)
: FlxObjBase(dolog), set_name(set_name), set1_name(set1_name), set2_name(set2_name), ov_name(ov_name), only_obsv(only_obsv), M(M), evtype(evtype),
  helpV(NULL), helpM(NULL)
{
  
}

FlxObjRBRV_mvn_post::~FlxObjRBRV_mvn_post()
{
  delete set_name;
  delete set1_name;
  delete set2_name;
  if (helpM) delete helpM;
  if (helpV) delete helpV;
}

void FlxObjRBRV_mvn_post::task()
{
  const std::string name = set_name->eval_word(true);
  // obtain the two relevant sets
    const std::string set1n = set1_name->eval_word(true);
    const std::string set2n = set2_name->eval_word(true);
    RBRV_set_base* set1 = data->rbrv_box.get_set(set1n,true);
    RBRV_set_base* set2 = data->rbrv_box.get_set(set2n,true);
    if (set1->get_NOX()!=set1->get_NOX_only_this()
      || set2->get_NOX()!=set2->get_NOX_only_this()
      || set1->get_NRV()!=set1->get_NRV_only_this()
      || set2->get_NRV()!=set2->get_NRV_only_this()
      || set1->get_NOX()!=set1->get_NRV()
      || set2->get_NOX()!=set2->get_NRV()
      || set1->get_NOX()!=set2->get_NRV()
    ) {
      std::ostringstream ssV;
      ssV << "Invalid sets '" << set1n << "' and '" << set2n << "'.";
      throw FlxException("FlxObjRBRV_mvn_post::task_1", ssV.str() );
    }
    tuint Ndim = set1->get_NRV();
    flxVec ov_vec(data->ConstMtxBox.get_Vec(ov_name,Ndim,true),Ndim);
  // get the covariance matrix
    RBRV_set_base* sb = data->rbrv_box.get_set(name,false);
    RBRV_set_MVN* ts = NULL;
    flxVec* mu = NULL;
    FlxMtxSym* CovM = NULL;
    if (sb) {
      ts = dynamic_cast<RBRV_set_MVN*>(sb);
    } 
    bool only_obsv_tmp = only_obsv;
    if (helpM) {
      if (helpV==NULL) throw FlxException_Crude("FlxObjRBRV_mvn_post::task_2" );
      if (helpM->ncols()!=Ndim || helpV->get_N()!=Ndim) {
        throw FlxException_Crude("FlxObjRBRV_mvn_post::task_3" );
      }
    } else {
      if (ts==NULL && only_obsv_tmp==true) {
        only_obsv_tmp = false;
      }
      helpV = new flxVec(Ndim);
      helpM = new FlxMtxSym(Ndim);
    }
    if (ts) {
      CovM = ts->get_CovM();
      if (only_obsv_tmp==false) {
        CovM->set_zeroMtx();
        helpM->set_zeroMtx();
      }
      mu = ts->get_mu();
      mu->set_zero();
    } else {
      CovM = new FlxMtxSym(Ndim);
      mu = new flxVec(Ndim);
    }
    try {
      if (only_obsv_tmp==false || ts==NULL) {
        set1->add_covMTX(*CovM);        // modelling error
        set2->add_covMTX(*helpM);        // measurement error
        CovM->Invert();
        helpM->Invert();
      }
      helpM->MultMv(ov_vec,*helpV);
      if (only_obsv_tmp==false || ts==NULL) {
        CovM->operator+=(*helpM);
        CovM->Invert();
      }
      CovM->MultMv(*helpV,*mu);
    } catch (FlxException &e) {
      FLXMSG("FlxObjRBRV_mvn_conv::task_5",1);
      if (ts==NULL) {
        delete CovM;
        delete mu;
      }
      throw;
    }
  if (ts) {        // update the set
    if (only_obsv_tmp==false) {
      ts->update_EVP();
    }
  } else {        // generate the set
    try {
      ts = new RBRV_set_MVN(false,Ndim,M,name,false,mu,CovM,evtype);
      data->rbrv_box.register_set(ts);
      GlobalVar.slog(4) << "rbrv_mvn: created new set '" << name << "' (as a posterior distribution)." << std::endl;
    } catch (FlxException& e) {
      FLXMSG("FlxObjRBRV_mvn_conv::task_6",1);
      if (ts) delete ts;
      throw;
    }
  }
}

FlxObjReadRBRV_mvn::FlxObjReadRBRV_mvn()
{
  AllDefParaBox->insert(new FlxOptionalParaFun(0.,"rbrv_mvn::m"));
    ParaBox.insert("m", "rbrv_mvn::m" ); 
  AllDefParaBox->insert(new FlxOptionalParaFun(2.,"rbrv_mvn::evtype"));
    ParaBox.insert("evtype", "rbrv_mvn::evtype" ); 
  AllDefParaBox->insert(new FlxOptionalParaFlxString("conv","rbrv_mvn::meth",true));
    ParaBox.insert("meth", "rbrv_mvn::meth" );
  AllDefParaBox->insert(new FlxOptionalParaFlxString("","rbrv_mvn::obsv",false));
    ParaBox.insert("obsv", "rbrv_mvn::obsv" );
  AllDefParaBox->insert(new FlxOptionalParaBool(false,"rbrv_mvn::only_obsv"));
    ParaBox.insert("only_obsv", "rbrv_mvn::only_obsv" ); 
}

FlxObjBase* FlxObjReadRBRV_mvn::read() {
  FlxString* set_name = new FlxString(false,false);
  FlxString* set1_name = NULL;
  FlxString* set2_name = NULL;
  try {
    reader->getChar('{');
    reader->getWord("set1");
    reader->getChar('=');
    set1_name = new FlxString(false,false);
    reader->getChar(';');
    reader->getWord("set2");
    reader->getChar('=');
    set2_name = new FlxString(false,false);
    reader->getChar('}');
    read_optionalPara(false);
    // get M
      const tuint M = get_optPara_tuint_from_FlxFunction("m",true,false);
    // get evtype
      const tuint evtype = get_optPara_tuint_from_FlxFunction("evtype",false,false);
    // get meth
      std::string meth = get_optPara_word_from_FlxString("meth",true);
    if (meth=="conv") {
      return new FlxObjRBRV_mvn_conv(get_doLog(),set_name, set1_name, set2_name, M, evtype);
    } else if (meth=="post") {
      std::string ov_name = get_optPara_word_from_FlxString("obsv",true);
      return new FlxObjRBRV_mvn_post(get_doLog(),set_name, set1_name, set2_name, ov_name, get_optPara_bool("only_obsv"), M, evtype);
    } else {
      std::ostringstream ssV;
      ssV << "Unknown method ID '" << meth << "'.";
      throw FlxException("FlxObjReadRBRV_mvn::read_3", ssV.str(), reader->getCurrentPos() );
    }
  } catch (FlxException& e) {
    FLXMSG("FlxObjReadRBRV_mvn::read_4",1);
    delete set_name;
    if (set1_name) delete set1_name;
    if (set2_name) delete set2_name;
    throw;
  }  
}

FlxObjRBRV_mvn_cond_obsv::FlxObjRBRV_mvn_cond_obsv(const bool dolog, FlxString* set_name, FlxString* tv_name)
: FlxObjBase(dolog), set_name(set_name), tv_name(tv_name)
{
  
}

FlxObjRBRV_mvn_cond_obsv::~FlxObjRBRV_mvn_cond_obsv()
{
  delete set_name;
  delete tv_name;
}

void FlxObjRBRV_mvn_cond_obsv::task()
{
  const std::string sn = set_name->eval_word(true);
  const std::string vn = tv_name->eval_word(true);
  RBRV_set_base* sb = data->rbrv_box.get_set(sn,true);
  RBRV_set_MVN_cond* ts = dynamic_cast<RBRV_set_MVN_cond*>(sb);
  if (ts==NULL) {
    std::ostringstream ssV;
    ssV << "The set '" << sn << "' is of wrong type.";
    throw FlxException("FlxObjRBRV_mvn_cond_obsv::task", ssV.str() );
  }
  tuint N = ts->get_Nobsv();
  flxVec tv(data->ConstMtxBox.get_Vec(vn,N,true),N);
  ts->set_x_obsv(tv);
}

FlxObjBase* FlxObjReadRBRV_mvn_cond_obsv::read() {
  FlxString* set_name = new FlxString(false,false);
  FlxString* tv_name = NULL;
  try {
    reader->getChar('(');
    tv_name = new FlxString(false,false);
    reader->getChar(')');
    read_optionalPara(false);
    return new FlxObjRBRV_mvn_cond_obsv(get_doLog(),set_name,tv_name);
  } catch (FlxException& e) {
    FLXMSG("FlxObjReadRBRV_mvn_cond_obsv::read",1);
    delete set_name;
    if (tv_name) delete tv_name;
    throw;
  }
}

FlxObjRBRV_psd::FlxObjRBRV_psd(const bool dolog, FlxString* set_name, std::vector< FlxString* > set_parents, FlxFunction* Nfun, FlxFunction* psd_fun, FlxFunction* lbfun, FlxFunction* ubfun)
: FlxObjBase(dolog), set_name(set_name), set_parents(set_parents), Nfun(Nfun), psd_fun(psd_fun), lbfun(lbfun), ubfun(ubfun)
{

}

FlxObjRBRV_psd::~FlxObjRBRV_psd()
{
  delete set_name;
  for (tuint i=0;i<set_parents.size();++i) {
    delete set_parents[i];
  }
  delete Nfun;
  delete psd_fun;
  delete lbfun;
  delete ubfun;
}

void FlxObjRBRV_psd::task()
{
  const std::string name = set_name->eval_word(true);
  const std::string family = name + "::";
  const tuint N = Nfun->cast2tuint(false);
  const tdouble lb = lbfun->cast2positive_or0(false);
  const tdouble ub = ubfun->cast2positive(false);
  if (ub<=lb) {
    throw FlxException("FlxObjRBRV_psd::task","Lower bound must be smaller than upper bound.");
  }
  RBRV_set_baseDPtr parents = nullptr;
  RBRV_entry_read_base::generate_set_base(data->rbrv_box,name,set_parents,parents);
  const tuint Nparents = set_parents.size();
  RBRV_set_psd* ts = NULL;
  try {
    ts = new RBRV_set_psd(false,name,N,new FlxFunction(*psd_fun),lb,ub,Nparents,parents,data->ConstantBox.getRef("gx"));
    parents = NULL;
    data->rbrv_box.register_set(ts);
    GlobalVar.slog(4) << "rbrv_psd: created new set '" << name << "'." << std::endl;
  } catch (FlxException& e) {
    FLXMSG("FlxObjRBRV_psd::task",1);
    if (parents) delete [] parents;
    if (ts) delete ts;
    throw;
  }
}

FlxObjBase* FlxObjReadRBRV_psd::read()
{
  FlxString* set_name = new FlxString(false,false);
  std::vector<FlxString*> set_parents;
  FlxFunction* Nfun = NULL;
  FlxFunction* psd_fun = NULL;
  FlxFunction* lbfun = NULL;
  FlxFunction* ubfun = NULL;
  try {
    RBRV_entry_read_base::read_parents(set_parents,false);
    reader->getChar('{');
    // N
      reader->getWord("n",false);
      reader->getChar('=');
      Nfun = new FlxFunction(funReader->read(false));
      reader->getChar(';');
    // psd_fun
      reader->getWord("psd",false);
      reader->getChar('=');
      psd_fun = new FlxFunction(funReader->read(false));
      reader->getChar(';');
    // bounds
      reader->getChar('[');
      lbfun = new FlxFunction(funReader->read(false));
      reader->getChar(',');
      ubfun = new FlxFunction(funReader->read(false));
      reader->getChar(']');
    reader->getChar('}');
    read_optionalPara(false);
    return new FlxObjRBRV_psd(get_doLog(), set_name, set_parents, Nfun, psd_fun, lbfun, ubfun );
  } catch (FlxException& e) {
    FLXMSG("FlxObjReadRBRV_psd::read",1);
    delete set_name;
    for (tuint i=0;i<set_parents.size();++i) {
      delete set_parents[i];
    }
    if (Nfun) delete Nfun;
    if (psd_fun) delete psd_fun;
    if (lbfun) delete lbfun;
    if (ubfun) delete ubfun;
    throw;
  }
}

FlxObjRBRV_sphere::FlxObjRBRV_sphere(const bool dolog, FlxString* set_name, std::vector< FlxString* > set_parents, FlxFunction* Nfun, FlxFunction* r)
: FlxObjBase(dolog), set_name(set_name), set_parents(set_parents), Nfun(Nfun), r(r)
{
  
}

FlxObjRBRV_sphere::~FlxObjRBRV_sphere()
{
  delete set_name;
  for (tuint i=0;i<set_parents.size();++i) {
    delete set_parents[i];
  }
  delete Nfun;
  delete r;
}

void FlxObjRBRV_sphere::task()
{
  const std::string name = set_name->eval_word(true);
  RBRV_set_baseDPtr parents = nullptr;
  RBRV_entry_read_base::generate_set_base(data->rbrv_box,name,set_parents,parents);
  const tuint Nparents = set_parents.size();
  RBRV_set_sphere* ts = NULL;
  try {
    const tuint Ndim = Nfun->cast2tuint(false);
    const std::string family = name + "::";
    ts = new RBRV_set_sphere(false,Ndim,name,false,Nparents,parents,new FlxFunction(*r));
    parents = NULL;
    data->rbrv_box.register_set(ts);
    GlobalVar.slog(4) << "rbrv_noise: created new set '" << name << "'." << std::endl;
  } catch (FlxException& e) {
    FLXMSG("FlxObjRBRV_noise::task",1);
    if (parents) delete [] parents;
    if (ts) delete ts;
    throw;
  }
}

FlxObjBase* FlxObjReadRBRV_sphere::read()
{
  FlxString* set_name = new FlxString(false,false);
  FlxFunction* Nfun = NULL;
  FlxFunction* r = NULL;
  std::vector<FlxString*> set_parents;
  try {
    RBRV_entry_read_base::read_parents(set_parents,false);
    reader->getChar('{');
    Nfun = new FlxFunction(funReader->read(false));
    reader->getChar(';');
    r = new FlxFunction(funReader,false);
    reader->getChar('}');
    read_optionalPara(false);
    return new FlxObjRBRV_sphere(get_doLog(), set_name, set_parents, Nfun, r );
  } catch (FlxException& e) {
    FLXMSG("FlxObjReadRBRV_sphere::read",1);
    delete set_name;
    if (Nfun) delete Nfun;
    for (tuint i=0;i<set_parents.size();++i) {
      delete set_parents[i];
    }
    if (r) delete r;
    throw;
  }  
}



FlxObjRBRV_vfset::FlxObjRBRV_vfset(const bool dolog, FlxString* set_name, std::vector< FlxString* > set_parents, FlxFunction* Nfun, FlxMtxConstFun* vecfun, FlxFunction* Ntrials, const tuint type)
: FlxObjBase(dolog), type(type), set_name(set_name), set_parents(set_parents), Nfun(Nfun), vecfun(vecfun), Ntrials(Ntrials)
{
  
}

FlxObjRBRV_vfset::~FlxObjRBRV_vfset()
{
  delete set_name;
  for (tuint i=0;i<set_parents.size();++i) {
    delete set_parents[i];
  }
  delete Nfun;
  delete vecfun;
  if (Ntrials) delete Ntrials;
}

void FlxObjRBRV_vfset::task()
{
  const std::string name = set_name->eval_word(true);
  RBRV_set_baseDPtr parents = nullptr;
  RBRV_entry_read_base::generate_set_base(data->rbrv_box,name,set_parents,parents);
  const tuint Nparents = set_parents.size();
  RBRV_set_parents* ts = nullptr;
  FlxMtxFun_MtxConst* vfun = nullptr;
  try {
    const tuint Ndim = Nfun->cast2tuint(false);
    const std::string family = name + "::";
      vfun = new FlxMtxFun_MtxConst(Ndim,*vecfun);
    switch (type) {
      case 0:
        ts = new RBRV_vfset(false,name,false,Ndim,vfun,Nparents,parents);
        break;
      case 1:
        ts = new RBRV_dirichlet(false,name,false,Ndim,vfun,Nparents,parents);
        break;
      case 2:
      {
        const tuint Ntri = Ntrials->cast2tuint(false);
        ts = new RBRV_multinomial(false,name,false,Ndim,vfun,Nparents,parents, Ntri);
        break;
      }
      default:
        throw FlxException_Crude("FlxObjRBRV_vfset::task_01");
    }
    vfun =  nullptr;
    parents = nullptr;
    data->rbrv_box.register_set(ts);
    GlobalVar.slog(4) << "rbrv_noise: created new set '" << name << "'." << std::endl;
  } catch (FlxException& e) {
    FLXMSG("FlxObjRBRV_vfset::task_02",1);
    if (vfun) delete vfun;
    if (parents) delete [] parents;
    if (ts) delete ts;
    throw;
  }
}

FlxObjBase* FlxObjReadRBRV_vfset::read()
{
  FlxString* set_name = new FlxString(false,false);
  FlxFunction* Nfun = NULL;
  FlxFunction* Ntrials = NULL;
  FlxMtxConstFun* vecfun = NULL;
  std::vector<FlxString*> set_parents;
  try {
    RBRV_entry_read_base::read_parents(set_parents,false);
    reader->getChar('{');
    Nfun = new FlxFunction(funReader->read(false));
    reader->getChar(';');
    vecfun = new FlxMtxConstFun(true);
    if (type==2) {
      reader->getChar(';');
      Ntrials = new FlxFunction(funReader->read(false));
    }
    reader->getChar('}');
    read_optionalPara(false);
    return new FlxObjRBRV_vfset(get_doLog(), set_name, set_parents, Nfun, vecfun, Ntrials, type );
  } catch (FlxException& e) {
    FLXMSG("FlxObjReadRBRV_rvset::read",1);
    delete set_name;
    if (Nfun) delete Nfun;
    for (tuint i=0;i<set_parents.size();++i) {
      delete set_parents[i];
    }
    if (vecfun) delete vecfun;
    if (Ntrials) delete Ntrials;
    throw;
  }  
}

void FlxObjRBRV_print::task()
{
  if (rbstr) {
    // create the constructor
      const std::string setstr = rbstr->eval(true);
      const std::vector<std::string> set_str_vec = parse_strseq_as_vec(setstr);
      RBRV_constructor* constr = new RBRV_constructor(set_str_vec,data->rbrv_box);
    try {
      sout() << "RBRV-sets: " << setstr << std::endl;
      constr->print_info(sout());
    } catch (FlxException& e) {
      FLXMSG("FlxObjRBRV_print::task_1",1);
      GlobalVar.alert.alert("FlxObjRBRV_print::task_2","There was a problem printing information about the sets '" + setstr + "':\n" + e.what());
    }
    delete constr;
  } else {
    sout() << "List of all RBRV-sets:" << std::endl;
    data->rbrv_box.print_sets(sout(),"  ");
  }
}

FlxObjBase* FlxObjReadRBRV_print::read() {
  FlxString* rbstr = NULL;
  reader->getChar('(');
  if (reader->whatIsNextChar()!=')') {
    rbstr = new FlxString(false,false);
  }
  try {
    reader->getChar(')');
    read_optionalPara(false);
    return new FlxObjRBRV_print(get_doLog(),get_stream(),get_verbose(),rbstr);
  } catch (FlxException& e) {
    FLXMSG("FlxObjReadRBRV_print::read",1);
    if (rbstr) delete rbstr;
    throw;
  }
}

void FlxObjRBRV_info::task()
{
  throw FlxException_NotImplemented("FlxObjRBRV_info::task");
}

FlxObjBase* FlxObjReadRBRV_info::read()
{
  RBRV_entry_RV_base* rep = NULL;
  try {
    reader->getChar('(');
    rep = RBRV_entry_read_base::read_gen_entry(false);
    reader->getChar(')');
    read_optionalPara(false);
    return new FlxObjRBRV_info(get_doLog(),get_stream(),rep);
  } catch (FlxException& e) {
    if (rep) delete rep;
    throw;
  }
}

FlxObjRBRV_vec_get::~FlxObjRBRV_vec_get()
{
  delete MtxConstStr; 
  delete rbstr;
  if (constr) delete constr;
}

void FlxObjRBRV_vec_get::task()
{
  // make sure that a constructor exists
    if (NOX==0) {
      #if FLX_DEBUG
        if (rbrvSet!=NULL || constr!=NULL) throw FlxException_Crude("FlxObjRBRV_vec_get::task_1");
      #endif
      if (only_this) {
        const std::string rbrvName = rbstr->eval_word(true);
        rbrvSet = data->rbrv_box.get_set(rbrvName,true);
        NOX = rbrvSet->get_NOX_only_this();
        NRV = rbrvSet->get_NRV_only_this();
        if ( (gType==y&&NRV==0) || NOX==0 ) {
          std::ostringstream ssV;
          ssV << "The set '" << rbrvName << "' does not contain any random variables.";
          throw FlxException("FlxObjRBRV_vec_get::task_2", ssV.str() );
        }
      } else {
        const std::string setstr = rbstr->eval(true);
        const std::vector<std::string> set_str_vec = parse_strseq_as_vec(setstr);
        constr = new RBRV_constructor(set_str_vec,data->rbrv_box);
        NOX = constr->get_NOX();
        NRV = constr->get_NRV();
        if ( (gType==y&&NRV==0) || NOX==0) {
          std::ostringstream ssV;
          ssV << "No random variables are contained in: " << setstr;
          throw FlxException("FlxObjRBRV_vec_get::task_3", ssV.str() );
        }
      }
      vecName = MtxConstStr->eval();
    }
  
  // get the vector
    tdouble* vp = data->ConstMtxBox.get_Vec(vecName,(gType==y)?NRV:NOX);
    if (only_this) {
      switch(gType) {
        case x:
          rbrvSet->get_x_only_this(vp);
          break;
        case y:
          rbrvSet->get_y_only_this(vp);
          break;
        case mean:
          rbrvSet->get_mean_only_this(vp);
          break;
        case sd:
          rbrvSet->get_sd_only_this(vp);
          break;
      }
    } else {
      switch(gType) {
        case x:
          constr->get_x_Vec(vp);
          break;
        case y:
          constr->get_y_Vec(vp);
          break;
        case mean:
          constr->get_mean_Vec(vp);
          break;
        case sd:
          constr->get_sd_Vec(vp);
          break;
      }
    }
}

FlxObjReadRBRV_vec_get::FlxObjReadRBRV_vec_get(): FlxObjReadBase()
{
  AllDefParaBox->insert(new FlxOptionalParaBool(false,"rbrv_vec::only_this"));
    ParaBox.insert("only_this", "rbrv_vec::only_this" );
}

FlxObjBase* FlxObjReadRBRV_vec_get::read()
{
  // get the return type
    FlxObjRBRV_vec_get::RBRVvecGetType gType;
    const std::string gIDstr = reader->getWord(true,false);
    if (gIDstr=="x") gType=FlxObjRBRV_vec_get::x;
    else if (gIDstr=="y") gType=FlxObjRBRV_vec_get::y;
    else if (gIDstr=="mean") gType=FlxObjRBRV_vec_get::mean;
    else if (gIDstr=="sd") gType=FlxObjRBRV_vec_get::sd;
    else {
      std::ostringstream ssV;
      ssV << "Unknown type-ID: " << gIDstr;
      throw FlxException("FlxObjReadRBRV_vec_get::read_1", ssV.str(), reader->getCurrentPos() );
    }
    reader->getChar(':');
  // read name of vector to store
    FlxMtxConstFun* MtxConstStr = new FlxMtxConstFun(false);
  // read rbrv-sets and finalize
    FlxString* rbstr = NULL;
    try {
      reader->getChar('=');
      rbstr = new FlxString(false,false);
      read_optionalPara(false);
      return new FlxObjRBRV_vec_get(get_doLog(),MtxConstStr,rbstr,get_optPara_bool("only_this"),gType);
    } catch (FlxException &e) {
      FLXMSG("FlxObjReadRBRV_vec_get::read_2",1);
      delete MtxConstStr;
      if (rbstr) delete rbstr;
      throw;
    }
}


FlxObjRBRV_vec_set::~FlxObjRBRV_vec_set()
{
  delete MtxConstStr; 
  delete rbstr;
  if (constr) delete constr;
}

void FlxObjRBRV_vec_set::task()
{
  // make sure that a constructor exists
    if (NOX==0) {
      #if FLX_DEBUG
        if (rbrvSet!=NULL || constr!=NULL) throw FlxException_Crude("FlxObjRBRV_vec_set::task_1");
      #endif
      if (only_this) {
        const std::string rbrvName = rbstr->eval_word(true);
        rbrvSet = data->rbrv_box.get_set(rbrvName,true);
        NOX = rbrvSet->get_NOX_only_this();
        NRV = rbrvSet->get_NRV_only_this();
        if ( (sType==y&&NRV==0) || NOX==0 ) {
          std::ostringstream ssV;
          ssV << "The set '" << rbrvName << "' does not contain any random variables.";
          throw FlxException("FlxObjRBRV_vec_set::task_2", ssV.str() );
        }
      } else {
        const std::string setstr = rbstr->eval(true);
        const std::vector<std::string> set_str_vec = parse_strseq_as_vec(setstr);
        constr = new RBRV_constructor(set_str_vec,data->rbrv_box);
        NOX = constr->get_NOX();
        NRV = constr->get_NRV();
        if ( (sType==y&&NRV==0) || NOX==0) {
          std::ostringstream ssV;
          ssV << "No random variables are contained in: " << setstr;
          throw FlxException("FlxObjRBRV_vec_set::task_3", ssV.str() );
        }
      }
    }
    vecName = MtxConstStr->eval();
  
  // get the vector
    tdouble* vp = data->ConstMtxBox.get_Vec((sType==y)?NRV:NOX,vecName,true);
    if (only_this) {
      switch(sType) {
        case x:
          rbrvSet->set_x_only_this(vp);
          rbrvSet->transform_x2y();
          break;
        case y:
          rbrvSet->set_y_only_this(vp);
          rbrvSet->transform_y2x();
          break;
      }
    } else {
      switch(sType) {
        case x:
          {
            flxVec vpv(vp,NOX);
            constr->set_smp_x_transform(vpv);
          }
          break;
        case y:
          {
            flxVec vpv(vp,NRV);
            constr->set_smp(vpv);
          }
          break;
      }
    }
}

FlxObjReadRBRV_vec_set::FlxObjReadRBRV_vec_set(): FlxObjReadBase()
{
  AllDefParaBox->insert(new FlxOptionalParaBool(false,"rbrv_vec::only_this"));
    ParaBox.insert("only_this", "rbrv_vec::only_this" );
}

FlxObjBase* FlxObjReadRBRV_vec_set::read()
{
  // get the return type
    FlxObjRBRV_vec_set::RBRVvecSetType sType;
    const std::string sIDstr = reader->getWord(true,false);
    if (sIDstr=="x") sType=FlxObjRBRV_vec_set::x;
    else if (sIDstr=="y") sType=FlxObjRBRV_vec_set::y;
    else {
      std::ostringstream ssV;
      ssV << "Unknown type-ID: " << sIDstr;
      throw FlxException("FlxObjReadRBRV_vec_set::read", ssV.str(), reader->getCurrentPos() );
    }
    reader->getChar(':');
  // read rbrv-sets
    FlxString* rbstr = new FlxString(false,false);
  // read name of vector to store and finalize
     FlxMtxConstFun* MtxConstStr = NULL;
    try {
      reader->getChar('=');
      MtxConstStr = new FlxMtxConstFun(true);
      read_optionalPara(false);
      return new FlxObjRBRV_vec_set(get_doLog(),MtxConstStr,rbstr,get_optPara_bool("only_this"),sType);
    } catch (FlxException &e) {
      FLXMSG("FlxObjReadRBRV_vec_set::read",1);
      delete MtxConstStr;
      if (rbstr) delete rbstr;
      throw;
    }
}

const tdouble FunRBRV::calc()
{
  if (thenumber==NULL) {
    thenumber = data->rbrv_box.get_entry(rbrv_name,true);
  }
  if (is_log) {
    return thenumber->get_value_log();
  } else {
    return thenumber->get_value();
  }
}

const std::string FunRBRV::write()
{
  return "rbrv(" + rbrv_name + ")";
}

const bool FunRBRV::dependOn_Const(const tdouble*const theconst)
{
  throw FlxException_NotImplemented("FunRBRV::dependOn_Const");
}

FunRBRV_fast::~FunRBRV_fast()
{
  delete rbrv_name;
}

const tdouble FunRBRV_fast::calc()
{
  RBRV_entry* thenumber = data->rbrv_box.get_entry(rbrv_name->eval(true),true);
  if (is_log) {
    return thenumber->get_value_log();
  } else {
    return thenumber->get_value();
  }
}

const std::string FunRBRV_fast::write()
{
  return "rbrv(" + rbrv_name->write() + ")";
}

const bool FunRBRV_fast::dependOn_Const(const tdouble*const theconst)
{
  throw FlxException_NotImplemented("FunRBRV_fast::dependOn_Const");
}

FunBase* FunReadFunRBRV::read(bool errSerious)
{
  FlxString* strV = new FlxString(false,errSerious);
  bool is_log = false;
  try {
    std::string rbrv_name;
    if (strV->is_static()) {
      rbrv_name = strV->eval(true);
      delete strV; strV=NULL;
    }
    if (reader->whatIsNextChar()==',') {
      reader->getChar(',');
      const char c = reader->getChar();
      if (c=='L'||c=='l') {
        is_log = true;
      } else if (c=='N'||c=='n') {
      } else {
        std::ostringstream ssV;
        ssV << "Unknown identifier '" << c << "'. Expected 'n' or 'l'.";
        throw FlxException("FunReadFunRBRV::read_1", ssV.str() );
      }
    }
    if (strV) {
      return new FunRBRV_fast( strV, is_log );
    } else {
      return new FunRBRV( rbrv_name, is_log );
    }
  } catch (FlxException &e) {
    FLXMSG("FunReadFunRBRV::read_2",1);
    if (strV) delete strV;
    throw;
  }
}

FunRBRV_prob::~FunRBRV_prob()
{
  delete MtxConstStr; 
  delete rbstr;
}

const tdouble FunRBRV_prob::calc()
{
  if (N==0) {
    const std::string rbrvName = rbstr->eval_word(true);
    rbrvSet = data->rbrv_box.get_set(rbrvName,true);
    N = rbrvSet->get_NOX_only_this();
    if (N==0) {
      std::ostringstream ssV;
      ssV << "The set '" << rbrvName << "' does not contain any random variables.";
      throw FlxException("FunRBRV_prob::task_1", ssV.str() );
    }
    vecName = MtxConstStr->eval_word(true);
  }

  // get the vector
    tuint Nv = 0;
    tdouble* vp = data->ConstMtxBox.get_Vec(vecName,Nv);
    if (Nv!=N) {
      std::ostringstream ssV;
      ssV << "The dimension of the vector (" << Nv << ") does not match the number of random variables in the set (" << N << ").";
      throw FlxException("FunRBRV_prob::task_2", ssV.str() );
    }
  // set vp as x
    rbrvSet->set_is_valid(false);
    rbrvSet->set_x_only_this(vp);
  // get PDF of x
    return rbrvSet->get_pdf_x_eval_log();
}

const bool FunRBRV_prob::search_circref(FlxFunction* fcr)
{
  return MtxConstStr->search_circref(fcr) || rbstr->search_circref(fcr);
}

const bool FunRBRV_prob::dependOn_Const(const tdouble*const theconst)
{
  throw FlxException_NotImplemented("FunRBRV_prob::dependOn_Const");
}

const std::string FunRBRV_prob::write()
{
  return "rbrv_proc(" + rbstr->write() + "," + MtxConstStr->write() + ")";
}

FunBase* FunReadFunRBRV_prob::read(bool errSerious)
{
  FlxString* setName = new FlxString(false,errSerious);
  FlxString* vecName = NULL;
  try {
    reader->getChar(',');
    vecName = new FlxString(false,errSerious);
    return new FunRBRV_prob(vecName,setName);
  } catch (FlxException &e) {
    FLXMSG("FunReadFunRBRV_prob::read",1);
    delete setName;
    if (vecName) delete vecName;
    throw;
  }
}

// ------------------------------------------------------------------------------------------------


const std::string FunRBRV_rp::write()
{
  std::string str1 = write_v() +"(";
  str1 += get_name();
  str1 += ',';
  str1 += child_1->write();
  str1 += ")";
  return str1;
}

const tdouble FunRBRV_rp_psd::calc()
{
  return rp->eval_realization(child_1->calc());
}

FunBase* FunReadFunRBRV_rp::read(bool errSerious)
{
  const std::string setname = reader->getWord(true,errSerious);
  RBRV_set_base* sb = data->rbrv_box.get_set(setname,true);
  reader->getChar(',',errSerious);
  RBRV_set_psd* sb_psd = dynamic_cast<RBRV_set_psd*>(sb);
  if (sb_psd) {
    return new FunRBRV_rp_psd( sb_psd, read_parameters(1,errSerious) );
  }
  std::ostringstream ssV;
  ssV << "The specified rbrv-set '" << setname << "' is not a random process (with specified power spectral density function).";
  throw FlxException_NeglectInInteractive("FunReadFunRBRV_rp_1", ssV.str(), reader->getCurrentPos() );
}

// ------------------------------------------------------------------------------------------------

const tdouble FunPDF::calc()
{
  rep->eval_para();
  return rep->calc_pdf_x(fun->calc(),true);
}

const bool FunPDF::dependOn_Const(const tdouble*const thenumber)
{
  throw FlxException_NotImplemented("FunPDF::dependOn_Const");
}

const bool FunPDF::search_circref(FlxFunction* fcr)
{
  if (fun) {
    return fun->search_circref(fcr) || rep->search_circref(fcr);
  } else {
    return rep->search_circref(fcr);
  }
}

const std::string FunPDF::write()
{
  return write_v() + "(...)";
}

FunBase* FunReadFunPDF::read( bool errSerious )
{
  FunBase* fun = NULL;
  try {
    if (((methID<3||methID>=7)&&methID!=9&&methID<=11)||methID==14) {
      fun = FunctionList->read(errSerious);
      reader->getChar(',');
    }
    RBRV_entry_RV_base* rep = RBRV_entry_read_base::read_gen_entry(errSerious);
    switch(methID) {
      case 0:
        return new FunPDF( fun, rep, true ); 
      case 1:
        return new FunCDF( fun, rep, true ); 
      case 2:
        return new FunCDF_inv( fun, rep, true ); 
      case 3:
        return new FunEntropy( fun, rep, true );
      case 4:
        return new FunRBRV_mean( fun, rep, true );
      case 5:
        return new FunRBRV_sd( fun, rep, true );
      case 6:
        return new FunRBRV_coeffofvar( fun, rep, true );
      case 7:
        return new FunRBRV_y2x( fun, rep, true );
      case 8:
        return new FunPDF_log( fun, rep, true );
      case 9:
        return new FunRndSample( fun, rep, true );
      case 10:
        return new FunRBRV_x2y( fun, rep, true );
      case 11:
        return new FunHPD( fun, rep, true );
      case 12:
        return new FunRBRV_median( fun, rep, true );
      case 13:
        return new FunRBRV_mode( fun, rep, true );
      case 14:
        return new FunSF( fun, rep, true );
    }
    throw FlxException_Crude("FunReadFunPDF::read");
  } catch (FlxException& e) {
    if (fun) delete fun;
    throw;
  }
}

// ------------------------------------------------------------------------------------------------

const tdouble FunPDF_log::calc()
{
  rep->eval_para();
  return rep->calc_pdf_x_log(fun->calc(),true);
}

// ------------------------------------------------------------------------------------------------

const tdouble FunCDF::calc()
{
  rep->eval_para();
  const tdouble x = fun->calc();
  const tdouble res = rep->calc_cdf_x(x,true);
  return res;
}

// ------------------------------------------------------------------------------------------------

const tdouble FunSF::calc()
{
  rep->eval_para();
  const tdouble x = fun->calc();
  const tdouble res = rep->calc_sf_x(x,true);
  return res;
}

// ------------------------------------------------------------------------------------------------

const tdouble FunCDF_inv::calc()
{
  rep->eval_para();
  const tdouble p = fun->calc();
  return rep->calc_icdf_x(p);
}

// ------------------------------------------------------------------------------------------------

const tdouble FunHPD::calc()
{
  rep->eval_para();
  return rep->get_HPD( fun->calc() );
}

// ------------------------------------------------------------------------------------------------

const tdouble FunRBRV_y2x::calc()
{
  rep->eval_para();
  const tdouble y = fun->calc();
  const tdouble res = rep->transform_y2x(y);
  #if FLX_DEBUG
    if (std::isnan(res)) {
      throw FlxException_Crude("FunRBRV_y2x::calc");
    }
  #endif
  return res;
}

// ------------------------------------------------------------------------------------------------

const tdouble FunRBRV_x2y::calc()
{
  rep->eval_para();
  const tdouble x = fun->calc();
  const tdouble res = rep->transform_x2y(x);
  #if FLX_DEBUG
    if (std::isnan(res)) {
      throw FlxException_Crude("FunRBRV_x2y::calc");
    }
  #endif
  return res;
}

// ------------------------------------------------------------------------------------------------

const tdouble FunEntropy::calc()
{
  rep->eval_para();
  return rep->calc_entropy();
}

// ------------------------------------------------------------------------------------------------

const tdouble FunRndSample::calc()
{
  rep->eval_para();
  const tdouble y = RndCreator->gen_smp();
  return rep->transform_y2x(y);
}

// ------------------------------------------------------------------------------------------------

const tdouble FunRBRV_mean::calc()
{
  rep->eval_para();
  return rep->get_mean_current_config();
}

// ------------------------------------------------------------------------------------------------

const tdouble FunRBRV_sd::calc()
{
  rep->eval_para();
  return rep->get_sd_current_config();
}

// ------------------------------------------------------------------------------------------------

const tdouble FunRBRV_coeffofvar::calc()
{
  rep->eval_para();
  return rep->get_sd_current_config()/rep->get_mean_current_config();
}

// ------------------------------------------------------------------------------------------------

const tdouble FunRBRV_median::calc()
{
  rep->eval_para();
  return rep->get_median_current_config();
}

// ------------------------------------------------------------------------------------------------

const tdouble FunRBRV_mode::calc()
{
  rep->eval_para();
  return rep->get_mode_current_config();
}

// ------------------------------------------------------------------------------------------------

FunExpectation_1d::~FunExpectation_1d()
{
  delete fun; 
  if (rbrv_name) delete rbrv_name; 
  delete ni; 
  delete ns; 
  delete rd; 
  delete lb; 
  delete ub;
}

const tdouble FunExpectation_1d::calc()
{
  if (thenumber==NULL) {
    const std::string rn = rbrv_name->eval(true);
    delete rbrv_name; rbrv_name = NULL;
    RBRV_entry* tn = data->rbrv_box.get_entry(rn,true);
    thenumber = dynamic_cast<RBRV_entry_RV_base*>(tn);
    if (tn==NULL) {
      std::ostringstream ssV;
      ssV << "'" << rn << "' cannot be sampled from directly.";
      throw FlxException("FunExpectation_1d::calc",ssV.str());
    }
  }
  calc_expectation_numerical_1D en(fun);
  const tulong NI = tulong_from(ni->calc(),"ni",false,false,ni);
  const tulong NS = tulong_from(ns->calc(),"ns",false,false,ns);
  const tdouble RD = rd->calc();
  const tdouble LB = lb->calc();
  const tdouble UB = ub->calc();
  if (LB>=UB) throw FlxException("FunExpectation_1d::calc","'lb' must not exceed 'ub'.");
  return en.calc_expectation(NI,NS,RD,*thenumber,LB,UB);  
}

const bool FunExpectation_1d::search_circref(FlxFunction* fcr)
{
   return fun->search_circref(fcr) 
          || ni->search_circref(fcr) 
          || ns->search_circref(fcr) 
          || rd->search_circref(fcr) 
          || lb->search_circref(fcr) 
          || ub->search_circref(fcr); 
}

const bool FunExpectation_1d::dependOn_Const(const tdouble*const theconst)
{
  return fun->dependOn_Const(theconst)
          || ni->dependOn_Const(theconst)
          || ns->dependOn_Const(theconst)
          || rd->dependOn_Const(theconst)
          || lb->dependOn_Const(theconst)
          || ub->dependOn_Const(theconst);
}

const bool FunExpectation_1d::optimize(FunBasePtr& optf, const Fun_OptimizeInfo& foi)
{
  child_optimize(fun,foi);
  if (is_number(fun) ) {
    optf = fun;
    fun = new FunDummy();
    return true;
  }
  child_optimize(ni,foi);
  child_optimize(ns,foi);
  child_optimize(rd,foi);
  child_optimize(lb,foi);
  child_optimize(ub,foi);
  return false;
}

FunBase* FunReadFunExpectation_1d::read(bool errSerious)
{
  FunBase* fun = FunctionList->read(errSerious);
  FlxString* strV = NULL;
  FunBase* ni=NULL;
  FunBase* ns=NULL;
  FunBase* rd=NULL;
  FunBase* lb=NULL;
  FunBase* ub=NULL;
  try {
    reader->getChar(',',false);
    strV = new FlxString(false,errSerious);
    while (reader->whatIsNextChar()==',') {
      reader->getChar(',',false);
      const std::string key = reader->getWord(true,false);
      if (key=="ni") {
        if (ni) { throw FlxException("FunReadFunExpectation_1d::read_1","Parameter 'ni' was already specified.", reader->getCurrentPos() ); }
        reader->getChar('=',false);
        ni = FunctionList->read(errSerious);
      } else if (key=="ns") {
        if (ns) { throw FlxException("FunReadFunExpectation_1d::read_2","Parameter 'ns' was already specified.", reader->getCurrentPos() ); }
        reader->getChar('=',false);
        ns = FunctionList->read(errSerious);
      } else if (key=="rd") {
        if (rd) { throw FlxException("FunReadFunExpectation_1d::read_3","Parameter 'rd' was already specified.", reader->getCurrentPos() ); }
        reader->getChar('=',false);
        rd = FunctionList->read(errSerious);
      } else if (key=="lb") {
        if (lb) { throw FlxException("FunReadFunExpectation_1d::read_4","Parameter 'lb' was already specified.", reader->getCurrentPos() ); }
        reader->getChar('=',false);
        lb = FunctionList->read(errSerious);
      } else if (key=="ub") {
        if (ub) { throw FlxException("FunReadFunExpectation_1d::read_5","Parameter 'ub' was already specified.", reader->getCurrentPos() ); }
        reader->getChar('=',false);
        ub = FunctionList->read(errSerious);
      } else {
        throw FlxException("FunReadFunExpectation_1d::read_6","Unexpected keyword '"+key+"'.",reader->getCurrentPos());
      }
    }
    if (*(data->ConstantBox.get("leak_check"))>ZERO) {
      if (ni==NULL) ni = new FunNumber(1e2);
      if (ns==NULL) ns = new FunNumber(1e4);
    } else {
      if (ni==NULL) ni = new FunNumber(1e4);
      if (ns==NULL) ns = new FunNumber(3e5);
    }
    if (rd==NULL) rd = new FunNumber(0.9);
    if (lb==NULL) lb = new FunNumber(-38.);
    if (ub==NULL) ub = new FunNumber(38.);
    return new FunExpectation_1d(fun,strV,ni,ns,rd,lb,ub);
  } catch (FlxException &e) {
    FLXMSG("FunReadFunExpectation_1d::read",1);
    delete fun;
    if (strV) delete strV;
    if (ni) delete ni;
    if (ns) delete ns;
    if (rd) delete rd;
    if (lb) delete lb;
    if (ub) delete ub;
    throw;
  }
}

FunExpectation_mci::~FunExpectation_mci()
{
  delete fun;
  if (rndBox) delete rndBox;
  if (rbrv_sets) delete rbrv_sets;
  delete ni;
  delete nsi;
  delete nsr;
  delete rd;
  delete bound;
}

const tdouble FunExpectation_mci::calc()
{
  if (rndBox==NULL) {
    const std::vector<std::string> set_str_vec = parse_strseq_as_vec(rbrv_sets->eval(true));
    rndBox = new RBRV_constructor(set_str_vec,data->rbrv_box);
    delete rbrv_sets; rbrv_sets = NULL;
  }
  calc_expectation_numerical_MCI en;
  const tulong NI = tulong_from(ni->calc(),"ni",false,false,ni);
  const tulong NSI = tulong_from(nsi->calc(),"nsi",false,false,nsi);
  const tulong NSR = tulong_from(nsr->calc(),"nsr",false,false,nsr);
  const tdouble RD = rd->calc();
  const tdouble B = bound->calc();
  return en.calc_expectation(fun,*rndBox,NI,NSI,NSR,RD,B);  
}

const bool FunExpectation_mci::search_circref(FlxFunction* fcr)
{
  return fun->search_circref(fcr) 
          || ni->search_circref(fcr) 
          || nsi->search_circref(fcr)           
          || nsr->search_circref(fcr) 
          || rd->search_circref(fcr) 
          || bound->search_circref(fcr); 
}

const bool FunExpectation_mci::dependOn_Const(const tdouble*const theconst)
{
  return fun->dependOn_Const(theconst)
          || ni->dependOn_Const(theconst)
          || nsi->dependOn_Const(theconst)
          || nsr->dependOn_Const(theconst)
          || rd->dependOn_Const(theconst)
          || bound->dependOn_Const(theconst);
}

const bool FunExpectation_mci::optimize(FunBasePtr& optf, const Fun_OptimizeInfo& foi)
{
  child_optimize(fun,foi);
  if (is_number(fun) ) {
    optf = fun;
    fun = new FunDummy();
    return true;
  }
  child_optimize(ni,foi);
  child_optimize(nsi,foi);
  child_optimize(nsr,foi);
  child_optimize(rd,foi);
  child_optimize(bound,foi);
  return false;
}

FunBase* FunReadFunExpectation_mci::read(bool errSerious)
{
  FunBase* fun = FunctionList->read(errSerious);
  FlxString* strV = NULL;
  FunBase* ni=NULL;
  FunBase* nsi=NULL;
  FunBase* nsr=NULL;
  FunBase* rd=NULL;
  FunBase* b=NULL;
  try {
    reader->getChar(',',false);
    strV = new FlxString(false,errSerious);
    while (reader->whatIsNextChar()==',') {
      reader->getChar(',',false);
      const std::string key = reader->getWord(true,false);
      if (key=="ni") {
        if (ni) { throw FlxException("FunReadFunExpectation_1d::read_1","Parameter 'ni' was already specified.", reader->getCurrentPos() ); }
        ni = FunctionList->read(errSerious);
      } else if (key=="nsi") {
        if (nsi) { throw FlxException("FunReadFunExpectation_1d::read_2","Parameter 'nsi' was already specified.", reader->getCurrentPos() ); }
        nsi = FunctionList->read(errSerious);
      } else if (key=="nsr") {
        if (nsr) { throw FlxException("FunReadFunExpectation_1d::read_3","Parameter 'nsr' was already specified.", reader->getCurrentPos() ); }
        nsr = FunctionList->read(errSerious);
      } else if (key=="rd") {
        if (rd) { throw FlxException("FunReadFunExpectation_1d::read_4","Parameter 'rd' was already specified.", reader->getCurrentPos() ); }
        rd = FunctionList->read(errSerious);
      } else if (key=="b") {
        if (b) { throw FlxException("FunReadFunExpectation_1d::read_5","Parameter 'b' was already specified.", reader->getCurrentPos() ); }
        b = FunctionList->read(errSerious);
      } else {
        throw FlxException("FunReadFunExpectation_1d::read_6","Unexpected keyword '"+key+"'.",reader->getCurrentPos());
      }
    }
    if (*(data->ConstantBox.get("leak_check"))>ZERO) {
      if (ni==NULL) ni = new FunNumber(100.);
      if (nsi==NULL) nsi = new FunNumber(1e2);
      if (nsr==NULL) nsr = new FunNumber(1e5);
    } else {
      if (ni==NULL) ni = new FunNumber(608.);
      if (nsi==NULL) nsi = new FunNumber(1e3);
      if (nsr==NULL) nsr = new FunNumber(1e7);
    }
    if (rd==NULL) rd = new FunNumber(0.2);
    if (b==NULL) b = new FunNumber(38.);
    return new FunExpectation_mci(fun,strV,ni,nsi,nsr,rd,b);
  } catch (FlxException &e) {
    FLXMSG("FunReadFunExpectation_mci::read",1);
    delete fun;
    if (strV) delete strV;
    if (ni) delete ni;
    if (nsi) delete nsi;
    if (nsr) delete nsr;
    if (rd) delete rd;
    if (b) delete b;
    throw;
  }
}

// ------------------------------------------------------------------------------------------------

FunRBRV_calc_R_for_rhoPrime::FunRBRV_calc_R_for_rhoPrime(RBRV_entry_RV_base* rv1, RBRV_entry_RV_base* rv2, FlxFunction* rhoF, const bool corr_approx)
: eps_x(data->ConstantBox.getRef("rbrv_corr_eps_x")), eps_y(data->ConstantBox.getRef("rbrv_corr_eps_y")),
  ubound(data->ConstantBox.getRef("rbrv_corr_eps_ubound")), intn(data->ConstantBox.getRef("rbrv_corr_eps_intn")),
  rv1(rv1), rv2(rv2), rhoF(rhoF), rho(ZERO), y1(ZERO), y2(ZERO), r(ZERO), m1(ZERO), m2(ZERO), s1(ZERO), s2(ZERO), corr_approx(corr_approx), last_num(false)
{
  // TODO there are cases with an analytical solution!!!
  FunBase *b2;
  // z1
    fun = new FunConst(&y1);
    fun = new FunRBRV_y2x(fun,rv1,false);
    fun = new FunSub(fun,new FunConst(&m1));
    fun = new FunMult_Div(fun,new FunConst(&s1));
  // z2
    b2 = new FunConst(&y2);
    b2 = new FunRBRV_y2x(b2,rv2,false);
    b2 = new FunSub(b2,new FunConst(&m2));
    b2 = new FunMult_Div(b2,new FunConst(&s2));
  // z1*z2
    fun = new FunMult(fun,b2);
  // bi-normal
    std::vector<FunBase*>* pv = new std::vector<FunBase*>(3);
    (*pv)[0] = new FunConst(&y1);
    (*pv)[1] = new FunConst(&y2);
    (*pv)[2] = new FunConst(&r);
    b2 = new FunPDFn2(pv);
  // function to integrate
    fun = new FunMult(fun,b2);
        
  const tdouble gauss = tdouble(2);
  fun = new FunInteg(fun,&y2,new FunNumber(-ubound),new FunNumber(ubound),new FunNumber(gauss),new FunNumber(intn));
  fun = new FunInteg(fun,&y1,new FunNumber(-ubound),new FunNumber(ubound),new FunNumber(gauss),new FunNumber(intn));
  fun = new FunSub(fun,new FunConst(&rho));
}

const tdouble FunRBRV_calc_R_for_rhoPrime::calc()
{
  return calc_(true);
}

const tdouble FunRBRV_calc_R_for_rhoPrime::calc_(const bool throwErr)
{
  rv1->eval_para();
  rv2->eval_para();
  last_num = false;
  // initial evaluations
    rho = rhoF->calc();
    if (rho<=-ONE||rho>=ONE) {
      std::ostringstream ssV;
      ssV << "'" << rv1->name << "' and '" << rv2->name << "' correlated with " << GlobalVar.Double2String(rho) << ".";
      throw FlxException("FunRBRV_calc_R_for_rhoPrime::calc_1", "Specified correlation is not valid.", ssV.str() );
    }
    r = rho;
  // try to find an approximation
    // Reference: Ditlevsen, page 128
    if (corr_approx) {
        std::pair<std::string,RBRV_entry_RV_base*> t1(rv1->get_type(),rv1);
        std::pair<std::string,RBRV_entry_RV_base*> t2(rv2->get_type(),rv2);
        // normal
          if (t2.first=="normal") {
            std::swap(t1,t2);
          }
          if ( t1.first=="normal" ) {
            if (t2.first=="normal" ) {
              return ONE*rho;
            } else if ( t2.first=="uniform" ) {
              return tdouble(1.023)*rho;
            } else if ( t2.first=="gumbel" ) {
              return tdouble(1.031)*rho;
            } else if ( t2.first=="logn" ) {
              RBRV_entry_RV_lognormal* rv2_ln = dynamic_cast<RBRV_entry_RV_lognormal*>(t2.second);
              tdouble V_j = rv2_ln->get_CoeffOfVar_withoutEpsilon();
              return (V_j/sqrt(log(ONE+pow2(V_j))))*rho;
            }
          }
        // logn
          if (t2.first=="logn") {
            std::swap(t1,t2);
          }
          if (t1.first=="logn") {
            if (t2.first=="logn") {
              RBRV_entry_RV_lognormal* rv_ln = dynamic_cast<RBRV_entry_RV_lognormal*>(t1.second);
              tdouble V_i = rv_ln->get_CoeffOfVar_withoutEpsilon();
              rv_ln = dynamic_cast<RBRV_entry_RV_lognormal*>(t2.second);
              tdouble V_j = rv_ln->get_CoeffOfVar_withoutEpsilon();
              if (fabs(rho)<=GlobalVar.TOL()) return ZERO;
              return (log(ONE+rho*V_i*V_j)/rho/sqrt(log(ONE+pow2(V_i))*log(ONE+pow2(V_j))))*rho;
            } else if ( t2.first=="uniform" ) {
              RBRV_entry_RV_lognormal* rv_ln = dynamic_cast<RBRV_entry_RV_lognormal*>(t1.second);
              tdouble V_i = rv_ln->get_CoeffOfVar_withoutEpsilon();
              return (1.019+0.010*pow2(rho)+0.014*V_i+0.249*pow2(V_i))*rho;
            } else if ( t2.first=="gumbel" ) {
              RBRV_entry_RV_lognormal* rv_ln = dynamic_cast<RBRV_entry_RV_lognormal*>(t1.second);
              tdouble V_i = rv_ln->get_CoeffOfVar_withoutEpsilon();
              return (1.029+0.001*rho+0.004*pow2(rho)+0.014*V_i+0.233*pow2(V_i)-0.197*V_i*rho)*rho;
            }
          } 
        // Gumbel
          if (t2.first=="gumbel") {
            std::swap(t1,t2);
          }
          if (t1.first=="gumbel") {
            if ( t2.first=="gumbel" ) {
              return (1.064-0.069*rho+0.005*pow2(rho))*rho;
            } else if ( t2.first=="uniform" ) {
              return (1.055+0.015*pow2(rho))*rho;
            }
          } 
        // Uniform
          if (t1.first=="uniform" && t2.first=="uniform") {
            return (1.047-0.047*pow2(rho))*rho;
          }
    }
  // compute a numerical soluation
    m1 = rv1->get_mean_current_config();
    s1 = rv1->get_sd_current_config();
    m2 = rv2->get_mean_current_config();
    s2 = rv2->get_sd_current_config();
    tdouble res;
    try {
      tdouble start = rho/2;
        if (start<=-ONE) start=-0.999;
      tdouble end = (rho+ONE)/2;
        if (end>=ONE) end = 0.999;
      res = FlxFun_RootSearch_RegulaFalsi(fun,&r,start,end,eps_x,eps_y,NULL);
    } catch (FlxException&e ) {
      if (throwErr) {
        std::ostringstream ssV;
        ssV << "'" << rv1->name << "' and '" << rv2->name << "' correlated with " << GlobalVar.Double2String(rho) << "." << std::endl << e.what();
        throw FlxException("FunRBRV_calc_R_for_rhoPrime::calc_2", "The correlation of the two underlying random variables could not be computed.", ssV.str() );
      } else {
        res = rho;
      }
    }
    if (res<=-ONE||res>=ONE) {
      std::ostringstream ssV;
      ssV << "The computed correlation '" << GlobalVar.Double2String(res) << "' of the underlying standard normal random variables is invalid." << std::endl;
      ssV << "'" << rv1->name << "' and '" << rv2->name << "' correlated with " << GlobalVar.Double2String(rho) << ".";
      throw FlxException("FunRBRV_calc_R_for_rhoPrime::calc_3", ssV.str() );
    }
    last_num = true;
    return res;
}

const std::string FunRBRV_calc_R_for_rhoPrime::write()
{
  throw FlxException_NotImplemented("FunRBRV_calc_R_for_rhoPrime::write");
}

const bool FunRBRV_calc_R_for_rhoPrime::dependOn_Const(const tdouble*const thenumber)
{
  throw FlxException_NotImplemented("FunRBRV_calc_R_for_rhoPrime::dependOn_Const");
}

const bool FunRBRV_calc_R_for_rhoPrime::search_circref(FlxFunction* fcr)
{
  throw FlxException_NotImplemented("FunRBRV_calc_R_for_rhoPrime::search_circref");
}

// ------------------------------------------------------------------------------------------------

RBRV_vfset::RBRV_vfset(const bool internal, const std::string& name, const bool noID, const tuint Nentries, FlxMtxFun_base* vecfun, const tuint Nparents, RBRV_set_base** const parents)
: RBRV_set_parents(internal,0,name,Nparents,parents,noID), Nentries(Nentries), x_of_set(Nentries), vecfun(vecfun)
{

}

void RBRV_vfset::transform_y2x()
{
  vecfun->eval();
  x_of_set = vecfun->get_res_vec();
  #if FLX_DEBUG
    valid = true;
  #endif
}

void RBRV_vfset::set_x(const tdouble*const x_vec)
{
  #if FLX_DEBUG
    if (valid) {
      throw FlxException_Crude("RBRV_vfset::set_x");
    }
  #endif
  flxVec tv(x_vec,Nentries);
  x_of_set = tv;
  #if FLX_DEBUG
    valid = true;
  #endif
}

void RBRV_vfset::get_x(tdouble*const x_vec)
{
  #if FLX_DEBUG
  if (valid) {
  #endif
  flxVec tv(x_vec,Nentries);
  tv = x_of_set;
  #if FLX_DEBUG
  } else {
    std::ostringstream ssV;
    ssV << "The set '" << name << "' does not have valid realizations.";
    throw FlxException("RBRV_vfset::get_x", ssV.str() );
  }
  #endif
}

void RBRV_vfset::get_mean(tdouble*const m_vec)
{
  throw FlxException_NotImplemented("RBRV_vfset::get_mean");
}

void RBRV_vfset::get_sd(tdouble*const s_vec)
{
  throw FlxException_NotImplemented("RBRV_vfset::get_sd");
}

void RBRV_vfset::print(std::ostream& sout, const std::string prelim, tuint& counter, const bool printID)
{
  sout << prelim << "- " << name << " (" << get_NRV_only_this() << "/" <<  get_NOX_only_this() << ")" << std::endl;
  sout << prelim << "  " << "vector function" << std::endl;
  counter += get_NOX_only_this();
}

// ------------------------------------------------------------------------------------------------

RBRV_dirichlet::RBRV_dirichlet(const bool internal, const std::string& name, const bool noID, const tuint Nentries, FlxMtxFun_base* vecfun, const tuint Nparents, RBRV_set_base** const parents, flxVec* avec, const tuint idim)
: RBRV_set_parents(internal,(idim==0)?Nentries:idim,name,Nparents,parents,noID), Nentries(Nentries), x_of_set(Nentries), alpha_vec(Nentries), vecfun(vecfun)
{
  if (avec) {
    if (vecfun || avec->get_N()!=Nentries) {
      if (vecfun) {
        delete vecfun;
        vecfun = nullptr;
      }
      throw FlxException_Crude("RBRV_dirichlet::RBRV_dirichlet_01");
    }
    alpha_vec = *avec;
    if (alpha_vec.get_min()<ZERO) {
      throw FlxException("RBRV_dirichlet::RBRV_dirichlet_02","Parameter value must not be smaller than zero.");
    }
  }
}

RBRV_dirichlet::~RBRV_dirichlet()
{
  if (vecfun) delete vecfun;
}

void RBRV_dirichlet::get_pars()
{
  // get alpha-values (parameters of Dirichlet distribution)
    if (vecfun) {
      vecfun->eval();
      alpha_vec = vecfun->get_res_vec();
      return;
    }
}

void RBRV_dirichlet::transform_y2x()
{
  get_pars();
  // use Gamma distribution for random number generation
    const tdouble*const y_vec = y_of_set.get_tmp_vptr_const();
    for (tuint i=0;i<Nentries;++i) {
      const tdouble y_val = y_vec[i];
      if (y_val<=ZERO) {
        x_of_set[i] = flxgamma_rl_inv(alpha_vec[i],rv_Phi(y_val));
      } else {
        x_of_set[i] = flxgamma_ru_inv(alpha_vec[i],rv_Phi(-y_val));
      }
    }
  // normalize vector
    x_of_set /= x_of_set.get_sum();
  #if FLX_DEBUG
    valid = true;
  #endif
}

void RBRV_dirichlet::set_x(const tdouble*const x_vec)
{
  #if FLX_DEBUG
    if (valid) {
      throw FlxException_Crude("RBRV_dirichlet::set_x_01");
    }
  #endif
  flxVec tv(x_vec,Nentries);
  x_of_set = tv;
  #if FLX_DEBUG
    valid = true;
  #endif
}

void RBRV_dirichlet::get_x(tdouble*const x_vec)
{
  #if FLX_DEBUG
  if (valid) {
  #endif
  flxVec tv(x_vec,Nentries);
  tv = x_of_set;
  #if FLX_DEBUG
  } else {
    std::ostringstream ssV;
    ssV << "The set '" << name << "' does not have valid realizations.";
    throw FlxException("RBRV_dirichlet::get_x", ssV.str() );
  }
  #endif
}
const bool RBRV_dirichlet::check_xVec(const tdouble* xp)
{
  flxVec tv(xp,Nentries);
  return (fabs(tv.get_sum()-ONE)<GlobalVar.TOL());
}

void RBRV_dirichlet::get_mean(tdouble*const m_vec)
{
  get_pars();
  const tdouble alpha0 = alpha_vec.get_sum();
  for (tuint i=0;i<Nentries;++i) m_vec[i] = alpha_vec[i]/alpha0;
}

void RBRV_dirichlet::get_sd(tdouble*const s_vec)
{
  get_pars();
  const tdouble alpha0 = alpha_vec.get_sum();
  for (tuint i=0;i<Nentries;++i) {
    const tdouble ai = alpha_vec[i]/alpha0;
    s_vec[i] = sqrt(ai*(ONE-ai)/(alpha0+ONE));
  }
}

void RBRV_dirichlet::print(std::ostream& sout, const std::string prelim, tuint& counter, const bool printID)
{
  get_pars();
  sout << prelim << "- " << name << " (" << get_NRV_only_this() << "/" <<  get_NOX_only_this() << ")" << std::endl;
  sout << prelim << "  " << "Dirichlet distribution" << std::endl;
  sout << prelim << "  " << "parameter vector: " << alpha_vec << std::endl;
  counter += get_NOX_only_this();
}

// ------------------------------------------------------------------------------------------------

RBRV_multinomial::RBRV_multinomial(const bool internal, const std::string& name, const bool noID, const tuint Nentries, FlxMtxFun_base* vecfun, const tuint Nparents, RBRV_set_base** const parents, const tuint Ntrials, flxVec* avec)
: RBRV_dirichlet(internal, name, noID, Nentries, vecfun, Nparents, parents, avec, Ntrials), Ntrials(Ntrials)
{
  if (avec) {
    alpha_vec /= alpha_vec.get_sum();
  }
}

RBRV_multinomial::~RBRV_multinomial()
{

}

void RBRV_multinomial::get_pars()
{
  RBRV_dirichlet::get_pars();
  if (vecfun) {
    alpha_vec /= alpha_vec.get_sum();
  }
}

void RBRV_multinomial::transform_y2x()
{
  get_pars();
  // for each trial, generate a realization of a categorical distribution with parameter vector alpha_vec
    x_of_set.set_zero();
    const tdouble*const y_vec = y_of_set.get_tmp_vptr_const();
    for (tuint i=0;i<Ntrials;++i) {
      const tuint rid = RndCreator->gen_smp_index2_help(rv_Phi(y_vec[i]),alpha_vec);
      x_of_set[rid] += ONE;
    }
  #if FLX_DEBUG
    valid = true;
  #endif
}

const bool RBRV_multinomial::check_xVec(const tdouble* xp)
{
  flxVec tv(xp,Nentries);
  return (tv.get_min()>=ZERO && fabs(tv.get_sum()-tdouble(Ntrials))/Ntrials<GlobalVar.TOL());
}

void RBRV_multinomial::get_mean(tdouble*const m_vec)
{
  get_pars();
  for (tuint i=0;i<Nentries;++i) m_vec[i] = alpha_vec[i]*Ntrials;
}

void RBRV_multinomial::get_sd(tdouble*const s_vec)
{
  get_pars();
  for (tuint i=0;i<Nentries;++i) {
    s_vec[i] = sqrt(alpha_vec[i]*(ONE-alpha_vec[i])*Ntrials);
  }
}

void RBRV_multinomial::print(std::ostream& sout, const std::string prelim, tuint& counter, const bool printID)
{
  get_pars();
  sout << prelim << "- " << name << " (" << get_NRV_only_this() << "/" <<  get_NOX_only_this() << ")" << std::endl;
  sout << prelim << "  " << "Multinomial distribution" << std::endl;
  sout << prelim << "  " << "Number of trials:    " << Ntrials << std::endl;
  sout << prelim << "  " << "event probabilities: " << alpha_vec << std::endl;
  counter += get_NOX_only_this();
}

// ------------------------------------------------------------------------------------------------
