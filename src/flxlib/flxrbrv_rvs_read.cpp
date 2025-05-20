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

#include "flxrbrv_rvs_read.h"
#include "flxobjrbrv.h"
#include "flxBayDA.h"

FlxObjRBRV_set_creator_box* RBRV_entry_read_base::rbrv_set_creator = NULL;


FlxObjRBRV_set_creator_box::FlxObjRBRV_set_creator_box()
{
  RBRV_entry_read_base::set_rbrv_set_creator(this);
}

FlxObjRBRV_set_creator_box::~FlxObjRBRV_set_creator_box()
{
  for (std::map<std::string,FlxObjRBRV_set_creator*>::iterator it=cmap.begin(); it!=cmap.end(); ++it) {
    FlxObjRBRV_set_creator* tp = it->second;
    delete tp;
  }
  cmap.clear();
}

void FlxObjRBRV_set_creator_box::create_new(const std::string& name, FlxObjRBRV_set_creator* ele)
{
  cmap.insert( std::pair<std::string,FlxObjRBRV_set_creator*>(name,ele) );
}

FlxObjRBRV_set_creator* FlxObjRBRV_set_creator_box::get_creator(const std::string& set_name, const bool throwErr)
{
  std::map<std::string,FlxObjRBRV_set_creator*>::const_iterator pos;
  pos = cmap.find(set_name);
  if ( pos == cmap.end() ) {
    if (throwErr) {
      std::ostringstream ssV;
      ssV << "No rbrv-set with the name (" << set_name << ") is in the construction-process.";
      throw FlxException("FlxObjRBRV_set_creator_box::get_creator", ssV.str() );
    } else {
      return NULL;
    }
  }
  return pos->second;
}

void FlxObjRBRV_set_creator_box::register_set(const std::string& name, RBRV_set_box& box)
{
  FlxObjRBRV_set_creator* crtr = get_creator(name);
  crtr->register_set(box,true);
  delete crtr;
  cmap.erase(name);
}


FlxObjRBRV_set_creator::FlxObjRBRV_set_creator( const std::string& set_name, RBRV_set_baseDPtr parents, const tuint Nparents, const bool allow_x2y)
: set_name(set_name), is_Nataf(false), is_Nataf_evalOnce(false), parents(parents), Nparents(Nparents), allow_x2y(allow_x2y), rID(0)
{

}

FlxObjRBRV_set_creator::FlxObjRBRV_set_creator( RBRV_set_box& box, const std::string& set_name, RBRV_set_baseDPtr parents, const tuint Nparents, const bool allow_x2y, std::vector< RBRV_entry_read_base* >& set_entriesV)
: set_name(set_name), is_Nataf(false), is_Nataf_evalOnce(false), parents(parents), Nparents(Nparents), allow_x2y(allow_x2y), rID(0)
{
  for (size_t i=0;i<set_entriesV.size();++i) {
    add_entry(box,set_entriesV[i]);
  }
}

FlxObjRBRV_set_creator::FlxObjRBRV_set_creator(const std::string& set_name, const bool eval_once)
: set_name(set_name), is_Nataf(true), is_Nataf_evalOnce(eval_once), parents(NULL), Nparents(0), allow_x2y(true), rID(0)
{

}

FlxObjRBRV_set_creator::~FlxObjRBRV_set_creator()
{
  if (parents) delete [] parents;
  for (tuint i=0;i<set_entries.size();++i) {
    delete set_entries[i];
  }
}

const bool FlxObjRBRV_set_creator::get_Nataf_evalOnce() const
{
  return (is_Nataf && is_Nataf_evalOnce );
}

void FlxObjRBRV_set_creator::add_entry(RBRV_set_box& box, RBRV_entry_RV_base* ep, FlxFunction* csVal, const std::string csNam, const bool csFix)
{
  try {
    for (tuint i=0;i<set_entries.size();++i) {
      if (ep->name==set_entries[i]->name) {
        throw FlxException("FlxObjRBRV_set_creator::add_entry_a01", "An entry with name '" + set_entries[i]->name + "' does already exist.");
      }
    }
    if (ep->get_iID()!=rID) {
      throw FlxException_Crude("FlxObjRBRV_set_creator::add_entry_a02");
    }
    ++rID;  // manually increase the running ID
    // register correlation (only for Rosenblatt transformation)
      if (csVal) {
        // consistency checks
            if (is_Nataf) {
                throw FlxException_NeglectInInteractive("FlxObjRBRV_set_creator::add_entry_a03", "Setting a correlation pair is not allowed for sets of random variables based on the Nataf transformation.");
            }
            // make sure that the current entry is a true random variable
                RBRV_entry_RV_base* rv_2 = dynamic_cast<RBRV_entry_RV_base*>(ep);
                if (rv_2==nullptr) {
                    throw FlxException("FlxObjRBRV_set_creator::add_entry_a04", "A correlation cannot be specified for'" + ep->name + "'." );
                }
        // search the random variable to correlate with
            RBRV_entry* rv_1b = nullptr;
            for (tuint i=0;i<set_entries.size();++i) {
                if (set_entries[i]->name==csNam) {
                    rv_1b = set_entries[i];
                    break;
                }
            }
            if (rv_1b==nullptr) {
                throw FlxException("FlxObjRBRV_set_creator::add_entry_a05", "An entry with name '" + csNam + "' was not found in the set." );
            }
        // make sure that it is also a true random variable
            RBRV_entry_RV_base* rv_1 = dynamic_cast<RBRV_entry_RV_base*>(rv_1b);
            if (rv_1==nullptr) {
              throw FlxException("FlxObjRBRV_set_creator::add_entry_a06", "A correlation cannot be specified for'" + rv_1b->name + "'." );
            }
        rv_2->set_corr(rv_1,csVal,csFix,true);
        delete csVal; csVal = nullptr;
      }
    box.register_entry(ep);
    set_entries.push_back(ep); ep=nullptr;
  } catch (FlxException &e) {
    if (ep) delete ep;
    if (csVal) delete csVal;
    throw;
  }
}

void FlxObjRBRV_set_creator::add_entry(RBRV_set_box& box, RBRV_entry_read_base* entry)
{
  const std::string family = set_name + "::";
  RBRV_entry* ep = NULL;
  try {
    if (is_Nataf && is_Nataf_evalOnce ) entry->activate_eval_once();
    ep = entry->generate_entry(family,rID);
    for (tuint i=0;i<set_entries.size();++i) {
      if (ep->name==set_entries[i]->name) {
        throw FlxException("FlxObjRBRV_set_creator::add_entry_b01", "An entry with name '" + set_entries[i]->name + "' does already exist.");
      }
    }
    box.register_entry(ep);
    set_entries.push_back(ep); ep=NULL;
    entry->generate_corr(set_entries,set_entries.size()-1,is_Nataf);
  } catch (FlxException &e) {
    if (ep) delete ep;
    throw;
  }
}

RBRV_entry_RV_base* FlxObjRBRV_set_creator::get_rv(const std::string rvn, const bool throwErr)
{
  for (tuint i=0;i<set_entries.size();++i) {
    RBRV_entry* te = set_entries[i];
    if (te->name==rvn) {
      RBRV_entry_RV_base* tr = dynamic_cast<RBRV_entry_RV_base*>(te);
      if (tr==NULL) {
        if (throwErr) {
          std::ostringstream ssV;
          ssV << "The entry '" << rvn << "' is not a basic random variable.";
          throw FlxException("FlxObjRBRV_set_creator::get_rv_1", ssV.str() );
        }
        return NULL;
      } else {
        return tr;
      }
    }
  }
  if (throwErr) {
    std::ostringstream ssV;
    ssV << "An entry with name '" << rvn << "' does not exist.";
    throw FlxException("FlxObjRBRV_set_creator::get_rv_2", ssV.str() );
  }
  return NULL;
}

const tuint FlxObjRBRV_set_creator::get_rvID(const std::string rvn)
{
  for (tuint i=0;i<set_entries.size();++i) {
    RBRV_entry* te = set_entries[i];
    if (te->name==rvn) {
      return i;
    }
  }
  throw FlxException_Crude("FlxObjRBRV_set_creator::get_rvID_2");
}

void FlxObjRBRV_set_creator::add_corr(const std::string& rv1, const std::string& rv2, const tdouble rho, const bool corr_approx, const bool rhogauss, const bool dolog)
{
  // consistency checks
    if (rv1==rv2) {
      std::ostringstream ssV;
      ssV << "The correlation of a random variable with itself is always one and cannot be changed.";
      throw FlxException("FlxObjRBRV_set_creator::add_corr_1", ssV.str() );
    }
    if (fabs(rho)>ONE) {
      std::ostringstream ssV;
      ssV << "The specified correlation value (" << GlobalVar.Double2String(rho) << ") must be in [-1,1].";
      throw FlxException("FlxObjRBRV_set_creator::add_corr_2", ssV.str() );
    }
  // compute equivalent correlation coefficient of the underlying Gaussian random variables
    tdouble rhoG = rho;
    if (rhogauss==false) {
      RBRV_entry_RV_base* rv1e = get_rv(rv1,true);
      RBRV_entry_RV_base* rv2e = get_rv(rv2,true);
      FunRBRV_calc_R_for_rhoPrime corrF( rv1e,rv2e,new FlxFunction(new FunNumber(rho)), corr_approx );
      rhoG = corrF.calc();
      if (dolog) {
        GlobalVar.slog(4) << "Nataf-set '" << set_name << "': correlation betwenn '" << rv1e->name << "' and '" << rv2e->name << "' is " << GlobalVar.Double2String(rho) << "  ("
          << GlobalVar.Double2String(rhoG) << " " << (corrF.last_was_numerical()?"numerical":"analytical") << ")" << std::endl;
      }
    }
  // memorize the correlation
    std::pair<std::string,std::string> tp;
    tp.first = (rv1<rv2)?rv1:rv2;
    tp.second = (rv1<rv2)?rv2:rv1;
    if (cormap.insert(std::pair<std::pair<std::string,std::string>,tdouble>(tp,rhoG)).second==false) {
      std::ostringstream ssV;
      ssV << "The correlation between '" << rv1 << "' and '" << rv2 << "' must be specified only once.";
      throw FlxException("FlxObjRBRV_set_creator::add_corr_3", ssV.str() );
    }
}

RBRV_set* FlxObjRBRV_set_creator::register_set_rbrv(RBRV_set_box& box, const bool doreg)
{
  if (is_Nataf) throw FlxException_Crude("FlxObjRBRV_set_creator::register_set_rbrv");
  const tuint Nentries = set_entries.size();
  RBRV_entry** entries = new RBRV_entry*[Nentries];
  RBRV_set* ts = NULL;
  for (tuint i=0;i<Nentries;++i) entries[i] = set_entries[i];
  set_entries.clear();  // memory management handed over to 'entries'
  try {
    // allocate the set
      ts = new RBRV_set(false,rID,set_name,false,Nentries,entries,Nparents,parents,allow_x2y);
      parents = NULL;
      entries = NULL;
    if (doreg) {
      box.register_set(ts);
      GlobalVar.slog(4) << "rbrv_set: created new set '" << set_name << "'." << std::endl;
    }
  } catch (FlxException &e) {
    FLXMSG("FlxObjRBRV_set_creator::register_set_rbrv",1);
    if (parents) delete [] parents;
    if (entries) {
      for (tuint i=0;i<Nentries;++i) {
        if (entries[i]) {
          delete entries[i];
        } else {
          break;
        }
      }
      delete [] entries;
    }
    if (ts) delete ts;
    throw;
  }
  return ts;
}

RBRV_set_Nataf* FlxObjRBRV_set_creator::register_set_Nataf(RBRV_set_box& box, const bool doreg)
{
  const tuint Nentries = set_entries.size();
  // handle correlations
    FlxMtxSparsLTri* L = NULL;
    if (cormap.size()>0) {
      // Assemble matrices
        FlxMtxSym rhoCorr(Nentries);
        for (std::map<std::pair<std::string,std::string>,tdouble>::iterator it=cormap.begin(); it!=cormap.end(); ++it) {
          const tuint id1 = get_rvID(it->first.first);
          const tuint id2 = get_rvID(it->first.second);
          rhoCorr.operator()(id1,id2) = it->second;
        }
        for (tuint i=0;i<Nentries;++i) {
          rhoCorr.operator()(i,i) = ONE;
        }
        FlxMtxSparsSym rhoPrime(rhoCorr);
      // CholeskyDecomposition
        FlxMtxLTri Lt(rhoPrime.nrows());
        Lt.CholeskyDec(rhoPrime);
        L = new FlxMtxSparsLTri(Lt);
    }
  RBRV_entry** entries = new RBRV_entry*[Nentries];
  RBRV_set_Nataf* ts = NULL;
  for (tuint i=0;i<Nentries;++i) entries[i] = set_entries[i];
  set_entries.clear();  // memory management handed over to 'entries'
  try {
    // allocate the set
      ts = new RBRV_set_Nataf(false,rID,set_name,false,Nentries,entries,L);
      parents = NULL;
      entries = NULL;
      L = NULL;
    if (doreg) {
      box.register_set(ts);
      GlobalVar.slog(4) << "rbrv_set: created new Nataf-set '" << set_name << "'." << std::endl;
    }
  } catch (FlxException &e) {
    FLXMSG("FlxObjRBRV_set_creator::register_set_Nataf",1);
    if (entries) {
      for (tuint i=0;i<Nentries;++i) {
        if (entries[i]) {
          delete entries[i];
        } else {
          break;
        }
      }
      delete [] entries;
    }
    if (ts) delete ts;
    if (L) delete L;
    throw;
  }
  return ts;
  
}

RBRV_set_base* FlxObjRBRV_set_creator::register_set(RBRV_set_box& box, const bool doreg)
{
  if (is_Nataf) {
    return register_set_Nataf(box,doreg);
  } else {
    return register_set_rbrv(box,doreg);
  }
}



void RBRV_entry_read_base::read_parents(std::vector< FlxString* >& set_parents, const bool errSerious)
{
  // read the parents
    if (reader->whatIsNextChar()=='(') {
      reader->getChar('(',errSerious);
      while (reader->whatIsNextChar()!=')') {
        set_parents.push_back(new FlxString(false,errSerious));
        // continue loop?
          if (reader->whatIsNextChar()==',') reader->getChar(',');
          else break;
      }
      reader->getChar(')',errSerious);
    }
}

RBRV_entry_RV_base* RBRV_entry_read_base::generate_entry_rv(const bool errSerious)
{
  const std::string family="dummy";
  tuint running_iID=0;
  RBRV_entry *rep_ = generate_entry(family,running_iID);
  // make sure it is a random variable
    RBRV_entry_RV_base* rep = dynamic_cast<RBRV_entry_RV_base*>(rep_);
    if (rep==NULL) {
      std::ostringstream ssV;
      ssV << "The specified RBRV has the wrong type.";
      delete rep_;
      FlxError(errSerious,"RBRV_entry_read_base::generate_entry_rv", ssV.str(), reader->getCurrentPos() );
    }
  return rep;
}

void RBRV_entry_read_base::read_corr(const bool errSerious)
{
  if (reader->whatIsNextChar()==':') {
    reader->getChar(':',errSerious);
    if (reader->whatIsNextChar()=='!') {
      corrFixed = true;
      reader->getChar('!');
    } else {
      corrFixed = false;
    }
    reader->getChar('(');
    corrName = new FlxString(false,errSerious);
    reader->getChar('=');
    corrVal = new FlxFunction(funReader,errSerious);
    reader->getChar(')');
  }
}

void RBRV_entry_read_base::generate_corr(std::vector< RBRV_entry* >& entries, const tuint IDinSet, const bool ensureEmpty)
{
  if (corrName==NULL) return;
  if (ensureEmpty) {
    std::ostringstream ssV;
    ssV << "This correlation statement is not allowed in a Nataf-set.";
    throw FlxException("RBRV_entry_read_base::generate_corr_0", ssV.str() );
  }
  #if FLX_DEBUG
    if (corrVal==NULL) throw FlxException_Crude("RBRV_entry_read_base::generate_corr_1");
  #endif
  // make sure that the current entry is a true random variable
    RBRV_entry_RV_base* rv_2 = dynamic_cast<RBRV_entry_RV_base*>(entries[IDinSet]);
    if (rv_2==NULL) {
      std::ostringstream ssV;
      ssV << "A correlation cannot be specified for'" << entries[IDinSet]->name << "'.";
      throw FlxException("RBRV_entry_read_base::generate_corr_2", ssV.str() );
    }
  // search the random variable to correlate with
    const std::string csNam = corrName->eval(true);
    RBRV_entry* rv_1b = NULL;
    for (tuint i=0;i<IDinSet;++i) {
      if (entries[i]->name==csNam) {
        rv_1b = entries[i];
        break;
      }
    }
    if (rv_1b==NULL) {
      std::ostringstream ssV;
      ssV << "An entry with name '" << csNam << "' was not found in the set.";
      throw FlxException("RBRV_entry_read_base::generate_corr_3", ssV.str() );
    }
  // make sure that it is also a true random variable
    RBRV_entry_RV_base* rv_1 = dynamic_cast<RBRV_entry_RV_base*>(rv_1b);
    if (rv_1==NULL) {
      std::ostringstream ssV;
      ssV << "A correlation cannot be specified for'" << entries[IDinSet]->name << "'.";
      throw FlxException("RBRV_entry_read_base::generate_corr_4", ssV.str() );
    }
  rv_2->set_corr(rv_1,corrVal,corrFixed,true);
}

void RBRV_entry_read_base::read(std::vector<RBRV_entry_read_base*>& set_entries, std::vector<FlxString*>& set_parents, const bool errSerious)
{
  read_parents(set_parents,errSerious);
  // read the entries
    reader->getChar('{',errSerious);
    while (true) {                                // at least one entry has to be defined
      set_entries.push_back(read_entry());
      set_entries[set_entries.size()-1]->read_corr(errSerious);
      // continue loop?
        if (reader->whatIsNextChar()==',') reader->getChar(',',errSerious);
        else break;
    }
    reader->getChar('}',errSerious);
}

RBRV_entry_read_base* RBRV_entry_read_base::read_entry(const bool readName, const bool readBrakets)
{
  const std::string sid = reader->getWord(true);
  if (readBrakets==false) {
    if (sid!="stdn") reader->getChar(',');
  }
  RBRV_entry_read_base* ptr = NULL;
  if (sid=="fun") {
    ptr = new RBRV_entry_read_fun(readName,readBrakets);
  } else if (sid=="stdn") {
    ptr = new RBRV_entry_read_stdn(readName,readBrakets);
  } else if (sid=="normal") {
    ptr = new RBRV_entry_read_normal(readName,readBrakets);
  } else if (sid=="logn") {
    ptr = new RBRV_entry_read_logn(readName,readBrakets);
  } else if (sid=="uniform") {
    ptr = new RBRV_entry_read_uniform(readName,readBrakets);
  } else if (sid=="gumbel") {
    ptr = new RBRV_entry_read_Gumbel(readName,readBrakets);
  } else if (sid=="normal_trunc") {
    ptr = new RBRV_entry_read_normal_trunc(readName,readBrakets);
  } else if (sid=="beta") {
    ptr = new RBRV_entry_read_beta(readName,readBrakets);
  } else if (sid=="exponential") {
    ptr = new RBRV_entry_read_exponential(readName,readBrakets);
  } else if (sid=="gamma") {
    ptr = new RBRV_entry_read_gamma(readName,readBrakets);
  } else if (sid=="poisson") {
    ptr = new RBRV_entry_read_Poisson(readName,readBrakets);
  } else if (sid=="binomial") {
    ptr = new RBRV_entry_read_Binomial(readName,readBrakets);
  } else if (sid=="cauchy") {
    ptr = new RBRV_entry_read_Cauchy(readName,readBrakets);
  } else if (sid=="weibull") {
    ptr = new RBRV_entry_read_Weibull(readName,readBrakets);
  } else if (sid=="chisquared") {
    ptr = new RBRV_entry_read_ChiSquared(true,readName,readBrakets);
  } else if (sid=="chi") {
    ptr = new RBRV_entry_read_ChiSquared(false,readName,readBrakets);
  } else if (sid=="studentst") {
    ptr = new RBRV_entry_read_StudentsT(readName,readBrakets);
  } else if (sid=="studentstgen") {
    ptr = new RBRV_entry_read_StudentsT_generalized(readName,readBrakets);
  } else if (sid=="laplace") {
    ptr = new RBRV_entry_read_Laplace(readName,readBrakets);
  } else if (sid=="usertransform") {
    ptr = new RBRV_entry_read_UserTransform(readName,readBrakets);
  } else if (sid=="truncated") {
    ptr = new RBRV_entry_read_Truncated(readName,readBrakets);
  } else if (sid=="maxmintransform") {
    ptr = new RBRV_entry_read_maxminTransform(readName,readBrakets);
  } else if (sid=="bayda") {
    ptr = new RBRV_entry_read_bayDA(readName,readBrakets);
  } else {
    std::ostringstream ssV;
    ssV << "Unknown keyword '" << sid << "'.";
    throw FlxException("RBRV_entry_read_base::read", ssV.str() );
  }
  if (readBrakets) {
    reader->getChar(')');
  }
  return ptr;
}

RBRV_entry_RV_base* RBRV_entry_read_base::read_gen_entry(bool errSerious)
{
  RBRV_entry_read_base* rbp = read_entry(false,false);
  RBRV_entry_RV_base* rep = rbp->generate_entry_rv(errSerious);
  delete rbp;
  return rep;
}

void RBRV_entry_read_base::generate_set_base_check_name(RBRV_set_box& box, const std::string& name)
{
  // check if a set with that name does already exist
    if (box.get_set(name,false)) {
      std::ostringstream ssV;
      ssV << "A rbrv-set with the same name (" << name << ") is already defined.";
      throw FlxException("RBRV_entry_read_base::generate_set_base_check_name_1", ssV.str() );
    }
  // check if a set with that name is waiting to be created
    if (rbrv_set_creator->get_creator(name,false)) {
      std::ostringstream ssV;
      ssV << "A rbrv-set with the same name (" << name << ") is already declared.";
      throw FlxException("RBRV_entry_read_base::generate_set_base_check_name_2", ssV.str() );
    }
}

void RBRV_entry_read_base::generate_set_base(RBRV_set_box& box, const std::string& name, std::vector< FlxString* > set_parents, RBRV_set_baseDPtr& parents)
{
  RBRV_entry_read_base::generate_set_base_check_name(box,name);
  const tuint Nparents = set_parents.size();
  parents = ((Nparents==0)?nullptr:(new RBRV_set_base*[Nparents]));
  try {
    // get the parents
      for (tuint i=0;i<Nparents;++i) {
        const std::string pn = set_parents[i]->eval_word(true);
        parents[i] = box.get_set(pn,true);
      }
  } catch (FlxException& e) {
    FLXMSG("RBRV_entry_read_base::generate_set_base_a01",1);
    if (parents) {
      delete [] parents;
      parents = nullptr;
    }
    throw;
  }
}

void RBRV_entry_read_base::generate_set_base(RBRV_set_box& box, std::vector<std::string> set_parents, RBRV_set_baseDPtr& parents)
{
  const tuint Nparents = set_parents.size();
  parents = ((Nparents==0)?nullptr:(new RBRV_set_base*[Nparents]));
  try {
    // get the parents
      for (tuint i=0;i<Nparents;++i) {
        parents[i] = box.get_set(set_parents[i],true);
      }
  } catch (FlxException& e) {
    FLXMSG("RBRV_entry_read_base::generate_set_base_b01",1);
    if (parents) {
      delete [] parents;
      parents = nullptr;
    }
    throw;
  }
}

RBRV_entry_read_base::RBRV_entry_read_base(const bool readName, const bool readBrakets, const bool hasParas)
:nameF(NULL), corrName(NULL), corrVal(NULL), corrFixed(false), eval_once(false)
{
  if (readBrakets) {
    reader->getChar('(');
  }
  if (readName) {
    try {
      nameF = new FlxString(false,true);
      if (hasParas) reader->getChar(',');
    } catch (FlxException& e) {
      FLXMSG("RBRV_entry_read_base::RBRV_entry_read_base",1);
      if (nameF) delete nameF;
    }
  } else {
    nameF = new FlxString(new FlxString_String("dummy",true),false);
  }
}

RBRV_entry_read_base::~RBRV_entry_read_base()
{
  if (nameF) delete nameF;
  if (corrName) delete corrName;
  if (corrVal) delete corrVal;
}

void RBRV_entry_read_base::read_eval_once()
{
  reader->getWord("eval_once");
  reader->getChar('=');
  eval_once = reader->getBool();
}

FlxFunction* RBRV_entry_read_base::read_opt_para(const std::string& pName)
{
  FlxFunction* res = NULL;
  if (reader->whatIsNextString(pName.length(),true)==pName) {
    reader->getWord(pName.c_str());
    reader->getChar('=');
    res = new FlxFunction(funReader);
  }
  return res;
}


RBRV_entry_read_fun::RBRV_entry_read_fun(const bool readName, const bool readBrakets)
: RBRV_entry_read_base(readName,readBrakets), fun(NULL)
{
  reader->getWord("fun");
  reader->getChar('=');
  fun = new FlxFunction(funReader);
}

RBRV_entry_read_fun::~RBRV_entry_read_fun()
{
  if (fun) delete fun;
}

RBRV_entry* RBRV_entry_read_fun::generate_entry(const std::string& family, tuint& running_iID)
{
  const std::string name = family + nameF->eval_word(true,false,true);
  return new RBRV_entry_fun(name,new FlxFunction(*fun));
}


RBRV_entry_read_normal::RBRV_entry_read_normal(const bool readName, const bool readBrakets)
: RBRV_entry_read_base(readName,readBrakets), pid(0), p1(NULL), p2(NULL), p3(NULL), p4(NULL)
{
  try {
    const std::string strV = reader->getWord(true);
    if (strV == "mu") {                        // mean, standard deviation
      reader->getChar('=');
      p1 = new FlxFunction(funReader);
      reader->getChar(',');
      reader->getWord("sd");
      reader->getChar('=');
      p2 = new FlxFunction(funReader);
      pid = 0;
    } else if (strV == "pr") {                // quantile values
      reader->getChar('(');
      p1 = new FlxFunction(funReader);
      reader->getChar(')');
      reader->getChar('=');
      p2 = new FlxFunction(funReader);
      reader->getChar(',');
      reader->getWord("pr",false);
      reader->getChar('(');
      p3 = new FlxFunction(funReader);
      reader->getChar(')');
      reader->getChar('=');
      p4 = new FlxFunction(funReader);
      pid = 1;
    } else if (strV == "cov") {                // C.o.V. and quantile value
      reader->getChar('=');
      p1 = new FlxFunction(funReader);                   // actually it is C.o.V.
      reader->getChar(',');
      reader->getWord("pr");
      reader->getChar('(');
      p2 = new FlxFunction(funReader);
      reader->getChar(')');
      reader->getChar('=');
      p3 = new FlxFunction(funReader);
      pid = 2;
    } else if (strV == "sd") {                // std.dev. and quantile value
      reader->getChar('=');
      p1 = new FlxFunction(funReader);  
      reader->getChar(',');
      reader->getWord("pr");
      reader->getChar('(');
      p2 = new FlxFunction(funReader);
      reader->getChar(')');
      reader->getChar('=');
      p3 = new FlxFunction(funReader);
      pid = 3;
    } else {
      std::ostringstream ssV;
      ssV << "Keyword '" << strV << "' not known.";
      throw FlxException("RBRV_entry_read_normal::RBRV_entry_read_normal_02", ssV.str() );
    }
    
    if (reader->whatIsNextChar()==',') {
      reader->getChar(',');
      read_eval_once();
    }
    
  } catch (FlxException& e) {
    FLXMSG("RBRV_entry_read_normal::RBRV_entry_read_normal_02",1);
    if (p1) delete p1;
    if (p2) delete p2;
    if (p3) delete p3;
    if (p4) delete p4;
    throw;
  }
}

RBRV_entry_read_normal::~RBRV_entry_read_normal()
{
  if (p1) delete p1;
  if (p2) delete p2;
  if (p3) delete p3;
  if (p4) delete p4;
}

RBRV_entry* RBRV_entry_read_normal::generate_entry(const std::string& family, tuint& running_iID)
{
  const std::string name = family + nameF->eval_word(true,false,true);
  return new RBRV_entry_RV_normal(name,running_iID++,pid,
    p1?(new FlxFunction(*p1)):NULL,
    p2?(new FlxFunction(*p2)):NULL,
    p3?(new FlxFunction(*p3)):NULL,
    p4?(new FlxFunction(*p4)):NULL,
    eval_once);
}

RBRV_entry_read_stdn::RBRV_entry_read_stdn(const bool readName, const bool readBrakets)
: RBRV_entry_read_base(readName, readBrakets,false)
{

}

RBRV_entry_read_stdn::~RBRV_entry_read_stdn()
{

}

RBRV_entry* RBRV_entry_read_stdn::generate_entry(const std::string& family, tuint& running_iID)
{
  const std::string name = family + nameF->eval_word(true,false,true);
  return new RBRV_entry_RV_stdN(name,running_iID++);
}

RBRV_entry_read_logn::RBRV_entry_read_logn(const bool readName, const bool readBrakets)
: RBRV_entry_read_base(readName,readBrakets), pid(0),p1(NULL),p2(NULL),p3(NULL),p4(NULL),epsilon(NULL)
{
  try {
    const std::string strV = reader->getWord(true);
    if (strV == "lambda") {
      reader->getChar('=');
      p1 = new FlxFunction(funReader);
      reader->getChar(',');
      reader->getWord("zeta");
      reader->getChar('=');
      p2 = new FlxFunction(funReader);
      pid = 0;
    } else if (strV == "mu") {
      reader->getChar('=');
      p1 = new FlxFunction(funReader);
      reader->getChar(',');
      reader->getWord("sd");
      reader->getChar('=');
      p2 = new FlxFunction(funReader);
      pid = 1;
    } else if (strV == "mode") {
      reader->getChar('=');
      p1 = new FlxFunction(funReader);
      reader->getChar(',');
      const std::string strV2 = reader->getWord(true);
      if (strV2=="sd") {
        pid = 2;
      } else if (strV2=="cov") {
        pid = 7;
      } else {
        std::ostringstream ssV;
        ssV << "Keyword '" << strV << "' not known.";
        throw FlxException("RBRV_entry_read_logn::RBRV_entry_read_logn_3", ssV.str() );
      }
      reader->getChar('=');
      p2 = new FlxFunction(funReader);
    } else if (strV == "median") {
      reader->getChar('=');
      p1 = new FlxFunction(funReader);
      reader->getChar(',');
      const std::string strV2 = reader->getWord(true);
      if (strV2=="sd") {
        pid = 3;
      } else if (strV2=="cov") {
        pid = 5;
      } else {
        std::ostringstream ssV;
        ssV << "Keyword '" << strV << "' not known.";
        throw FlxException("RBRV_entry_read_logn::RBRV_entry_read_logn_2", ssV.str() );
      }
      reader->getChar('=');
      p2 = new FlxFunction(funReader);
    } else if (strV == "pr") {
      reader->getChar('(');
      p1 = new FlxFunction(funReader);
      reader->getChar(')');
      reader->getChar('=');
      p2 = new FlxFunction(funReader);
      reader->getChar(',');
      reader->getWord("pr",false);
      reader->getChar('(');
      p3 = new FlxFunction(funReader);
      reader->getChar(')');
      reader->getChar('=');
      p4 = new FlxFunction(funReader);
      pid = 4;
    } else if (strV == "cov") {
      reader->getChar('=');
      p1 = new FlxFunction(funReader);
      reader->getChar(',');
      reader->getWord("pr",false);
      reader->getChar('(');
      p2 = new FlxFunction(funReader);
      reader->getChar(')');
      reader->getChar('=');
      p3 = new FlxFunction(funReader);
      pid = 6;
    } else {
      std::ostringstream ssV;
      ssV << "Keyword '" << strV << "' not known.";
      throw FlxException("RBRV_entry_read_logn::RBRV_entry_read_logn_1", ssV.str() );
    }
    
    if (reader->whatIsNextChar()==',') {
      reader->getChar(',');
      epsilon = read_opt_para("epsilon");
      bool readN = (epsilon==NULL);
      if (readN==false) {
        if (reader->whatIsNextChar()==',') {
          reader->getChar(',');
          readN=true;
        }
      }
      if (readN) read_eval_once();
    }
    if (epsilon==NULL) epsilon = new FlxFunction(new FunNumber(ZERO));
    
  } catch (FlxException& e) {
    FLXMSG("RBRV_entry_read_logn::RBRV_entry_read_logn_2",1);
    if (p1) delete p1; 
    if (p2) delete p2; 
    if (p3) delete p3;
    if (p4) delete p4;
    if (epsilon) delete epsilon;
    throw;
  }
}

RBRV_entry_read_logn::~RBRV_entry_read_logn()
{
  if (p1) delete p1;
  if (p2) delete p2;
  if (p3) delete p3;
  if (p4) delete p4;
  if (epsilon) delete epsilon;
}

RBRV_entry* RBRV_entry_read_logn::generate_entry(const std::string& family, tuint& running_iID)
{
  const std::string name = family + nameF->eval_word(true,false,true);
  return new RBRV_entry_RV_lognormal(name,running_iID++,pid,new FlxFunction(*p1),new FlxFunction(*p2),p3?(new FlxFunction(*p3)):NULL,p4?(new FlxFunction(*p4)):NULL,new FlxFunction(*epsilon),eval_once);
}

RBRV_entry_read_uniform::RBRV_entry_read_uniform(const bool readName, const bool readBrakets)
: RBRV_entry_read_base(readName,readBrakets), a(NULL), b(NULL)
{
  try {
    reader->getChar('a');
    reader->getChar('=');
    a = new FlxFunction(funReader);
    reader->getChar(',');
    reader->getChar('b');
    reader->getChar('=');
    b = new FlxFunction(funReader);
    
    if (reader->whatIsNextChar()==',') {
      reader->getChar(',');
      read_eval_once();
    }
    
  } catch (FlxException& e) {
    FLXMSG("RBRV_entry_read_uniform::RBRV_entry_read_uniform",1);
    if (a) delete a; 
    if (b) delete b; 
    throw;
  }
}

RBRV_entry_read_uniform::~RBRV_entry_read_uniform()
{
  if (a) delete a;
  if (b) delete b;
}

RBRV_entry* RBRV_entry_read_uniform::generate_entry(const std::string& family, tuint& running_iID)
{
  const std::string name = family + nameF->eval_word(true,false,true);
  return new RBRV_entry_RV_uniform(name,running_iID++,new FlxFunction(*a),new FlxFunction(*b),eval_once);
}


RBRV_entry_read_Gumbel::RBRV_entry_read_Gumbel(const bool readName, const bool readBrakets)
: RBRV_entry_read_base(readName,readBrakets), methID(-1),p1(NULL),p2(NULL),p3(NULL),p4(NULL)
{
  try {
    const std::string strV = reader->getWord(true);
    if (strV == "u") {
      reader->getChar('=');
      p1 = new FlxFunction(funReader);
      reader->getChar(',');
      reader->getWord("alpha");
      reader->getChar('=');
      p2 = new FlxFunction(funReader);
      methID = 0;
    } else if (strV == "mu") {
      reader->getChar('=');
      p1 = new FlxFunction(funReader);
      reader->getChar(',');
      reader->getWord("sd");
      reader->getChar('=');
      p2 = new FlxFunction(funReader);
      methID = 1;
    } else if (strV == "pr") {
      reader->getChar('(');
      p1 = new FlxFunction(funReader);
      reader->getChar(')');
      reader->getChar('=');
      p2 = new FlxFunction(funReader);
      reader->getChar(',');
      reader->getWord("pr",false);
      reader->getChar('(');
      p3 = new FlxFunction(funReader);
      reader->getChar(')');
      reader->getChar('=');
      p4 = new FlxFunction(funReader);
      methID = 2;
    } else {
      std::ostringstream ssV;
      ssV << "Keyword '" << strV << "' not known.";
      throw FlxException("RBRV_entry_read_Gumbel::RBRV_entry_read_Gumbel_1", ssV.str() );
    }
    
    if (reader->whatIsNextChar()==',') {
      reader->getChar(',');
      read_eval_once();
    }
    
  } catch (FlxException& e) {
    FLXMSG("RBRV_entry_read_Gumbel::RBRV_entry_read_Gumbel_2",1);
    if (p1) delete p1; 
    if (p2) delete p2; 
    if (p3) delete p3;
    if (p4) delete p4;
    throw;
  }
}

RBRV_entry_read_Gumbel::~RBRV_entry_read_Gumbel()
{
  if (p1) delete p1;
  if (p2) delete p2;
}

RBRV_entry* RBRV_entry_read_Gumbel::generate_entry(const std::string& family, tuint& running_iID)
{
  const std::string name = family + nameF->eval_word(true,false,true);
  return new RBRV_entry_RV_Gumbel(name,running_iID++,methID,new FlxFunction(*p1),new FlxFunction(*p2),(p3==NULL)?NULL:(new FlxFunction(*p3)),(p4==NULL)?NULL:(new FlxFunction(*p4)),eval_once);
}

RBRV_entry_read_normal_trunc::RBRV_entry_read_normal_trunc(const bool readName, const bool readBrakets)
: RBRV_entry_read_base(readName,readBrakets), m(NULL), s(NULL), a(NULL), b(NULL)
{
  try {
    reader->getChar('m');
    reader->getChar('=');
    m = new FlxFunction(funReader);
    reader->getChar(',');
    reader->getChar('s');
    reader->getChar('=');
    s = new FlxFunction(funReader);
    bool stop = false;
    if (reader->whatIsNextChar()==',') {
      reader->getChar(',');
      if (reader->whatIsNextChar()=='a') {
        reader->getChar('a');
        reader->getChar('=');
        a = new FlxFunction(funReader);
        if (reader->whatIsNextChar()==',') {
          reader->getChar(',');
        } else {
          stop = true;
        }
      }
      if (reader->whatIsNextChar()=='b') {
        reader->getChar('b');
        reader->getChar('=');
        b = new FlxFunction(funReader);
        if (reader->whatIsNextChar()==',') {
          reader->getChar(',');
        } else {
          stop = true;
        }
      }
      if (!stop) {
        read_eval_once();
      }
    }
  } catch (FlxException& e) {
    FLXMSG("RBRV_entry_read_normal_trunc::RBRV_entry_read_normal_trunc",1);
    if (m) delete m;
  if (s) delete s;
  if (a) delete a;
  if (b) delete b;
    throw;
  }
}

RBRV_entry_read_normal_trunc::~RBRV_entry_read_normal_trunc()
{
  if (m) delete m;
  if (s) delete s;
  if (a) delete a;
  if (b) delete b;
}

RBRV_entry* RBRV_entry_read_normal_trunc::generate_entry(const std::string& family, tuint& running_iID)
{
  const std::string name = family + nameF->eval_word(true,false,true);
  return new RBRV_entry_RV_normal_trunc(name,running_iID++,new FlxFunction(*m),new FlxFunction(*s),a?(new FlxFunction(*a)):NULL,b?(new FlxFunction(*b)):NULL,eval_once);
}


RBRV_entry_read_beta::RBRV_entry_read_beta(const bool readName, const bool readBrakets)
: RBRV_entry_read_base(readName,readBrakets), is_mean(false),p1(NULL),p2(NULL),a(NULL),b(NULL)
{
  try {
    const std::string strV = reader->getWord(true);
    if (strV == "alpha") {
      reader->getChar('=');
      p1 = new FlxFunction(funReader);
      reader->getChar(',');
      reader->getWord("beta");
      reader->getChar('=');
      p2 = new FlxFunction(funReader);
      is_mean = false;
    } else if (strV == "mu") {
      reader->getChar('=');
      p1 = new FlxFunction(funReader);
      reader->getChar(',');
      reader->getWord("sd");
      reader->getChar('=');
      p2 = new FlxFunction(funReader);
      is_mean = true;
    } else {
      std::ostringstream ssV;
      ssV << "Keyword '" << strV << "' not known.";
      throw FlxException("RBRV_entry_read_beta::RBRV_entry_read_beta_1", ssV.str() );
    }
    
    if (reader->whatIsNextChar()==',') {
      reader->getChar(',');
      bool readN = true;
      if (reader->whatIsNextChar()=='a') {
        readN = false;
        reader->getChar('a');
        reader->getChar('=');
        a = new FlxFunction(funReader);
        reader->getChar(',');
        reader->getChar('b');
        reader->getChar('=');
        b = new FlxFunction(funReader);
        if (reader->whatIsNextChar()==',') {
          reader->getChar(',');
          readN=true;
        }
      }
      if (readN) read_eval_once();
    }
    if (a==NULL) a = new FlxFunction(new FunNumber(ZERO));
    if (b==NULL) b = new FlxFunction(new FunNumber(ONE));
    
  } catch (FlxException& e) {
    FLXMSG("RBRV_entry_read_beta::RBRV_entry_read_beta_2",1);
    if (p1) delete p1; 
    if (p2) delete p2; 
    if (a) delete a;
    if (b) delete b;
    throw;
  }
}

RBRV_entry_read_beta::~RBRV_entry_read_beta()
{
  if (p1) delete p1;
  if (p2) delete p2;
  if (a) delete a;
  if (b) delete b;
}

RBRV_entry* RBRV_entry_read_beta::generate_entry(const std::string& family, tuint& running_iID)
{
  const std::string name = family + nameF->eval_word(true,false,true);
  return new RBRV_entry_RV_beta(name,running_iID++,is_mean,new FlxFunction(*p1),new FlxFunction(*p2),new FlxFunction(*a),new FlxFunction(*b),eval_once);
}


RBRV_entry_read_exponential::RBRV_entry_read_exponential(const bool readName, const bool readBrakets)
: RBRV_entry_read_base(readName,readBrakets), lambda(NULL), epsilon(NULL)
{
  try {
    reader->getWord("lambda");
    reader->getChar('=');
    lambda = new FlxFunction(funReader);
    
    if (reader->whatIsNextChar()==',') {
      reader->getChar(',');
      reader->getWord("epsilon");
      reader->getChar('=');
      epsilon = new FlxFunction(funReader);
    } else {
      epsilon = new FlxFunction(new FunNumber(ZERO));
    }
    
    if (reader->whatIsNextChar()==',') {
      reader->getChar(',');
      read_eval_once();
    }
    
  } catch (FlxException& e) {
    FLXMSG("RBRV_entry_read_exponential::RBRV_entry_read_exponential",1);
    if (lambda) delete lambda; 
    throw;
  }
}

RBRV_entry_read_exponential::~RBRV_entry_read_exponential()
{
  if (lambda) delete lambda;
  if (epsilon) delete epsilon;
}

RBRV_entry* RBRV_entry_read_exponential::generate_entry(const std::string& family, tuint& running_iID)
{
  const std::string name = family + nameF->eval_word(true,false,true);
  if (eval_once) {
    const tdouble lambdaV = lambda->cast2positive();
    const tdouble eps = epsilon->calc();
    return new RBRV_entry_RV_exponential(name,running_iID++,new FlxFunction(new FunNumber(lambdaV)),new FlxFunction(new FunNumber(eps)));
  } else {
    return new RBRV_entry_RV_exponential(name,running_iID++,new FlxFunction(*lambda),new FlxFunction(*epsilon));
  }
}


RBRV_entry_read_gamma::RBRV_entry_read_gamma(const bool readName, const bool readBrakets)
: RBRV_entry_read_base(readName,readBrakets), is_mean(false), p1(NULL), p2(NULL), epsilon(NULL)
{
    try {
    const std::string strV = reader->getWord(true);
    if (strV == "k") {
      reader->getChar('=');
      p1 = new FlxFunction(funReader);
      reader->getChar(',');
      reader->getWord("lambda");
      reader->getChar('=');
      p2 = new FlxFunction(funReader);
      is_mean = false;
    } else if (strV == "mu") {
      reader->getChar('=');
      p1 = new FlxFunction(funReader);
      reader->getChar(',');
      reader->getWord("sd");
      reader->getChar('=');
      p2 = new FlxFunction(funReader);
      is_mean = true;
    } else {
      std::ostringstream ssV;
      ssV << "Keyword '" << strV << "' not known.";
      throw FlxException("RBRV_entry_read_gamma::RBRV_entry_read_gamma_1", ssV.str() );
    }
    
    if (reader->whatIsNextChar()==',') {
      reader->getChar(',');
      epsilon = read_opt_para("epsilon");
      bool readN = (epsilon==NULL);
      if (readN==false) {
        if (reader->whatIsNextChar()==',') {
          reader->getChar(',');
          readN=true;
        }
      }
      if (readN) read_eval_once();
    }
    if (epsilon==NULL) epsilon = new FlxFunction(new FunNumber(ZERO));
    
  } catch (FlxException& e) {
    FLXMSG("RBRV_entry_read_gamma::RBRV_entry_read_gamma_2",1);
    if (p1) delete p1; 
    if (p2) delete p2; 
    throw;
  }
}

RBRV_entry_read_gamma::~RBRV_entry_read_gamma()
{
  if (p1) delete p1;
  if (p2) delete p2;
  if (epsilon) delete epsilon;
}

RBRV_entry* RBRV_entry_read_gamma::generate_entry(const std::string& family, tuint& running_iID)
{
  const std::string name = family + nameF->eval_word(true,false,true);
  return new RBRV_entry_RV_gamma(name,running_iID++,is_mean, new FlxFunction(*p1), new FlxFunction(*p2), new FlxFunction(*epsilon),eval_once);
}


RBRV_entry_read_Poisson::RBRV_entry_read_Poisson(const bool readName, const bool readBrakets)
: RBRV_entry_read_base(readName,readBrakets), mean(NULL)
{
  try {
    reader->getWord("mu");
    reader->getChar('=');
    mean = new FlxFunction(funReader);
    
    if (reader->whatIsNextChar()==',') {
      reader->getChar(',');
      read_eval_once();
    }
    
  } catch (FlxException& e) {
    FLXMSG("RBRV_entry_read_Poisson::RBRV_entry_read_Poisson",1);
    if (mean) delete mean; 
    throw;
  }
}

RBRV_entry_read_Poisson::~RBRV_entry_read_Poisson()
{
  if (mean) delete mean;
}

RBRV_entry* RBRV_entry_read_Poisson::generate_entry(const std::string& family, tuint& running_iID)
{
  const std::string name = family + nameF->eval_word(true,false,true);
  if (eval_once) {
    const tdouble meanV = mean->cast2positive();
    return new RBRV_entry_RV_Poisson(name,running_iID++,new FlxFunction(new FunNumber(meanV)));
  } else {
    return new RBRV_entry_RV_Poisson(name,running_iID++,new FlxFunction(*mean));
  }
}


RBRV_entry_read_Binomial::RBRV_entry_read_Binomial(const bool readName, const bool readBrakets)
: RBRV_entry_read_base(readName,readBrakets), p(NULL), N(NULL)
{
  try {
    reader->getWord("p");
      reader->getChar('=');
      p = new FlxFunction(funReader);
      reader->getChar(',');
    reader->getWord("n");
      reader->getChar('=');
      N = new FlxFunction(funReader);
      
    if (reader->whatIsNextChar()==',') {
      reader->getChar(',');
      read_eval_once();
    }
    
  } catch (FlxException& e) {
    FLXMSG("RBRV_entry_read_Binomial::RBRV_entry_read_Binomial",1);
    if (p) delete p; 
    if (N) delete N;
    throw;
  }
}

RBRV_entry_read_Binomial::~RBRV_entry_read_Binomial()
{
  if (p) delete p;
  if (N) delete N;
}

RBRV_entry* RBRV_entry_read_Binomial::generate_entry(const std::string& family, tuint& running_iID)
{
  const std::string name = family + nameF->eval_word(true,false,true);
  return new RBRV_entry_RV_Binomial(name,running_iID++,new FlxFunction(*p),new FlxFunction(*N),eval_once);
}


RBRV_entry_read_Cauchy::RBRV_entry_read_Cauchy(const bool readName, const bool readBrakets)
: RBRV_entry_read_base(readName, readBrakets), loc(NULL), scale(NULL)
{
  try {
    reader->getChar('l');
    reader->getChar('=');
    loc = new FlxFunction(funReader);
    reader->getChar(',');
    reader->getChar('s');
    reader->getChar('=');
    scale = new FlxFunction(funReader);
    
    if (reader->whatIsNextChar()==',') {
      reader->getChar(',');
      read_eval_once();
    }
    
  } catch (FlxException& e) {
    FLXMSG("RBRV_entry_read_Cauchy::RBRV_entry_read_Cauchy",1);
    if (loc) delete loc; 
    if (scale) delete scale; 
    throw;
  }
}

RBRV_entry_read_Cauchy::~RBRV_entry_read_Cauchy()
{
  delete loc;
  delete scale;
}

RBRV_entry* RBRV_entry_read_Cauchy::generate_entry(const std::string& family, tuint& running_iID)
{
  const std::string name = family + nameF->eval_word(true,false,true);
  if (eval_once) {
    const tdouble locv = loc->calc();
    const tdouble scalev = scale->cast2positive();
    return new RBRV_entry_RV_Cauchy(name,running_iID++,new FlxFunction(new FunNumber(locv)),new FlxFunction(new FunNumber(scalev)));
  } else {
    return new RBRV_entry_RV_Cauchy(name,running_iID++,new FlxFunction(*loc),new FlxFunction(*scale));
  }
}

RBRV_entry_read_Weibull::RBRV_entry_read_Weibull(const bool readName, const bool readBrakets)
: RBRV_entry_read_base(readName,readBrakets), is_mean(false), p1(NULL), p2(NULL), epsilon(NULL)
{
  try {
    const std::string strV = reader->getWord(true);
    if (strV == "k") {
      reader->getChar('=');
      p1 = new FlxFunction(funReader);
      reader->getChar(',');
      reader->getWord("lambda");
      reader->getChar('=');
      p2 = new FlxFunction(funReader);
      is_mean = false;
    } else if (strV == "mu") {
      reader->getChar('=');
      p1 = new FlxFunction(funReader);
      reader->getChar(',');
      reader->getWord("sd");
      reader->getChar('=');
      p2 = new FlxFunction(funReader);
      is_mean = true;
    } else {
      std::ostringstream ssV;
      ssV << "Keyword '" << strV << "' not known.";
      throw FlxException("RBRV_entry_read_Weibull::RBRV_entry_read_Weibull_1", ssV.str() );
    }
    
    if (reader->whatIsNextChar()==',') {
      reader->getChar(',');
      epsilon = read_opt_para("epsilon");
      bool readN = (epsilon==NULL);
      if (readN==false) {
        if (reader->whatIsNextChar()==',') {
          reader->getChar(',');
          readN=true;
        }
      }
      if (readN) read_eval_once();
    }
    if (epsilon==NULL) epsilon = new FlxFunction(new FunNumber(ZERO));
    
  } catch (FlxException& e) {
    FLXMSG("RBRV_entry_read_Weibull::RBRV_entry_read_Weibull_2",1);
    if (p1) delete p1; 
    if (p2) delete p2; 
    throw;
  }
}

RBRV_entry_read_Weibull::~RBRV_entry_read_Weibull()
{
  if (p1) delete p1;
  if (p2) delete p2;
  if (epsilon) delete epsilon;
}

RBRV_entry* RBRV_entry_read_Weibull::generate_entry(const std::string& family, tuint& running_iID)
{
  const std::string name = family + nameF->eval_word(true,false,true);
  return new RBRV_entry_RV_Weibull(name,running_iID++,is_mean, new FlxFunction(*p1), new FlxFunction(*p2), new FlxFunction(*epsilon),eval_once);
}

RBRV_entry_read_ChiSquared::RBRV_entry_read_ChiSquared(const bool isSquared, const bool readName, const bool readBrakets)
: RBRV_entry_read_base(readName, readBrakets), isSquared(isSquared), p1(NULL)
{
  reader->getWord("dof",false);
  reader->getChar('=',false);
  p1 = new FlxFunction(funReader);
  if (reader->whatIsNextChar()==',') {
    reader->getChar(',',false);
    read_eval_once();
  }
}

RBRV_entry_read_ChiSquared::~RBRV_entry_read_ChiSquared()
{
  if (p1) delete p1;
}

RBRV_entry* RBRV_entry_read_ChiSquared::generate_entry(const std::string& family, tuint& running_iID)
{
  const std::string name = family + nameF->eval_word(true,false,true);
  if (isSquared) {
    return new RBRV_entry_RV_ChiSquared(name,running_iID++, new FlxFunction(*p1),eval_once);
  } else {
    return new RBRV_entry_RV_Chi(name,running_iID++, new FlxFunction(*p1),eval_once);
  }
}

RBRV_entry_read_StudentsT::RBRV_entry_read_StudentsT(const bool readName, const bool readBrakets)
: RBRV_entry_read_base(readName, readBrakets), p1(NULL)
{
  reader->getWord("dof",false);
  reader->getChar('=',false);
  p1 = new FlxFunction(funReader);
  if (reader->whatIsNextChar()==',') {
    reader->getChar(',',false);
    read_eval_once();
  }
}

RBRV_entry_read_StudentsT::~RBRV_entry_read_StudentsT()
{
  if (p1) delete p1;
}

RBRV_entry* RBRV_entry_read_StudentsT::generate_entry(const std::string& family, tuint& running_iID)
{
  const std::string name = family + nameF->eval_word(true,false,true);
  return new RBRV_entry_RV_StudentsT(name,running_iID++, new FlxFunction(*p1),eval_once);
}

RBRV_entry_read_StudentsT_generalized::RBRV_entry_read_StudentsT_generalized(const bool readName, const bool readBrakets)
: RBRV_entry_read_base(readName, readBrakets), nu(NULL), loc(NULL), scale(NULL)
{
  try {
    reader->getWord("dof",false);
    reader->getChar('=',false);
    nu = new FlxFunction(funReader);
    reader->getChar(',',false);
    reader->getWord("loc",false);
    reader->getChar('=',false);
    loc = new FlxFunction(funReader);
    reader->getChar(',',false);
    reader->getWord("scale",false);
    reader->getChar('=',false);
    scale = new FlxFunction(funReader);
  } catch (FlxException& e) {
    if (nu) delete nu;
    if (loc) delete loc;
    if (scale) delete scale;
    throw;
  }
}

RBRV_entry_read_StudentsT_generalized::~RBRV_entry_read_StudentsT_generalized()
{
    if (nu) delete nu;
    if (loc) delete loc;
    if (scale) delete scale;
}

RBRV_entry* RBRV_entry_read_StudentsT_generalized::generate_entry(const std::string& family, tuint& running_iID)
{
  const std::string name = family + nameF->eval_word(true,false,true);
  return new RBRV_entry_RV_StudentsT_generalized(name,running_iID++, new FlxFunction(*nu), new FlxFunction(*loc), new FlxFunction(*scale));
}

RBRV_entry_read_Laplace::RBRV_entry_read_Laplace(const bool readName, const bool readBrakets)
: RBRV_entry_read_base(readName, readBrakets), loc(NULL), scale(NULL)
{
  try {
    reader->getWord("loc",false);
    reader->getChar('=',false);
    loc = new FlxFunction(funReader);
    reader->getChar(',',false);
    reader->getWord("scale",false);
    reader->getChar('=',false);
    scale = new FlxFunction(funReader);
  } catch (FlxException& e) {
    if (loc) delete loc;
    if (scale) delete scale;
  }
}

RBRV_entry_read_Laplace::~RBRV_entry_read_Laplace()
{
    if (loc) delete loc;
    if (scale) delete scale;
}

RBRV_entry* RBRV_entry_read_Laplace::generate_entry(const std::string& family, tuint& running_iID)
{
  const std::string name = family + nameF->eval_word(true,false,true);
  return new RBRV_entry_RV_Laplace(name,running_iID++, new FlxFunction(*loc), new FlxFunction(*scale));
}

RBRV_entry_read_UserTransform::RBRV_entry_read_UserTransform(const bool readName, const bool readBrakets)
: RBRV_entry_read_base(readName, readBrakets), is_z2x(true), t1(NULL), t2(NULL), dh(NULL), checkXf(NULL), rv_z(NULL)
{
  FunReadPara::set_NumbOfPara(1);
  try {
    const std::string strV = reader->getWord(true);
    if (strV == "z2x") {
      reader->getChar('=');
      t1 = new FlxFunction(funReader);
      reader->getChar(',');
      if (reader->whatIsNextString(3,true)=="x2z") {
        reader->getWord("x2z");
        reader->getChar('=');
        t2 = new FlxFunction(funReader);
        reader->getChar(',');
      }
      if (reader->whatIsNextString(4,true)=="dx2z") {
        reader->getWord("dx2z");
        reader->getChar('=');
        dh = new FlxFunction(funReader);
        reader->getChar(',');
      }
      if (reader->whatIsNextString(6,true)=="checkx") {
        reader->getWord("checkx");
        reader->getChar('=');
        checkXf = new FlxFunction(funReader);
        reader->getChar(',');
      }
      is_z2x = true;
    } else if (strV == "y2z") {
      reader->getChar('=');
      t1 = new FlxFunction(funReader);
      reader->getChar(',');
      if (reader->whatIsNextString(3,true)=="z2y") {
        reader->getWord("z2y");
        reader->getChar('=');
        t2 = new FlxFunction(funReader);
        reader->getChar(',');
      }
      if (reader->whatIsNextString(4,true)=="dz2y") {
        reader->getWord("dz2y");
        reader->getChar('=');
        dh = new FlxFunction(funReader);
        reader->getChar(',');
      }
      is_z2x = false;
    } else {
      std::ostringstream ssV;
      ssV << "Keyword '" << strV << "' not known.";
      throw FlxException("RBRV_entry_read_UserTransform::RBRV_entry_read_UserTransform_01", ssV.str() );
    }
    FunReadPara::set_NumbOfPara(0);
    rv_z = RBRV_entry_read_base::read_entry(false,false);
  } catch (FlxException& e) {
    FLXMSG("RBRV_entry_read_UserTransform::RBRV_entry_read_UserTransform_02",1);
    FunReadPara::set_NumbOfPara(0);
    if (t1) delete t1; 
    if (t2) delete t2; 
    if (dh) delete dh;
    if (checkXf) delete checkXf;
    if (rv_z) delete rv_z;
    throw;
  }
}

RBRV_entry_read_UserTransform::~RBRV_entry_read_UserTransform()
{
  if (t1) delete t1; 
  if (t2) delete t2; 
  if (dh) delete dh;
  if (checkXf) delete checkXf;
  if (rv_z) delete rv_z;
}

RBRV_entry* RBRV_entry_read_UserTransform::generate_entry(const std::string& family, tuint& running_iID)
{
  const std::string name = family + nameF->eval_word(true,false,true);
  RBRV_entry_RV_base* rv_ze = rv_z->generate_entry_rv();
  return new RBRV_entry_RV_UserTransform(name,running_iID++, is_z2x, new FlxFunction(*t1), t2?(new FlxFunction(*t2)):NULL, dh?(new FlxFunction(*dh)):NULL, checkXf?(new FlxFunction(*checkXf)):NULL, rv_ze );
}

RBRV_entry_read_Truncated::RBRV_entry_read_Truncated(const bool readName, const bool readBrakets)
: RBRV_entry_read_base(readName, readBrakets), a(NULL), b(NULL), rv_z(NULL)
{
  try {
    if (reader->whatIsNextString(5,true)=="lower") {
      reader->getWord("lower");
      reader->getChar('=');
      a = new FlxFunction(funReader);
      reader->getChar(',');
    }
    if (reader->whatIsNextString(5,true)=="upper") {
      reader->getWord("upper");
      reader->getChar('=');
      b = new FlxFunction(funReader);
      reader->getChar(',');
    }
    rv_z = RBRV_entry_read_base::read_entry(false,false);
  } catch (FlxException& e) {
    FLXMSG("RBRV_entry_read_Truncated::RBRV_entry_read_Truncated",1);
    FunReadPara::set_NumbOfPara(0);
    if (a) delete a;
    if (b) delete b;
    if (rv_z) delete rv_z;
    throw;
  }
}

RBRV_entry_read_Truncated::~RBRV_entry_read_Truncated()
{
  if (a) delete a;
  if (b) delete b;
  if (rv_z) delete rv_z;
}

RBRV_entry* RBRV_entry_read_Truncated::generate_entry(const std::string& family, tuint& running_iID)
{
  const std::string name = family + nameF->eval_word(true,false,true);
  RBRV_entry_RV_base* rv_ze = rv_z->generate_entry_rv();
  return new RBRV_entry_RV_Truncated(name,running_iID++, a?(new FlxFunction(*a)):NULL, b?(new FlxFunction(*b)):NULL, rv_ze );
}

RBRV_entry_read_maxminTransform::RBRV_entry_read_maxminTransform(const bool readName, const bool readBrakets)
: RBRV_entry_read_base(readName, readBrakets), is_max(false), n(NULL), rv_z(NULL)
{
  try {
    const std::string strV = reader->getWord(true);
    if (strV == "min") {
      is_max = false;
    } else if (strV == "max" ) {
      is_max = true;
    } else {
      std::ostringstream ssV;
      ssV << "Keyword '" << strV << "' not known.";
      throw FlxException("RBRV_entry_read_maxminTransform::RBRV_entry_read_maxminTransform_01", ssV.str() );
    }
    reader->getChar('=');
    n = new FlxFunction(funReader);
    reader->getChar(',');
    rv_z = RBRV_entry_read_base::read_entry(false,false);
  } catch (FlxException& e) {
    FLXMSG("RBRV_entry_read_maxminTransform::~RBRV_entry_read_maxminTransform_02",1);
    if (n) delete n; 
    if (rv_z) delete rv_z;
    throw;
  }
}

RBRV_entry_read_maxminTransform::~RBRV_entry_read_maxminTransform()
{
  if (n) delete n;
  if (rv_z) delete rv_z;
}

RBRV_entry* RBRV_entry_read_maxminTransform::generate_entry(const std::string& family, tuint& running_iID)
{
  const std::string name = family + nameF->eval_word(true,false,true);
  RBRV_entry_RV_base* rv_ze = rv_z->generate_entry_rv();
  return new RBRV_entry_RV_maxminTransform(name,running_iID++, is_max, new FlxFunction(*n), rv_ze );
}







