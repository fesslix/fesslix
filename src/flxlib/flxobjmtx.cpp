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

#include "flxobjmtx.h"
#include "flxmtxfun_fun.h"


void FlxCreateObjReaders_Mtx::createObjReaders(FlxObjectReadBox* objReadBox) {
  objReadBox->insert("mtxcoeff", new FlxObjReadMtxCoeff());  
  objReadBox->insert("mtxconst_print", new FlxObjReadMtxCalc());
  objReadBox->insert("mtxconst_2octave", new FlxObjReadTransformMtx2Octave());
  objReadBox->insert("mtxconst_free", new FlxObjReadMtxConst_free());
  objReadBox->insert("mtxconst_new", new FlxObjReadMtxConstNew());
  objReadBox->insert("mtxconst_op", new FlxObjReadMtxConstOP());
  objReadBox->insert("mtxconst_sub", new FlxObjReadMtxConstSub());
  objReadBox->insert("mtxconst_mult", new FlxObjReadMtxConstMult());
  objReadBox->insert("mtxconst_fromfile", new FlxObjReadMtxConstFromFile());
  objReadBox->insert("mtxconst_transpose", new FlxObjReadMtxConstTranspose());
}

void FlxCreateObjReaders_Mtx::createFunReaders(FlxData* dataBox)
{
  flxmtxfun_fun_insert(dataBox->FunBox);
}


// -------------------------- OBJECTS ----------------------------------------------------------------------


void FlxObjMtxCoeff::task()
{
  const std::string& cname_str = cname->eval();
  FlxSMtx& mtx = *(data->ConstMtxBox.get(cname_str,true));
  tuint i = i_fun->cast2tuintW0(false);
  tuint j = j_fun->cast2tuintW0(false);
  tdouble d = c_fun->calc();
  if (i>=mtx.get_nrows() || j>=mtx.get_ncols()) {
    std::ostringstream ssV;
    ssV << "Index of coefficient (" << i << "," << j << ") are not within the matrix '" << cname_str << "'.";
    throw FlxException_NeglectInInteractive("FlxObjMtxCoeff::task", ssV.str() );
  }
  mtx.insert(i,j,d);
}


void FlxObjMtxCalc::task()
{
  const std::string mn = fun->eval();
  if (!only_coefs) {
    sout() << mn << " = " << std::endl;
    sout() << "{";
  }
  FlxSMtx* mp = data->ConstMtxBox.get(mn,true);
  sout() << *mp;
  if (!only_coefs) {
    sout() << " }" << "(" << mp->get_nrows() << "," << mp->get_ncols() << ")";
  }
  sout() << std::endl;
}

void FlxObjTransformMtx2Octave::task()
{
  const std::string mtxName = fun->eval();
  FlxSMtx* mp = data->ConstMtxBox.get(mtxName);
  sout() << "[";
  const tuint Nrows = mp->get_nrows();
  const tuint Ncols = mp->get_ncols();
  for (tuint i=0;i<Nrows;++i) {
    if (i>0) sout() << ";";
    for (tuint j=0;j<Ncols;++j) {
      sout() << " " << GlobalVar.Double2String(mp->operator()(i,j));
    }
  }
  sout() << "]" << std::endl;
}

void FlxObjMtxConst_free::task()
{
  data->ConstMtxBox.freeC(funStr->eval());
}

FlxObjMtxConstSeq::~FlxObjMtxConstSeq()
{
  delete mcn;
  delete startF;
  delete funCond;
  delete funConst;
}

void FlxObjMtxConstSeq::task()
{
  const tdouble cv_prev = *cv;
  *cv = startF->calc();
  std::stack<tdouble,std::list<tdouble> > q;
  while ( funCond->calc() > ZERO ) {
    q.push(*cv);
    *cv = funConst->calc();
  }
  tdouble* t = data->ConstMtxBox.get_Vec(q.size(),mcn->eval());
  for (size_t i = q.size(); i > 0; --i) {
    t[i-1] = q.top();
    q.pop();
  }
  *cv = cv_prev;
}

FlxObjMtxConstNew::~FlxObjMtxConstNew()
{
  delete mcn;
  if (mtx_right) delete mtx_right;
  if (rows) delete rows;
  if (cols) delete cols;
  if (val) delete val;
}

void FlxObjMtxConstNew::task()
{
  const std::string& mln = mcn->eval();
  if (mtx_right) {
    if (rows) throw FlxException_Crude("FlxObjMtxConstNew::task_2"); 
    const std::string& mrn = mtx_right->eval();
    if (mln==mrn) {
      std::ostringstream ssV;
      ssV << "The left-hand side (" << mln << ") must be different from the right-hand side!";
      throw FlxException("FlxObjMtxConstNew::task", ssV.str() );
    }
    FlxSMtx* mrp = data->ConstMtxBox.get(mrn,true);
    FlxSMtx* mlp = data->ConstMtxBox.get(mln,mrp->get_nrows(),mrp->get_ncols(),false);
    *mlp = *mrp;
  } else {
    if (rows==NULL) throw FlxException_Crude("FlxObjMtxConstNew::task_3"); 
    const tuint nr = rows->cast2tuint(false);
    const tuint nc = cols?(cols->cast2tuint(false)):1;
    const tdouble v = val?(val->calc()):ZERO;
    FlxSMtx* mlp = data->ConstMtxBox.get(mln,nr,nc,false);
    *mlp = v;
  }
}

FlxObjMtxConstNewU::~FlxObjMtxConstNewU()
{
  delete mcn;
  for (size_t i = 0; i < vecV.size(); ++i) {
    delete vecV[i];
  }
}

void FlxObjMtxConstNewU::task()
{
  tdouble* dp = data->ConstMtxBox.get_Mtx(nrows,ncols,mcn->eval());
  for (size_t i = 0; i < vecV.size(); ++i) {
    dp[i] = vecV[i]->calc();
  }
}

FlxObjMtxConstOP::~FlxObjMtxConstOP()
{
  delete mcn;
  if (f) delete f;
  if (s) delete s;
}

void FlxObjMtxConstOP::task()
{
  tdouble dprev = ZERO;
  if (dp) dprev = *dp;
  try {
    if ( f==NULL&&s==NULL ) throw FlxException_Crude("FlxObjMtxConstOP::task_1");
    const std::string& str_m1 = mcn->eval();
    FlxSMtx* mtx = data->ConstMtxBox.get(str_m1,true);
    tdouble* mtxp = mtx->get_internalPtr(true);
    const size_t N = mtx->get_Ncoeff();
    flxVec mtxv(mtxp,N);
    FlxSMtx* mtx2 = NULL;
    if (s) {
      const tdouble d2 = f?(f->calc()):ONE;
      const std::string& str_m2 = s->eval();
      mtx2 = data->ConstMtxBox.get(str_m2,true);
      if (mtx->get_ncols()!=mtx2->get_ncols() || mtx->get_nrows()!=mtx2->get_nrows()) {
        std::ostringstream ssV;
        ssV << "The matrices '" << str_m1 << "' and '" << str_m2 << "' do not have the same size: "
            << "rows(" << mtx->get_nrows() << "&" << mtx2->get_nrows() << ") "
            << "columns(" << mtx->get_ncols() << "&" << mtx2->get_ncols() << ")";
        throw FlxException("FlxObjMtxConstOP::task_2","Matrices have not the same size.",ssV.str());
      }
      tdouble* mtx2p = mtx2->get_internalPtr(true);
      flxVec mtx2v(mtx2p,N);
      switch (c) {
        case '+':
          if (f) {
            mtxv.add(mtx2v,d2);
          } else {
            mtxv += mtx2v;
          }
          break;
        case '-':
          if (f) {
            mtxv.add(mtx2v,-d2);
          } else {
            mtxv -= mtx2v;
          }
          break;
        case '*':
          for (size_t i=0;i<N;++i) {
            mtxp[i] *= mtx2p[i];
          }
          break;
        case '/':
          for (size_t i=0;i<N;++i) {
            mtxp[i] /= mtx2p[i];
          }
          break;
        case '^':
          for (size_t i=0;i<N;++i) {
            mtxp[i] = pow(mtxp[i],mtx2p[i]);
          }
          break;
        case ':':
          mtxv = mtx2v;
          break;
        default:
          throw FlxException_NotImplemented("FlxObjMtxConstOP::task_6");
      }
    } else {
      const tdouble d2 = (c=='(')?ZERO:f->calc();
      switch (c) {
        case '+':
          mtxv += d2;
          break;
        case '-':
          mtxv -= d2;
          break;
        case '*':
          mtxv *= d2;
          break;
        case '/':
          mtxv /= d2;
          break;
        case '^':
          {
            for (size_t i=0;i<N;++i) {
              mtxp[i] = pow(mtxp[i],d2);
            }
          }
          break;
        case ':':
          mtxv = d2;
          break;
        case '(':
          {
            for (size_t i=0;i<N;++i) {
              *dp = mtxp[i];
              mtxp[i] = f->calc();
            }
          }
          break;
        default:
          throw FlxException_NotImplemented("FlxObjMtxConstOP::task_5");
      }
    }
  } catch (FlxException& e) {
    if (dp) *dp = dprev;
    throw;
  }
  if (dp) *dp = dprev;
}

FlxObjMtxConstSub::~FlxObjMtxConstSub()
{
  delete mcn_target;
  delete mcn_from;
  for (tuint i=0;i<hv.size();++i) {
    if (hv[i]) delete hv[i];
  }
}

void FlxObjMtxConstSub::task()
{
  const std::string& mcn_from_str = mcn_from->eval();
  const std::string& mcn_target_str = mcn_target->eval();
  if (mcn_from_str==mcn_target_str) {
    std::ostringstream ssV;
    ssV << "You have to specify two different names - and not twice '" << mcn_from_str << "'.";
    throw FlxException("FlxObjMtxConstSub::task_a1", ssV.str() );
  }
  const std::string& m_full = (extract?mcn_target_str:mcn_from_str);
  const std::string& m_sub = (extract?mcn_from_str:mcn_target_str);
  FlxSMtx* mtx = data->ConstMtxBox.get(m_sub,true);
  switch (methID) {
    case cols:
    {
      tuint Nr = mtx->get_nrows();
      tuint Nc_mtx = mtx->get_ncols();
      tuint Nc = hv.size();
      // evaluate column numbers
        std::vector<tuint> hv_(Nc);
        for (tuint i=0;i<Nc;++i) {
          hv_[i] = hv[i]->cast2tuint(false);
        }
      // make sure that all given columns are valid
        for (tuint i=0;i<Nc;++i) {
          if (hv_[i]>Nc_mtx) {
            std::ostringstream ssV;
            ssV << "Matrix (" << m_sub << ") has only " << Nc_mtx << " columns. ";
            ssV << "You referred to column " << hv_[i] << " of this matrix - which is not valid.";
            throw FlxException("FlxObjMtxConstSub::task_b1", ssV.str() );
          }
        }
      tdouble* tp = data->ConstMtxBox.get_Mtx(m_full,Nr,Nc,!extract);        // target pointer
      tdouble* mp = mtx->get_internalPtr(false);        // from pointer
      if (mp) {
        if (extract) {
          for (tuint i=0;i<Nc;++i) {        // loop over selected columns
            for (tuint j=0;j<Nr;++j) {        // loop over rows of matrix
              tp[j*Nc+i] = mp[j*Nc_mtx+(hv_[i]-1)];
            }
          }
        } else {
          for (tuint i=0;i<Nc;++i) {        // loop over selected columns
            for (tuint j=0;j<Nr;++j) {        // loop over rows of matrix
              mp[j*Nc_mtx+(hv_[i]-1)] = tp[j*Nc+i];
            }
          }
        }
      } else {        // try a less efficient approach
        if (extract) {
          for (tuint i=0;i<Nc;++i) {        // loop over selected columns
            for (tuint j=0;j<Nr;++j) {        // loop over rows of matrix
              tp[j*Nc+i] = mtx->operator()(j,hv_[i]-1);
            }
          }
        } else {
          for (tuint i=0;i<Nc;++i) {        // loop over selected columns
            for (tuint j=0;j<Nr;++j) {        // loop over rows of matrix
              mtx->insert(j,hv_[i]-1,tp[j*Nc+i]);
            }
          }
        }
      }
      break;
    }
    case rows:
    {
      tuint Nr_mtx = mtx->get_nrows();
      tuint Nc = mtx->get_ncols();
      tuint Nr = hv.size();
      // evaluate row numbers
        std::vector<tuint> hv_(Nr);
        for (tuint i=0;i<Nr;++i) {
          hv_[i] = hv[i]->cast2tuint(false);
        }
      // make sure that all given rows are valid
        for (tuint i=0;i<Nr;++i) {
          if (hv_[i]>Nr_mtx) {
            std::ostringstream ssV;
            ssV << "Matrix (" << m_sub << ") has only " << Nr_mtx << " rows. ";
            ssV << "You referred to row " << hv_[i] << " of this matrix - which is not valid.";
            throw FlxException("FlxObjMtxConstSub::task_c1", ssV.str() );
          }
        }
      tdouble* tp = data->ConstMtxBox.get_Mtx(m_full,Nr,Nc,!extract);        // target pointer
      tdouble* mp = mtx->get_internalPtr(false);        // from pointer
      if (mp) {
        if (extract) {
          for (tuint i=0;i<Nr;++i) {        // loop over selected rows
            for (tuint j=0;j<Nc;++j) {        // loop over columns of matrix
              tp[i*Nc+j] = mp[(hv_[i]-1)*Nc+j];
            }
          }
        } else {
          for (tuint i=0;i<Nr;++i) {        // loop over selected rows
            for (tuint j=0;j<Nc;++j) {        // loop over columns of matrix
              mp[(hv_[i]-1)*Nc+j] = tp[i*Nc+j];
            }
          }
        }
      } else {        // try a less efficient approach
        if (extract) {
          for (tuint i=0;i<Nc;++i) {        // loop over selected rows
            for (tuint j=0;j<Nr;++j) {        // loop over columns of matrix
              tp[i*Nc+j] = mtx->operator()((hv_[i]-1),j);
            }
          }
        } else {
          for (tuint i=0;i<Nc;++i) {        // loop over selected rows
            for (tuint j=0;j<Nr;++j) {        // loop over columns of matrix
              mtx->insert(hv_[i]-1,j,tp[i*Nc+j]);
            }
          }
        }
      }
      break;
    }
    case seq:
    {
      tuint Nr_mtx = mtx->get_nrows();
      tuint Nc_mtx = mtx->get_ncols();
      const tuint rs = (hv[0]==NULL)?1:(hv[0]->cast2tuint(false));
      const tuint re = (hv[1]==NULL)?Nr_mtx:(hv[1]->cast2tuint(false));
      const tuint cs = (hv[2]==NULL)?1:(hv[2]->cast2tuint(false));
      const tuint ce = (hv[3]==NULL)?Nc_mtx:(hv[3]->cast2tuint(false));
      if (rs==0||cs==0||rs>re||cs>ce||re>Nr_mtx||ce>Nc_mtx) {
        std::ostringstream ssV;
          ssV << "Inconsistency observed.";
          throw FlxException("FlxObjMtxConstSub::task_c1", ssV.str() );
      }
      tuint Nr = (1+re)-rs;
      tuint Nc = (1+ce)-cs;
      tdouble* tp = data->ConstMtxBox.get_Mtx(m_full,Nr,Nc,!extract);        // target pointer
      tdouble* mp = mtx->get_internalPtr(false);        // from pointer
      if (mp) {
        if (extract) {
          for (tuint j=0;j<Nr;++j) {        // loop over rows of matrix
            for (tuint i=0;i<Nc;++i) {        // loop over selected columns
              tp[j*Nc+i] = mp[(j+(rs-1))*Nc_mtx+i+(cs-1)];
            }
          }
        } else {
          for (tuint i=0;i<Nc;++i) {        // loop over selected columns
            for (tuint j=0;j<Nr;++j) {        // loop over rows of matrix
              mp[(j+(rs-1))*Nc_mtx+i+(cs-1)] = tp[j*Nc+i];
            }
          }
        }
      } else {        // try a less efficient approach
        if (extract) {
          for (tuint j=0;j<Nr;++j) {        // loop over rows of matrix
            for (tuint i=0;i<Nc;++i) {        // loop over selected columns
              tp[j*Nc+i] = mtx->operator()(j+(rs-1),i+(cs-1));
            }
          }
        } else {
          for (tuint i=0;i<Nc;++i) {        // loop over selected columns
            for (tuint j=0;j<Nr;++j) {        // loop over rows of matrix
              mtx->insert(j+(rs-1),i+(cs-1),tp[j*Nc+i]);
            }
          }
        }
      }
      break;
    }
    default:
      throw FlxException_Crude("FlxObjMtxConstSub::task_z1");
  };
}

FlxObjMtxConstMult::~FlxObjMtxConstMult()
{
  delete mcn;
  delete mn1;
  delete mn2;
}

void FlxObjMtxConstMult::task()
{
  const std::string& mcn_str = mcn->eval();
  const std::string& mn1_str = mn1->eval();
  const std::string& mn2_str = mn2->eval();
  if (mcn_str==mn1_str || mcn_str==mn2_str) {
    std::ostringstream ssV;
    ssV << "The matrix on the left-hand side (" << mcn_str << ") must not appear on the right-hand side!";
    throw FlxException("FlxObjMtxConstNew::task", ssV.str() );
  }
  FlxSMtx* m1 = data->ConstMtxBox.get(mn1_str,true);
  FlxSMtx* m2 = data->ConstMtxBox.get(mn2_str,true);
  // check size of matrices
    if (m1->get_ncols()!=m2->get_nrows()) {
      std::ostringstream ssV;
      ssV << "Matrices can not be multiplied.";
      throw FlxException("FlxObjMtxConstMult::task_1", ssV.str() );
    }
  // request left-hand side
    FlxSMtx* mres = data->ConstMtxBox.get(mcn_str,m1->get_nrows(),m2->get_ncols(),false);
  mres->mult(*m1,*m2);
}

FlxObjMtxConstFromFile::~FlxObjMtxConstFromFile()
{
  delete mcn_target;
  delete cols;
  delete strV;
}

void FlxObjMtxConstFromFile::task()
{
  // read numbers from file
    const tuint c = cols->cast2tuint();
    FlxIstream& istrm = data->IstreamBox.get(strV->eval_word(true));
    std::vector<tdouble> iVec;
    tdouble d;
    while (istrm.get_value(d,true)) {
      iVec.push_back(d);
    }
    const size_t S = iVec.size();
    if (S%c!=0) {
      throw FlxException("FunMtxFromFile::calc","Size mismatch");
    }
  // assign numbers to matrix
    const tuint r = S/c;
    tdouble* dptr = data->ConstMtxBox.get(mcn_target->eval(),r,c,false)->get_internalPtr(true);
    for (size_t i=0;i<S;++i) {
      dptr[i] = iVec[i];
    }
}

void FlxObjMtxConstTranspose::task()
{
  const std::string target_name = mcn_target->eval();
  FlxSMtx* mtx_old = data->ConstMtxBox.get(target_name,true);
  const tuint Nrows = mtx_old->get_nrows();
  const tuint Ncols = mtx_old->get_ncols();
  FlxSMtx* mtx_new = new FlxSMtx(Ncols,Nrows,ZERO); 
  for (tuint i=0;i<Nrows;++i) {
    for (tuint j=0;j<Ncols;++j) {
      mtx_new->insert(j,i,mtx_old->operator()(i,j));
    }
  }
  data->ConstMtxBox.insert(target_name,mtx_new);
}


// -------------------------- READERS ----------------------------------------------------------------------

FlxObjBase* FlxObjReadMtxCoeff::read()
{
  FlxMtxConstFun* cname = new FlxMtxConstFun(false);
  FlxFunction* i_fun = NULL;
  FlxFunction* j_fun = NULL;
  FlxFunction* c_fun=NULL;
  try {
    reader->getChar('(',false);
    i_fun = new FlxFunction(funReader,false);
    j_fun=NULL;
    if (reader->whatIsNextChar()==',') {
      reader->getChar(',',false);
      j_fun = new FlxFunction(funReader,false);
    } else {
      j_fun = new FlxFunction(new FunNumber(ZERO));
    }
    reader->getChar(')',false);
    reader->getChar('=',false);
    c_fun = new FlxFunction(funReader,false);
    read_optionalPara(false);
    return new FlxObjMtxCoeff(get_doLog(),cname,i_fun,j_fun,c_fun);
  } catch (FlxException &e) {
    FLXMSG("FlxObjReadMtxCoeff::read",1);
    delete cname;
    if (i_fun) delete i_fun;
    if (j_fun) delete j_fun;
    if (c_fun) delete c_fun;
    throw;
  }
}

FlxObjReadMtxCalc::FlxObjReadMtxCalc(): FlxObjReadOutputBase()
{
  // only_coefs
    AllDefParaBox->insert(new FlxOptionalParaBool(false,"mtxconst_print::only_coefs"));
    ParaBox.insert("only_coefs", "mtxconst_print::only_coefs" );
}

FlxObjBase* FlxObjReadMtxCalc::read()
{
  FlxMtxConstFun* f1 = new FlxMtxConstFun(true);
  try {
    read_optionalPara(false);
    return new FlxObjMtxCalc(get_doLog(), f1, get_stream(), get_optPara_bool("only_coefs") );
  } catch ( FlxException &e) {
    FLXMSG("FlxObjReadMtxCalc::read",1);
    if (f1) delete f1;
    throw;
  }
}

FlxObjBase* FlxObjReadTransformMtx2Octave::read()
{
  FlxMtxConstFun* funStr = new FlxMtxConstFun(true);
  try {
    read_optionalPara(false);
    return new FlxObjTransformMtx2Octave(get_doLog(), funStr, get_stream() );
  } catch ( FlxException &e) {
    FLXMSG("FlxObjReadTransformMtx2Octave::read",1);
    delete funStr;
    throw;
  }
}

FlxObjBase* FlxObjReadMtxConst_free::read()
{
  FlxMtxConstFun* funStr = new FlxMtxConstFun(false);
  try {
    read_optionalPara(false);
    return new FlxObjMtxConst_free(get_doLog(), funStr );
  } catch ( FlxException &e) {
    FLXMSG("FlxObjReadMtxConst_free::read",1);
    if (funStr) delete funStr;
    throw;
  }
}

void FlxObjReadMtxConstNew::read_mtx(std::vector< FlxFunction* >& vecV, tuint& nrows, tuint& ncols)
{
  nrows = 1;
  ncols = 1;
  try {
    reader->getChar('{',false);
    // obtain number of columns
      FlxFunction* fp = new FlxFunction(funReader,false);
      vecV.push_back(fp);
      while ( reader->whatIsNextChar() == ',' ) {
        reader->getChar(',',false);
        ncols++;
        fp = new FlxFunction(funReader,false);
        vecV.push_back(fp);
      };
    // obtain number of rows
      while (reader->whatIsNextChar() == ';' ) {
        reader->getChar(';',false);
        nrows++;
        fp = new FlxFunction(funReader,false);
        vecV.push_back(fp);
        for (tuint i = 1; i < ncols; ++i) {
          reader->getChar(',',false);
          fp = new FlxFunction(funReader,false);
          vecV.push_back(fp);
        };
      }
    reader->getChar('}',false);
  } catch (FlxException &e) {
    for (tuint i = 0; i < vecV.size(); ++i) {
      delete vecV[i];
    }
    throw;
  }
}

void FlxObjReadMtxConstNew::read_mtx_Matlab(std::vector< FlxFunction* >& vecV, tuint& nrows, tuint& ncols)
{
  nrows = 1;
  ncols = 1;
  try {
    reader->getChar('[',false);
    // obtain number of columns
      FlxFunction* fp = new FlxFunction(funReader,false);
      vecV.push_back(fp);
      while ( true ) {
        // check 'extended' loop condition
          const char nc = reader->whatIsNextChar();
          if (nc == ',') {                // the comma-separator is optional in the Matlab syntax
            reader->getChar(',',false);
          } else if ( nc == ';' || nc == ']') {
            break;
          }
        ncols++;
        fp = new FlxFunction(funReader,false);
        vecV.push_back(fp);
      };
    // obtain number of rows
      while (reader->whatIsNextChar() == ';' ) {
        reader->getChar(';',false);
        nrows++;
        fp = new FlxFunction(funReader,false);
        vecV.push_back(fp);
        for (tuint i = 1; i < ncols; ++i) {
          if (reader->whatIsNextChar()==',') {
            reader->getChar(',',false);
          }
          fp = new FlxFunction(funReader,false);
          vecV.push_back(fp);
        };
      }
    reader->getChar(']',false);
  } catch (FlxException &e) {
    for (tuint i = 0; i < vecV.size(); ++i) {
      delete vecV[i];
    }
    throw;
  }
}

void FlxObjReadMtxConstNew::read_seq(tdoublePtr& cv, FlxFunctionPtr& startF, FlxFunctionPtr& funCondL, FlxFunctionPtr& funConst)
{
  cv = NULL; startF = NULL; funCondL = NULL; funConst = NULL;
  try {
    reader->getChar('(',false);
    if (reader->getNextType() != ReadStream::STRING) {
      std::ostringstream ssV;
      ssV << "Name of the 'seq' variable to use expected.";
      throw FlxException_NeglectInInteractive("FlxObjReadMtxConstNew::read_seq_1", ssV.str(), reader->getCurrentPos());
    }
    cv = data->ConstantBox.get(reader->getWord(true,false),true);
    reader->getChar('=',false);
    startF = new FlxFunction(funReader,false);
    reader->getChar(',',false);
    funCondL = new FlxFunction(funReader,false);
    reader->getChar(',',false);
    funConst = new FlxFunction(funReader,false);
    reader->getChar(')',false);
  } catch (FlxException &e) {
    FLXMSG("FlxObjReadMtxConstNew::read_seq_2",1);
    if (startF) delete startF;
    if (funCondL) delete funCondL;
    if (funConst) delete funConst;
    throw;
  }
}

FlxObjBase* FlxObjReadMtxConstNew::read()
{
  FlxMtxConstFun* mcn = new FlxMtxConstFun(false);
  try {
    if (reader->whatIsNextChar()=='=') {
      reader->getChar('=',false);
      if (reader->whatIsNextChar()=='{') {
        std::vector<FlxFunction*> vecV;
        tuint nrows, ncols;
        read_mtx(vecV,nrows,ncols);
        return new FlxObjMtxConstNewU(get_doLog(),mcn,vecV,nrows,ncols);
      } else if (reader->whatIsNextString(3,true)=="seq") {
        reader->getWord("seq");
        tdouble* cv; FlxFunction* startF; FlxFunction* funCondL; FlxFunction* funConst;
        read_seq(cv, startF, funCondL, funConst);
        return new FlxObjMtxConstSeq(get_doLog(),mcn,cv,startF,funCondL,funConst);
      } else if (reader->whatIsNextString(6,true)=="matlab") {
        reader->getWord("matlab");
        reader->getChar(':',false);
        std::vector<FlxFunction*> vecV;
        tuint nrows, ncols;
        read_mtx_Matlab(vecV,nrows,ncols);
        return new FlxObjMtxConstNewU(get_doLog(),mcn,vecV,nrows,ncols);
      } else {
        std::ostringstream ssV;
        ssV << "Unexpected character.";
        throw FlxException_NeglectInInteractive("FlxObjReadMtxConstNew::read_3", ssV.str(), reader->getCurrentPos());
      }
    } else {
      FlxMtxConstFun* mtx_right = NULL;
      FlxFunction* rows = NULL;
      FlxFunction* cols = NULL;
      FlxFunction* val = NULL;
      FlxObjMtxConstNew* obj = NULL;
      try {
        reader->getChar('(',false);
        if (reader->whatIsNextChar()=='{') {
          reader->getChar('{',false);
          mtx_right = new FlxMtxConstFun(false);
          reader->getChar('}',false);
        } else {
          rows = new FlxFunction(funReader,false);
          if (reader->whatIsNextChar()==',') {
            reader->getChar(',',false);
            cols = new FlxFunction(funReader,false);
            if (reader->whatIsNextChar()==',') {
              reader->getChar(',',false);
              val = new FlxFunction(funReader,false);
            }
          }
        }
        reader->getChar(')',false);
        read_optionalPara(false);
        return new FlxObjMtxConstNew(get_doLog(),mcn,mtx_right,rows,cols,val);
      } catch (FlxException& e) {
        if (mtx_right) delete mtx_right;
        if (rows) delete rows;
        if (cols) delete cols;
        if (val) delete val;
        if (obj) delete obj;
        throw;
      }
    }
  } catch (FlxException& e) {
      delete mcn;
      throw;
    }
}

FlxObjBase* FlxObjReadMtxConstOP::read()
{
  FlxMtxConstFun* mcn = new FlxMtxConstFun(false);
  FlxMtxConstFun *s = NULL;
  FlxFunction *f = NULL;
  tdouble* dp = NULL;
  try {
    const char c = reader->getChar();
    if (c=='(') {
      const std::string cn = reader->getWord(true,false);
      dp = data->ConstantBox.get(cn,true);
      reader->getChar(')',false);
    } else {
      if (c!='+'&&c!='-'&&c!='*'&&c!='/'&&c!='^'&&c!=':')  {
        std::ostringstream ssV;
        ssV << "Unknown operator '" << c << "'.";
        throw FlxException_NeglectInInteractive("FlxObjReadMtxConstOP::read_1", ssV.str(), reader->getCurrentPos());
      }
    }
    reader->getChar('=',false);
    if (reader->whatIsNextChar()=='{') {
      reader->getChar('{',false);
      s = new FlxMtxConstFun(true);
      reader->getChar('}',false);
      if (c=='+'||c=='-') {
        if (reader->whatIsNextChar()=='*') {
          reader->getChar();
          f = new FlxFunction(funReader,false);
        }
      };
    } else {
      f = new FlxFunction(funReader,false);
    }
    read_optionalPara(false);
    return new FlxObjMtxConstOP(get_doLog(),mcn,c,f,s,dp);
  } catch (FlxException& e) {
    delete mcn;
    if (f) delete f;
    if (s) delete s;
    throw;
  }
}

void FlxObjReadMtxConstSub::read_subInfo(FlxObjMtxConstSub::Meth& meth, std::vector< FlxFunction* >& hv)
{
  try {
    reader->getChar('(',false);
    // extract the method id
      const std::string meth_str = reader->getWord(true);
      if (meth_str=="col") {
        meth = FlxObjMtxConstSub::cols;
      } else if (meth_str=="row") {
        meth = FlxObjMtxConstSub::rows;
      } else if (meth_str=="seq") {
        meth = FlxObjMtxConstSub::seq;
      } else {
        std::ostringstream ssV;
        ssV << "Unknown method-ID '" << meth_str << "'.";
        throw FlxException_NeglectInInteractive("FlxObjReadMtxConstSub::read_1", ssV.str(), reader->getCurrentPos());
      }
    reader->getChar(':',false);
    if (meth==FlxObjMtxConstSub::seq) {
      if (reader->whatIsNextChar()==':') {
        hv.push_back(NULL);
      } else {
        hv.push_back(new FlxFunction(funReader,false));
      }
      reader->getChar(':');
      if (reader->whatIsNextChar()==','||reader->whatIsNextChar()==')') {
        hv.push_back(NULL);
      } else {
        hv.push_back(new FlxFunction(funReader,false));
      }
      if (reader->whatIsNextChar()==')') {
        hv.push_back(NULL);
        hv.push_back(NULL);
      } else {
        reader->getChar(',');
        if (reader->whatIsNextChar()==':') {
          hv.push_back(NULL);
        } else {
          hv.push_back(new FlxFunction(funReader,false));
        }
        reader->getChar(':');
        if (reader->whatIsNextChar()==')') {
          hv.push_back(NULL);
        } else {
          hv.push_back(new FlxFunction(funReader,false));
        }
      }
    } else {
      while (reader->whatIsNextChar()!=')') {
        if (hv.empty()==false) reader->getChar(',',false);
        hv.push_back(new FlxFunction(funReader,false));
      }
    }
    reader->getChar(')',false);
    if (hv.empty()) {
      std::ostringstream ssV;
      ssV << "Expected entry 'FlxFunction' and not ')'.";
      throw FlxException_NeglectInInteractive("FlxObjReadMtxConstSub::read_subInfo_2", ssV.str(), reader->getCurrentPos());
    }
  } catch (FlxException& e) {
    for (tuint i=0;i<hv.size();++i) {
      if (hv[i]) delete hv[i];
    }
    hv.clear();
    throw;
  }
}

FlxObjBase* FlxObjReadMtxConstSub::read()
{
  FlxMtxConstFun* mcn_target = new FlxMtxConstFun(false);
  FlxMtxConstFun *mcn_from = NULL;
  std::vector<FlxFunction*> hv;
  bool extract = true;
  try {
    FlxObjMtxConstSub::Meth meth;
    if (reader->whatIsNextChar()=='(') {
      extract = false;
      read_subInfo(meth,hv);
    }
    reader->getChar('=',false);
    mcn_from = new FlxMtxConstFun(false);
    if (extract) {
      read_subInfo(meth,hv);
    }
    read_optionalPara(false);
    return new FlxObjMtxConstSub(get_doLog(),mcn_target,mcn_from,meth,hv,extract);
  } catch (FlxException& e) {
    delete mcn_target;
    if (mcn_from) delete mcn_from;
    for (tuint i=0;i<hv.size();++i) {
      if (hv[i]) delete hv[i];
    }
    throw;
  }
}

FlxObjBase* FlxObjReadMtxConstMult::read()
{
  FlxMtxConstFun* mcn_target = new FlxMtxConstFun(false);
  FlxMtxConstFun* mn1 = NULL;
  FlxMtxConstFun* mn2 = NULL;
  try {
    reader->getChar('=');
    mn1 = new FlxMtxConstFun(false);
    reader->getChar('*');
    mn2 = new FlxMtxConstFun(false);
    read_optionalPara(false);
    return new FlxObjMtxConstMult(get_doLog(),mcn_target,mn1,mn2);
  } catch (FlxException& e) {
    delete mcn_target;
    if (mn1) delete mn1;
    if (mn2) delete mn2;
    throw;
  }
}

FlxObjBase* FlxObjReadMtxConstFromFile::read()
{
  FlxMtxConstFun* mcn_target = new FlxMtxConstFun(false);
  FlxFunction* colN = NULL;
  FlxString* istrm = NULL;
  try {
    reader->getChar('(');
    reader->getWord("col",false);
    reader->getChar('=');
    colN = new FlxFunction(funReader,false);
    reader->getChar(')');
    reader->getChar('=');
    istrm = new FlxString(false,false);
    read_optionalPara(false);
    return new FlxObjMtxConstFromFile(get_doLog(),mcn_target,colN,istrm);
  } catch (FlxException& e) {
    delete mcn_target;
    if (colN) delete colN;
    if (istrm) delete istrm;
    throw;
  }
}

FlxObjBase* FlxObjReadMtxConstTranspose::read()
{
  FlxMtxConstFun* mcn_target = new FlxMtxConstFun(true);
  try {
    read_optionalPara(false);
    return new FlxObjMtxConstTranspose(get_doLog(),mcn_target);
  } catch (FlxException& e) {
    delete mcn_target;
    throw;
  }
}

