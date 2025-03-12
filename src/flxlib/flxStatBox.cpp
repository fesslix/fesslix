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


#include "flxStatBox.h"

#include "flxMtx_Eigen.h"


FlxStatBox::FlxStatBox(tuint N, tuint M)
: N(N), M(M), Nc(0), tp(new tdouble[size_t(N)*size_t(M)])
{
  if (N==0 || M==0) {
    std::ostringstream ssV;
    ssV << "Zero size is not allowed.";
    throw FlxException("FlxStatBox::FlxStatBox_1",ssV.str());
  }
}

FlxStatBox::~FlxStatBox()
{
  delete [] tp;
}

void FlxStatBox::add(const flxVec& sv)
{
  if (sv.get_N()!=M) {
    std::ostringstream ssV;
    ssV << "Input vector has wrong dimension: " << sv.get_N() << " and not " << M << ".";
    throw FlxException("FlxStatBox::add_1",ssV.str());
  }
  if (Nc>=N) {
    std::ostringstream ssV;
    ssV << "No more elements can be added.";
    throw FlxException("FlxStatBox::add_2",ssV.str());
  }
  const tdouble* svp = sv.get_tmp_vptr_const();
  for (tuint i=0;i<M;++i) {
    tp[size_t(N)*i+size_t(Nc)] = svp[i];
  }
  ++Nc;
}

void FlxStatBox::compute_ExpSd(flxVec& eV, flxVec& sV)
{
  #if FLX_DEBUG
    if (eV.get_N()!=M) {
      throw FlxException_Crude("FlxStatBox::compute_ExpSd_1");
    }
    if (sV.get_N()!=M) {
      throw FlxException_Crude("FlxStatBox::compute_ExpSd_2");
    }
  #endif
  pdouble pp; tdouble tt;
  for (tuint i=0;i<M;++i) {
    pp = ZERO;
    size_t st = size_t(N)*i;
    for (tuint j=0;j<Nc;++j) {
        pp += tp[st+j];
    }
    tt = pp.cast2double()/Nc;
    pp = ZERO;
    for (tuint j=0;j<Nc;++j) {
        pp += pow2(tp[st+j]-tt);
    }
    eV[i] = tt;
    sV[i] = sqrt(pp.cast2double()/Nc);
  }
}

void FlxStatBox::compute_mean(flxVec& eV)
{
  #if FLX_DEBUG
    if (eV.get_N()!=M) {
      throw FlxException_Crude("FlxStatBox::compute_mean");
    }
  #endif
  pdouble pp; tdouble tt;
  for (tuint i=0;i<M;++i) {
    pp = ZERO;
    size_t st = size_t(N)*i;
    for (tuint j=0;j<Nc;++j) {
        pp += tp[st+j];
    }
    tt = pp.cast2double()/Nc;
    eV[i] = tt;
  }
}


void FlxStatBox::plot_confidence(std::ostream& ostrm, flxVec& clevels, const int fixW, flxVec* initCol)
{
  if (Nc==0) {
    std::ostringstream ssV;
    ssV << "No samples!";
    throw FlxException("FlxStatBox::plot_confidence_0",ssV.str());
  }
  if (initCol && initCol->get_N()!=M) {
    std::ostringstream ssV;
    ssV << "The user-defined vector does not have the correct dimension.";
    throw FlxException("FlxStatBox::plot_confidence_1",ssV.str());
  }
  // sort clevels
    tdouble* clp = clevels.get_tmp_vptr();
    const tuint clN = clevels.get_N();
    if (clN==0) {
      std::ostringstream ssV;
      ssV << "An empty vector is not allowed as input.";
      throw FlxException("FlxStatBox::plot_confidence_2",ssV.str());
    }
    std::sort(clp,clp+clN);
    if (clevels[0]<ZERO || clevels[clN-1]>ONE) {
      std::ostringstream ssV;
      ssV << "At least on coefficient of the vector is not a probability.";
      throw FlxException("FlxStatBox::plot_confidence_3",ssV.str());
    }
    const tdouble clh = (ONE)/(2*Nc);
  // plot headings
    {
    ostrm << '#';
    if (fixW<0) {
      if (initCol) ostrm << "User\t";
      ostrm << "Mean\tStd.Dev\t";
    } else {
      const tuint tsl = GlobalVar.Double2String(ONE,false,-1,fixW).length();
      // init col
        if (initCol) {
          std::string tstr = "User";
          if (tstr.length()>=tsl-1) {
            tstr = tstr.substr(0,tsl-1);
          } else {
            for (tuint i=tstr.length();i<tsl-1;++i) tstr += ' ';
          }
          ostrm << tstr;
          }
      // mean
        const tuint tsl2 = (initCol?tsl:tsl-1);
        std::string tstr = "Mean";
        if (tstr.length()>=tsl2) {
          tstr = tstr.substr(0,tsl2);
        } else {
          for (tuint i=tstr.length();i<tsl2;++i) tstr += ' ';
        }
        ostrm << tstr;
      // standard deviation
        tstr = "Std.Dev.";
        if (tstr.length()>=tsl) {
          tstr = tstr.substr(0,tsl);
        } else {
          for (tuint i=tstr.length();i<tsl;++i) tstr += ' ';
        }
        ostrm << tstr << ((fixW<0)?"\t":" ");
    }
    for (tuint i=0;i<clN;++i) {
      ostrm << GlobalVar.Double2String(clp[i],false,-1,fixW) << ((fixW<0)?"\t":" ");
    }
    ostrm << std::endl;
    }
  // get expectations and standard deviations
    flxVec eV(M);
    flxVec sV(M);
    compute_ExpSd(eV,sV);
  // do the actual printing
    flxVec sorVv(Nc);
    tdouble* sorVp = sorVv.get_tmp_vptr();
    for (tuint i=0;i<M;++i) {
      // user-defined column
        if (initCol) {
          ostrm << GlobalVar.Double2String(initCol->operator[](i),false,-1,fixW) << ((fixW<0)?"\t":" ");
        }
      // print expectation and standard deviation
        ostrm 
        << GlobalVar.Double2String(eV[i],false,-1,fixW) << ((fixW<0)?"\t":" ")
        << GlobalVar.Double2String(sV[i],false,-1,fixW) << ((fixW<0)?"\t":" ");
      // compute the confidence levels
        sorVv = flxVec(tp+size_t(N)*i,Nc);
        std::sort(sorVp,sorVp+Nc);
        for (tuint j=0;j<clN;++j) {
          tdouble dt;
          if (clp[j] <= clh) {
            dt = sorVp[0];
          } else if (clp[j] >= (ONE-clh)) {
            dt = sorVp[Nc-1];
          } else {
            const tuint id = (clp[j]-clh)/(2*clh);
            const tdouble pd = clp[j]-(id+ONE/2)/Nc;
            #if FLX_DEBUG
              if (pd<ZERO||pd>2*clh) throw FlxException_Crude("FlxStatBox::plot_confidence_4");
            #endif
            dt = sorVp[id]+(sorVp[id+1]-sorVp[id])*(pd/(2*clh));
          }
          ostrm << GlobalVar.Double2String(dt,false,-1,fixW) << ((fixW<0)?"\t":" ");
        }
      // end the line ;)
        ostrm << std::endl;
    } 
}

void FlxStatBox::get_smpl_eigen(const tdouble p, const tuint pNmax, flxVec& sde, FlxMtxSym& covmtx, flxVec& hvN1, flxVec& hvN2, flxVec& hvM, iVec& ivN, std::vector< flxVec >& Eigenvectors, FlxMtx& Tinv, const tdouble p_single, const tuint pNmax_single)
{
  if (Nc==0) {
    std::ostringstream ssV;
    ssV << "The statBox does not contain any entries.";
    throw FlxException_NeglectInInteractive("FlxStatBox::get_smpl_eigen_0", ssV.str() );
  }
  tdouble* hvN1p = hvN1.get_tmp_vptr();
  tdouble* hvN2p = hvN2.get_tmp_vptr();
  tdouble* hvMp = hvM.get_tmp_vptr();
  tdouble* sdep = sde.get_tmp_vptr();
  tuint* ivNp = &(ivN[0]);                        // contains the Np IDs that are considered
  // determine the relevant samples
    for (tuint i=0;i<Nc;++i) {                // compute distances to y-vec
      hvM = sde;
      for (tuint j=0;j<M;++j) {
        hvMp[j] -= tp[i+j*size_t(N)];
      }
      hvN1p[i] = hvM.get_Norm2_NOroot();
    }
    hvN2 = hvN1;
    std::sort(hvN2p,hvN2p+Nc);
    tuint Np = Nc*p;        // the number of samples to consider
    if (pNmax>0 && Np>pNmax) Np = pNmax;
    // find the threshold value
      tuint Np_single_max = Np*p_single;        // the maximum number of the same sample appearing multiple times
      if (pNmax_single>0 && Np_single_max>pNmax_single) Np_single_max = pNmax_single;
      tdouble thr = -ONE;
      tuint Np_count = 0;
      tuint single_count = 0;
      tuint lai = 0;
      for (tuint i=0;i<Nc;++i) {
        if (hvN2p[i]==thr) {        // repeated sample
          ++single_count;
          if (single_count>Np_single_max) {        // sample appears to often ... get rid of one of them
            #if FLX_DEBUG
              tuint lai_prev = lai;
            #endif
            for (tuint j=lai;j<Nc;++j) {
              if (hvN1p[j]==thr) {
                hvN1p[j] = 9999999.123;
                lai = j+1;
                break;
              }
            }
            #if FLX_DEBUG
              if (lai_prev == lai) throw FlxException_Crude("FlxStatBox::get_smpl_eigen_1");
            #endif
          } else {                                // repetition still not critical
            ++Np_count;
          }
        } else {                // sample not repeated (so far)
          thr = hvN2p[i];
          ++Np_count;
          single_count = 1;
          lai = 0;
        }
        // the stopping criterion
          if (Np_count>=Np) {
            if (i+1<Nc) {
              thr = (hvN2p[i]+hvN2p[i+1])/2;        // if (hvN1p[i]<=thr) this index is in the selected set
            }
            break;
          }
      }
      if (Np_count!=Np) {
        GlobalVar.alert.alert("FlxStatBox::get_smpl_eigen_2","Too many multiple occurrences.");
        Np = Np_count;
      }
    // count the actual number of samples that are considered
      Np = 0;
      for (tuint i=0;i<Nc;++i) {
        if (hvN1p[i]<=thr) {
          ivNp[Np++] = i;
        }
      }
  // compute the mean and standard deviation
    pdouble pp; tdouble tt; tdouble sc = ZERO;
    for (tuint i=0;i<M;++i) {
      pp = ZERO;
      const size_t st = size_t(N)*i;
      for (tuint j=0;j<Np;++j) {
        pp += tp[st + ivNp[j]];
      }
      tt = pp.cast2double()/Np;
      pp = ZERO;
      for (tuint j=0;j<Np;++j) {
          pp += pow2(tp[st + ivNp[j]]-tt);
      }
      sdep[i] = tt;                                        // this is the ith mean
      tt = pp.cast2double()/Np;
      sc += tt;
      covmtx(i,i) = tt;                                        // this is the ith variance
    }
    if (sc<GlobalVar.TOL()) {
      GlobalVar.alert.alert("FlxStatBox::get_smpl_eigen_3","There seems to be a convergence problem.");
    }
  // compute the covariances
    for (tuint i=0;i<M;++i) {
      for (tuint j=0;j<i;++j) {
        pp = ZERO;
        const size_t st1 = size_t(N)*j;
        const size_t st2 = size_t(N)*i;
        for (tuint k=0;k<Np;++k) {
          pp += (tp[st1 + ivNp[k]]-sdep[j]) * (tp[st2 + ivNp[k]]-sdep[i]);
        }
        covmtx(i,j) = pp.cast2double()/Np;
      }
    }
  // solve the eigenvalue problem
    MtxEigenValue(covmtx,M,sde,Eigenvectors,2);
    for (tuint i=0;i<M;++i) {                // normalize the eigenvectors
      flxVec& et = Eigenvectors[i];
      const tdouble etl = et.get_Norm2();
      et /= etl;
    }
  // assemble the transformation matrix
    for (tuint i=0;i<M;++i) {
      flxVec& et = Eigenvectors[i];
      for (tuint j=0;j<M;++j) {
        Tinv(j,i) = et[j];
      }
    }
  // compute the standard deviations in eigen-space
    for (tuint i=0;i<M;++i) {
      if (sde[i]<ZERO) sde[i] = ZERO;        // to prevent -nan 
      else sde[i] = sqrt(sde[i]);
    }
}

void FlxStatBox::get_smpl_eigen(FlxMtxSym& covmtx, flxVec& evs, std::vector< flxVec >& Eigenvectors, FlxMtx& Tinv)
{
  if (Nc==0) {
    std::ostringstream ssV;
    ssV << "The statBox does not contain any entries.";
    throw FlxException_NeglectInInteractive("FlxStatBox::get_smpl_eigen_100", ssV.str() );
  }
  tdouble* sdep = evs.get_tmp_vptr();
  // compute the mean and standard deviation
    pdouble pp; tdouble tt; tdouble sc = ZERO;
    for (tuint i=0;i<M;++i) {
      pp = ZERO;
      const size_t st = size_t(N)*i;
      for (tuint j=0;j<Nc;++j) {
        pp += tp[st + j];
      }
      tt = pp.cast2double()/Nc;
      pp = ZERO;
      for (tuint j=0;j<Nc;++j) {
          pp += pow2(tp[st + j]-tt);
      }
      sdep[i] = tt;                                        // this is the ith mean
      tt = pp.cast2double()/Nc;
      sc += tt;
      covmtx(i,i) = tt;                                        // this is the ith variance
    }
    if (sc<GlobalVar.TOL()) {
      GlobalVar.alert.alert("FlxStatBox::get_smpl_eigen_101","There seems to be a convergence problem.");
    }
  // compute the covariances
    for (tuint i=0;i<M;++i) {
      for (tuint j=0;j<i;++j) {
        pp = ZERO;
        const size_t st1 = size_t(N)*j;
        const size_t st2 = size_t(N)*i;
        for (tuint k=0;k<Nc;++k) {
          pp += (tp[st1 + k]-sdep[j]) * (tp[st2 + k]-sdep[i]);
        }
        covmtx(i,j) = pp.cast2double()/Nc;
      }
    }
  // solve the eigenvalue problem
    MtxEigenValue(covmtx,M,evs,Eigenvectors,2);
    for (tuint i=0;i<M;++i) {                // normalize the eigenvectors
      flxVec& et = Eigenvectors[i];
      const tdouble etl = et.get_Norm2();
      et /= etl;
    }
  // assemble the transformation matrix
    for (tuint i=0;i<M;++i) {
      flxVec& et = Eigenvectors[i];
      for (tuint j=0;j<M;++j) {
        Tinv(j,i) = et[j];
      }
    }
  // compute the standard deviations in eigen-space
    for (tuint i=0;i<M;++i) {
      if (evs[i]<ZERO) evs[i] = ZERO;        // to prevent -nan 
      else evs[i] = sqrt(evs[i]);
    }
}



