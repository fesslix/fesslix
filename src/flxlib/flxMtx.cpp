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

#include "flxMtx.h"

#include <algorithm>
#include <execution>


//=================================================================
// Do we want to activate '__cpp_lib_parallel_algorithm'
//=================================================================
// this section needs to be placed after "#include <execution>"
#ifdef __WINDOWS__
    // do not activate it under Windows, as compilation run into problems
    #define FLX_parallel_algorithm 0
#else
    // otherwise, thrust the flag '__cpp_lib_parallel_algorithm'
    #ifdef __cpp_lib_parallel_algorithm
        #define FLX_parallel_algorithm 1
    #else
        #define FLX_parallel_algorithm 0
    #endif
#endif




void MtxProd_BTKB_mtx(const FlxMtx& B, const FlxMtxSym& E, FlxMtxSym& K)
{
  #if FLX_DEBUG
    if (K.ncols() != B.ncols() || B.nrows() != E.nrows() ) {
      std::ostringstream ssV;
      ssV << "Matrices can not be multiplied.";
      throw FlxException("MtxProd_BTKB_1", ssV.str() );
    }
  #endif
  const tdouble* Bp = B.get_mtx_flxVec().get_tmp_vptr_const();
  const tdouble* Ep = &E.mtx[0];
  tdouble* Kp = &K.mtx[0];
  const tnlong a = K.nrows();
  const tnlong b = E.nrows();
  tVec h(b); 
  tdouble* hp = &h[0];
  #if FLX_KAHAN_MTX_FULL
    pdouble t;
  #else
    tdouble t;
  #endif
  tnlong n, m, k, j;
  for (tnlong i = 0; i < a; ++i) {
    h=ZERO;
    for (n=0; n < b; ++n) {
      t = ZERO;
      for (m=0; m < b; ++m) {
        t+=Bp[m*a+i]*((m<=n)?(Ep[(n*n+n)/2+m]):(Ep[(m*m+m)/2+n]));
      }
      #if FLX_KAHAN_MTX_FULL
        hp[n]=t.cast2double();
      #else
        hp[n]=t;
      #endif
    }
    for (j=i; j<a; ++j) {
      k=(j*j+j)/2+i;
      t = ZERO;
      for (n=0; n<b; ++n) {
        t+=hp[n]*Bp[n*a+j];
      }
      #if FLX_KAHAN_MTX_FULL
        Kp[k]=t.cast2double();
      #else
        Kp[k]=t;
      #endif
    }
  }
}

FlxMtxSym MtxProd_BTKB(const FlxMtx& B, const FlxMtxSym& E)
{
  FlxMtxSym K(B.ncols());
  MtxProd_BTKB_mtx(B,E,K);
  return K;
}

void MtxProd_BTKB(const FlxMtxTransformation& B, const FlxMtxSym& E, FlxMtxSym& K)
{
  #if FLX_DEBUG
    if (K.ncols() != B.ncols() || B.nrows() != E.nrows() ) {
      std::ostringstream ssV;
      ssV << "Matrices can not be multiplied.";
      throw FlxException("MtxProd_BTKB_2", ssV.str() );
    }
  #endif
  const size_t S = B.Ttm.size();
  tnlong c = 0; tnlong r;
  const FlxMtx *M, *N;
  size_t Mrows,Mcols,Nrows;
  flxVec& Km = K.mtx;
  tnlong l,j,k,i;        // loop variables
  tnlong m,n; 
  flxVec h(K.nrows());
  // search maximum dim of g and f
    m=0;
    for (i=0;i<S;++i) {
      if (B.Ttm[i]->nrows() > m) m = B.Ttm[i]->nrows();
    }
    tVec g(m);
    flxVec f(m);
  for (i=0; i<S;++i) {
    M = B.Ttm[i];
    const flxVec& Mm = M->get_mtx_flxVec();
    Mrows = M->nrows(); Mcols = M->ncols();
    flxVec g(Mrows);
    for (j=0;j<Mcols;++j) {
      g.slice(Mm.get_tmp_vptr_const()+j,Mrows);
      E.MultMv_slice(g,h,c,c+Mcols-1);
      r = 0;
      for (k=0;k<S;++k) {
        N = B.Ttm[k];
        Nrows = N->nrows();
        if (r+Nrows-1>=c+j) { 
          const flxVec hhelp(h.get_tmp_vptr_const(),Nrows);
          N->TransposeMmultVec(hhelp,f);
          n = (r>c+j)?(r):(c+j);
          for (l=r+Nrows;l>n;--l) {
            m = l - 1;
            Km[(m*m+m)/2+c+j]=f[l-r-1];
          }
        }
        r+=Nrows;
      }
    }
    c += Mcols;
  }
}

FlxMtxSym MtxProd_BTKB(const FlxMtxTransformation& B, const FlxMtxSym& E)
{
  FlxMtxSym K(B.ncols());
  MtxProd_BTKB(B,E,K);
  return K;
}

void MtxProd_BTKB_ltri(const FlxMtxLTri& B, const FlxMtxSym& E, FlxMtxSym& K)
{
  #if FLX_DEBUG
    if (K.ncols() != B.ncols() || B.nrows() != E.nrows() ) {
      std::ostringstream ssV;
      ssV << "Matrices can not be multiplied.";
      throw FlxException("MtxProd_BTKB_3", ssV.str() );
    }
  #endif
  const tdouble* Bp = &B.mtx[0];
  const tdouble* Ep = &E.mtx[0];
  tdouble* Kp = &K.mtx[0];
  const tnlong a = K.nrows();
  const tnlong b = E.nrows();
  #if FLX_KAHAN_MTX_LTRI
    pdouble t;
  #else
    tdouble t;
  #endif
  tVec h(b);
  tdouble* hp = &h[0];
  tnlong i,n,j,k,m;
  for (i = 0; i < a; ++i) {
    for (n=0; n < b; ++n) {
      t = ZERO;
      for (m=i; m < b; ++m) {
        if (n>=m) {
          t+=Bp[(m*m+m)/2+i]*Ep[(n*n+n)/2+m];
        } else {
          t+=Bp[(m*m+m)/2+i]*Ep[(m*m+m)/2+n];
        }
      }
      #if FLX_KAHAN_MTX_LTRI
        hp[n]=t.cast2double();
      #else
        hp[n]=t;
      #endif
    }
    for (j=i; j<a; ++j) {
      k=(j*j+j)/2+i;
      t = ZERO;
      for (n=j; n< b; ++n) {
        t+=hp[n]*Bp[(n*n+n)/2+j];
      }
      #if FLX_KAHAN_MTX_LTRI
        Kp[k]=t.cast2double();
      #else
        Kp[k]=t;
      #endif
    }
  } 
}

FlxMtxSym MtxProd_BTKB(const FlxMtxLTri& B, const FlxMtxSym& E)
{
  FlxMtxSym K(B.ncols());
  MtxProd_BTKB_ltri(B,E,K);
  return K;
}


void MtxProd_BTKB(const FlxMtx_base& B, const FlxMtx_baseS& E, FlxMtxSym& K)
{
  #if FLX_DEBUG
    if (K.ncols() != B.ncols() || B.nrows() != E.nrows() ) {
      std::ostringstream ssV;
      ssV << "Matrices can not be multiplied.";
      throw FlxException("MtxProd_BTKB_4", ssV.str() );
    }
  #endif
  tdouble* Kp = &K.mtx[0];
  const tnlong a = K.nrows();
  const tnlong b = E.nrows();
  #if FLX_KAHAN_MTX_SYM
    pdouble t;
  #else
    tdouble t;
  #endif
  tVec h(b);
  tdouble* hp = &h[0];
  tnlong i,n,m,j,k;
  for (i = 0; i < a; i++) {
    for (n=0; n < b; n++) {
      t = ZERO;
      for (m=0; m < b; m++) {
        t+=B(m,i)*E(n,m);
      }
      #if FLX_KAHAN_MTX_SYM
        hp[n] = t.cast2double();
      #else
        hp[n] = t;
      #endif
    }
    for (j=i; j<a; j++) {
      k=(j*j+j)/2+i;
      t = ZERO;
      for (n=0; n<b; n++) {
        t+=hp[n]*B(n,j);
      }
      #if FLX_KAHAN_MTX_SYM
        Kp[k] = t.cast2double();
      #else
        Kp[k] = t;
      #endif
    }
  }
}

FlxMtxSym MtxProd_BTKB(const FlxMtx_base& B, const FlxMtx_baseS& E)
{
  FlxMtxSym K(B.ncols());
  MtxProd_BTKB(B,E,K);
  return K;
}

void MtxProd_BTKB_1D(const FlxMtx& B, const FlxMtxSym& E, FlxMtxSym& K)
{
  #if FLX_DEBUG
    if (K.ncols()!=6 || B.nrows()!=3 || B.ncols()!=3 || E.nrows()!=2 ) {
      std::ostringstream ssV;
      ssV << "Matrices can not be multiplied.";
      throw FlxException("MtxProd_BTKB_1D_1", ssV.str() );
    }
  #endif
  tdouble t1,t2,t3;
  const flxVec& Emtx = E.get_mtx_flxVec();
  const flxVec& Bmtx = B.get_mtx_flxVec();
  tdoublePtr Kmtx = K.get_VecPointer();
  for (tuint i=0; i<3; ++i) {
    t1 = Emtx[0] * Bmtx[i];
    t2 = Emtx[1] * Bmtx[i];
    t3 = Emtx[2] * Bmtx[i];
    for (tuint j=i; j<3; ++j) {
      Kmtx[(j*j+j)/2+i] = t1 * Bmtx[j];
    }
    for (tuint j=0; j<3; ++j) {
      Kmtx[(j*j+7*j)/2+i+6] = t2 * Bmtx[j];
    }
    for (tuint j=i; j<3; ++j) {
      Kmtx[(j*j+7*j)/2+i+9] = t3 * Bmtx[j];
    }
  }
}

void MtxProd_BTKB_1D_part(const FlxMtx& B, const FlxMtxSym& E, FlxMtxSym& K, bool tranformFirst)
{
  #if FLX_DEBUG
    if (K.ncols()!=4 || B.nrows()!=3 || B.ncols()!=3 || E.nrows()!=2 ) {
      std::ostringstream ssV;
      ssV << "Matrices can not be multiplied.";
      throw FlxException("MtxProd_BTKB_1D_part_1", ssV.str() );
    }
  #endif
  tdouble t1,t2,t3;
  const flxVec& Emtx = E.get_mtx_flxVec();        // the element stiffness matrix (local)
  const flxVec& Bmtx = B.get_mtx_flxVec();        // the transformation matrix
  tdoublePtr Kmtx = K.get_VecPointer();        // the element stiffness matrix (global)
  if (tranformFirst) {
    Kmtx[9] = Emtx[2];
    for (tuint i=0; i<3; ++i) {
      t1 = Emtx[0] * Bmtx[i];
      t2 = Emtx[1] * Bmtx[i];
      for (tuint j=i; j<3; ++j) {
        Kmtx[(j*j+j)/2+i] = t1 * Bmtx[j];
      }
      Kmtx[6+i] = t2;
    }
  } else {
    Kmtx[0] = Emtx[0];
    for (tuint i=0; i<3; ++i) {
      t2 = Emtx[1] * Bmtx[i];
      t3 = Emtx[2] * Bmtx[i];
      for (tuint j=i; j<3; ++j) {
        Kmtx[(j*(j+3))/2+i+2] = t3 * Bmtx[j];
      }
      Kmtx[(i*(i+3))/2+1] = t2;
    }
  }
}

FlxMtxSym VecDyadProd(const tVec& v)
{
  const size_t N = v.size();
  FlxMtxSym B(N);
  tdoublePtr Bmtx = B.get_VecPointer();
  for (size_t i = 0; i < N; ++i) {
    for (size_t j = 0; j <= i; ++j) {
      Bmtx[(i*i+i)/2+j]=v[i]*v[j];
    }
  }
  return B;
}


std::ostream& operator<<(std::ostream& os, const FlxMtx_base& M)
{
  M.output_Mtx(os);
  return os;
}

void FlxMtx_base::MultMtx(const FlxMtx_base& rhs, FlxMtx& res)
{
  const FlxMtx_base& lhs = *this;
  #if FLX_DEBUG
    if (lhs.ncols() != rhs.nrows()) throw FlxException_Crude("FlxMtx_base::MultMtx_1");
    if (lhs.nrows()!=res.nrows() || rhs.ncols()!=res.ncols()) throw FlxException_Crude("FlxMtx_base::MultMtx_2");
  #endif
  const tnlong a = lhs.nrows();
  const tnlong b = rhs.ncols();
  const tnlong c = lhs.ncols();
  tnlong k; tnlong j; tnlong i;
  for (i = 0; i<a;i++) {
    for (j = 0; j < b; j++) {
      for (k = 0; k < c; k++) {
        res(i,j)+=lhs(i,k)*rhs(k,j);
      }
    }
  }
}


const flxVec FlxMtx_base::operator*(const flxVec& v ) const
{
  #if FLX_DEBUG
    if (v.get_N() != ncols()) {
      std::ostringstream ssV;
      ssV << "Vector/Matrix sizes do not match.";
      throw FlxException("FlxMtxSparsSym::operator*_2", ssV.str() );
    }
  #endif
  flxVec w(nrows());
  MultMv(v,w);
  return w;
}

const flxpVec FlxMtx_base::operator*(const flxpVec& v) const
{
  #if FLX_DEBUG
    if (v.get_N() != ncols()) {
      std::ostringstream ssV;
      ssV << "Vector/Matrix sizes do not match.";
      throw FlxException("FlxMtxSparsSym::operator*_3", ssV.str() );
    }
  #endif
  flxpVec w(nrows());
  MultMv(v,w);
  return w;
}

void FlxMtx_base::MultMv(const flxpVec& v, flxpVec& w) const
{
  throw FlxException_NotImplemented("FlxMtx_base::MultMv");
}

FlxMtx_base* FlxMtx_base::copy()
{
  throw FlxException_NotImplemented("FlxMtx_base::copy");
}

void FlxMtx_base::output_Mtx(std::ostream& sout) const
{
  tnlong N = nrows();
  tnlong M = ncols();
  for (tnlong i = 0; i<N;i++) {
    for (tnlong j = 0; j<M; j++) {
      sout << GlobalVar.Double2String( operator()(i,j) ) << "\t";
    }
    sout << std::endl;
  }
}

void FlxMtx_base::output_OctaveMtx(std::ostream& sout, bool checkTOLv, bool doendl) const
{
  tnlong N = nrows();
  tnlong M = ncols();
  sout << "[ ";
  for (tnlong i = 0; i<N;i++) {
    for (tnlong j = 0; j<M; j++) {
      sout << GlobalVar.Double2String(operator()(i,j),checkTOLv);
      if (j < M-1) sout << " ";
    }
    if (i < N-1) sout << "; ";
    else sout << "]";
    if (doendl) sout << std::endl;
  }
}

const tdouble FlxMtx_base::max() const
{
  const tnlong N = nrows();
  const tnlong M = ncols();
  tdouble d = operator()(0,0);
  tdouble dt;
  for (tnlong i = 0; i<N;i++) {
    for (tnlong j = 0; j<M; j++) {
      dt = operator()(i,j);
      if ( dt > d) {
        d = dt;
      }
    }
  }
  return d;
}

const tdouble FlxMtx_base::min() const
{
  const tnlong N = nrows();
  const tnlong M = ncols();
  tdouble d = operator()(0,0);
  tdouble dt;
  for (tnlong i = 0; i<N;i++) {
    for (tnlong j = 0; j<M; j++) {
      dt = operator()(i,j);
      if (dt < d) {
        d = dt;
      }
    }
  }
  return d;
}

const size_t FlxMtx_base::maxID() const
{
  const tnlong N = nrows();
  const tnlong M = ncols();
  tdouble d = operator()(0,0);
  tdouble dt;
  size_t c = 0;
  for (tnlong i = 0; i<N;i++) {
    for (tnlong j = 0; j<M; j++) {
      dt = operator()(i,j);
      if ( dt > d) {
        d = dt;
        c = i*M+j;
      }
    }
  }
  return c;
}

const size_t FlxMtx_base::minID() const
{
  const tnlong N = nrows();
  const tnlong M = ncols();
  tdouble d = operator()(0,0);
  tdouble dt;
  size_t c = 0;
  for (tnlong i = 0; i<N;i++) {
    for (tnlong j = 0; j<M; j++) {
      dt = operator()(i,j);
      if ( dt < d) {
        d = dt;
        c = i*M+j;
      }
    }
  }
  return c;
}

tdouble* FlxMtx_base::get_internalPtr(const bool throwErr)
{
  if (throwErr) {
    throw FlxException_NotImplemented("FlxMtx_base::get_internalPtr");
  } else {
    return NULL;
  }
}

void FlxMtx_baseS::assembleMinv(int i)
{
  if (Minv!=NULL) {
    delete Minv;
    Minv = NULL;
  }
  if (i == 0) {
    Minv = new FlxMtxIdentity(this->nrows());
  } else if ( i == 1) {
    FlxMtxDiag* d = new FlxMtxDiag(*this);
    Minv = d->get_Inverse();
    delete d;
  } else if ( i == 3) {
    FlxMtxSparsSym* S = dynamic_cast<FlxMtxSparsSym*>(this);
    if (S != NULL) {
      FlxMtxLTri L(S->nrows());
      L.CholeskyDec(*S);
      L.Invert();
      FlxMtxSym* MinvS = new FlxMtxSym(S->nrows());
      MinvS->assign_LTL(L);
      Minv = MinvS;
      return;
    } 
    FlxMtxSym* SF = dynamic_cast<FlxMtxSym*>(this);
    if (SF != NULL) {
      FlxMtxLTri L(SF->nrows());
      L.CholeskyDec(*SF);
      L.Invert();
      FlxMtxSym* MinvS = new FlxMtxSym(SF->nrows());
      MinvS->assign_LTL(L);
      Minv = MinvS;
      return;
    } 
    std::ostringstream ssV;
    ssV << "Matrix is not a symmetric sparse matrix nor an symmetric dense matrix.";
    throw FlxException("FlxMtx_baseS::assembleMinv_1", ssV.str() );
  } else if ( i==4 || i == 5 || i == 6 ) {
    FlxMtxSparsSym* S = dynamic_cast<FlxMtxSparsSym*>(this);
    if (S != NULL) {
      if (i==4) {
        Minv = new FlxMtxSparsSymLU(*S);
      } else if (i==5) {
        Minv = new FlxMtxSparsSymILU(*S);
      } else {
        Minv = new FlxMtxSparsSymILU(*S,true);
      }
      return;
    }
    std::ostringstream ssV;
    ssV << "Matrix is not a symmetric sparse matrix.";
    throw FlxException("FlxMtx_baseS::assembleMinv_2", ssV.str() );
  } else {
    Minv = NULL;
  }
}

FlxMtx_baseS::~FlxMtx_baseS()
{
  if (Minv) delete Minv;
}

void FlxMtx_baseS::preconding(const flxVec& r, flxVec& z, const int i) const
{
  throw FlxException_NotImplemented("FlxMtx_baseS::preconding_1");
}

void FlxMtx_baseS::preconding(const flxpVec& r, flxpVec& z, const int i) const
{
  throw FlxException_NotImplemented("FlxMtx_baseS::preconding_2");
}

const bool FlxMtx_baseS::solve_CG(flxVec& x, const flxVec& f, tdouble& eps, tnlong& iter, const tuint pcn, const bool startZero)
{
  if (ncols() != nrows() || nrows() != f.get_N() || f.get_N() != x.get_N() ) {
    std::ostringstream ssV;
    ssV << "ncols=" << ncols() << "; nrows=" << nrows() << "; f_size=" << f.get_N() << "; x_size=" << x.get_N() << ";";
    throw FlxException("FlxMtx_base::solve_CG_1", "Vector/Matrix sizes do not match.", ssV.str() );
  }
  
  // check for trivial solution
  if (f.get_NormMax()==ZERO) {
    x.set_zero();
    iter = 0;
    eps = ZERO;
    return true;
  }
  
  assembleMinv(pcn);
  
  if (startZero) { x.set_zero();; }
  
  const tnlong maxiter = iter;
  
  #if FLX_KAHAN_CG_STAGE
    flxpVec X(x);
    flxpVec tttp(f.get_N());
  #endif
  flxVec ttt(f.get_N());
  MultMv(x,ttt);
  #if FLX_KAHAN_CG_STAGE_4
    flxpVec r(f);
    r -= ttt;
    flxpVec z(f.get_N());
  #else
    flxVec r(f);
    r -= ttt;
    #if FLX_KAHAN_CG_STAGE_3 || FLX_KAHAN_CG_STAGE_2
      flxpVec R(r);
    #endif
    flxVec z(f.get_N());
  #endif
  if (Minv==NULL) {
    preconding(r, z, pcn);
  } else {
    Minv->MultMv(r,z);
  }
  #if FLX_KAHAN_CG_STAGE_4
    flxpVec p(z);
    pdouble zr = r*z;
  #else
    flxVec p(z);
    #if FLX_KAHAN_CG_STAGE_3
      flxpVec P(p);
    #endif
    tdouble zr = r*z;
  #endif
  const double stp = eps*f.get_Norm2();
 
  #if FLX_KAHAN_CG_STAGE_4
    flxVec rr(r);
  #else
    flxVec& rr = r;
  #endif
  if (rr.get_NormMax() <= GlobalVar.TOL() ) {
    eps = ZERO; 
    iter = 0;
    return true;
  }

  #if FLX_KAHAN_CG_STAGE_4
    pdouble pap,alpha,zrold,beta;
    flxpVec mp(p.get_N());
  #else
    tdouble pap,alpha,zrold,beta;
    flxVec mp(p.get_N());
  #endif
  for (iter = 0; iter < maxiter; iter++) {
    MultMv(p,mp);        //     mp = (*this)*p;
    pap = p*mp;
    alpha = zr/pap;
    #if FLX_KAHAN_CG_STAGE_4
      X.add(p,alpha);
      r.add(mp,-1.*alpha);
      rr = r;
    #else
      #if FLX_KAHAN_CG_STAGE_3 || FLX_KAHAN_CG_STAGE_2
        X.add(p,alpha);
        R.add(mp,-1.*alpha);
        r = R;
      #else
        #if FLX_KAHAN_CG_STAGE_1
          X.add(p,alpha);
        #else
          x.add(p,alpha);
        #endif
        r.add(mp,-1.*alpha);
      #endif
    #endif
    if (rr.get_Norm2() <= stp) {
      break;
    }
    if (Minv==NULL) {
      preconding(r, z, pcn);
    } else {
      Minv->MultMv(r,z);
    }
    zrold = zr;
    zr = r*z;
    beta = zr/zrold;
    #if FLX_KAHAN_CG_STAGE_3
      P *= beta;
      P += z;
      p = P;
    #else
      p*=beta;
      p+=z;
    #endif
  }

  #if FLX_KAHAN_CG_STAGE_4
    rr = r;
  #endif
  eps = rr.get_Norm2();
  #if FLX_KAHAN_CG_STAGE
    x = X;
    MultMv(X,tttp);
    ttt = tttp;
  #else
    MultMv(x,ttt);
  #endif
  rr = f;
  rr -= ttt;
  
  // check for difference ...
    tdouble de = rr.get_Norm2();
    if (de > stp && eps > ZERO) {
      std::ostringstream ssV;
      ssV << "Drifting error (" << GlobalVar.Double2String(fabs(de-eps)/eps) << "). ";
      ssV << "Residuum (" << de << ") not within bound (" << stp << "). (iter=" << iter << "; unknowns=" << x.get_N() << ")";
      GlobalVar.alert.alert("FlxMtx_base::solve_CG_3",ssV.str());
    }
    #if FLX_DEBUG
      if (eps>ZERO) {
        GlobalVar.slog(4) << "  Drifting Error: " << GlobalVar.Double2String(fabs(de-eps)/eps) << std::endl;
      }
    #endif
    eps = de;
  x.check_TOL();;
  if (iter==maxiter) return false;
  else return true;
}

FlxMtxDiag* FlxMtx_baseS::get_diag() const
{
  return new FlxMtxDiag(*this);
}

FlxMtx_baseS* FlxMtx_baseS::get_Inverse()
{
  throw FlxException_NotImplemented("FlxMtx_baseS::get_Inverse");
}

void FlxMtx_baseS::Invert()
{
  throw FlxException_NotImplemented("FlxMtx_baseS::Invert");
}

tdouble* FlxMtx_baseS::get_VecPointer()
{
  throw FlxException_NotImplemented("FlxMtx_baseS::get_VecPointer");
}

void FlxMtx_baseS::get_VecPointer_full(tdouble* vec) const
{
  const tnlong msize = nrows();
  for (tnlong i=0;i<msize;++i) {
    for (tnlong j=0;j<i;++j) {
      vec[j*msize+i] = vec[i*msize+j] = this->operator()(i,j);
    }
    vec[i*msize+i] = this->operator()(i,i);
  }
}


void FlxMtx_baseS::output_Mtx(std::ostream& sout) const
{
  for (tnlong i =0;i<nrows();i++) {
    for (tnlong j = 0; j <= i; j++) {
      sout << GlobalVar.Double2String( operator()(i,j), false,-1,0 ) << " ";
    }
    sout << std::endl;
  }
}

FlxMtx& FlxMtx::operator+=(const FlxMtx& m)
{
  #if FLX_DEBUG
    if (m.ncols() != ncols() || m.nrows() != nrows() ) {
      std::ostringstream ssV;
      ssV << "Matrices can not be added.";
      throw FlxException("FlxMtxSym::operator+=", ssV.str() );
    }
  #endif
  mtx+=m.mtx;
  return *this;
}

// FlxMtx FlxMtx::operator*(const FlxMtx_base& m) const
// {
//   #if FLX_DEBUG
//     if (m.nrows() != ncols()) {
//       std::ostringstream ssV;
//       ssV << "Matrices can not be multiplied.";
//       throw FlxException("FlxMtx::operator*", ssV.str() );
//     }
//   #endif
//   FlxMtx r(nrows(),m.ncols());
//   for (tnlong i = 0; i < nrows(); i++) { // rows
//     for (tnlong j = 0; j < r.ncols(); j++) { // colums
//       r(i,j)=ZERO;
//       for (tnlong k = 0; k < ncols(); k++) {  // Scalar-Product
//         r(i,j)+=operator()(i,k)*r(k,j);
//       }
//     }
//   }
// }

tdouble& FlxMtx::operator()(const tnlong& i, const tnlong& j)
{
  #if FLX_DEBUG
    if (i < 0 || j < 0 || i >= nrows() || j >= ncols()) {
      std::ostringstream ssV;
      ssV << "Index (" << i << ", " << j << ") not within matrix.";
      throw FlxException("FlxMtx::operator()_1", ssV.str() );
    }
  #endif
  return mtx[i*csize+j];
}

const tdouble FlxMtx::operator()(const tnlong& i, const tnlong& j) const
{
  #if FLX_DEBUG
    if (i < 0 || j < 0 || i >= nrows() || j >= ncols()) {
      std::ostringstream ssV;
      ssV << "Index (" << i << "," << j << ") not within matrix (" << nrows() << "," << ncols() << ").";
      throw FlxException("FlxMtx::operator()_2", ssV.str() );
    }
  #endif
  return mtx[i*csize+j];
}

void FlxMtx::MultMv(const flxVec& v, flxVec& w) const
{
  #if FLX_KAHAN_MTX_FULL
    pdouble t;
  #else
    tdouble t;
  #endif
  tnlong n = 0;
  const tdouble* mtxp = mtx.get_tmp_vptr_const();
  const tdouble* vp = v.get_tmp_vptr_const();
  tdouble* wp = w.get_tmp_vptr();
  for (tnlong i = 0; i<rsize;++i) {
    t = ZERO;
    for (tnlong j=0;j<csize;++j) {
      t+=mtxp[n]*vp[j];
      ++n;
    }
    #if FLX_KAHAN_MTX_FULL
      wp[i] = t.cast2double();
    #else
      wp[i] = t;
    #endif
  }
}

FlxMtx::FlxMtx(const FlxMtx_base& mB)
: rsize(mB.nrows()), csize(mB.ncols()), mtx(rsize*csize)
{
  tnlong c=0;
  for (tnlong i=0;i<rsize;++i) {
    for (tnlong j=0;j<csize;++j) {
      mtx[c++] = mB.operator()(i,j);
    }
  }
}

FlxMtx::FlxMtx(const FlxMtxSym* symM)
: rsize(symM->nrows()), csize(symM->ncols()), mtx(rsize*csize)
{
  symM->get_VecPointer_full(mtx.get_tmp_vptr());
}

FlxMtx* FlxMtx::copy()
{
  return new FlxMtx(*this);
}

void FlxMtx::add_value(const tnlong& i, const tnlong& j, const tdouble& v)
{
  #if FLX_DEBUG
    if (i < 0 || j < 0 || i >= nrows() || j >= ncols()) {
      std::ostringstream ssV;
      ssV << "Index not within matrix.";
      throw FlxException("FlxMtx::add_value", ssV.str() );
    }
  #endif
  mtx[i*ncols()+j]+=v;
}

void FlxMtx::TransposeMmultM(FlxMtxSym& res) const
{
  #if FLX_DEBUG
    if (res.nrows() != ncols()) {
      std::ostringstream ssV;
      ssV << "ERROR.";
      throw FlxException("FlxMtx::TransposeMmultM", ssV.str() );
    }
  #endif
  const tnlong c = ncols();
  const tnlong r = nrows();
  #if FLX_KAHAN_MTX_FULL
    pdouble d;
  #else
    tdouble d;
  #endif
  const tdouble* mtxp = mtx.get_tmp_vptr_const();
  tnlong k; tnlong j; tnlong i;
  for (i = 0; i < c; ++i) {
    for (j = 0; j <= i; ++j) {
      d = ZERO;
      for (k = 0; k < r; ++k) {
        d+= mtxp[k*ncols()+i]*mtxp[k*c+j];
      }
      #if FLX_KAHAN_MTX_FULL
        res.operator()(i,j) = d.cast2double();
      #else
        res.operator()(i,j) = d;
      #endif
    }
  }
}

void FlxMtx::TransposeMmultVec(const flxVec& v, flxVec& w) const
{
  #if FLX_DEBUG
    if (v.get_N() != nrows() || w.get_N() != ncols() ) {
      std::ostringstream ssV;
      ssV << "ERROR.";
      throw FlxException("FlxMtx::TransposeMmultM", ssV.str() );
    }
  #endif
  #if FLX_KAHAN_MTX_FULL
    pdouble d;
  #else
    tdouble d;
  #endif
  const tdouble* mtxp = mtx.get_tmp_vptr_const();
  const tdouble* vp = v.get_tmp_vptr_const();
  tdouble* wp = w.get_tmp_vptr();
  tnlong k,i;
  for (i = 0; i < csize; ++i) {
    d = ZERO;
    for (k = 0; k < rsize; ++k) {
      d+= mtxp[k*csize+i]*vp[k];
    }
    #if FLX_KAHAN_MTX_FULL
      wp[i] = d.cast2double();
    #else
      wp[i] = d;
    #endif
  }
}

tdouble* FlxMtx::get_internalPtr(const bool throwErr)
{
  if (rsize*csize==0) throw FlxException_Crude("FlxMtx::get_internalPtr");
  return mtx.get_tmp_vptr();
}

FlxMtxTransformation::FlxMtxTransformation(const std::vector< FlxMtx* >& Ttm): Ttm(Ttm)
{
  rsize=0;
  for (size_t i=0; i < Ttm.size(); ++i) {
    #if FLX_DEBUG
      if (Ttm[i]->ncols() != Ttm[i]->nrows() ) {
        std::ostringstream ssV;
        ssV << "ERROR.";
        throw FlxException("FlxMtxTransformation::FlxMtxTransformation_1", ssV.str() );
      }
    #endif
    rsize+= Ttm[i]->ncols();
  }
}

FlxMtxTransformation::~FlxMtxTransformation()
{
  for (size_t i=0; i<Ttm.size(); ++i) {
    if (Ttm[i]!=NULL) {
      for (size_t j=i+1;j<Ttm.size();++j) {
        if (Ttm[i]==Ttm[j]) {
          Ttm[j] = NULL;
        }
      }
      delete Ttm[i];
      Ttm[i] = NULL;
    }
  }
}

FlxMtxTransformation& FlxMtxTransformation::operator*=(const tdouble& s)
{
  throw FlxException_NotImplemented("FlxMtxTransformation::operator*=");
}

void FlxMtxTransformation::add_value(const tnlong& i, const tnlong& j, const tdouble& v)
{
  throw FlxException_NotImplemented("FlxMtxTransformation::add_value");
}

const tdouble FlxMtxTransformation::operator()(const tnlong& n, const tnlong& m) const
{
  tnlong i = n;
  tnlong j = m;
  #if FLX_DEBUG
    if (i < 0 || j < 0 || i >= nrows() || j >= ncols()) {
      std::ostringstream ssV;
      ssV << "Index (" << i << "," << j << ") not within matrix (" << nrows() << "," << ncols() << ").";
      throw FlxException("FlxMtxTransformation::operator()_1", ssV.str() );
    }
  #endif
  tnlong c=0;
  for (size_t ii=0;ii<Ttm.size();++ii) {
    if (i<c+Ttm[ii]->nrows()) {
      if (j<c+Ttm[ii]->nrows() && j>=c) {
        i-=c;
        j-=c;
        return Ttm[ii]->operator()(i,j);
      } else {
        return ZERO;
      }
    }
    c+=Ttm[ii]->nrows();
  }
  std::ostringstream ssV;
  ssV << "ERROR";
  throw FlxException("FlxMtxTransformation::operator()_2", ssV.str() );
}

void FlxMtxTransformation::MultMv(const flxVec& v, flxVec& w) const
{
  #if FLX_DEBUG
    if (v.get_N()!=w.get_N() || v.get_N()!=rsize ) {
      std::ostringstream ssV;
      ssV << "ERROR";
      throw FlxException("FlxMtxTransformation::MultMv_1", ssV.str() );
    }
  #endif
  tnlong c = 0;
  tnlong s;
  const std::size_t Ttms = Ttm.size();
  for (std::size_t i = 0; i < Ttms; ++i) {
    s = Ttm[i]->nrows();
    flxVec wtmp(&(w[c]),s);
    const flxVec vtmp(&(v[c]),s);
    wtmp = (*Ttm[i])*vtmp;
    c+=s;
  }
}

inline void FlxMtxSym::add_value(const tnlong& i, const tnlong& j, const tdouble& v)
{
  #if FLX_DEBUG
    if (i < 0 || j < 0 || i >= nrows() || j >= ncols()) {
      std::ostringstream ssV;
      ssV << "Index not within matrix.";
      throw FlxException("FlxMtxSym::add_value", ssV.str() );
    }
  #endif
  tnlong iV = i; tnlong jV=j;
  if (jV>iV) std::swap(iV,jV);
  mtx[(iV*iV+iV)/2+jV]+=v;
}

void FlxMtxSym::add_mtx(const FlxMtx_baseS& m, const tdouble f)
{
  tdouble* mtxp = &(mtx[0]);
  tnlong c = 0;
  for (tnlong i=0;i<msize;++i) {
    for (tnlong j=0;j<i;++j) {
      mtxp[c++] += m(i,j)*f;
    }
    mtxp[c++] += m(i,i)*f;
  }
}

FlxMtxSym& FlxMtxSym::operator+=(const FlxMtxSym& m)
{
  #if FLX_DEBUG
    if (m.ncols() != ncols() || m.nrows() != nrows() ) {
      std::ostringstream ssV;
      ssV << "Matrices can not be added.";
      throw FlxException("FlxMtxSym::operator+=", ssV.str() );
    }
  #endif
  mtx+=m.mtx;
  return *this;
}

tdouble& FlxMtxSym::operator()(const tnlong& i, const tnlong& j)
{
  #if FLX_DEBUG
    if (i < 0 || j < 0 || i >= nrows() || j >= ncols()) {
      std::ostringstream ssV;
      ssV << "Index not within matrix.";
      throw FlxException("FlxMtxSym::operator()_1", ssV.str() );
    }
  #endif
  tnlong iV = i; tnlong jV=j;
  if (jV>iV) std::swap(iV,jV);
  return mtx[(iV*iV+iV)/2+jV];
}

void FlxMtxSym::get_VecPointer_full(tdouble* vec) const
{
  const tdouble* mtxp = &mtx[0];
  tnlong i,j;
  for (i=0;i<msize;++i) {
    for (j=0;j<i;++j) {
      vec[j*msize+i] = vec[i*msize+j] = mtxp[(i*i+i)/2+j];
    }
    vec[i*msize+i] = mtxp[(i*i+i)/2+i];
  }
}

void FlxMtxSym::Invert()
{
  FlxMtxLTri L(msize);
  Invert(L);
}

void FlxMtxSym::Invert(FlxMtxLTri& Ltri)
{
  Ltri.CholeskyDec(*this);
  Ltri.Invert();
  assign_LTL(Ltri);
}

void FlxMtxSym::assign_LTL(const FlxMtxLTri& L)
{
  const tnlong N = nrows();
  const tdouble* LTLmtxp = &L.mtx[0];
  tdouble* mtxp = &mtx[0];

  auto inner_fun = [&](const tnlong i, const tnlong j, tdouble &dp) {
        #if FLX_KAHAN_MTX_SYM
          pdouble sum;
        #else
          tdouble sum = ZERO;
        #endif
        for (tnlong k = i; k < N; ++k) {
          const tnlong m = (k*k+k)/2;
          sum += LTLmtxp[m+j]*LTLmtxp[m+i];
        }
        #if FLX_KAHAN_MTX_SYM
          dp = sum.cast2double();
        #else
          dp = sum;
        #endif
  };
  auto outer_fun = [&](const tnlong n, tdouble &dp) {
    const tnlong i = (tnlong)((sqrt(ONE+8*n)-ONE)/2);
    const tnlong N_i = (i*i+i)/2;
    const tnlong j = n-N_i;
    inner_fun(i,j,dp);
  };

  #if FLX_parallel_algorithm
  if (N>1000) {
    std::for_each(
      std::execution::par_unseq,
      mtx.begin(),
      mtx.end(),
      [&](auto& val){
            const tnlong n_ = &val-mtxp;
            outer_fun(n_,val);
      });
  } else {
  #endif
    tnlong c = 0;
    for (tnlong i = 0; i < N; ++i) {
      for (tnlong j = 0; j <= i; ++j) {
        inner_fun(i,j,mtxp[c++]);
      }
    }
  #if FLX_parallel_algorithm
  }
  #endif
}

FlxMtxSym::FlxMtxSym(const FlxMtx_baseS& S)
: msize(S.nrows()), mtx((msize*msize+msize)/2)
{
  tnlong i,j;
  for (i = 0; i < msize; ++i) {
    for ( j = 0; j <= i; ++j) {
      mtx[(i*i+i)/2+j] = S.operator()(i,j);
      #if FLX_DEBUG
        if ( fabs(mtx[(i*i+i)/2+j]) <= GlobalVar.TOL() ) mtx[(i*i+i)/2+j] = ZERO;
      #endif
    }
  }
}

const tdouble FlxMtxSym::operator()(const tnlong& i, const tnlong& j) const
{
  #if FLX_DEBUG
    if (i < 0 || j < 0 || i >= nrows() || j >= ncols()) {
      std::ostringstream ssV;
      ssV << "Index not within matrix.";
      throw FlxException("FlxMtxSym::operator()_2", ssV.str() );
    }
  #endif
  tnlong iV = i; tnlong jV=j;
  if (jV>iV) std::swap(iV,jV);
  return mtx[(iV*iV+iV)/2+jV];
}

const bool FlxMtxSym::isPositiveDefinite(tuint& r, const bool fixIt)
{
  return isPositiveDefinite_ext(mtx,msize,r,fixIt);
}

const bool FlxMtxSym::isPositiveDefinite_ext(flxVec& mtxe, tuint rows, tuint& r, const bool fixIt)
{
  const tdouble* mtxp = mtxe.get_tmp_vptr_const();
  const tdouble GVT = GlobalVar.TOL();
  for (size_t i = 0; i < rows; ++i) {
    if ( mtxp[(i*i+i)/2+i] <= GVT ) {
      bool stop = false;
      if (fixIt) {
        stop = true;
        for (tuint j=0;j<i && stop;++j) {
          if (mtxp[(i*i+i)/2+j] > GVT) {
            stop = false;
          }
        }
        for (tuint j=i+1;j<rows;++j) {
          if (mtxp[(j*j+j)/2+i] > GVT) {
            stop = false;
          }
        }
        if (stop) {
          mtxe[(i*i+i)/2+i] = ONE;
          std::ostringstream ssV;
          ssV << "Fixed row " << i+1 << ".";
          GlobalVar.alert.alert("FlxMtxSym::isPositiveDefinite",ssV.str());
        }
      }
      if (!stop) {
        r = i+1;
        return false;
      }
    }
  }
  return true;
}

inline void FlxMtxSym::preconding(const flxVec& r, flxVec& z, const int precn) const
{
  if (precn == 2) {
    const tdouble omega=1.2;
    tdouble sum;
    #if FLX_KAHAN_MTX_SYM
      pdouble t;
    #else
      tdouble& t = sum;
    #endif
    z.set_zero();
    for (tnlong i=0;i<msize;++i) {
      t = ZERO;
      for (tnlong j = 0;j<i;++j) t += operator()(i,j)*z[j];
      #if FLX_KAHAN_MTX_SYM
        sum = t.cast2double();
      #endif
      z[i] = (r[i]-omega*sum)/operator()(i,i);
    }
    for (tnlong i = msize;i>0;--i) {
      t = ZERO;
      for (tnlong j = i; j<msize;++j) {
        t+=operator()(i-1,j)*z[j];
      }
      #if FLX_KAHAN_MTX_SYM
        sum = t.cast2double();
      #endif
      z[i-1]-=omega*sum/operator()(i-1,i-1);
    }
  } else {
    FlxMtx_baseS::preconding(r,z,precn);
  }
}

void FlxMtxSym::MultMv(const flxVec& v, flxVec& w) const
{
  #if FLX_KAHAN_MTX_SYM
    flxpVec W(nrows());
    pdouble* Wp = W.get_tmp_vptr();
  #else
    flxVec& W = w;
    tdouble* Wp = W.get_tmp_vptr();
  #endif
  const tdouble* vp = v.get_tmp_vptr_const();
  const tdouble* mtxp = &mtx[0];
  tnlong n = 0;
  for (tnlong i = 0; i<msize;++i) {
    Wp[i] = ZERO;
    for (tnlong j=0;j<i;++j) {
      Wp[i]+=mtxp[n]*vp[j];
      Wp[j]+=mtxp[n]*vp[i];
      ++n;
    }
    Wp[i]+=mtxp[n]*vp[i];
    ++n;
  }
  #if FLX_KAHAN_MTX_SYM
    w = W;
  #endif
}

void FlxMtxSym::MultMv_slice(const flxVec& v, flxVec& w, const tuint& li, const tuint& ui) const
{
  #if FLX_DEBUG
    if (w.get_N() != nrows() || li > ui || ui >= nrows() || v.get_N() != ui-li+1 ) {
      std::ostringstream ssV;
      ssV << "ERROR";
      throw FlxException("FlxMtxSym::MultMv_slice", ssV.str() );
    }
  #endif
  size_t j, i;
  size_t n, m;
  #if FLX_KAHAN_MTX_SYM
    flxpVec W(nrows());
    pdouble* Wp = W.get_tmp_vptr();
  #else
    flxVec& W = w;
    tdouble* Wp = W.get_tmp_vptr();
  #endif
  const tdouble* vp = v.get_tmp_vptr_const();
  const tdouble* mtxp = &(mtx[0]);
  for (i = 0; i<msize;++i) {
    W[i] = ZERO;
    // lower part of the matrix
    n = (i*i+i)/2+li;
    m = ((ui<i)?ui:i);
    for (j=li;j<=m;++j) {
      Wp[i]+=mtxp[n]*vp[j-li];
      ++n;
    }
    // uper part of the matrix
    m = ((li<=i)?i+1:li);
    for (j=m;j<=ui;++j) {
      Wp[i]+=mtxp[(j*j+j)/2+i]*vp[j-li];
    }
  }
  #if FLX_KAHAN_MTX_SYM
    tVec_assemble(w,W);
  #endif
}

void FlxMtxSym::MultMv(const flxpVec& v, flxpVec& w) const
{
  const pdouble* vp = v.get_tmp_vptr_const();
  pdouble* wp = w.get_tmp_vptr();
  const tdouble* mtxp = &mtx[0];
  tnlong n = 0;
  tnlong i, j;
  for (i = 0; i<msize;++i) {
    wp[i] = ZERO;
    for (j=0;j<i;++j) {
      wp[i]+=vp[j]*mtxp[n];
      wp[j]+=vp[i]*mtxp[n];
      ++n;
    }
    wp[i]+=vp[i]*mtxp[n];
    ++n;
  }
}

FlxMtxSym* FlxMtxSym::copy()
{
  return new FlxMtxSym(*this);
}


FlxMtxSymBand* FlxMtxSymBand::copy()
{
  return new FlxMtxSymBand(*this);
}

void FlxMtxSymBand::MultMv(const flxVec& v, flxVec& w) const
{
  #if FLX_KAHAN_MTX_SYM
    flxpVec W(nrows());
    pdouble* Wp = W.get_tmp_vptr();
  #else
    flxVec& W = w;
    tdouble* Wp = W.get_tmp_vptr();
  #endif
  const tdouble* vp = v.get_tmp_vptr_const();
  const tdouble* mtxp = &mtx[0];
  tnlong n = 0;
  for (tnlong i = 0; i<msize;++i) {
    Wp[i] = ZERO;
    const tnlong jl = ((bsize>=i)?0:(i-bsize));
    const tnlong ju = ((i+bsize+1>=msize)?msize:(i+bsize+1));
    for (tnlong j=jl;j<ju;++j) {
      Wp[i] += mtxp[n++]*vp[j];
    }
  }
  #if FLX_KAHAN_MTX_SYM
    w = W;
  #endif
}

void FlxMtxSymBand::assign_LTL(const FlxMtxLTriBand& L)
{
  const tnlong N = nrows();
  tdouble sum;
  #if FLX_KAHAN_MTX_SYM
    pdouble t;
  #else
    tdouble& t = sum;
  #endif
  const tdouble* LTLmtxp = &L.mtx[0];
  tdouble* mtxp = &mtx[0];
  tnlong c = 0;
  for (tnlong i = 0; i < N; ++i) {
    const tnlong jl = ((bsize>=i)?0:(i-bsize));
    for (tnlong j = jl; j <= i; ++j) {
      t = ZERO;
      tnlong ku = j+bsize+1; if (ku>N) ku=N;
      for (tnlong k = i; k < ku; ++k) {
        tnlong m = L.countUp2RDiag(k);
        t += LTLmtxp[m-(k-j)]*LTLmtxp[m-(k-i)];
      }
      #if FLX_KAHAN_MTX_SYM
        sum = t.cast2double();
      #endif
      mtxp[c++] = sum;
    }
  }
}

void FlxMtxSymBand::Invert()
{
  throw FlxException_NotImplemented("FlxMtxSymBand::Invert");
}

FlxMtxSymBand& FlxMtxSymBand::operator+=(const FlxMtxSymBand& m)
{
  #if FLX_DEBUG
    if (msize != m.msize || bsize != m.bsize ) {
      std::ostringstream ssV;
      ssV << "Matrices can not be added.";
      throw FlxException("FlxMtxSymBand::operator+=_1", ssV.str() );
    }
  #endif
  mtx+=m.mtx;
  return *this;
}

FlxMtxSymBand& FlxMtxSymBand::operator+=(const FlxMtxDiag& m)
{
  #if FLX_DEBUG
    if (msize != m.nrows() ) {
      std::ostringstream ssV;
      ssV << "Matrices can not be added.";
      throw FlxException("FlxMtxSymBand::operator+=_2", ssV.str() );
    }
  #endif
  tnlong c = 0;
  for (tnlong i=0;i<msize;++i) {
    c += (bsize<=i)?bsize:i;
    mtx[c++]+=m(i,i);
    c += (bsize<=(msize-i-1))?bsize:(msize-i-1);
  }
  return *this;
}

void FlxMtxSymBand::add_mtx(const FlxMtxSymBand& m, const tdouble& s)
{
  #if FLX_DEBUG
    if (msize != m.msize || bsize != m.bsize ) {
      std::ostringstream ssV;
      ssV << "Matrices can not be added.";
      throw FlxException("FlxMtxSymBand::operator+=", ssV.str() );
    }
  #endif
  tdouble* tp = &(mtx[0]);
  const tdouble* const mp = &(m.mtx[0]);
  const size_t N = mtx.size();
  for (tnlong i=0;i<N;++i) {
    tp[i] += s*mp[i];
  }
}

void FlxMtxSymBand::add_mtx(const FlxMtxDiag& m, const tdouble& s)
{
  #if FLX_DEBUG
    if (msize != m.nrows() ) {
      std::ostringstream ssV;
      ssV << "Matrices can not be added.";
      throw FlxException("FlxMtxSymBand::operator+=_2", ssV.str() );
    }
  #endif
  tnlong c = 0;
  for (tnlong i=0;i<msize;++i) {
    c += (bsize<=i)?bsize:i;
    mtx[c++]+=m(i,i)*s;
    c += (bsize<=(msize-i-1))?bsize:(msize-i-1);
  }
}

void FlxMtxSymBand::preconding(const flxVec& r, flxVec& z, const int precn) const
{
  throw FlxException_NotImplemented("FlxMtxSymBand::preconding");
}

const tnlong FlxMtxSymBand::countUp2Row(const tnlong R) const
{
  tnlong j = 0;
  for (tnlong i=0;i<R;++i) {
    ++j;
    j += (bsize<=i)?bsize:i;
    j += (bsize<=(msize-i-1))?bsize:(msize-i-1);
  }
  return j;
}

const tdouble FlxMtxSymBand::operator()(const tnlong& i, const tnlong& j) const
{
  #if FLX_DEBUG
    if (i < 0 || j < 0 || i >= nrows() || j >= ncols()) {
      std::ostringstream ssV;
      ssV << "Index not within matrix.";
      throw FlxException("FlxMtxSymBand::operator()_2", ssV.str() );
    }
  #endif
  tnlong iV = i; tnlong jV=j;
  if (jV>iV) std::swap(iV,jV);
  if ((iV-jV) > bsize) return ZERO;
  return mtx[countUp2Row(jV)+((bsize<=jV)?bsize:jV)+(iV-jV)];
}

inline void FlxMtxSymBand::add_value(const tnlong& i, const tnlong& j, const tdouble& v)
{
  #if FLX_DEBUG
    if (i < 0 || j < 0 || i >= nrows() || j >= ncols()) {
      std::ostringstream ssV;
      ssV << "Index not within matrix.";
      throw FlxException("FlxMtxSymBand::add_value_1", ssV.str() );
    }
  #endif
  tnlong iV = i; tnlong jV=j;
  if (jV>iV) std::swap(iV,jV);
  if (iV-jV > bsize) {
    std::ostringstream ssV;
    ssV << "Index not within writeable region of band-matrix.";
    throw FlxException("FlxMtxSymBand::add_value_2", ssV.str() );
  }
  mtx[countUp2Row(jV)+((bsize<=jV)?bsize:jV)+(iV-jV)] += v;
  if (iV==jV) return;
  mtx[countUp2Row(iV)+((bsize<=iV)?(bsize-(iV-jV)):jV)] += v;
}

const bool FlxMtxSymBand::isPositiveDefinite(tuint& r, const bool fixIt)
{
  throw FlxException_NotImplemented("FlxMtxSymBand::isPositiveDefinite");
}


void FlxMtxDiag::add_value(const tnlong& ii, const tnlong& jj, const tdouble& v)
{
  if (v==ZERO) return;
  tnlong i=ii;
  tnlong j=jj;
  #if FLX_DEBUG
    if (i < 0 || j < 0 || i >= nrows() || j >= nrows()) {
      std::ostringstream ssV;
      ssV << "Index not within matrix.";
      throw FlxException("FlxMtxDiag::add_value_1", ssV.str() );
    }
  #endif
  if (i == j) {
    mtx[i]+=v;
  } else {
    std::ostringstream ssV;
    ssV << "Not allowed to add value '" << v << "' at this position (" << i << ", " << j << ") - diagonal matrix.";
    throw FlxException("FlxMtxDiag::add_value_2", ssV.str() );
  }
}

const bool FlxMtxDiag::isPositiveDefinite(tuint& r, const bool fixIt)
{
  if (mtx.min() <= GlobalVar.TOL() ) return false;
  else return true;
}

FlxMtx_baseS* FlxMtxDiag::get_Inverse()
{
  FlxMtxDiag* Inv = new FlxMtxDiag(msize);
  tdouble* Invmtxp = &Inv->mtx[0];
  const tdouble* mtxp = &mtx[0];
  for (tnlong i = 0; i < msize; ++i) {
    Invmtxp[i] = 1/mtxp[i];
  }
  return Inv;
}

FlxMtxDiag::FlxMtxDiag(const FlxMtx_baseS& mb)
: msize(mb.ncols()),mtx(mb.ncols())
{
  for (tnlong i = 0; i < msize; ++i) {
    mtx[i] = mb.operator()(i,i);
  }
}

FlxMtx_baseS* FlxMtxIdentity::get_Inverse()
{
  return new FlxMtxIdentity(msize);
}

FlxMtxIdentity& FlxMtxIdentity::operator*=(const tdouble& s)
{
  std::ostringstream ssV;
  ssV << "Invalid operation.";
  throw FlxException("FlxMtxIdentity::operator*=", ssV.str() );
}

FlxMtxIdentity& FlxMtxIdentity::operator+=(const FlxMtx_base& m)
{
  std::ostringstream ssV;
  ssV << "Invalid operation.";
  throw FlxException("FlxMtxIdentity::operator+=", ssV.str() );
}

void FlxMtxIdentity::add_value(const tnlong& i, const tnlong& j, const tdouble& v)
{
  std::ostringstream ssV;
  ssV << "Invalid operation.";
  throw FlxException("FlxMtxIdentity::add_value", ssV.str() );
}

void FlxMtxSparsSym::set_value(const tnlong& ii, const tnlong& jj, const tdouble& v)
{
  tnlong i=ii;
  tnlong j=jj;
  #if FLX_DEBUG
    if (i < 0 || j < 0 || i >= nrows() || j >= nrows()) {
      std::ostringstream ssV;
      ssV << "Index not within matrix.";
      throw FlxException("FlxMtxSparsSym::set_value_1", ssV.str() );
    }
  #endif
  if (i == j) {
    sa[i]=v;
    return;
  }
  if (j > i) std::swap(i, j);
  for (tnlong k = ija[i];k<ija[i+1];k++) {
    if (ija[k]==j) {
      sa[k]=v;
      return;
    }
    else if (ija[k]>j) {
      break;
    }
  }
  std::ostringstream ssV;
  ssV << "Not allowed to set value '" << v << "' at this position (" << i << ", " << j << ") - sparse matrix.";
  throw FlxException("FlxMtxSparsSym::set_value_2", ssV.str() );
}

void FlxMtxSparsSym::add_value(const tnlong& ii, const tnlong& jj, const tdouble& v)
{
  if (v==ZERO) return;
  tnlong i=ii;
  tnlong j=jj;
  #if FLX_DEBUG
    if (i < 0 || j < 0 || i >= nrows() || j >= nrows()) {
      std::ostringstream ssV;
      ssV << "Index not within matrix.";
      throw FlxException("FlxMtxSparsSym::add_value_1", ssV.str() );
    }
  #endif
  if (i == j) {
    sa[i]+=v;
    return;
  }
  if (j > i) std::swap(i, j);
  for (tnlong k = ija[i];k<ija[i+1];k++) {
    if (ija[k]==j) {
      sa[k]+=v;
      return;
    }
    else if (ija[k]>j) {
      break;
    }
  }
  std::ostringstream ssV;
  ssV << "Not allowed to add value '" << v << "' at this position (" << i << ", " << j << ") - sparse matrix.";
  throw FlxException("FlxMtxSparsSym::add_value_2", ssV.str() );
}

void FlxMtxDiag::MultMv(const flxVec& v, flxVec& w) const
{
  const tdouble* mtxp = &mtx[0];
  const tdouble* vp = v.get_tmp_vptr_const();
  tdouble* wp = w.get_tmp_vptr();
  for (tnlong i = 0; i < msize; ++i) {
    wp[i] = vp[i]*mtxp[i];
  }
}

void FlxMtxDiag::MultMv(const flxpVec& v, flxpVec& w) const
{
  const tdouble* mtxp = &mtx[0];
  const pdouble* vp = v.get_tmp_vptr_const();
  pdouble* wp = w.get_tmp_vptr();
  for (tnlong i = 0; i < msize; ++i) {
    wp[i] = vp[i]*mtxp[i];
  }
}

void FlxMtxDiag::MultInv(const flxVec& v, flxVec& w)
{
  const tdouble* mtxp = &mtx[0];
  const tdouble* vp = v.get_tmp_vptr_const();
  tdouble* wp = w.get_tmp_vptr();
  for (tnlong i = 0; i < msize; ++i) {
    wp[i] = vp[i]/mtxp[i];
  }
}

void FlxMtxSparsSym::MultMv(const flxVec& v, flxVec& w) const
{
  const tnlong N = nrows();
  #if FLX_KAHAN_MTX_SSYM
    flxpVec W(N);
    pdouble* Wp = W.get_tmp_vptr();
  #else
    flxVec& W = w;
    tdouble* Wp = W.get_tmp_vptr();
  #endif
  const tdouble* vp = v.get_tmp_vptr_const();
  for (tnlong i=0;i<N;++i) {
    Wp[i] = sa[i]*vp[i];
  }
  for (tnlong i=0;i<N;++i) {
    for (tnlong k=ija[i];k<ija[i+1];++k) {
      Wp[i] += sa[k]*vp[ija[k]];
      Wp[ija[k]]+=sa[k]*vp[i];
    }
  }
  #if FLX_KAHAN_MTX_SSYM
    w = W;
  #endif
}

void FlxMtxSparsSym::MultMv(const flxpVec& v, flxpVec& w) const
{
  pdouble* wp = w.get_tmp_vptr();
  const pdouble* vp = v.get_tmp_vptr_const();
  const tnlong N = nrows();
  tnlong k,i;
  for (i=0;i<N;++i) {
    wp[i] = vp[i]*sa[i];
  }
  for (i=0;i<N;++i) {
    for (k=ija[i];k<ija[i+1];++k) {
      wp[i] += vp[ija[k]]*sa[k];
      wp[ija[k]]+=vp[i]*sa[k];
    }
  }
}

FlxMtx_base& FlxMtxSparsSym::operator*=(const tdouble& s)
{
  const tnlong N=ija[nrows()];
  for (tnlong i=0;i<N;++i) sa[i]*=s;
  return *this;
}

FlxMtxSparsSym::FlxMtxSparsSym(const FlxMtxSym& Mtx)
{
  const tnlong n = Mtx.nrows();
  const flxVec& mtx = Mtx.get_mtx_flxVec();
  const tdouble* mtxp = mtx.get_tmp_vptr_const();
  const tnlong mtxs = mtx.get_N();
  tnlong nmax = n+1;
  // find largest entry
    tdouble mtx_max = ONE;
    if (mtxs>0) {        // only if size is larger than 0
      mtx_max = fabs(mtxp[0]);
      for (size_t i=1;i<mtxs;++i) {
        if (fabs(mtxp[i])>mtx_max) mtx_max = fabs(mtxp[i]);
      }
    }
  const tdouble TOL = GlobalVar.TOL()*mtx_max;
  tnlong i, j;
  for (i = 0; i < mtxs; ++i) {
    if (fabs(mtxp[i])>TOL ) ++nmax;
  }
  for (i=1; i<=n; ++i) {
    if (fabs(mtxp[(i*i+i)/2-1])>TOL ) --nmax;
  }

  sa = new tdouble[nmax];
  ija = new tnlong[nmax];
  ija[0] = n+1;
  tnlong k=n;
  tnlong c=0;
  for (i=0;i<n;++i) {
    for (j=0;j<i;++j) {
      if (fabs(mtxp[c])>TOL ) {
        ++k;
        sa[k]=mtxp[c];
        ija[k]=j;
      }
      ++c;
    }
    ija[i+1]=k+1;
    sa[i] = mtxp[c++];
  }
  
}

FlxMtxSparsSym::~FlxMtxSparsSym()
{
  delete[] sa;
  delete[] ija;
}

const tdouble FlxMtxSparsSym::operator()(const tnlong& ii, const tnlong& jj) const
{
  tnlong i=ii;
  tnlong j=jj;
  #if FLX_DEBUG
    if (i < 0 || j < 0 || i >= nrows() || j >= nrows()) {
      std::ostringstream ssV;
      ssV << "Index not within matrix.";
      throw FlxException("FlxMtxSparsSym::operator()_2", ssV.str() );
    }
  #endif
  if (i == j) return sa[i];
  if (j > i) std::swap(i, j);
  for (tnlong k = ija[i];k<ija[i+1];k++) {
    if (ija[k]==j) {
      return sa[k];
    }
    else if (ija[k]>j) {
      return ZERO;
    }
  }
  return ZERO;
}

const bool FlxMtxSparsSym::isPositiveDefinite(tuint& r, const bool fixIt)
{
  const tdouble GVT = GlobalVar.TOL();
  for (tnlong i = nrows(); i > 0; --i) {
    if ( sa[i-1] <= GVT ) { 
      bool stop = false;
      if (fixIt) {
        stop = true;
        if (ija[i-1]!=ija[i]) stop = false;
        for (tuint j=i;j<nrows();++j) {
          if (operator()(j,i-1) > GVT) {
            stop = false;
          }
        }
        if (stop) {
          sa[i-1] = ONE;
          std::ostringstream ssV;
          ssV << "Fixed row " << i+1 << ".";
          GlobalVar.alert.alert("FlxMtxSparsSym::isPositiveDefinite",ssV.str());
        }
      }
      if (!stop) {
        r = i;
        return false;
      }
    }
  }
  return true;
}

void FlxMtxSparsSym::output_Mtx(std::ostream& sout) const
{
  for (tnlong i =0;i<nrows();i++) {
    tnlong k=0;
    for (tnlong j = ija[i]; j<ija[i+1];j++) {
      while (k<ija[j]) {
        sout << GlobalVar.Double2String(0.,false,-1,0) << " ";
        k++;
      }
      sout << GlobalVar.Double2String(sa[j],false,-1,0) << " ";
      k++;
    }
    while (k<i) {
      sout << GlobalVar.Double2String(0.,false,-1,0) << " ";
      k++;
    }
    sout << GlobalVar.Double2String(sa[i],false,-1,0);
    sout << std::endl;
  }
}

void FlxMtxSparsSym::preconding(const flxVec& r, flxVec& z, const int precn) const
{
  tdouble sum = 0;
  #if FLX_KAHAN_MTX_SYM
    pdouble t;
  #else
    tdouble& t = sum;
  #endif
  if (precn == 2) {
    const tnlong N = nrows();
    tnlong i, j;
    const tdouble omega=1.2;
    for (i=0;i<N;++i) {
      t = 0;
      for (j = ija[i];j<ija[i+1];++j) {
        t += sa[j]*z[ija[j]];
      }
      #if FLX_KAHAN_MTX_SYM
        sum = t.cast2double();
      #endif
      z[i] = (r[i]-omega*sum)/sa[i];
    }
    for (i = N;i>=1;i--) {
      t = 0;
      for (j = i; j<N;++j) {
        t+=operator()(j,i-1)*z[j];
      }
      #if FLX_KAHAN_MTX_SYM
        sum = t.cast2double();
      #endif
      z[i-1]-=omega*sum/sa[i-1];
    }
  } else {
    FlxMtx_baseS::preconding(r,z,precn);
  }
}

FlxMtxLTri::FlxMtxLTri(const FlxMtxSparsLTri& Msp)
: msize(Msp.nrows()), mtx((msize*msize+msize)/2), is_ldl(false)
{
  const tnlong n = nrows();
  const tnlong* ija = Msp.ija;
  const tdouble* sa = Msp.sa;
  for (tnlong i=0;i<n;++i) {
    const tnlong m = (i*i+i)/2;
    for (tnlong j=ija[i];j<ija[i+1];++j) {
      mtx[m+ija[j]] = sa[j];
    }
    mtx[m+i] = sa[i];
  }
}

FlxMtxLTri* FlxMtxLTri::copy()
{
  return new FlxMtxLTri(*this);
}

void FlxMtxLTri::MultMv(const flxVec& v, flxVec& w) const
{
  if (is_ldl) throw FlxException_NotImplemented("FlxMtxLTri::MultMv");
  size_t n = mtx.get_N();
  #if FLX_KAHAN_MTX_LTRI
    pdouble t;
  #else
    tdouble t;
  #endif
  const tdouble* vp = v.get_tmp_vptr_const();
  tdouble* wp = w.get_tmp_vptr();
  const tdouble* mtxp = &mtx[0];
  for (tnlong i = msize; i>0;--i) {
    t = ZERO;
    for (tnlong j=i;j>0;--j) {
      t+=mtxp[n-1]*vp[j-1];
      --n;
    }
    #if FLX_KAHAN_MTX_LTRI
      wp[i-1]=t.cast2double();
    #else
      wp[i-1]=t;
    #endif
  }
}

void FlxMtxLTri::MultInv(const flxVec& v, flxVec& w, bool ldl_sqrt)
{
  const tnlong n = nrows();
  tdouble sum;
  #if FLX_KAHAN_MTX_SLTRI
    pdouble t;
  #else
    tdouble& t=sum;
  #endif
  const tdouble* vp = v.get_tmp_vptr_const();
  tdouble* wp = w.get_tmp_vptr();
  const tdouble* mtxp = &mtx[0];
  tnlong i, j;
  tnlong c = 0;
  if (is_ldl && ldl_sqrt) {
    for (i=0;i<n;++i) {
      t=ZERO;
      for (j=0;j<i;++j) {
        t += mtxp[c++]*sqrt(mtxp[(j*j+j)/2+j])*wp[j];
      }
      #if FLX_KAHAN_MTX_SLTRI
        sum = t.cast2double();
      #endif
      wp[i]=(vp[i]-sum)/sqrt(mtxp[c++]);
    }
  } else {
    for (i=0;i<n;++i) {
      t=ZERO;
      for (j=0;j<i;++j) {
        t += mtxp[c++]*wp[j];
      }
      #if FLX_KAHAN_MTX_SLTRI
        sum = t.cast2double();
      #endif
      wp[i]=(vp[i]-sum);
      if (is_ldl) {
        c++;
      } else {
        wp[i]/=mtxp[c++];
      };
    }
  }
}

void FlxMtxLTri::TransMultInv(const flxVec& v, flxVec& w)
{
  tdouble sum;
  #if FLX_KAHAN_MTX_SLTRI
    pdouble t;
  #else
    tdouble& t=sum;
  #endif
  const tdouble* vp = v.get_tmp_vptr_const();
  tdouble* wp = w.get_tmp_vptr();
  const tdouble* mtxp = &mtx[0];
  for (tnlong i = msize; i>0;--i) {
    t=ZERO;
    for (tnlong j=i+1;j<=msize;++j) {
      t += mtxp[(j*j-j)/2+i-1]*wp[j-1];
    }
    #if FLX_KAHAN_MTX_SLTRI
      sum = t.cast2double();
    #endif
    if (is_ldl) {
      wp[i-1]=vp[i-1]/mtxp[(i*i-i)/2+i-1]-sum;
    } else {
      wp[i-1]=(vp[i-1]-sum)/mtxp[(i*i-i)/2+i-1];
    }
  }
}


const tdouble FlxMtxLTri::det_log() const
{
  const tdouble* mtxp = &mtx[0];
  tdouble r = ZERO;
  for (tnlong i=0;i<msize;++i) {
    r += log(mtxp[(i*i+i)/2+i]);
  }
  if (is_ldl) {
    r /= 2;
  }
  return r;
}

FlxMtxLTri& FlxMtxLTri::operator*=(const tdouble& s) {
  if (is_ldl) throw FlxException_NotImplemented("FlxMtxLTri::operator*=");
  mtx*=s;
  return *this;
}

const tdouble FlxMtxLTri::operator()(const tnlong& i, const tnlong& j) const
{
  if (is_ldl) {
    throw FlxException_NotImplemented("FlxMtxLTri::operator()");
  }
  #if FLX_DEBUG
    if (i < 0 || j < 0 || i >= nrows() || j >= ncols()) {
      std::ostringstream ssV;
      ssV << "Index not within matrix.";
      throw FlxException("FlxMtxLTri::operator()_1", ssV.str() );
    }
  #endif
  if (j>i) return ZERO;
  else return mtx[(i*i+i)/2+j];
}

void FlxMtxLTri::add_value(const tnlong& i, const tnlong& j, const tdouble& v)
{
  if (is_ldl) throw FlxException_NotImplemented("FlxMtxLTri::add_value_0");
  #if FLX_DEBUG
    if (i < 0 || j < 0 || i >= nrows() || j >= ncols()) {
      std::ostringstream ssV;
      ssV << "Index not within matrix.";
      throw FlxException("FlxMtxLTri::add_value_1", ssV.str() );
    }
  #endif
  if (j<=i) mtx[(i*i+i)/2+j]+=v;
  else {
    std::ostringstream ssV;
      ssV << "Not allowed to add value at this position. (" << i << ", " << j << ")";
      throw FlxException("FlxMtxLTri::add_value_2", ssV.str() );
  }
}

FlxMtxLTri& FlxMtxLTri::CholeskyDec(FlxMtxSym& Sm, const bool do_ldl)
{
  #if FLX_DEBUG
    if (Sm.nrows()!=nrows()) throw FlxException_Crude("FlxMtxLTri::CholeskyDec_0");
  #endif
  return CholeskyDec(Sm.get_mtx_flxVec(),do_ldl);
}

FlxMtxLTri& FlxMtxLTri::CholeskyDec(const flxVec& mtxV, const bool do_ldl)
{
  is_ldl = do_ldl;
  mtx = mtxV;
  tdouble* mtxp = &mtx[0];
  const tnlong N = nrows();
  tdouble sum;
  #if FLX_KAHAN_MTX_LTRI
    pdouble t;
  #else
    tdouble& t=sum;
  #endif
  if (do_ldl) {
    for (tnlong i = 0; i < N; ++i) {
      for (tnlong j = 0; j < i; ++j) {
        t = mtxp[(i*i+i)/2+j];
        for (tnlong k=0;k<j;k++) {
          t-=mtxp[(i*i+i)/2+k]*mtxp[(j*j+j)/2+k]*mtxp[(k*k+k)/2+k];
        }
        #if FLX_KAHAN_MTX_LTRI
          sum = t.cast2double();
        #endif
        mtxp[(i*i+i)/2+j] = sum/mtxp[(j*j+j)/2+j];
      }
      t = mtxp[(i*i+i)/2+i];
      for (tnlong k=0;k<i;++k) {
        const tdouble d = mtxp[(i*i+i)/2+k];
        t-=pow2(d)*mtxp[(k*k+k)/2+k];
      }
      #if FLX_KAHAN_MTX_LTRI
        sum = t.cast2double();
      #endif
      if (sum<=ZERO) {
        std::ostringstream ssV;
        ssV << "Matrix is not positiv definite.";
        std::ostringstream ssV2;
        ssV2 << "  line = " << i+1 << std::endl;
        ssV2 << "  sum  = " << sum << std::endl;
        //ssV2 << "  min  = " << this->min() << std::endl;
        //ssV2 << "  max  = " << this->max() << std::endl;
        #if FLX_DEBUG
          //ssV2 << "Correlation matrix:" << std::endl;
          //this->output_Mtx(ssV2);
        #endif
          throw FlxException("FlxMtxLTri::CholeskyDec_1", ssV.str(), ssV2.str() );
      } else {
        // GlobalVar.slogcout(1) << " Frieda " << i << " " << sum << std::endl;
        mtxp[(i*i+i)/2+i] = sum;
      }
    }
  } else {
    for (tnlong i = 0; i < N; ++i) {
      for (tnlong j = 0; j < i; ++j) {
        t = mtxp[(i*i+i)/2+j];
        for (tnlong k=0;k<j;k++) {
          t-=mtxp[(i*i+i)/2+k]*mtxp[(j*j+j)/2+k];
        }
        // t = ZERO;
        // for (tnlong k=0;k<j;k++) {
        //   t+=mtxp[(i*i+i)/2+k]*mtxp[(j*j+j)/2+k];
        // }
        // t -= mtxp[(i*i+i)/2+j];
        // t *= -ONE;
        #if FLX_KAHAN_MTX_LTRI
          sum = t.cast2double();
        #endif
        mtxp[(i*i+i)/2+j] = sum/mtxp[(j*j+j)/2+j];
      }
      t = mtxp[(i*i+i)/2+i];
      for (tnlong k=0;k<i;++k) {
        const tdouble d = mtxp[(i*i+i)/2+k];
        t-=pow2(d);
      }
      // t = ZERO;
      // for (tnlong k=0;k<i;++k) {
      //   const tdouble d = mtxp[(i*i+i)/2+k];
      //   t+=pow2(d);
      // }
      // t -= mtxp[(i*i+i)/2+i];
      // t *= -ONE;
      #if FLX_KAHAN_MTX_LTRI
        sum = t.cast2double();
      #endif
      if (sum<=ZERO) {
        std::ostringstream ssV;
        ssV << "Matrix is not positiv definite.";
        std::ostringstream ssV2;
        ssV2 << "  line = " << i+1 << std::endl;
        ssV2 << "  sum  = " << sum << std::endl;
        ssV2 << "  min  = " << this->min() << std::endl;
        ssV2 << "  max  = " << this->max() << std::endl;
        #if FLX_DEBUG
          //ssV2 << "Correlation matrix:" << std::endl;
          //this->output_Mtx(ssV2);
        #endif
        throw FlxException("FlxMtxLTri::CholeskyDec_1", ssV.str(), ssV2.str() );
      } else {
        // GlobalVar.slogcout(1) << " Frieda " << i << " " << sum << std::endl;
        mtxp[(i*i+i)/2+i] = sqrt(sum);
      }
    }
  }
  return *this;
}

FlxMtxLTri& FlxMtxLTri::CholeskyDec(FlxMtxSparsSym& Sm)
{
  is_ldl = false;
  const tnlong N = nrows();
  tdouble sum, d, d1;
  #if FLX_KAHAN_MTX_LTRI
    pdouble t;
    pdouble sum2;
  #else
    tdouble& t=sum;
    tdouble sum2;
  #endif
  tdouble* mtxp = &mtx[0];
  tnlong j, k, l, ii, jj;
  for (tnlong i = 0; i < N; ++i) {
    l = Sm.ija[i];
    ii = (i*i+i)/2;
    sum2=ZERO;
    for (j = 0; j < i; ++j) {
      jj = (j*j+j)/2;
      if (l<Sm.ija[i+1] && j==Sm.ija[l] ) {
        t = Sm.sa[l];
        ++l;
      } else {
        t = ZERO;
      }
      for (k=0;k<j;++k) {
        d=mtxp[ii+k];
        d1=mtxp[jj+k];
        t-=d*d1;
      }
      #if FLX_KAHAN_MTX_LTRI
        sum = t.cast2double();
      #endif
      sum/=mtxp[jj+j];
      mtxp[ii+j] = sum;
      sum2+=pow2(sum);
    }
    // calculate diagonal entries
    #if FLX_KAHAN_MTX_LTRI
      sum = Sm.sa[i] - sum2.cast2double();
    #else
      sum = Sm.sa[i] - sum2;
    #endif
    if (sum<=ZERO) {
      std::ostringstream ssV;
      ssV << "Matrix is not positiv definite.";
      std::ostringstream ssV2;
      ssV2 << "  sum = " << sum << std::endl;
      #if FLX_DEBUG
        ssV2 << "Correlation matrix:" << std::endl;
        this->output_Mtx(ssV2);
      #endif
      throw FlxException("FlxMtxLTri::CholeskyDec_2", ssV.str(), ssV2.str() );
    } else {
      mtxp[ii+i] = sqrt(sum);
    }
  }
  return *this;
}


FlxMtxLTri& FlxMtxLTri::Invert()
{
  if (is_ldl) throw FlxException_NotImplemented("FlxMtxLTri::Invert");
  const tnlong n = this->ncols();
  tdouble* mtxp = &mtx[0];
  tnlong i,j, m;
  size_t k = 0;
  tdouble a;
  #if FLX_KAHAN_MTX_LTRI
    pdouble t;
  #else 
    tdouble& t=a;
  #endif
  for (i = 0; i < n; ++i) {                // rows
    for (j = 0; j < i; ++j) {                // columns
      t = ZERO;
      for (m = j; m < i; ++m) {                // summation
        t+=mtxp[(i*i+i)/2+m]*mtxp[(m*m+m)/2+j];
      }
      #if FLX_KAHAN_MTX_LTRI
        a = t.cast2double();
      #endif
      mtxp[k]=-a/mtxp[(i*i+i)/2+i];
      ++k;
    }
    mtxp[k]=1/(mtxp[k]);
    ++k;
  }
  
  return *this;  
}

void FlxMtxLTri::TransMultVec(flxVec& v)
{
  if (is_ldl) throw FlxException_NotImplemented("FlxMtxLTri::TransMultVec");
  #if FLX_KAHAN_MTX_LTRI
    pdouble a;
  #else
    tdouble a;
  #endif
  const tuint vs = v.get_N();
  const tdouble* mtxp = &mtx[0];
  tdouble* vp = &v[0];
  for (tuint i = 0; i < vs; i++) {
    a = ZERO;
    for (tuint j = i; j < vs; j++) {
      a+=mtxp[(j*j+j)/2+i]*vp[j];
    }
    #if FLX_KAHAN_MTX_LTRI
      vp[i]=a.cast2double();
    #else
      vp[i]=a;
    #endif
  }
}

const size_t FlxMtxLTri::count_nan() const
{
  return mtx.count_nan();
}


FlxMtxLTriBand* FlxMtxLTriBand::copy()
{
  return new FlxMtxLTriBand(*this);
}

const tnlong FlxMtxLTriBand::countUp2RDiag(const tnlong R) const
{
  if (R==0) return 0;
  const tnlong bs = ((bsize<=R-1)?bsize:R-1);
  return (bs+1)*R - (bs*bs+bs)/2 + ((bsize<=R)?bsize:R);
}

FlxMtxLTriBand& FlxMtxLTriBand::CholeskyDec(FlxMtxSymBand& Sm)
{
  const tnlong N = nrows();
  #if FLX_DEBUG
    if (N!=Sm.nrows() || bsize!=Sm.get_bsize()) {
      throw FlxException_Crude("FlxMtxLTriBand::CholeskyDec_0");
    }
  #endif
  // assign Sm to this
    tdouble* mtxp = &mtx[0];
    const tdouble* Smp = &((Sm.get_mtx_tVec())[0]);
    tnlong cS = 0;
    tnlong c = 0;
    for (tnlong i=0;i<N;++i) {
      const tnlong jl = ((bsize>=i)?0:(i-bsize));
      for (tnlong j=jl;j<=i;++j) {
        mtxp[c++] = Smp[cS++];
      }
      cS += (bsize<=(msize-i-1))?bsize:(msize-i-1);
    }
  tdouble sum, d;
  #if FLX_KAHAN_MTX_LTRI
    pdouble t;
  #else
    tdouble& t=sum;
  #endif
  c = 0;
  for (tnlong i = 0; i < N; ++i) {
    const tnlong jl = ((bsize>=i)?0:(i-bsize));
    for (tnlong j = jl; j < i; ++j) {
      const tnlong c2 = countUp2RDiag(j);
      t = mtxp[c];
      for (tnlong k=j-jl;k>=1;--k) {
        t -= mtxp[c-k]*mtxp[c2-k];
      }
      #if FLX_KAHAN_MTX_LTRI
        sum = t.cast2double();
      #endif
      mtxp[c++] = sum/mtxp[c2];
    }
    t = mtxp[c];
    for (tnlong k=i-jl;k>=1;--k) {
      d = mtxp[c-k];
      t-=pow2(d);
    }
    #if FLX_KAHAN_MTX_LTRI
      sum = t.cast2double();
    #endif
    if (sum<=ZERO) {
      std::ostringstream ssV;
      ssV << "Matrix is not positiv definite. (sum = " << sum << ") [line: " << i+1 << "]";
      #if FLX_DEBUG
        std::ostringstream ssV2;
        ssV2 << "Matrix:" << std::endl;
        this->output_Mtx(ssV2);
        throw FlxException("FlxMtxLTriBand::CholeskyDec_1dbg", ssV.str(), ssV2.str() );
      #else
        throw FlxException("FlxMtxLTriBand::CholeskyDec_1", ssV.str() );
      #endif
    } else {
      mtxp[c++] = sqrt(sum);
    }
  }
  return *this;
}

FlxMtxLTriBand& FlxMtxLTriBand::Invert()
{
  throw FlxException_NotImplemented("FlxMtxLTriBand::Invert");
}

void FlxMtxLTriBand::MultMv(const flxVec& v, flxVec& w) const
{
  throw FlxException_NotImplemented("FlxMtxLTriBand::MultMv");
}

void FlxMtxLTriBand::MultInv(const flxVec& v, flxVec& w)
{
  const tnlong n = nrows();
  tdouble sum;
  #if FLX_KAHAN_MTX_SLTRI
    pdouble t;
  #else
    tdouble& t=sum;
  #endif
  const tdouble* vp = v.get_tmp_vptr_const();
  tdouble* wp = w.get_tmp_vptr();
  const tdouble* mtxp = &mtx[0];
  tnlong c = 0;
  for (tnlong i=0;i<n;++i) {
    t=ZERO;
    const tnlong jl = ((bsize>=i)?0:(i-bsize));
    for (tnlong j=jl;j<i;++j) {
      t += mtxp[c++]*wp[j];
    }
    #if FLX_KAHAN_MTX_SLTRI
      sum = t.cast2double();
    #endif
    wp[i]=(vp[i]-sum)/mtxp[c++];
  }
}

void FlxMtxLTriBand::MultInvTrans(const flxVec& v, flxVec& w)
{
  const tnlong n = nrows();
  tdouble sum;
  #if FLX_KAHAN_MTX_SLTRI
    pdouble t;
  #else
    tdouble& t=sum;
  #endif
  const tdouble* vp = v.get_tmp_vptr_const();
  tdouble* wp = w.get_tmp_vptr();
  const tdouble* mtxp = &mtx[0];
  tnlong c = mtx.size();
  for (tnlong i=0;i<n;++i) {
    const tnlong ii = n-1-i;
    t=ZERO;
    const tnlong ju = n-((bsize>i)?0:(i-bsize));
    for (tnlong j=ii+1;j<ju;++j) {
      t += mtxp[--c]*wp[j];
    }
    #if FLX_KAHAN_MTX_SLTRI
      sum = t.cast2double();
    #endif
    wp[ii]=(vp[ii]-sum)/mtxp[--c];
  }
}

void FlxMtxLTriBand::add_value(const tnlong& i, const tnlong& j, const tdouble& v)
{
  throw FlxException_NotImplemented("FlxMtxLTriBand::add_value");
}

void FlxMtxLTriBand::TransMultVec(flxVec& v)
{
  throw FlxException_NotImplemented("FlxMtxLTriBand::TransMultVec");
}


const tdouble FlxMtxLTriBand::operator()(const tnlong& i, const tnlong& j) const
{
  #if FLX_DEBUG
    if (i < 0 || j < 0 || i >= nrows() || j >= ncols()) {
      std::ostringstream ssV;
      ssV << "Index not within matrix.";
      throw FlxException("FlxMtxLTriBand::operator()_1", ssV.str() );
    }
  #endif
  if ((i-j) > bsize) return ZERO;
  return mtx[countUp2RDiag(i)-(i-j)];
}

FlxMtxSpars::FlxMtxSpars(const tnlong nmax)
:sa(new tdouble[nmax]),ija(new tnlong[nmax])
{

}

FlxMtxSpars::FlxMtxSpars(const FlxMtxSpars& mtx): FlxMtx_base(mtx)
{
  const tnlong nmax = mtx.ija[mtx.ija[0]-1];
  sa = new tdouble[nmax];
  ija = new tnlong[nmax];
  const tdouble* mtxsap = &mtx.sa[0];
  const tnlong* mtxijap = &mtx.ija[0];
  for (tnlong i=0;i<nmax;++i) {
    sa[i] = mtxsap[i];
    ija[i] = mtxijap[i];
  }
}

FlxMtxSpars::FlxMtxSpars(const FlxMtx& mtx): FlxMtx_base(mtx)
{
  const tnlong n = mtx.nrows();
  const tnlong m = mtx.ncols();
  tnlong nmax = n+1;
  const tdouble* mtxmtxp = mtx.mtx.get_tmp_vptr_const();
  const tdouble GVT = GlobalVar.TOL();
  tnlong j,i;
  for (i = 0; i<n;++i) {
    for (j = 0; j<m;++j) {
      if (fabs(mtxmtxp[i*m+j])>GVT ) nmax++;
    }
  }  
  sa = new tdouble[nmax];
  ija = new tnlong[nmax];
  ija[0] = n+1;
  tnlong k=n;
  tdouble dt;
  for (i=0;i<n;++i) {
    for (j=0;j<m;++j) {
      dt = mtxmtxp[i*m+j];
      if (fabs(dt)>GVT ) {
        ++k;
        sa[k]=dt;
        ija[k]=j;
      }
    }
    ija[i+1]=k+1;
  }
}

FlxMtxSpars::~FlxMtxSpars()
{
  delete[] sa;
  delete[] ija;
}

FlxMtxSpars* FlxMtxSpars::copy()
{
  return new FlxMtxSpars(*this);
}

FlxMtxSpars& FlxMtxSpars::operator*=(const tdouble& s)
{
  const tnlong N=ija[nrows()];
  for ( tnlong i=0;i<N;i++) sa[i]*=s;
  return *this;
}

const tdouble FlxMtxSpars::operator()(const tnlong& i, const tnlong& j) const
{
  #if FLX_DEBUG
    if (i < 0 || j < 0 || i >= nrows() || j >= nrows()) {
      std::ostringstream ssV;
      ssV << "Index not within matrix.";
      throw FlxException("FlxMtxSpars::operator()_2", ssV.str() );
    }
  #endif
  for ( tnlong k = ija[i];k<ija[i+1];++k) {
    if (ija[k]==j) {
      return sa[k];
    }
    else if (ija[k]>j) {
      return ZERO;
    }
  }
  return ZERO;
}

void FlxMtxSpars::add_value(const tnlong& i, const tnlong& j, const tdouble& v)
{
  if (v==ZERO) return;
  #if FLX_DEBUG
    if (i < 0 || j < 0 || i >= nrows() || j >= nrows()) {
      std::ostringstream ssV;
      ssV << "Index not within matrix.";
      throw FlxException("FlxMtxSpars::add_value_1", ssV.str() );
    }
    if (j > i) {
      std::ostringstream ssV;
      ssV << "Not allowed to add value '" << v << "' at this position (" << i << ", " << j << ") - sparse matrix.";
      throw FlxException("FlxMtxSpars::add_value_2", ssV.str() );
    }
  #endif
  
  for (tnlong k = ija[i];k<ija[i+1];++k) {
    if (ija[k]==j) {
      sa[k]+=v;
      return;
    }
    else if (ija[k]>j) {
      break;
    }
  }
  std::ostringstream ssV;
  ssV << "Not allowed to add value '" << v << "' at this position (" << i << ", " << j << ") - sparse matrix.";
  throw FlxException("FlxMtxSpars::add_value_3", ssV.str() );
}

void FlxMtxSpars::MultMv(const flxVec& v, flxVec& w) const
{
  const tnlong N = nrows();
  #if FLX_KAHAN_MTX_SLTRI
    pdouble t;
  #else
    tdouble t;
  #endif
  const tdouble* vp = v.get_tmp_vptr_const();
  tdouble* wp = w.get_tmp_vptr();
  tnlong i, k;
  for (i=0;i<N;++i) {
    t=ZERO;
    for (k=ija[i];k<ija[i+1];++k) {
      t+= sa[k]*vp[ija[k]];
    }
    #if FLX_KAHAN_MTX_SLTRI
      wp[i]=t.cast2double();
    #else
      wp[i]=t;
    #endif
  }
}

const tdouble FlxMtxSpars::MultRowV(const flxVec& v, const tnlong rownumber) const
{
  #if FLX_KAHAN_MTX_SLTRI
    pdouble t = ZERO;
  #else
    tdouble t = ZERO;
  #endif
  const tdouble* vp = v.get_tmp_vptr_const();
  for (tnlong k=ija[rownumber];k<ija[rownumber+1];++k) {
    t+= sa[k]*vp[ija[k]];
  }
  #if FLX_KAHAN_MTX_SLTRI
    return t.cast2double();
  #else
    return t;
  #endif
}

void FlxMtxSpars::TransMultVec(const flxVec& v, flxVec& w)
{
  const tnlong N = nrows();
  w.set_zero();
  tdouble* wp = w.get_tmp_vptr();
  const tdouble* vp = v.get_tmp_vptr_const();
  tnlong i, k;
  for (i=0;i<N;++i) {
    for (k=ija[i];k<ija[i+1];++k) {
      wp[ija[k]]+= sa[k]*vp[i];
    }
  }
}

FlxMtxSparsLTri::FlxMtxSparsLTri(FlxMtxDiag& mtx)
{
  const tuint N = mtx.nrows();
  sa = new tdouble[N+1];
  ija = new tnlong[N+1];
  
  tnlong i;
  for (i=0;i<=N;++i) {
    ija[i] = N+1;
  }
  for (i=0;i<N;++i) {
    sa[i] = mtx.operator()(i,i);
  }
}

FlxMtxSparsLTri& FlxMtxSparsLTri::operator=(const FlxMtxDiag& mtx)
{
  const tuint N = mtx.nrows();
  for (tnlong i=0;i<=N;++i) {
    ija[i] = N+1;
  }
  for (tnlong i=0;i<N;++i) {
    sa[i] = mtx.operator()(i,i);
  }
  return *this;
}

FlxMtxSparsLTri::FlxMtxSparsLTri(const FlxMtxLTri& mtx)
{
  const tnlong n = mtx.nrows();
  tnlong nmax = n+1;
  const tdouble* mtxmtxp = &mtx.mtx[0];
  const tdouble GVT = GlobalVar.TOL();
  tnlong j,i;
  for (i = 1; i<n;i++) {
    for (j = 0; j<i;j++) {
      if (fabs(mtxmtxp[(i*i+i)/2+j])>GVT ) nmax++;
    }
  }  
  sa = new tdouble[nmax];
  ija = new tnlong[nmax];
  for (i = 0;i<n;i++) sa[i] = mtx(i,i);
  ija[0] = n+1;
  tnlong k=n;
  tdouble dt;
  for (i=0;i<n;i++) {
    for (j=0;j<i;j++) {
      dt = mtx(i,j);
      if (fabs(dt)>GVT ) {
        k++;
        sa[k]=dt;
        ija[k]=j;
      }
    }
    ija[i+1]=k+1;
  }
}

FlxMtxSparsLTri::FlxMtxSparsLTri(const FlxMtxSparsLTri& mtx)
:FlxMtxSpars(mtx)
{

}

void FlxMtxSparsLTri::CholeskyDec(tdoublePtr& saP, tnlong*& ijaP, FlxMtxSparsSym& mtx)
{
  tnlong N = mtx.nrows();
  // reserve memory in advance (on the safe side)
    // TODO: maybe a smarter solution is possible???
    tVec saT(N*N+1);
    tdouble* saTp = &saT[0];
    std::valarray<tnlong>ijaT(saT.size());
    tnlong* ijaTp = &ijaT[0];
  const tdouble* mtxsap = &mtx.sa[0];
  const tnlong* mtxijap = &mtx.ija[0];
  tnlong i,j,k,k2,l,pos;
  tdouble sum;
  const tnlong& lmax = mtx.ija[N];
  
  // copy initial setting
    for (i=0;i<N;++i) {
      saTp[i] = mtxsap[i];
    }
    pos=N+1;
    ijaTp[0] = pos;

  #if FLX_KAHAN_MTX_SLTRI
    pdouble t;
    pdouble sum2;
  #else
    tdouble& t = sum;
    tdouble sum2;
  #endif
  for (tnlong i = 0; i < N; ++i) {
    sum2=ZERO;
    l = mtxijap[i];
    for (j = 0; j < i; ++j) {
      if (l<lmax && j==mtxijap[l] && l<mtxijap[i+1]) {
        t = mtxsap[l];
        ++l;
      } else {
        t = ZERO;
      }
      k2 = ijaTp[j];
      for (k=ijaTp[i];k<pos;++k) {
        while (ijaTp[k2]<ijaTp[k] && k2 < ijaTp[j+1]) ++k2;
        if (k2>=ijaTp[j+1]) break;
        if (ijaTp[k]==ijaTp[k2]) {
          t-=saTp[k]*saTp[k2];
        }
      }
      #if FLX_KAHAN_MTX_SLTRI
        sum = t.cast2double();
      #endif
      if (sum!=ZERO) {
        sum/=saTp[j];
        saTp[pos] = sum;
        sum2+=pow2(sum);
        ijaTp[pos]=j;
        ++pos;
      }
    }
    // calculate diagonal entry
      #if FLX_KAHAN_MTX_SLTRI
        saTp[i] = saTp[i]-sum2.cast2double();
      #else
        saTp[i] = saTp[i]-sum2;
      #endif
      if (saTp[i]<=ZERO) {
        std::ostringstream ssV;
        ssV << "Matrix is not positiv definite. (sum = " << saT[i] << " in row " << i+1 << ")";
        throw FlxException("FlxMtxSparsLTri::CholeskyDec", ssV.str() );
      }
      saTp[i] = sqrt(saTp[i]);
    ijaTp[i+1] = pos;
  }
  
  saP = new tdouble[pos];
  ijaP = new tnlong[pos];
  for (i=0;i<pos;++i) {
    saP[i] = saTp[i];
    ijaP[i] = ijaTp[i];
  }
  #if FLX_DEBUG
    if (N>1) {
      GlobalVar.slog(4) << "  Sparseness: ";
      GlobalVar.slog(4) << GlobalVar.Double2String((pos-N-1)*2.0/(N*N-N)) << " (original: ";
      GlobalVar.slog(4) << GlobalVar.Double2String((mtx.ija[mtx.ija[0]-1]-mtx.ija[0]-1)*2.0/(mtx.ija[0]*mtx.ija[0]-mtx.ija[0])) << ")";
      GlobalVar.slog(4) << std::endl;
    }
  #endif
}

FlxMtxSparsLTri::FlxMtxSparsLTri(FlxMtxSparsSym& Sm)
{
  CholeskyDec(sa,ija,Sm);
}

FlxMtxSparsLTri* FlxMtxSparsLTri::copy()
{
  return new FlxMtxSparsLTri(*this);
}

void FlxMtxSparsLTri::MultInv(const flxVec& v, flxVec& w)
{
  const tnlong n = nrows();
  tdouble sum;
  #if FLX_KAHAN_MTX_SLTRI
    pdouble t;
  #else
    tdouble& t=sum;
  #endif
  const tdouble* vp = v.get_tmp_vptr_const();
  tdouble* wp = w.get_tmp_vptr();
  tnlong i, j;
  for (i=0;i<n;++i) {
    t=ZERO;
    for (j=ija[i];j<ija[i+1];++j) {
      t+=sa[j]*wp[ija[j]];
    }
    #if FLX_KAHAN_MTX_SLTRI
      sum = t.cast2double();
    #endif
    wp[i]=(vp[i]-sum)/sa[i];
  }
}

void FlxMtxSparsLTri::MultMv(const flxVec& v, flxVec& w) const
{
  const tnlong N = nrows();
  #if FLX_KAHAN_MTX_SLTRI
    pdouble t;
  #else
    tdouble t;
  #endif
  const tdouble* vp = v.get_tmp_vptr_const();
  tdouble* wp = w.get_tmp_vptr();
  for (tnlong i=N;i>0;--i) {
    t=sa[i-1]*vp[i-1];
    for (tnlong k=ija[i];k>ija[i-1];--k) {
      t+= sa[k-1]*vp[ija[k-1]];
    }
    #if FLX_KAHAN_MTX_SLTRI
      wp[i-1]=t.cast2double();
    #else
      wp[i-1]=t;
    #endif
  }
}

const tdouble FlxMtxSparsLTri::MultRowV(const flxVec& v, const tnlong rownumber) const
{
  throw FlxException_NotImplemented("FlxMtxSparsLTri::MultRowV_1");
}

void FlxMtxSparsLTri::TransMultVec(flxVec& v)
{
  tnlong N = nrows();
  tnlong i, k;
  #if FLX_KAHAN_MTX_SLTRI
    flxpVec w(v);
    pdouble* wp = w.get_tmp_vptr();
  #else
    flxVec& w = v;
    tdouble* wp = w.get_tmp_vptr();
  #endif
  for (i=0;i<N;++i) {
    for (k=ija[i];k<ija[i+1];++k) {
      wp[ija[k]] += wp[i]*sa[k];
    }
    wp[i]=wp[i]*sa[i];
  }
  #if FLX_KAHAN_MTX_SLTRI
    v = w;
  #endif
}

void FlxMtxSparsLTri::TransMultVec(const flxVec& v, flxVec& w)
{
  throw FlxException_NotImplemented("FlxMtxSparsLTri::TransMultVec");
}

// FlxMtxSparsLTri* FlxMtxSparsLTri::MultM_Mdiag(const FlxMtxDiag& D) const
// {
//   #if FLX_DEBUG
//     if (nrows()!=D.nrows()) {
//       throw FlxException("FlxMtxSparsLTri::MultM_Mdiag", "ERROR" );
//     }
//   #endif
//   const tnlong nmax = ija[ija[0]-1];
//   FlxMtxSparsLTri* res = new FlxMtxSparsLTri(nmax);
//   
//   for (register tnlong i=0;i<nmax;++i) {
//     res->ija[i] = ija[i];
//   }
//   const tnlong N = nrows();
//   register tnlong i, k;
//   register tdouble t;
//   for (i=0;i<N;++i) {
//     t = D.operator()(i,i);
//     res->sa[i] = sa[i]*t;
//     for (k=ija[i];k<ija[i+1];++k) {
//       res->sa[k] = sa[k]*t;
//     }
//   }
//   
//   return res;
// }

void FlxMtxSparsLTri::MultM_Mdiag(const FlxMtxDiag& D, FlxMtxSparsLTri& res) const
{
  const tnlong nmax = get_nmax();
  #if FLX_DEBUG
    if (nrows()!=D.nrows()) {
      throw FlxException_Crude("FlxMtxSparsLTri::MultM_Mdiag_1");
    }
  #endif
  for (tnlong i=0;i<nmax;++i) {
    res.ija[i] = ija[i];
  }
  const tnlong N = nrows();
  tnlong i, k;
  tdouble t;
  for (i=0;i<N;++i) {
    t = D.operator()(i,i);
    res.sa[i] = sa[i]*t;
    for (k=ija[i];k<ija[i+1];++k) {
      res.sa[k] = sa[k]*t;
    }
  }
}

const tdouble FlxMtxSparsLTri::det_log() const
{
  tdouble r = ZERO;
  const tnlong msize = ija[0]-1;
  for (tnlong i=0;i<msize;++i) {
    r += log(sa[i]);
  }
  return r;
}

const tdouble FlxMtxSparsLTri::operator()(const tnlong& i, const tnlong& j) const
{
  if (i == j) return sa[i];
  if (j > i) return ZERO;
  return FlxMtxSpars::operator()(i,j);
}

void FlxMtxSparsLTri::add_value(const tnlong& i, const tnlong& j, const tdouble& v)
{
  if (v==ZERO) return;
  if (i == j) {
    sa[i]+=v;
    return;
  }
  FlxMtxSpars::add_value(i,j,v);
}

FlxMtxSparsSymLU::FlxMtxSparsSymLU(FlxMtxSparsSym& mtx)
{
  sa = NULL; 
  ija = NULL;
  FlxMtxSparsLTri::CholeskyDec(sa,ija,mtx);
}

FlxMtxSparsSymILU::FlxMtxSparsSymILU(const FlxMtxSparsSym& mtx, bool mod0diagentry)
{
  const tnlong nmax = mtx.ija[mtx.nrows()];
  sa = new tdouble[nmax];
  ija = new tnlong[nmax];
  const tdouble* mtxsap = &mtx.sa[0];
  const tnlong* mtxijap = &mtx.ija[0];
  for (tnlong i = 0; i < nmax; ++i) {
    sa[i] = mtxsap[i];
    ija[i] = mtxijap[i];
  }
  
  const tnlong N = nrows();
  tdouble sum;
  tnlong i, j, k;
  
  #if FLX_KAHAN_MTX_SSILU
    pdouble t;
  #else
    tdouble& t = sum;
  #endif
  
  for (i = 0; i < N; ++i) {
    for (j=ija[i];j<ija[i+1];++j) {
      t = sa[j];
      for (k=ija[i];k<j;++k) {
        t-=sa[k]*operator()(ija[j],ija[k]);  // TODO more efficient version possible ???
      }
      #if FLX_KAHAN_MTX_SSILU
        sum = t.cast2double();
      #endif
      sa[j] = sum/(sa[ija[j]]);
    }
    t = sa[i];
    for (k=ija[i];k<ija[i+1];++k) {
      t-=pow2(sa[k]);
    }
    #if FLX_KAHAN_MTX_SSILU
      sum = t.cast2double();
    #endif
    if (sum<=ZERO) {
      if (!mod0diagentry) {
        std::ostringstream ssV;
        ssV << "Matrix is not positiv definite. (sum=" << sum << " in i=" << i << ")"; // <<  std::endl << *this;
        throw FlxException("FlxMtxSparsSymILU::FlxMtxSparsSymILU_1", ssV.str() );
      } else {
        sa[i] = sqrt(sa[i]);
        GlobalVar.alert.alert("FlxMtxSparsSymILU::FlxMtxSparsSymILU_2","ILU: modified diagonal entry.");
      }
    } else {
      sa[i] = sqrt(sum);
    }
  }
}

void FlxMtxSparsSymILU::MultMv(const flxVec& v, flxVec& w) const
{
  // calculate L*y=f
  const tnlong N = nrows();
  #if FLX_KAHAN_MTX_SSILU
    flxpVec W(N);
    pdouble d;
    pdouble* Wp = W.get_tmp_vptr();
  #else 
    flxVec& W = w;
    tdouble* Wp = W.get_tmp_vptr();
    tdouble d;
  #endif
  const tdouble* vp = v.get_tmp_vptr_const();
  const tdouble* wp = w.get_tmp_vptr_const();
  
  for (tnlong i = 0; i < N; ++i) {
    d=vp[i];
    for (tnlong j=ija[i];j<ija[i+1];++j) {
      d-=sa[j]*wp[ija[j]];
    }
    Wp[i] = d/sa[i];
  }
  //calculate L^T*x=y
  for (tnlong i = N; i > 0; --i) {
    Wp[i-1] = Wp[i-1]/sa[i-1];
    for (tnlong j=ija[i-1];j<ija[i];++j) {
      Wp[ija[j]]-=Wp[i-1]*sa[j];
    }
  }
  #if FLX_KAHAN_MTX_SSILU
    w = W;
  #endif
}

void FlxMtxSparsSymILU::MultMv(const flxpVec& v, flxpVec& w) const
{
  // calculate L*y=f
  const tnlong N = nrows();
  tnlong i,j;
  pdouble d;
  const pdouble* vp = w.get_tmp_vptr_const();
  pdouble* wp = w.get_tmp_vptr();
  for (i = 0; i < N; ++i) {
    d=vp[i];
    for (j=ija[i];j<ija[i+1];++j) {
      d-=wp[ija[j]]*sa[j];
    }
    wp[i] = d/sa[i];
  }
  //calculate L^T*x=y
  for (i = N; i > 0; --i) {
    wp[i-1] = wp[i-1]/sa[i-1];
    for (j=ija[i-1];j<ija[i];++j) {
      wp[ija[j]]-=wp[i-1]*sa[j];
    }
  }
}

FlxMtxSparsSymILU::~FlxMtxSparsSymILU()
{
  delete [] sa;
  delete [] ija;
}

const tdouble FlxMtxSparsSymILU::operator()(const tnlong& i, const tnlong& j) const
{
  #if FLX_DEBUG
    if (i < 0 || j < 0 || i >= nrows() || j >= nrows()) {
      std::ostringstream ssV;
      ssV << "Index not within matrix.";
      throw FlxException("FlxMtxSparsSym::operator()_2", ssV.str() );
    }
  #endif
  if (i == j) return sa[i];
  if (j > i) return ZERO;
  for (tnlong k = ija[i];k<ija[i+1];k++) {
    if (ija[k]==j) {
      return sa[k];
    }
    else if (ija[k]>j) {
      return ZERO;
    }
  }
  return ZERO;
}

FlxMtx_base& FlxMtxSparsSymILU::operator*=(const tdouble& s)
{
  throw FlxException_NotImplemented("FlxMtxSparsSymILU::operator*=");
}

void FlxMtxSparsSymILU::add_value(const tnlong& i, const tnlong& j, const tdouble& v)
{
  throw FlxException_NotImplemented("FlxMtxSparsSymILU::add_value");
}

const tdouble FlxMtxSparsSFEMSym::operator()(const tnlong& i, const tnlong& j) const
{
  tnlong ii = i - (i%Kdim);
  tnlong in = ii/Kdim;
  ii = i%Kdim;
  tnlong jj = j - (j%Kdim);
  tnlong jn = jj/Kdim;
  jj = j%Kdim;
  
  if (in == jn) {
    return sa[in]*(sb[in]->operator()(ii,jj));
  }
  if (jn > in) {
    std::swap(in, jn);
    std::swap(ii, jj);
  }
  for (tnlong k = ija[in];k<ija[in+1];k++) {
    if (ija[k]==jn) {
      if (sa[k]==ZERO) {
        return ZERO;
      } else {
        return sa[k]*(sb[k]->operator()(ii,jj));
      }
    }
    else if (ija[k]>jn) {
      return ZERO;
    }
  }
  return ZERO;
}

const bool FlxMtxSparsSFEMSym::isPositiveDefinite(tuint& r, const bool fixIt)
{
  return sb[0]->isPositiveDefinite(r,fixIt);
}


void FlxMtxSparsSFEMSym::add_value(const tnlong& i, const tnlong& j, const tdouble& v)
{
  throw FlxException_NotImplemented("FlxMtxSparsSFEMSym::add_value");
}

const bool FlxMtxSparsSFEMSym::solve_CG(flxVec& x, const flxVec& f, tdouble& eps, tnlong& iter, const tuint pcn, const bool startZero)
{
  if (!startZero) {
    flxVec xx(Kdim);
    flxVec ff(f.get_tmp_vptr_const(),Kdim);
    tnlong iterV = iter;
    tdouble epsV = eps;
    sb[0]->solve_CG(xx,ff,epsV,iterV,1,true);
    flxVec xtmp(x.get_tmp_vptr(),Kdim);
    xtmp = xx;
  }
  return FlxMtx_baseS::solve_CG(x, f, eps, iter, pcn, startZero);
}

// FlxMtxDiag* FlxMtxSparsSFEMSym::get_diag() const
// {
//   FlxMtxDiag* d = new FlxMtxDiag(nrows());
//   FlxMtxDiag* di = sb[0]->get_diag();
//   
//   tnlong N = (ija[0]-1);
//   for (tnlong i = 0; i < N; ++i) {
//     for (tnlong j = 0; j < Kdim; ++j) {
//       d[i*Kdim+j] = di[j];
//     }
//   }
//   delete di;
//   return d;
// }

void FlxMtxSparsSFEMSym::output_Mtx(std::ostream& sout) const
{
  tnlong N = ija[0]-1;
  
  for (tnlong i =0;i<N;++i) {
    tnlong k=0;
    for (tnlong j = ija[i]; j<ija[i+1];++j) {
      while (k<ija[j]) {
        sout << GlobalVar.Double2String(0.,false,-1,0) << " ";
        k++;
      }
      sout << GlobalVar.Double2String(sa[j]) << "*K"<< (box.find(sb[j]))->second << " ";
      k++;
    }
    while (k<i) {
      sout << GlobalVar.Double2String(0.,false,-1,0) << " ";
      k++;
    }
    sout << GlobalVar.Double2String(sa[i])  << "*K"<< (box.find(sb[i]))->second;
    sout << std::endl;
  }
  
  for (std::map<FlxMtx_baseS*, tuint>::const_iterator it = box.begin(); it != box.end(); ++it) {
    sout << "K" << it->second << "=";
    it->first->output_Mtx(sout);
  }
}

void FlxMtxSparsSFEMSym::assembleMinv(int i)
{
  if (i == 3) {
    tnlong N = (ija[0]-1);
    tVec c0ii(N);
    for (tnlong j = 0; j < N; ++j) {
      c0ii[j] = sa[j];
    }
    FlxMtxSparsSym* S = dynamic_cast<FlxMtxSparsSym*>(sb[0]);
    if (S == NULL) {
      std::ostringstream ssV;
      ssV << "Matrix is not a symmetric sparse matrix.";
      throw FlxException("FlxMtxSparsSFEMSym::assembleMinv_1", ssV.str() );
    }
    Minv = new FlxMtxPrecnInvSFEMSym(*S,c0ii);
  } else if (i == 4 || i==5 || i==6) {
    tnlong N = (ija[0]-1);
    tVec c0ii(N);
    for (tnlong j = 0; j < N; ++j) {
      c0ii[j] = sa[j];
    }
    FlxMtxSparsSym* S = dynamic_cast<FlxMtxSparsSym*>(sb[0]);
    if (S == NULL) {
      std::ostringstream ssV;
      ssV << "Matrix is not a symmetric sparse matrix.";
      throw FlxException("FlxMtxSparsSFEMSym::assembleMinv_1", ssV.str() );
    }
    if (i == 4 ) {
      Minv = new FlxMtxPrecnILUSFEMSym(*S,c0ii,true);
    } else if (i == 5) {
      Minv = new FlxMtxPrecnILUSFEMSym(*S,c0ii);
    } else {
      Minv = new FlxMtxPrecnILUSFEMSym(*S,c0ii,false,true);
    }
  } else {
    FlxMtx_baseS::assembleMinv(i);
  }
}

FlxMtxSparsSFEMSym::~FlxMtxSparsSFEMSym()
{
  delete[] sa;
  delete[] sb;
  delete[] ija;
}

FlxMtxSparsSFEMSym::FlxMtxSparsSFEMSym(std::valarray< tdouble >& fM, std::valarray< int >& KM, FlxMtx_baseS** StfMtxV, tulong P)
:Kdim(StfMtxV[0]->ncols())
{
  // assign box
    for ( int i = 0; i <= KM.max(); ++i) {
      std::pair<FlxMtx_baseS*, tuint> Element(StfMtxV[i],i);
      box.insert(Element);
    }
  
  tulong nmax = P+1;
  for (tulong i = 1; i<P;++i) {
    for (tulong j = 0; j<i;++j) {
      if (fM[(i*i+i)/2+j]!=ZERO) nmax++;
    }
  }  
  sa = new tdouble[nmax];
  sb = new FlxMtx_baseS*[nmax];
  ija = new tnlong[nmax];
  for (tnlong i = 0;i<P;i++) {
    sa[i] = fM[(i*i+i)/2+i];
    sb[i] = StfMtxV[KM[(i*i+i)/2+i]];
  }
  ija[0] = P+1;
  tnlong k=P;
  for (tnlong i=0;i<P;i++) {
    for (tnlong j=0;j<i;j++) {
      if (fM[(i*i+i)/2+j]!=ZERO) {
        k++;
        sa[k]=fM[(i*i+i)/2+j];
        sb[k] = StfMtxV[KM[(i*i+i)/2+j]];
        ija[k]=j;
      }
    }
    ija[i+1]=k+1;
  }
}

FlxMtx_base& FlxMtxSparsSFEMSym::operator*=(const tdouble& s)
{
  const tnlong N=ija[ija[0]-1];
  for (tnlong i=0;i<N;i++) sa[i]*=s;
  return *this;
}

void FlxMtxSparsSFEMSym::MultMv(const flxVec& v, flxVec& w) const
{
  tnlong N = ija[0]-1;
  #if FLX_KAHAN_MTX_SSFEM
    flxpVec W(N);
  #else
    flxVec& W = w;
  #endif
  flxVec t(Kdim);
  for (tnlong i=0;i<N;++i) {
    const flxVec vtmp(&(v[i*Kdim]),Kdim);
    #if FLX_KAHAN_MTX_SSFEM
      flxpVec wtmp(&(W[i*Kdim]),Kdim);
      (*sb[i]).MultMv(vtmp,t);
      t *= sa[i];
      wtmp = t;
    #else
      flxVec wtmp(&(W[i*Kdim]),Kdim);
      (*sb[i]).MultMv(vtmp,wtmp);
      wtmp *= sa[i];
    #endif
  }
  for (tnlong i=0;i<N;++i) {
    const flxVec vtmp_i(&(v[i*Kdim]),Kdim);
    #if FLX_KAHAN_MTX_SSFEM
      flxpVec wtmp_i(&(W[i*Kdim]),Kdim);
    #else
      flxVec wtmp_i(&(W[i*Kdim]),Kdim);
    #endif
    for (tnlong k=ija[i];k<ija[i+1];++k) {
      const flxVec vtmp_k(&(v[ija[k]*Kdim]),Kdim);
      #if FLX_KAHAN_MTX_SSFEM
        flxpVec wtmp_k(&(W[ija[k]*Kdim]),Kdim);
      #else
        flxVec wtmp_k(&(W[ija[k]*Kdim]),Kdim);
      #endif
      (*sb[k]).MultMv(vtmp_k,t);
      t *= sa[k];
      wtmp_i += t;
      (*sb[k]).MultMv(vtmp_i,t);
      t *= sa[k];
      wtmp_k += t;
    }
  }
  #if FLX_KAHAN_MTX_SSFEM
    w = W;
  #endif
}

void FlxMtxSparsSFEMSym::MultMv(const flxpVec& v, flxpVec& w) const
{
  const tnlong N = ija[0]-1;
  for (tnlong i=0;i<N;++i) {
    flxpVec wtmp(&(w[i*Kdim]),Kdim);
    const flxpVec vtmp(&(v[i*Kdim]),Kdim);
    (*sb[i]).MultMv(vtmp,wtmp);
    wtmp *= sa[i];
  }
  flxVec t(Kdim);
  for (tnlong i=0;i<N;++i) {
    flxpVec wtmp_i(&(w[i*Kdim]),Kdim);
    const flxpVec vtmp_i(&(v[i*Kdim]),Kdim);
    for (tnlong k=ija[i];k<ija[i+1];++k) {
      flxpVec wtmp_k(&(w[ija[k]*Kdim]),Kdim);
      const flxpVec vtmp_k(&(v[ija[k]*Kdim]),Kdim);
      (*sb[k]).MultMv(vtmp_k,t);
      t *= sa[k];
      wtmp_i += t;
      (*sb[k]).MultMv(vtmp_i,t);
      t *= sa[k];
      wtmp_k += t;
    }
  }
}

void FlxMtxPrecnInvSFEMSym::add_value(const tnlong& i, const tnlong& j, const tdouble& v)
{
  throw FlxException_NotImplemented("FlxMtxPrecnInvSFEMSym::add_value");
}

const tdouble FlxMtxPrecnInvSFEMSym::operator()(const tnlong& i, const tnlong& j) const
{
  throw FlxException_NotImplemented("FlxMtxPrecnInvSFEMSym::operator()");
}

FlxMtx_base& FlxMtxPrecnInvSFEMSym::operator*=(const tdouble& s)
{
  throw FlxException_NotImplemented("FlxMtxPrecnInvSFEMSym::operator*=");
}

FlxMtxPrecnInvSFEMSym::FlxMtxPrecnInvSFEMSym(FlxMtxSparsSym& K, const tVec& c0ii): c0iiInv(c0ii), Kdim(K.nrows())
{
  const tnlong N=c0ii.size();
  tdouble* c0iiInvp = &c0iiInv[0];
  for (tnlong i = 0; i < N; ++i) {
    c0iiInvp[i] = 1/c0iiInvp[i];
  }
  
  FlxMtxLTri L(Kdim);
  
//   FlxMtxSym S(K);
//   L.CholeskyDec(S);
  
  L.CholeskyDec(K);
  L.Invert();
  Kinv = new FlxMtxSym(Kdim);
  Kinv->assign_LTL(L);
}

FlxMtxPrecnInvSFEMSym::~FlxMtxPrecnInvSFEMSym()
{
  delete Kinv;
}

void FlxMtxPrecnInvSFEMSym::MultMv(const flxVec& v, flxVec& w) const
{
  const tnlong N = c0iiInv.size();
  for (tnlong i=0;i<N;i++) {
    flxVec wtmp(&(w[i*Kdim]),Kdim);
    const flxVec vtmp(&(v[i*Kdim]),Kdim);
    (*Kinv).MultMv(vtmp,wtmp);
    wtmp *= c0iiInv[i];
  }
}

void FlxMtxPrecnILUSFEMSym::add_value(const tnlong& i, const tnlong& j, const tdouble& v)
{
  throw FlxException_NotImplemented("FlxMtxPrecnILUSFEMSym::add_value");
}

FlxMtx_base& FlxMtxPrecnILUSFEMSym::operator*=(const tdouble& s)
{
  throw FlxException_NotImplemented("FlxMtxPrecnILUSFEMSym::operator*=");
}

const tdouble FlxMtxPrecnILUSFEMSym::operator()(const tnlong& i, const tnlong& j) const
{
  throw FlxException_NotImplemented("FlxMtxPrecnILUSFEMSym::operator()");
}

FlxMtxPrecnILUSFEMSym::~FlxMtxPrecnILUSFEMSym()
{
  delete ILU;
}

void FlxMtxPrecnILUSFEMSym::MultMv(const flxVec& v, flxVec& w) const
{
  const tnlong N = c0iiInv.size();
  for (tnlong i=0;i<N;i++) {
    flxVec wtmp(&(w[i*Kdim]),Kdim);
    const flxVec vtmp(&(v[i*Kdim]),Kdim);
    (*ILU).MultMv(vtmp,wtmp);
    wtmp *= c0iiInv[i];
  }
}

void FlxMtxPrecnILUSFEMSym::MultMv(const flxpVec& v, flxpVec& w) const
{
  const tnlong N = c0iiInv.size();
  for (tnlong i=0;i<N;i++) {
    flxpVec wtmp(&(w[i*Kdim]),Kdim);
    const flxpVec vtmp(&(v[i*Kdim]),Kdim);
    (*ILU).MultMv(vtmp,wtmp);
    wtmp *= c0iiInv[i];
  }
}

FlxMtxPrecnILUSFEMSym::FlxMtxPrecnILUSFEMSym(FlxMtxSparsSym& K, const tVec& c0ii, bool FullDecomp, bool mod0diagentry)
: c0iiInv(c0ii), Kdim(K.nrows())
{        
  const tnlong N=c0ii.size();
  tdouble* c0iiInvp = &c0iiInv[0];
  for (tnlong i = 0; i < N; ++i) {
    c0iiInvp[i] = 1/c0iiInvp[i];
  }
  
  if (FullDecomp) {
    ILU = new FlxMtxSparsSymLU(K);
  } else {
    ILU = new FlxMtxSparsSymILU(K,mod0diagentry);
  }
}

