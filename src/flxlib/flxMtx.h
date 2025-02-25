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

#ifndef fesslix_FlxMtx_H
#define fesslix_FlxMtx_H

#include "flxVec.h"

#include <algorithm>
#include <map>
#include <vector>

class FlxMtx;
class FlxMtxSym;


// -------------------------- MATRIX ----------------------------------------------------------------------

/* Matrix Types defined in this library
 * ====================================
 * 
 * FlxMtx_base (*)
 * 
 *   FlxMtx                        full nxm-Matrix
 *   FlxMtxTransformation        a transformating matrix consisting of (nxn)-matrices of type FlxMtx
 *   FlxMtxLTri                        full nxn lower triangular matrix
 *   FlxMtxSparsLTri                sparse nxn lower triangular matrix
 * 
 *   FlxMtx_baseS (*)
 *     FlxMtxSparsSymILU        incomplete LU factorization of a sparse symmetric matrix
 *     FlxMtxSparsSymLU                Cholesky decomposition of a (sparse) symmetric matrix
 *     FlxMtxSym                full symmetric nxn matrix
 *     FlxMtxDiag                diagonal matrix
 *     FlxMtxIdentity                identity matrix
 *     FlxMtxSparsSym                sparse symmetric nxn matrix
 *     FlxMtxSparsSFEMSym        sparse symmetric matrix needed for SFEM problems (coefficients are of type FlxMtxSparsSym)
 *     FlxMtxPrecnInvSFEMSym                -> preconding
 *     FlxMtxPrecnILUSFEMSym                -> preconding
 *     
 */


/**
* @brief The base class of matrices in Fesslix
*/
class FLXLIB_EXPORT FlxMtx_base {
  public:
    virtual ~FlxMtx_base() {};
    virtual FlxMtx_base* copy();
    virtual const flxVec operator*(const flxVec& v) const;
    virtual const flxpVec operator*(const flxpVec& v) const;
    virtual FlxMtx_base& operator*=(const tdouble& s) = 0;
    virtual const tnlong nrows() const = 0;
    virtual const tnlong ncols() const = 0;
    /**
    * @brief Matrix-vector product: w=M*v;
    */
    virtual void MultMv(const flxVec& v, flxVec& w) const = 0;
    /**
    * @brief Matrix-vector product: w=M*v;
    */
    virtual void MultMv(const flxpVec& v, flxpVec& w) const;
    /**
    * @brief add value v to coefficient (i,j)  -  has to be within range of matrix (not checked!!!)
    */
    virtual void add_value(const tnlong& i, const tnlong& j, const tdouble& v)=0;
    /**
    * @brief get coefficient (i,j)  -  has to be within range of matrix (not checked!!!)
    */
    virtual const tdouble operator()(const tnlong& i, const tnlong& j) const = 0;
    /**
    * @brief sets all coefficients in the matrix to zero
    */
    void set_zeroMtx() {operator*=(ZERO);};
    /**
    * @brief res = *this * rhs
    */
    void MultMtx(const FlxMtx_base& rhs, FlxMtx& res);
    
    /**
    * @brief output-function of the matrix
    */
    virtual void output_Mtx(std::ostream& sout) const;
    /**
    * @brief output matrix in Octave format
    */
    virtual void output_OctaveMtx(std::ostream& sout, bool checkTOLv=false, bool doendl=true) const;
    /**
    * @brief returns the largest coefficient
    */
    virtual const tdouble max() const;
    /**
    * @brief returns the smallest coefficient
    */
    virtual const tdouble min() const;
    /**
    * @brief returns the index of the largest coefficient
    */
    virtual const size_t maxID() const;
    /**
    * @brief returns the index of the smallest coefficient
    */
    virtual const size_t minID() const;
    /**
    * @brief returns a pointer to the internal storage of the matrix
    * @note implement this only, if the matrix is a full matrix and the internal storage is of size Nc*Nr
    */
    virtual tdouble* get_internalPtr(const bool throwErr);
};

FLXLIB_EXPORT std::ostream& operator<<(std::ostream& os, const FlxMtx_base& M);

class FlxMtxDiag;

/**
* @brief the base class for symmetric matrices in Fesslix
*/
class FLXLIB_EXPORT FlxMtx_baseS : public FlxMtx_base {
  protected:
    FlxMtx_baseS* Minv;
    virtual void assembleMinv(int i);
    virtual void preconding(const flxVec& r, flxVec& z, const int i=0) const;
    virtual void preconding(const flxpVec& r, flxpVec& z, const int i=0) const;
  public:
    FlxMtx_baseS() : Minv(NULL) {}
    virtual ~FlxMtx_baseS();
    /**
    * @brief Conjugate Gradient algorithm
    */
    virtual const bool solve_CG(flxVec& x, const flxVec& f, tdouble &eps, tnlong& iter, const tuint pcn, const bool startZero);
    /**
    * @brief returns the diagnol matrix of this matrix
    */
    virtual FlxMtxDiag* get_diag() const;
    /**
    * @brief returns the inverse of this matrix
    */
    virtual FlxMtx_baseS* get_Inverse();
    virtual void Invert();
    /**
    * @brief checks if the matrix is positive definite
    */
    virtual const bool isPositiveDefinite(tuint& r, const bool fixIt=false) = 0;
    /**
    * @brief returns a pointer to the internal storage
    *         values are stored below the diagonal !!!
    */
    virtual tdouble* get_VecPointer();
    /**
    * @brief assembles vec (dim=N*N) with the elements of this matrix (returns a full matrix)
    * @param vec NxN-vector (size is NOT checked)
    */
    virtual void get_VecPointer_full(tdouble* vec) const;
    
    virtual void output_Mtx(std::ostream& sout) const;
};

class FlxMtxTransformation;
class FlxMtxSpars;

/**
* @brief stores a full matrix
*/
class FLXLIB_EXPORT FlxMtx : public FlxMtx_base {
  private:
    tnlong rsize;
    tnlong csize;
    flxVec mtx;
  public:
    const tnlong nrows() const {return rsize;};
    const tnlong ncols() const {return csize;};

    FlxMtx(const tnlong rsize,const tnlong csize) : rsize(rsize),csize(csize),mtx(rsize*csize) {};
    FlxMtx(const tnlong rsize,const tnlong csize, const tdouble* dptr) : rsize(rsize),csize(csize),mtx(dptr,rsize*csize) {};
    FlxMtx(const FlxMtx_base& mB);
    FlxMtx(const FlxMtxSym* symM);
    FlxMtx* copy();
    /**
    * @brief Matrix-vector product: w=M*v;
    */
    virtual void MultMv(const flxVec& v, flxVec& w) const;
    FlxMtx& operator+=(const FlxMtx& m);
    FlxMtx& operator*=(const tdouble& s) {mtx*=s; return *this; };
//     FlxMtx operator*(const FlxMtx_base& m) const;
    tdouble& operator()(const tnlong& i, const tnlong& j);
    const tdouble operator()(const tnlong& i, const tnlong& j) const;
    void add_value(const tnlong& i, const tnlong& j, const tdouble& v);
    /**
    * @brief res = Transpose(M).M
    */
    void TransposeMmultM(FlxMtxSym& res) const;
    /**
    * @brief w = Transpose(M).v
    */
    void TransposeMmultVec(const flxVec& v, flxVec& w) const;
    /**
    * @brief returns the matrix-vector (constant!)
    */
    const flxVec& get_mtx_flxVec() const { return mtx; }
    tdouble* get_internalPtr(const bool throwErr);
    
    friend class FlxMtxSpars;
    friend void MtxProd_BTKB(const FlxMtx&, const FlxMtxSym&, FlxMtxSym&);
    friend void MtxProd_BTKB(const FlxMtxTransformation&, const FlxMtxSym&, FlxMtxSym&);
};

/**
* @brief stores a transformation matrix
*/
class FLXLIB_EXPORT FlxMtxTransformation : public FlxMtx_base {
  private:
    tnlong rsize;
    std::vector<FlxMtx*> Ttm;
  public:
    const tnlong nrows() const {return rsize;};
    const tnlong ncols() const {return rsize;};

    FlxMtxTransformation(const std::vector<FlxMtx*> &Ttm);
    ~FlxMtxTransformation();
//     FlxMtx* copy();
    FlxMtxTransformation& operator*=(const tdouble& s);
    const tdouble operator()(const tnlong& i, const tnlong& j) const;
    void add_value(const tnlong& i, const tnlong& j, const tdouble& v);
    /**
    * @brief Matrix-vector product: w=M*v;
    */
    virtual void MultMv(const flxVec& v, flxVec& w) const;
    
    friend void MtxProd_BTKB(const FlxMtxTransformation&, const FlxMtxSym&, FlxMtxSym&);
};


class FlxMtxSparsSym;
class FlxMtxSparsLTri;

/**
* @brief stores a lower triangular matrix
*/
class FLXLIB_EXPORT FlxMtxLTri : public FlxMtx_base {
  protected:
    tnlong msize;
    flxVec mtx;
    bool is_ldl;
  public:
    const tnlong nrows() const {return msize;};
    const tnlong ncols() const {return msize;};

    FlxMtxLTri(const tnlong& msize) :msize(msize),mtx((msize*msize+msize)/2), is_ldl(false) {};
    FlxMtxLTri(const FlxMtxSparsLTri& Msp);
    FlxMtxLTri* copy();
    /**
    * @brief Matrix-vector product: w=M*v;
    * NOTE: for this special case, w and v can refer to the same vector
    */
    virtual void MultMv(const flxVec& v, flxVec& w) const;
    /**
    * @brief Calculates w=L^-1*v -> inverse matrix is not explicitly computed !!!   
    * NOTE: for this special case, w and v can refer to the same vector
    * a.k.a. forward substitution
    * for is_ldl, this solves Lw=v (unless ldl_sqrt=true, then L*sqrt(D) w = v is solved
    */
    void MultInv(const flxVec& v,flxVec& w, bool ldl_sqrt=false);
    /**
    * @brief Calculates w=Transpose(*this)^-1*v -> inverse matrix is not explicitly computed !!!
    * NOTE: for this special case, w and v can refer to the same vector
    * a.k.a. back substitution  
    * for is_ldl, this solves D(L^T)w=v
    */
    void TransMultInv(const flxVec& v,flxVec& w);
    /**
    * @brief Calculates the determinant of the matrix (log-transform)
    */
    const tdouble det_log() const;
//     FlxMtx& operator+=(const FlxMtxLTri& m);
    FlxMtxLTri& operator*=(const tdouble& s);
    const tdouble operator()(const tnlong& i, const tnlong& j) const;
    void add_value(const tnlong& i, const tnlong& j, const tdouble& v);
    /**
    * @brief Assembles this matrix by doing a CholeskyDecomposition of Sm
    */
    FlxMtxLTri& CholeskyDec(FlxMtxSym& Sm, const bool do_ldl=false);
    FlxMtxLTri& CholeskyDec(const flxVec& mtxV, const bool do_ldl=false);        // mtxV -> vector of a FlxMtxSym
    FlxMtxLTri& CholeskyDec(FlxMtxSparsSym& Sm);
    /**
    * @brief Invert the Matrix
    */
    FlxMtxLTri& Invert();
    /**
    * @brief solves the following problem: v=Transpose(*this)*v;
    */
    void TransMultVec(flxVec& v);

    const size_t count_nan() const;
    
    friend FLXLIB_EXPORT void MtxProd_BTKB_ltri(const FlxMtxLTri& B, const FlxMtxSym& E, FlxMtxSym& K);
    friend class FlxMtxSym;
    friend class FlxMtxSparsLTri;
};

class FlxMtxSymBand;

/**
* @brief stores a lower triangular matrix
*/
class FLXLIB_EXPORT FlxMtxLTriBand : public FlxMtx_base {
  protected:
    tnlong msize;
    tnlong bsize;
    tVec mtx;
  public:
    const tnlong countUp2RDiag(const tnlong R) const;  // count up to diagonal entry of row R ( R=0...msize-1)
    const tnlong nrows() const {return msize;};
    const tnlong ncols() const {return msize;};
    const tuint get_bsize() const { return bsize; }
    
    FlxMtxLTriBand(const tnlong& msize,const tnlong& bsize) :msize(msize),bsize(bsize),mtx((bsize+1)*msize-(bsize*bsize+bsize)/2) {};
    FlxMtxLTriBand* copy();
    /**
    * @brief Matrix-vector product: w=M*v;
    * NOTE: for this special case, w and v can refer to the same vector
    */
    virtual void MultMv(const flxVec& v, flxVec& w) const;
    /**
    * @brief Calculates v=L^-1*v -> inverse matrix is not explicitly computed !!!
    */
    void MultInv(const flxVec& v,flxVec& w);
    void MultInvTrans(const flxVec& v,flxVec& w);
    FlxMtxLTriBand& operator*=(const tdouble& s) {mtx*=s; return *this; };
    const tdouble operator()(const tnlong& i, const tnlong& j) const;
    void add_value(const tnlong& i, const tnlong& j, const tdouble& v);
    /**
    * @brief Assembles this matrix by doing a CholeskyDecomposition of Sm
    */
    FlxMtxLTriBand& CholeskyDec(FlxMtxSymBand& Sm);
    /**
    * @brief Invert the Matrix
    */
    FlxMtxLTriBand& Invert();
    /**
    * @brief solves the following problem: v=Transpose(*this)*v;
    */
    void TransMultVec(flxVec& v);

    friend class FlxMtxSymBand;
};


// Matrix - Example:
//    0: 1: 2: 3:
// 0: 1
// 1:    2
// 2: 3     4
// 3:    8     6
// 
//      0  1  2  3  4  5  6  7  8  9  10
// ija  5  6  7  9  11 0  1  0  2  1  3
// sa   X  X  X  X  X  1  2  3  4  8  6

// length of the arrays: ija[ija[0]-1]
// dimension: ija[0]-1
/**
* @brief stores a sparse matrix (nxm)
*/
class FLXLIB_EXPORT FlxMtxSpars : public FlxMtx_base {
  protected:
    tdouble* sa;
    tnlong* ija;
    FlxMtxSpars() : sa(NULL),ija(NULL) {}
    /**
    * @brief allocates just the memory ...
    */
    FlxMtxSpars(const tnlong nmax);
  public:
    const tnlong nrows() const {return ija[0]-1;};
    const tnlong ncols() const {return ija[0]-1;};

    /**
    * @brief creates a sparse matrix out of a full matrix
    */
    FlxMtxSpars(const FlxMtx& mtx);
    /**
    * @brief creates a copy of the matrix
    */
    FlxMtxSpars(const FlxMtxSpars& mtx);
    virtual ~FlxMtxSpars();
    FlxMtxSpars* copy();
    /**
    * @brief Matrix-vector product: w=M*v;
    * NOTE: does not check dimensions of v and w
    */
    virtual void MultMv(const flxVec& v, flxVec& w) const;
    /**
    * @brief multiplies a certain row with a vector
    * @param rownumber number of the corresponding row (please make sure index is valid)
    */
    virtual const tdouble MultRowV(const flxVec& v, const tnlong rownumber) const;
    /**
    * @brief solves the following problem: w=Transpose(*this)*v;
    */
    virtual void TransMultVec(const flxVec& v, flxVec& w);
    FlxMtxSpars& operator*=(const tdouble& s);
    const tdouble operator()(const tnlong& i, const tnlong& j) const;
    void add_value(const tnlong& i, const tnlong& j, const tdouble& v);

};


// Matrix - Example:
//    0: 1: 2: 3:
// 0: 1
// 1:    2
// 2: 3     4
// 3:    7     6
// 
//      0  1  2  3  4  5  6 
// ija  5  5  5  6  7  0  1
// sa   1  2  4  6     3  7
// length of the arrays: ija[ija[0]-1]
// dimension: ija[0]-1
/**
* @brief stores a sparse lower triangular matrix
*/
class FLXLIB_EXPORT FlxMtxSparsLTri : public FlxMtxSpars {
  public:
    /**
    * @brief allocates just the memory ...
    */
    FlxMtxSparsLTri(const tnlong nmax) : FlxMtxSpars(nmax) {}
    /**
    * @brief creates a sparse lower triangular matrix out of the lower triangular matrix
    */
    FlxMtxSparsLTri(const FlxMtxLTri& mtx);
    /**
    * @brief creates a copy of the matrix
    */
    FlxMtxSparsLTri(const FlxMtxSparsLTri& mtx);
    /**
    * @brief Assembles this matrix by doing a CholeskyDecomposition of Sm
    */
    FlxMtxSparsLTri(FlxMtxSparsSym& Sm);
    /**
    * @brief ... creates a diagonal matrix
    */
    FlxMtxSparsLTri(FlxMtxDiag& mtx);
    /**
    * NOTE: YOU have to make sure that nmax is set to mtx.nrows()+1
    */
    FlxMtxSparsLTri& operator=(const FlxMtxDiag& mtx);
    FlxMtxSparsLTri* copy();
    
    const tnlong get_nmax() const { return ija[ija[0]-1]; }
    
    /**
    * @brief Matrix-vector product: w=M*v;
    * NOTE: for this special case, w and v can refer to the same vector
    */
    virtual void MultMv(const flxVec& v, flxVec& w) const;
    virtual const tdouble MultRowV(const flxVec& v, const tnlong rownumber) const;
    /**
    * @brief Calculates v=L^-1*v -> inverse matrix is not explicitly computed !!!
    */
    void MultInv(const flxVec& v,flxVec& w);
    const tdouble operator()(const tnlong& i, const tnlong& j) const;
    void add_value(const tnlong& i, const tnlong& j, const tdouble& v);
    /**
    * @brief solves the following problem: v=Transpose(*this)*v;
    */
    void TransMultVec(flxVec& v);
    virtual void TransMultVec(const flxVec& v, flxVec& w);
    /**
    * @brief solves the following problem: return=D*(*this);
    * NOTE: does not check if nmax of both matrices is the same !!!
    */
//     FlxMtxSparsLTri* MultM_Mdiag(const FlxMtxDiag& D) const;
    void MultM_Mdiag(const FlxMtxDiag& D, FlxMtxSparsLTri& res) const;
    
    /**
    * @brief Calculates the determinant of the matrix (log-transform)
    */
    const tdouble det_log() const;
    
    friend class FlxMtxLTri;
    
    /**
    * @brief Doing a sparse Cholesky decomposition into saP and ijaP
    */
    static void CholeskyDec(tdoublePtr& saP, tnlong*& ijaP, FlxMtxSparsSym& mtx);
};

/**
* @brief sparse symmetric matrix - incomplete LU factorization of a sparse symmetric matrix
*/
class FLXLIB_EXPORT FlxMtxSparsSymILU : public FlxMtx_baseS {
  protected:
    tdouble* sa;
    tnlong* ija;
    FlxMtxSparsSymILU() {};
  public:
    const tnlong nrows() const {return ija[0]-1;};
    const tnlong ncols() const {return ija[0]-1;};
    
    FlxMtxSparsSymILU(const FlxMtxSparsSym& mtx, bool mod0diagentry=false);
    virtual ~FlxMtxSparsSymILU();
    /**
    * @brief NOTE: this is not a simple Matrix-Vector multiplication ... this solves a linear system of equations (THIS)*w=v
    */
    virtual void MultMv(const flxVec& v, flxVec& w) const;
    virtual void MultMv(const flxpVec& v, flxpVec& w) const;
    FlxMtx_base& operator*=(const tdouble& s);
    const tdouble operator()(const tnlong& i, const tnlong& j) const;
    void add_value(const tnlong& i, const tnlong& j, const tdouble& v);
    const bool isPositiveDefinite(tuint& r, const bool fixIt=false) { return true;};
};

/**
* @brief sparse symmetric matrix - Cholesky decomposition of a (sparse) symmetric matrix
*/
class FLXLIB_EXPORT FlxMtxSparsSymLU : public FlxMtxSparsSymILU {
  public:
    FlxMtxSparsSymLU(FlxMtxSparsSym& mtx);
    const tdouble Percentage_of_NonZero() const { return (ija[ija[0]-1]+ija[0])/(pow2(ija[0]-1.)-ija[0]+1)*2;}
};

/**
* @brief stores a full symmetric matrix
*/
class FLXLIB_EXPORT FlxMtxSym : public FlxMtx_baseS {
  private:
    tnlong msize;
    flxVec mtx;
    void preconding(const flxVec& r, flxVec& z, const int precn = 0) const;
  public:
    const tnlong nrows() const {return msize;};
    const tnlong ncols() const {return msize;};
    FlxMtxSym(const tnlong& msizeV) :msize(msizeV),mtx((msizeV*msizeV+msizeV)/2) {};
    FlxMtxSym(const FlxMtx_baseS& S);
    FlxMtxSym* copy();
    virtual void MultMv(const flxVec& v, flxVec& w) const;
    /**
    * @brief Multiplies a slice vector v (a lot of zeros in the beginning and the end) with this matrix
    * @param li lower index to start
    * @param ui upper index to stop
    * @param v the slice vector (has to be defined between li and ui)
    * @param w the result vector (has to be of the same size as the dimension of the matrix)
    */
    void MultMv_slice(const flxVec& v, flxVec& w, const tuint& li, const tuint& ui) const;
    void MultMv(const flxpVec& v, flxpVec& w) const;
    void add_mtx(const FlxMtx_baseS& m, const tdouble f);
    FlxMtxSym& operator+=(const FlxMtxSym& m);
    FlxMtxSym& operator*=(const tdouble& s) {mtx*=s; return *this; };
    tdouble& operator()(const tnlong& i, const tnlong& j);
    const tdouble operator()(const tnlong& i, const tnlong& j) const;
    void add_value(const tnlong& i, const tnlong& j, const tdouble& v);
    const bool isPositiveDefinite(tuint& r, const bool fixIt=false);
    /**
    * @brief returns a pointer to the internal storage
    *         values are stored below the diagonal !!!
    */
    tdouble* get_VecPointer() { return &mtx[0];}
    /**
    * @brief assembles vec (dim=N*N) with the elements of this matrix (returns a full matrix)
    * @param vec NxN-vector (size is NOT checked)
    */
    void get_VecPointer_full(tdouble* vec) const;
    /**
    * @brief returns the matrix-vector (constant!)
    */
    const flxVec& get_mtx_flxVec() const { return mtx; }
    /**
    * @brief assigns *this = L^T L
    */
    void assign_LTL(const FlxMtxLTri& L);
    /**
    * @brief computes the inverse of the matrix
    */
    void Invert();
    void Invert(FlxMtxLTri& Ltri);
    
    static const bool isPositiveDefinite_ext(flxVec& mtxe, tuint rows, tuint& r, const bool fixIt = false);
    friend FLXLIB_EXPORT void MtxProd_BTKB_mtx(const FlxMtx& B, const FlxMtxSym& E, FlxMtxSym& K);
    friend FLXLIB_EXPORT void MtxProd_BTKB_ltri(const FlxMtxLTri& B, const FlxMtxSym& E, FlxMtxSym& K);
    friend void MtxProd_BTKB(const FlxMtx_base& B, const FlxMtx_baseS& E, FlxMtxSym& K);
    friend void MtxProd_BTKB(const FlxMtxTransformation&, const FlxMtxSym&, FlxMtxSym&);
};

/**
* @brief stores a band-symmetric matrix
* the band-width (bsize) is defined as the additional number of entries to one side of the diagonal
* bsize must not be larger than msize-1 !!!
*/
class FLXLIB_EXPORT FlxMtxSymBand : public FlxMtx_baseS {
  private:
    tnlong msize;
    tuint bsize;
    tVec mtx;
    void preconding(const flxVec& r, flxVec& z, const int precn = 0) const;
    const tnlong countUp2Row(const tnlong R) const; // count entries up to row R ( R=0...msize); but excluding row R
  public:
    const tnlong nrows() const {return msize;}
    const tnlong ncols() const {return msize;}
    const tuint get_bsize() const { return bsize; }
    FlxMtxSymBand(const tnlong msizeV,const tuint bsizeV) : msize(msizeV),bsize(bsizeV),mtx((bsizeV*2+1)*msizeV-(bsizeV*bsizeV+bsizeV)) {};
    FlxMtxSymBand* copy();
    virtual void MultMv(const flxVec& v, flxVec& w) const;
    //void MultMv(const flxpVec& v, flxpVec& w) const;
    FlxMtxSymBand& operator+=(const FlxMtxSymBand& m);
    FlxMtxSymBand& operator+=(const FlxMtxDiag& m);
    FlxMtxSymBand& operator*=(const tdouble& s) {mtx*=s; return *this; };
    const tdouble operator()(const tnlong& i, const tnlong& j) const;
    void add_value(const tnlong& i, const tnlong& j, const tdouble& v);
    void add_mtx(const FlxMtxSymBand& m,const tdouble& s);
    void add_mtx(const FlxMtxDiag& m,const tdouble& s);
    const bool isPositiveDefinite(tuint& r, const bool fixIt=false);
    /**
    * @brief returns the matrix-vector (constant!)
    */
    const tVec& get_mtx_tVec() const { return mtx; }
    tdouble* get_VecPointer() { return &mtx[0];}
    /**
    * @brief assigns *this = L^T L
    */
    void assign_LTL(const FlxMtxLTriBand& L);
    /**
    * @brief computes the inverse of the matrix
    */
    void Invert();
};

/**
* @brief stores a diagonal matrix
*/
class FLXLIB_EXPORT FlxMtxDiag : public FlxMtx_baseS {
  private:
    tnlong msize;
    tVec mtx;
  public:
    const tnlong nrows() const {return msize;};
    const tnlong ncols() const {return msize;};

    FlxMtxDiag(const tnlong& msizeV,const tdouble defVal = ZERO) :msize(msizeV),mtx(defVal,msizeV) {};
    FlxMtxDiag(const FlxMtx_baseS& mb);
    virtual void MultMv(const flxVec& v, flxVec& w) const;
    void MultInv(const flxVec& v,flxVec& w);
    void MultMv(const flxpVec& v, flxpVec& w) const;
    FlxMtxDiag& operator+=(const FlxMtxDiag& m) { mtx+=m.mtx; return *this; };
    FlxMtxDiag& operator*=(const tdouble& s) {mtx*=s; return *this; };
    const tdouble operator()(const tnlong& i, const tnlong& j) const { if (i==j) return mtx[i]; else return ZERO; }
    tdouble& operator()(const tnlong&i) { return mtx[i]; }
    void add_value(const tnlong& i, const tnlong& j, const tdouble& v);
    const bool isPositiveDefinite(tuint& r, const bool fixIt=false);
    /**
    * @brief returns the matrix-vector (constant!)
    */
    const tVec& get_mtx_tVec() const { return mtx; }
    
    virtual FlxMtx_baseS* get_Inverse();
};

/**
* @brief the identity matrix
*/
class FLXLIB_EXPORT FlxMtxIdentity : public FlxMtx_baseS {
  private:
    tnlong msize;
  public:
    const tnlong nrows() const {return msize;};
    const tnlong ncols() const {return msize;};

    FlxMtxIdentity(const tnlong& msizeV) :msize(msizeV) {};
    virtual void MultMv(const flxVec& v, flxVec& w) const { w = v; }
    void MultMv(const flxpVec& v, flxpVec& w) const { w = v; }
    FlxMtxIdentity& operator+=(const FlxMtx_base& m);
    FlxMtxIdentity& operator*=(const tdouble& s);
    const tdouble operator()(const tnlong& i, const tnlong& j) const { if (i==j) return 1.0; else return 0.0; }
    void add_value(const tnlong& i, const tnlong& j, const tdouble& v);
    const bool isPositiveDefinite(tuint& r, const bool fixIt=false) { return true; }
    
    virtual FlxMtx_baseS* get_Inverse();
};

/**
* @brief stores a sparse symmetric matrix
*/
class FLXLIB_EXPORT FlxMtxSparsSym : public FlxMtx_baseS {
  private:
    tdouble* sa;
    tnlong* ija;
    void preconding(const flxVec& r, flxVec& z, const int precn = 0) const;
  public:
    const tnlong nrows() const {return ija[0]-1;};
    const tnlong ncols() const {return ija[0]-1;};

    FlxMtxSparsSym(const FlxMtxSym& Mtx);
    virtual ~FlxMtxSparsSym();
    virtual void MultMv(const flxVec& v, flxVec& w) const;
    void MultMv(const flxpVec& v, flxpVec& w) const;
    FlxMtx_base& operator*=(const tdouble& s);
    const tdouble operator()(const tnlong& i, const tnlong& j) const;
    void add_value(const tnlong& i, const tnlong& j, const tdouble& v);
    void set_value(const tnlong& i, const tnlong& j, const tdouble& v);
    const bool isPositiveDefinite(tuint& r, const bool fixIt=false);
    
    virtual void output_Mtx(std::ostream& sout) const;
    const tdouble Percentage_of_NonZero() const { return (ija[ija[0]-1]+ija[0])/(pow2(ija[0]-1.)-ija[0]+1)*2;}
    
    friend class FlxMtxLTri;
    friend class FlxMtxSparsSymILU;
    friend class FlxMtxSparsLTri;
};

class FLXLIB_EXPORT FlxMtxSparsSFEMSym : public FlxMtx_baseS {
  private:
    /**
    * @brief Dimension of the stiffness matrix
    */
    tnlong Kdim;
    tdouble* sa;
    FlxMtx_baseS** sb;
    tnlong* ija;
    std::map<FlxMtx_baseS*, tuint> box;
    virtual void assembleMinv(int i);
  public:
    const tnlong nrows() const {return (ija[0]-1)*Kdim;};
    const tnlong ncols() const {return nrows();};

    /**
    * @brief Assembles this matrix
    * fM: the multiply coefficients
    * KM: the matrix coefficients
    * StfMtxV: the indexed stiffness matrices as pointers
    * P: the size of the homogeneous chaos
    */
    FlxMtxSparsSFEMSym(std::valarray< tdouble >& fM, std::valarray< int >& KM, FlxMtx_baseS** StfMtxV, tulong P);
    virtual ~FlxMtxSparsSFEMSym();
    virtual void MultMv(const flxVec& v, flxVec& w) const;
    void MultMv(const flxpVec& v, flxpVec& w) const;
    FlxMtx_base& operator*=(const tdouble& s);
    const tdouble operator()(const tnlong& i, const tnlong& j) const;
    void add_value(const tnlong& i, const tnlong& j, const tdouble& v);
    const bool isPositiveDefinite(tuint& r, const bool fixIt=false);
    
    virtual const bool solve_CG(flxVec& x, const flxVec& f, tdouble &eps, tnlong& iter, const tuint pcn, const bool startZero);
    
    virtual void output_Mtx(std::ostream& sout) const;
};

class FLXLIB_EXPORT FlxMtxPrecnInvSFEMSym : public FlxMtx_baseS {
  private:
    tVec c0iiInv;
    tnlong Kdim;
    FlxMtxSym* Kinv;
  public:
    const tnlong nrows() const {return c0iiInv.size()*Kdim; }
    const tnlong ncols() const {return nrows();};

    FlxMtxPrecnInvSFEMSym(FlxMtxSparsSym& K, const tVec& c0ii);
    ~FlxMtxPrecnInvSFEMSym();
    virtual void MultMv(const flxVec& v, flxVec& w) const;
    FlxMtx_base& operator*=(const tdouble& s);
    const tdouble operator()(const tnlong& i, const tnlong& j) const;
    void add_value(const tnlong& i, const tnlong& j, const tdouble& v);
    const bool isPositiveDefinite(tuint& r, const bool fixIt=false) { return true; }
};

class FLXLIB_EXPORT FlxMtxPrecnILUSFEMSym : public FlxMtx_baseS {
  private:
    tVec c0iiInv;
    tnlong Kdim;
    FlxMtxSparsSymILU* ILU;
  public:
    const tnlong nrows() const {return c0iiInv.size()*Kdim; }
    const tnlong ncols() const {return nrows();};

    FlxMtxPrecnILUSFEMSym(FlxMtxSparsSym& K, const tVec& c0ii, bool FullDecomp=false, bool mod0diagentry=false);
    ~FlxMtxPrecnILUSFEMSym();
    virtual void MultMv(const flxVec& v, flxVec& w) const;
    void MultMv(const flxpVec& v, flxpVec& w) const;
    FlxMtx_base& operator*=(const tdouble& s);
    const tdouble operator()(const tnlong& i, const tnlong& j) const;
    void add_value(const tnlong& i, const tnlong& j, const tdouble& v);
    const bool isPositiveDefinite(tuint& r, const bool fixIt=false) { return true; }
};


/**
* @brief Matrix multiplication of type K=Trans(B)*E*B
*/
FLXLIB_EXPORT void MtxProd_BTKB_mtx(const FlxMtx& B, const FlxMtxSym& E, FlxMtxSym& K);
FlxMtxSym MtxProd_BTKB(const FlxMtx& B, const FlxMtxSym& E);
FLXLIB_EXPORT void MtxProd_BTKB_ltri(const FlxMtxLTri& B, const FlxMtxSym& E, FlxMtxSym& K);
FlxMtxSym MtxProd_BTKB(const FlxMtxLTri& B, const FlxMtxSym& E);
void MtxProd_BTKB(const FlxMtx_base& B, const FlxMtx_baseS& E, FlxMtxSym& K);
FlxMtxSym MtxProd_BTKB(const FlxMtx_base& B, const FlxMtx_baseS& E);
void MtxProd_BTKB(const FlxMtxTransformation& B, const FlxMtxSym& E, FlxMtxSym& K);
FlxMtxSym MtxProd_BTKB(const FlxMtxTransformation& B, const FlxMtxSym& E);
void MtxProd_BTKB_1D(const FlxMtx& B, const FlxMtxSym& E, FlxMtxSym& K);
/**
* @brief Matrix multiplication (for a transformation where just a part of K has to be transformed)
* @param B the transformation matrix (3x3-matrix)
* @param E the element stiffness matrix in local coordinates (2x2-matrix)
* @param K the element stiffenss matrix in global coordinates (4x4-matrix)
* @param tranformFirst if true: the first DOF should be transfomed; if false: the second DOF should be transfomed
*/
FLXLIB_EXPORT void MtxProd_BTKB_1D_part(const FlxMtx& B, const FlxMtxSym& E, FlxMtxSym& K, bool tranformFirst);

/**
* @brief Calculates v*Transpose(v)
*/
FLXLIB_EXPORT FlxMtxSym VecDyadProd(const tVec& v);


#endif // fesslix_FlxMtx_H


