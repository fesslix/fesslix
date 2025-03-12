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

#include "flxMtx_Eigen_GSL.h"

#if FLX_USE_GSL

#include <gsl/gsl_math.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>

void MtxEigenValue_GSL(FlxMtx_baseS& Amtx, const int M, flxVec& EigenValues, std::vector< flxVec >& Eigenvectors)
{
  const std::size_t n = Amtx.nrows();
  // get the matrix to solve
    tdouble* data = new tdouble[n*n];
    Amtx.get_VecPointer_full(data);
    gsl_matrix_view m = gsl_matrix_view_array (data, n, n);
  // allocate some memory
    gsl_vector *eval = gsl_vector_alloc (n);
    gsl_matrix *evec = gsl_matrix_alloc (n, n);
    gsl_eigen_symmv_workspace * w = gsl_eigen_symmv_alloc (n);
  // solve the problem
    gsl_eigen_symmv (&m.matrix, eval, evec, w);
    gsl_eigen_symmv_free (w);
    gsl_eigen_symmv_sort (eval, evec, GSL_EIGEN_SORT_VAL_DESC);
  
  // assign the results
    for (int i = 0; i < M; ++i)
      {
	EigenValues[i] = gsl_vector_get (eval, i);
	gsl_vector_view evec_i = gsl_matrix_column (evec, i);
	for (tuint j=0;j<n;++j) {
	  Eigenvectors[i][j] = gsl_vector_get(&evec_i.vector,j);
	}
      }

  // free the memory
    gsl_vector_free (eval);
    gsl_matrix_free (evec);
    delete [] data;

}

const bool ensure_complex_is_real_GSL(gsl_complex cn) 
{  
  const tdouble ca = gsl_complex_abs(cn);
  const tdouble ia = fabs(GSL_IMAG(cn));
  return (ia/ca <= GlobalVar.TOL());	
}


void MtxEigenValue_GSL(FlxMtx& Amtx, const int M, flxVec& EigenValues, std::vector< flxVec >& Eigenvectors)
{
  const std::size_t n = Amtx.nrows();
  // get the matrix to solve
    tdouble* data = &(Amtx.operator()(0,0));
    gsl_matrix_view m = gsl_matrix_view_array (data, n, n);
  // allocate some memory
    gsl_vector_complex *eval = gsl_vector_complex_alloc (n);
    gsl_matrix_complex *evec = gsl_matrix_complex_alloc (n, n);
    gsl_eigen_nonsymmv_workspace * w = gsl_eigen_nonsymmv_alloc (n);
  // solve the problem
    gsl_eigen_nonsymmv (&m.matrix, eval, evec, w);
    gsl_eigen_nonsymmv_free (w);
    gsl_eigen_nonsymmv_sort (eval, evec, GSL_EIGEN_SORT_ABS_DESC);
  
  // assign the results
    for (int i = 0; i < M; ++i)
      {
	gsl_complex eval_i = gsl_vector_complex_get (eval, i);
	gsl_vector_complex_view evec_i = gsl_matrix_complex_column (evec, i);
	
	// make sure the eigenvalues are real
	if (!ensure_complex_is_real_GSL(eval_i)) {
	  throw FlxException_Crude("MtxEigenValue_GSL_1");
	}
	EigenValues[i] = GSL_REAL(eval_i);
	
	// make sure the eigenvectors are real
	  for (tuint j=0;j<n;++j) {
	    gsl_complex z = gsl_vector_complex_get(&evec_i.vector, j);
	    if (!ensure_complex_is_real_GSL(z)) {
	      throw FlxException_Crude("MtxEigenValue_GSL_2");
	    }
	    Eigenvectors[i][j] = GSL_REAL(z);
	  }
      }

  // free the memory
    gsl_vector_complex_free(eval);
    gsl_matrix_complex_free(evec);

}

void MtxEigenValue_GSL(FlxMtx_baseS& Amtx, FlxMtx_baseS& Bmtx, const int M, flxVec& EigenValues, std::vector< flxVec >& Eigenvectors, const tuint sortID)
{
  const size_t n = Amtx.nrows();
  // get the matrix to solve
    tdouble* data = new tdouble[n*n];
    tdouble* dataB = new tdouble[n*n];
    Amtx.get_VecPointer_full(data);
    Bmtx.get_VecPointer_full(dataB);
    gsl_matrix_view m = gsl_matrix_view_array (data, n, n);
    gsl_matrix_view b = gsl_matrix_view_array (dataB, n, n);
  // allocate some memory
    gsl_vector *eval = gsl_vector_alloc (n);
    gsl_matrix *evec = gsl_matrix_alloc (n, n);
    gsl_eigen_gensymmv_workspace * w = gsl_eigen_gensymmv_alloc(n);
  // solve the problem
    gsl_eigen_gensymmv (&m.matrix,&b.matrix, eval, evec, w);
    gsl_eigen_gensymmv_free (w);
    if (sortID==1) {
      gsl_eigen_gensymmv_sort (eval, evec, GSL_EIGEN_SORT_VAL_DESC);
    } else {
      gsl_eigen_gensymmv_sort (eval, evec, GSL_EIGEN_SORT_VAL_ASC);
    }
  
  // assign the results
    for (int i = 0; i < M; ++i)
      {
	EigenValues[i] = gsl_vector_get (eval, i);
	gsl_vector_view evec_i = gsl_matrix_column (evec, i);
	for (tuint j=0;j<n;++j) {
	  Eigenvectors[i][j] = gsl_vector_get(&evec_i.vector,j);
	}
      }

  // free the memory
    gsl_vector_free (eval);
    gsl_matrix_free (evec);
    delete [] data;
    delete [] dataB;
}

void MtxEigenValue_GSLstab(FlxMtx_baseS& Amtx, const int M, flxVec& EigenValues, std::vector< flxVec >& Eigenvectors)
{
  throw FlxException_NotImplemented("MtxEigenValue_GSLstab");
}

void MtxEigenValue_GSLstab(FlxMtx_baseS& Amtx, FlxMtx_baseS& Bmtx, const int M, flxVec& EigenValues, std::vector< flxVec >& Eigenvectors)
{
  GlobalVar.slog(4) << std::endl;
  GlobalVar.slog(4) << "Solving eigenvalue problem - stabilized version" << std::endl;
  GlobalVar.slog(4) << "-----------------------------------------------" << std::endl;
  GlobalVar.slog(4) << "  Problem: Bx = lMx" << std::endl << std::endl;
  const size_t n = Amtx.nrows();
  // get the matrix to solve
    tdouble* data = new tdouble[n*n];
    tdouble* dataB = new tdouble[n*n];
    Amtx.get_VecPointer_full(data);
    Bmtx.get_VecPointer_full(dataB);
    gsl_matrix_view m = gsl_matrix_view_array (data, n, n);
    gsl_matrix_view b = gsl_matrix_view_array (dataB, n, n);
  // allocate some memory
    gsl_vector *eval = gsl_vector_alloc (n);
    gsl_matrix *evec = gsl_matrix_alloc (n, n);
  // stabilize matrix B
    const tdouble CONDMIN = 1e-7;
    gsl_eigen_symmv_workspace * wB = gsl_eigen_symmv_alloc(n);
    gsl_eigen_symmv(&b.matrix,eval, evec, wB);
    gsl_eigen_symmv_free(wB);
    gsl_eigen_symmv_sort (eval,evec,GSL_EIGEN_SORT_VAL_DESC);
    tdouble maxev = gsl_vector_get (eval, 0);
    tdouble minev = maxev;
    tdouble negev = 0.0;
    tuint corn = 0;
    for (tuint i = 1; i < n; ++i) {
      const tdouble cev = gsl_vector_get (eval, i);
      if (cev<0.0) {
	gsl_vector_set(eval,i,CONDMIN*maxev);
	++corn;
	negev = cev;
      } else {
	minev = cev;
	if (cev/maxev < CONDMIN) {
	  gsl_vector_set(eval,i,CONDMIN*maxev);
	  ++corn;
	}
      }
    }
    GlobalVar.slog(4) << "  Matrix M" << std::endl;
    GlobalVar.slog(4) << "  --------" << std::endl;
    GlobalVar.slog(4) << "  condition number:      " << GlobalVar.Double2String(maxev/minev) << std::endl;
    GlobalVar.slog(4) << "  corrected eigenvalues: " << corn << std::endl;
    GlobalVar.slog(4) << "  largest positive ev.:  " << maxev << std::endl;
    GlobalVar.slog(4) << "  smallest positive ev.: " << minev << std::endl;
    if (negev<0.0) {
    GlobalVar.slog(4) << "  smallest negative ev.: " << negev << std::endl;
    }
    if (corn>0) {
    GlobalVar.slog(4) << "  new condition number:  " << GlobalVar.Double2String(1.0/CONDMIN) << std::endl;
    }
    for (tuint k=0;k<n;++k) {
      for (tuint j=0;j<n;++j) {
	tdouble t = 0.0;
	for (tuint i=0;i<n;++i) {
	  t += gsl_matrix_get(evec,k,n-i-1)*gsl_matrix_get(evec,j,n-i-1)*gsl_vector_get (eval,n-i-1);
	}
	gsl_matrix_set(&b.matrix,k,j,t);
      }
    }

  // solve the problem
    gsl_eigen_gensymmv_workspace * w = gsl_eigen_gensymmv_alloc(n);
    gsl_eigen_gensymmv (&m.matrix,&b.matrix, eval, evec, w);
    gsl_eigen_gensymmv_free (w);
    gsl_eigen_gensymmv_sort (eval, evec, GSL_EIGEN_SORT_VAL_DESC);
  
  // assign the results
    for (int i = 0; i < M; ++i)
      {
	EigenValues[i] = gsl_vector_get (eval, i);
	gsl_vector_view evec_i = gsl_matrix_column (evec, i);
	for (tuint j=0;j<n;++j) {
	  Eigenvectors[i][j] = gsl_vector_get(&evec_i.vector,j);
	}
      }

  // free the memory
    gsl_vector_free (eval);
    gsl_matrix_free (evec);
    delete [] data;
    delete [] dataB;
  GlobalVar.slog(4) << std::endl;
}

#else

void MtxEigenValue_GSL(FlxMtx& Amtx, const int M, flxVec& EigenValues, std::vector< flxVec >& Eigenvectors)
{
  std::ostringstream ssV;
  ssV << "Fesslix has been compiled without GSL-support.";
  throw FlxException("MtxEigenValue_GSL_76", ssV.str() );
}

void MtxEigenValue_GSL(FlxMtx_baseS& Amtx, const int M, flxVec& EigenValues, std::vector< flxVec >& Eigenvectors)
{
  std::ostringstream ssV;
  ssV << "Fesslix has been compiled without GSL-support.";
  throw FlxException("MtxEigenValue_GSL_77", ssV.str() );
}

void MtxEigenValue_GSL(FlxMtx_baseS& Amtx, FlxMtx_baseS& Bmtx, const int M, flxVec& EigenValues, std::vector< flxVec >& Eigenvectors, const tuint sortID)
{
  std::ostringstream ssV;
  ssV << "Fesslix has been compiled without GSL-support.";
  throw FlxException("MtxEigenValue_GSL_78", ssV.str() );
}

void MtxEigenValue_GSLstab(FlxMtx_baseS& Amtx, const int M, flxVec& EigenValues, std::vector<flxVec>& Eigenvectors)
{
  std::ostringstream ssV;
  ssV << "Fesslix has been compiled without GSL-support.";
  throw FlxException("MtxEigenValue_GSLstab_77", ssV.str() );
}

void MtxEigenValue_GSLstab(FlxMtx_baseS& Amtx, FlxMtx_baseS& Bmtx, const int M, flxVec& EigenValues, std::vector<flxVec>& Eigenvectors)
{
  std::ostringstream ssV;
  ssV << "Fesslix has been compiled without GSL-support.";
  throw FlxException("MtxEigenValue_GSLstab_78", ssV.str() );
}

#endif

