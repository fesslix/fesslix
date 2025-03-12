/* Fesslix - Stochastic Analysis
 * Copyright (C) 2010-2018 Wolfgang Betz
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

#include "flxMtx_Eigen_ARP.h"

#if FLX_USE_ARPACK

#include <arpack++/ardsmat.h>
#include <arpack++/ardssym.h>
#include <arpack++/ardgsym.h>
#include <arpack++/blas1c.h>
#include <arpack++/lapackc.h>


template<class ARMATRIX, class ARFLOAT>
void Solution(ARMATRIX &A, ARluSymStdEig<ARFLOAT> &Prob)
/*
  Prints eigenvalues and eigenvectors of symmetric eigen-problems
  on standard "cout" stream.
*/
{

  int     i, n, nconv, mode;
  ARFLOAT *Ax;
  ARFLOAT *ResNorm;

  n     = Prob.GetN();
  nconv = Prob.ConvergedEigenvalues();
  mode  = Prob.GetMode();

  *(GlobalVar.get_cout()) << std::endl << std::endl << "Testing ARPACK++ class ARluSymStdEig \n";
  *(GlobalVar.get_cout()) << "Real symmetric eigenvalue problem: A*x - lambda*x" << std::endl;
  switch (mode) {
  case 1:
    *(GlobalVar.get_cout()) << "Regular mode" << std::endl;
    break;
  case 3:
    *(GlobalVar.get_cout()) << "Shift and invert mode" << std::endl;
  }
  *(GlobalVar.get_cout()) << std::endl;

  *(GlobalVar.get_cout()) << "Dimension of the system            : " << n              << std::endl;
  *(GlobalVar.get_cout()) << "Number of 'requested' eigenvalues  : " << Prob.GetNev()  << std::endl;
  *(GlobalVar.get_cout()) << "Number of 'converged' eigenvalues  : " << nconv          << std::endl;
  *(GlobalVar.get_cout()) << "Number of Arnoldi vectors generated: " << Prob.GetNcv()  << std::endl;
  *(GlobalVar.get_cout()) << "Number of iterations taken         : " << Prob.GetIter() << std::endl;
  *(GlobalVar.get_cout()) << std::endl;

  if (Prob.EigenvaluesFound()) {

    // Printing eigenvalues.

    *(GlobalVar.get_cout()) << "Eigenvalues:" << std::endl;
    for (i=0; i<nconv; i++) {
      *(GlobalVar.get_cout()) << "  lambda[" << (i+1) << "]: " << Prob.Eigenvalue(i) << std::endl;
    }
    *(GlobalVar.get_cout()) << std::endl;
  }

  if (Prob.EigenvectorsFound()) {

    // Printing the residual norm || A*x - lambda*x ||
    // for the nconv accurately computed eigenvectors.

    Ax      = new ARFLOAT[n];
    ResNorm = new ARFLOAT[nconv+1];

    for (i=0; i<nconv; i++) {
      A.MultMv(Prob.RawEigenvector(i), Ax);
      axpy(n, -Prob.Eigenvalue(i), Prob.RawEigenvector(i), 1, Ax, 1);
      ResNorm[i] = nrm2(n, Ax, 1)/fabs(Prob.Eigenvalue(i));
    }

    for (i=0; i<nconv; i++) {
      *(GlobalVar.get_cout()) << "||A*x(" << (i+1) << ") - lambda(" << (i+1);
      *(GlobalVar.get_cout()) << ")*x(" << (i+1) << ")||: " << ResNorm[i] << "\n";
    }
    *(GlobalVar.get_cout()) << "\n";

    delete[] Ax;
    delete[] ResNorm;

  }

} // Solution


void MtxEigenValue_ARP(FlxMtx_baseS& Amtx, const int M, flxVec& EigenValues, std::vector< flxVec >& Eigenvectors)
{
  // Defining variables;
  const int     n = Amtx.ncols();   	// Dimension of the problem.
  tdouble* A= Amtx.get_VecPointer();   	// Pointer to an array that stores the lower triangular elements of A.
  
  ARdsSymMatrix<tdouble> matrix(n, A,'U');
  ARluSymStdEig<tdouble> dprob(M, matrix,(char*)"LM");

  // Finding eigenvalues and eigenvectors.
  dprob.FindEigenvectors();
  
  if (M != dprob.ConvergedEigenvalues() ) {
    std::ostringstream ssV;
    ssV << "Not all eigenvalues converged.";
    throw FlxException("MtxEigenValue_ARP_2", ssV.str() );
  }
  if (!dprob.EigenvaluesFound()) {
    std::ostringstream ssV;
    ssV << "Eigenvalues could not be obtained.";
    throw FlxException("MtxEigenValue_ARP_3", ssV.str() );
  }
  if (!dprob.EigenvectorsFound()) {
    std::ostringstream ssV;
    ssV << "Eigenvectors could not be obtained.";
    throw FlxException("MtxEigenValue_ARP_4", ssV.str() );
  }

  for (int i=0; i<M; i++) {
    EigenValues[M-1-i] = dprob.Eigenvalue(i);
  }
  
  for (int i=0; i<M; i++) {
    tdouble* dp = dprob.RawEigenvector(i);
    for (int j=0;j<n;j++) {
      Eigenvectors[M-1-i][j]=dp[j];
    }
  }
  
//   *(GlobalVar.get_cout()) << "Dimension of the system            : " << n              << std::endl;
//   *(GlobalVar.get_cout())<< "Number of 'requested' eigenvalues  : " << dprob.GetNev()  << std::endl;
//   *(GlobalVar.get_cout()) << "Number of 'converged' eigenvalues  : " << dprob.ConvergedEigenvalues()          << std::endl;
//   *(GlobalVar.get_cout()) << "Number of Arnoldi vectors generated: " << dprob.GetNcv()  << std::endl;
//   *(GlobalVar.get_cout()) << "Number of iterations taken         : " << dprob.GetIter() << std::endl;
//   *(GlobalVar.get_cout()) << std::endl;
  
  // Printing solution.
//   Solution(matrix, dprob);

}

void MtxEigenValue_ARP(FlxMtx_baseS& Amtx, FlxMtx_baseS& Bmtx, const int M, flxVec& EigenValues, std::vector< flxVec >& Eigenvectors)
{
  // Defining variables;
  const int n = Amtx.ncols();   	// Dimension of the problem.
  tdouble* A= Amtx.get_VecPointer();   	// Pointer to an array that stores the lower triangular elements of A.
  tdouble* B=Bmtx.get_VecPointer();
  const char uplo='U';
  
  ARdsSymMatrix<tdouble> matrixA(n, A,uplo);
  ARdsSymMatrix<tdouble> matrixB(n, B,uplo);
  ARluSymGenEig<tdouble> dprob(M, matrixA, matrixB,(char*)"LM");

  // Finding eigenvalues and eigenvectors.
  dprob.FindEigenvectors();
  
  if (M != dprob.ConvergedEigenvalues() ) {
    std::ostringstream ssV;
    ssV << "Not all eigenvalues converged.";
    throw FlxException("MtxEigenValue_ARP_302", ssV.str() );
  }
  if (!dprob.EigenvaluesFound()) {
    std::ostringstream ssV;
    ssV << "Eigenvalues could not be obtained.";
    throw FlxException("MtxEigenValue_ARP_303", ssV.str() );
  }
  if (!dprob.EigenvectorsFound()) {
    std::ostringstream ssV;
    ssV << "Eigenvectors could not be obtained.";
    throw FlxException("MtxEigenValue_ARP_304", ssV.str() );
  }

  for (int i=0; i<M; i++) {
    EigenValues[M-1-i] = dprob.Eigenvalue(i);
  }
  
  for (int i=0; i<M; i++) {
    tdouble* dp = dprob.RawEigenvector(i);
    for (int j=0;j<n;j++) {
      Eigenvectors[M-1-i][j]=dp[j];
    }
  }
}

#else

void MtxEigenValue_ARP(FlxMtx_baseS& Amtx, const int M, flxVec& EigenValues, std::vector< flxVec >& Eigenvectors)
{
  std::ostringstream ssV;
  ssV << "Fesslix has been compiled without ARPACK-support.";
  throw FlxException("MtxEigenValue_ARP_77", ssV.str() );
}

void MtxEigenValue_ARP(FlxMtx_baseS& Amtx, FlxMtx_baseS& Bmtx, const int M, flxVec& EigenValues, std::vector< flxVec >& Eigenvectors)
{
  std::ostringstream ssV;
  ssV << "Fesslix has been compiled without ARPACK-support.";
  throw FlxException("MtxEigenValue_ARP_78", ssV.str() );
}

#endif

