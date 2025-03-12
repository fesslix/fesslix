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

#include "flxMtx_Eigen.h"



void MtxEigenValue(FlxMtx_baseS& Amtx, const tuint M, flxVec& EigenValues, std::vector< flxVec >& Eigenvectors, const int Mode)
{
  const tuint N = Amtx.nrows();
  if (N < M) {
    std::ostringstream ssV;
    ssV << "Cannot compute more Eigenvalues (" << M << ") than number of DOFs in the system(" <<  N << ").";
    throw FlxException("MtxEigenValue_G_0.0", ssV.str());
  }
  
  #if FLX_DEBUG
    if (EigenValues.get_N() != M) {
      std::ostringstream ssV;
      ssV << "Problem with size of vector. (" << EigenValues.get_N() << "; " << M << ")";
      throw FlxException("MtxEigenValue_G_0.1", ssV.str() );
    }
    for (tuint i=0;i<M;i++) {
      if (Eigenvectors[i].get_N()!=N) {
	std::ostringstream ssV;
	ssV << "Problem with size of vector.";
	throw FlxException("MtxEigenValue_G_0.2", ssV.str() );
      }
    }
  #endif
  
  if (Mode == 1) {
    MtxEigenValue_ARP(Amtx,M,EigenValues,Eigenvectors);
  } else if (Mode == 2) {
    MtxEigenValue_GSL(Amtx,M,EigenValues,Eigenvectors);
  } else if (Mode == 3) {
    MtxEigenValue_GSLstab(Amtx,M,EigenValues,Eigenvectors);
  } else {
    throw FlxException_Crude("MtxEigenValue_G_0.3");
  }
  
  EV_orientation(M,Eigenvectors);
  
}

void MtxEigenValue(FlxMtx& Amtx, const tuint M, flxVec& EigenValues, std::vector< flxVec >& Eigenvectors, const int Mode)
{
  const tuint N = Amtx.nrows();
  #if FLX_DEBUG
    if (N!=Amtx.ncols()) {
      throw FlxException_Crude("MtxEigenValue_G_2.1");
    }
  #endif
  if (N < M) {
    std::ostringstream ssV;
    ssV << "Cannot compute more Eigenvalues (" << M << ") than number of DOFs in the system(" <<  N << ").";
    throw FlxException("MtxEigenValue_G_2.2", ssV.str());
  }
  
  #if FLX_DEBUG
    if (EigenValues.get_N() != M) {
      std::ostringstream ssV;
      ssV << "Problem with size of vector. (" << EigenValues.get_N() << "; " << M << ")";
      throw FlxException("MtxEigenValue_G_2.1", ssV.str() );
    }
    for (tuint i=0;i<M;i++) {
      if (Eigenvectors[i].get_N()!=N) {
	std::ostringstream ssV;
	ssV << "Problem with size of vector.";
	throw FlxException("MtxEigenValue_G_2.2", ssV.str() );
      }
    }
  #endif
  
  if (Mode == 1) {
    throw FlxException_NotImplemented("MtxEigenValue_G_2.3");
  } else if (Mode == 2) {
    MtxEigenValue_GSL(Amtx,M,EigenValues,Eigenvectors);
  } else {
    throw FlxException_Crude("MtxEigenValue_G_2.3");
  }
  
  EV_orientation(M,Eigenvectors);
  
}

void MtxEigenValue(FlxMtx_baseS& Amtx, FlxMtx_baseS& Bmtx, const tuint M, flxVec& EigenValues, std::vector< flxVec >& Eigenvectors, const int Mode)
{
  const tuint N = Amtx.nrows();
  if (N < M) {
    std::ostringstream ssV;
    ssV << "Cannot compute more Eigenvalues (" << M << ") than number of DOFs in the system (" <<  N << ").";
    throw FlxException("MtxEigenValue_G_1.0", ssV.str());
  }
  
  #if FLX_DEBUG
    if (EigenValues.get_N() != M) {
      std::ostringstream ssV;
      ssV << "Problem with size of vector. (" << EigenValues.get_N() << "; " << M << ")";
      throw FlxException("MtxEigenValue_G_1.1", ssV.str() );
    }
    for (tuint i=0;i<M;i++) {
      if (Eigenvectors[i].get_N()!=Amtx.ncols()) {
	std::ostringstream ssV;
	ssV << "Problem with size of vector.";
	throw FlxException("MtxEigenValue_G_1.2", ssV.str() );
      }
    }
    if (Amtx.ncols()!=Bmtx.ncols()) {
      std::ostringstream ssV;
      ssV << "Matrix sizes do not match";
      throw FlxException("MtxEigenValue_G_1.3", ssV.str() );
    }
  #endif
  
  if (Mode == 1) {
    MtxEigenValue_ARP(Amtx,Bmtx,M,EigenValues,Eigenvectors);
  } else if (Mode == 2) {
    MtxEigenValue_GSL(Amtx,Bmtx,M,EigenValues,Eigenvectors,1);
  } else if (Mode == 3) {
    MtxEigenValue_GSLstab(Amtx,Bmtx,M,EigenValues,Eigenvectors);
  } else {
    throw FlxException_Crude("MtxEigenValue_G_1.4");
  }
  
  EV_orientation(M,Eigenvectors);
}

void EV_orientation(const tuint M, std::vector< flxVec >& Eigenvectors)
{
  // define orientation of eigenvectors
  for (tuint i=0;i<M;++i) {
    flxVec& EV = Eigenvectors[i];
    const tdouble mN = EV.get_absMean();
    const tdouble r = mN*0.1;
    const tuint N = EV.get_N();
    const tdouble* EVp = EV.get_tmp_vptr_const();
    #if FLX_DEBUG
    bool b=false;
    #endif
    for (tuint j=0;j<N;++j) {
      if (fabs(EVp[j])>=r) {
	if (EVp[j]<0.) {
	  EV *= -1;
	}
	#if FLX_DEBUG
	b=true;
	#endif
	break;
      }
    }
    #if FLX_DEBUG
    if (b==false) {
      throw FlxException_Crude("EV_orientation");
    }
    #endif
  }
}


void flxeigen_logInfo(std::ostream& lout)
{
  lout << " FlxEigen: " << std::endl;
  lout << "   compiled with the options ..." << std::endl;
  // FLXDEBUG
    lout << "     FLX_DEBUG                     ";
    #if FLX_DEBUG
      lout << "ON";
    #else 
      lout << "OFF";
    #endif
    lout << std::endl;
  // FLXDEBUG_COUT
    lout << "     FLX_DEBUG_COUT                ";
    #if FLX_DEBUG_COUT
      lout << "ON";
    #else 
      lout << "OFF";
    #endif
    lout << std::endl;
  // FLX_USE_ARPACK
    lout << "     FLX_USE_ARPACK                ";
    #if FLX_USE_ARPACK
      lout << "ON";
    #else 
      lout << "OFF";
    #endif  
    lout << std::endl;
  // FLX_USE_GSL
    lout << "     FLX_USE_GSL                   ";
    #if FLX_USE_GSL 
      lout << "ON";
    #else 
      lout << "OFF";
    #endif  
    lout << std::endl;
}
