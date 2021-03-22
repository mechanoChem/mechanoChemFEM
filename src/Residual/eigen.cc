#include "../../include/Residual.h"

template <class T, int dim>
T Residual<T,dim>::det(const dealii::Table<2, T >& b)
{
  if (dim == 1){
    return b[0][0];
  }
  else if (dim == 2){
    return b[0][0]*b[1][1] - b[1][0]*b[0][1];
  }
  else { // dim == 3
    return b[0][0]*(b[1][1]*b[2][2] - b[2][1]*b[1][2]) -
      b[0][1]*(b[1][0]*b[2][2] - b[2][0]*b[1][2]) +
      b[0][2]*(b[1][0]*b[2][1] - b[2][0]*b[1][1]);
  }
}

template <class T, int dim>
void Residual<T,dim>::spec_decomp_sym3x3(const dealii::Table<2, T >& be_tmp, dealii::Table<1, T >& lambda_sqr, dealii::Table<3, T >& nn)
{
  // Perform the spectral decomposition on a real symmetric 3x3 matrix, where be = \sum_{A=1}^3 \lambda_A^2 m_a, and m_a = n_a \outer_product n_a
  
  // Eigenvalue algorithm: Oliver Smith, "Eigenvalues of a Symmetric 3 x 3 Matrix," https://dl.acm.org/doi/pdf/10.1145/355578.366316
  nn.reinit(TableIndices<3>(3,3,3));
  lambda_sqr.reinit(TableIndices<1>(3));

  // To help the automatic differentiation, it might be necessary to enforce symmetry in this way:
  dealii::Table<2, T > be(3,3);
  dealii::Table<2,double> pert_tmp(3,3), pert(3,3);
  /*
  for (int i=0; i<3; ++i){
    for (int j=0; j<3; ++j){
      pert_tmp[i][j] = 1.e-8*(static_cast <double> (rand())/(static_cast <double>(RAND_MAX))-0.5);
    }
  }
  */
  for (int i=0; i<3; ++i){
    for (int j=0; j<3; ++j){
      be[i][j] = 0.5*(be_tmp[i][j] + be_tmp[j][i]) + 1.e-8*(static_cast <double> (rand())/(static_cast <double>(RAND_MAX))-0.5)*(i==j);
      //pert[i][j] = 0.5*(pert_tmp[i][j] + pert_tmp[j][i]);
    }
  }
  /*
  if (be[0][1]*be[0][1] + be[0][2]*be[0][2] + be[1][2]*be[1][2] < 1e-12){ 
    for (int i=0; i<3; ++i){
      for (int j=0; j<3; ++j){
	be[i][j] += pert[i][j];
      }
    }
  }
  */
  // Following a tip from the Wikipedia algorith, and also to help the automatic differentiation, treat diagonal as special case:
  //*
  if(be[0][1]*be[0][1] + be[0][2]*be[0][2] + be[1][2]*be[1][2] < 1e-15){
    lambda_sqr[0] = be[0][0];
    lambda_sqr[1] = be[1][1];
    lambda_sqr[2] = be[2][2];
  }
  else{
    // */
    T m = 1./3.*(be[0][0] + be[1][1] + be[2][2]); // 3m = trace(be)
    dealii::Table<2, T > B(3,3);
    for (int i=0; i<3; ++i){
      for (int j=0; j<3; ++j){
	B[i][j] = be[i][j] - m*(i==j);
      }
    }
    T q = 0.5*det(B);
    /*
      if (std::abs(q) < 1.e-13){
      //std::cout << "Perturb...\n";
      for (int i=0; i<3; ++i){
      for (int j=0; j<3; ++j){
      be[i][j] += pert[i][j];
      }
      }
      m = 1./3.*(be[0][0] + be[1][1] + be[2][2]); // 3m = trace(be)
      for (int i=0; i<3; ++i){
      for (int j=0; j<3; ++j){
      B[i][j] = be[i][j] - m*(i==j);
      }
      }
      q = 0.5*det(B);
      }
    */
    T p = 1./6.*(B[0][0]*B[0][0] + B[1][1]*B[1][1] + B[2][2]*B[2][2] +
		 2.*B[0][1]*B[0][1] + 2.*B[0][2]*B[0][2] + 2.*B[1][2]*B[1][2]);
    T phi;
    phi = 1./3.*std::acos(q/std::sqrt(p*p*p));
    //T phi = 1./3.*std::atan2(std::sqrt(p*p*p - q*q),q);

    lambda_sqr[0] = m + 2.*std::sqrt(p)*std::cos(phi);
    lambda_sqr[1] = m - std::sqrt(p)*(std::cos(phi) + std::sqrt(3.)*std::sin(phi));
    lambda_sqr[2] = m - std::sqrt(p)*(std::cos(phi) - std::sqrt(3.)*std::sin(phi));
  }
  // Eigenvector algorithm: Simo & Hughes, Computational Inelasticity, chapter: "Nonlinear Continuum Mechanics and Phenomenological Plasticity Models"
  if (false){//be[0][1]*be[0][1] + be[0][2]*be[0][2] + be[1][2]*be[1][2] == 0){
    for (int i=0 ;i<3; ++i){
      for (int j=0; j<3; ++j){
	for (int k=0; k<3; ++k){
	  nn[i][j][k] = (i==j)*(i==k);
	}
      }
    }
  }
  else{
    for (int i1=0;i1<3; ++i1){
      int i2 = (i1+1)%3;
      int i3 = (i1+2)%3;
      for (int j=0; j<3; ++j){
	for (int k=0; k<3; ++k){
	  nn[i1][j][k] = 0.;
	  for (int l=0; l<3; ++l){
	    nn[i1][j][k] += (be[j][l] - lambda_sqr[i2]*(j==l))*(be[l][k] - lambda_sqr[i3]*(l==k))/
	      ((lambda_sqr[i1] - lambda_sqr[i2])*(lambda_sqr[i1] - lambda_sqr[i3]));
	  }
	}
      }
    }
  }

  // 3 distinct eigenvalues
  /*
  for (int i1=0;i1<3; ++i1){
    int i2 = (i1+1)%3;
    int i3 = (i1+2)%3;
    for (int j=0; j<3; ++j){
      for (int k=0; k<3; ++k){
	nn[i1][j][k] = 0.;
	for (int l=0; l<3; ++l){
	  nn[i1][j][k] += (be[j][l] - lambda_sqr[i2]*(j==l) + 1.e-8)*(be[l][k] - lambda_sqr[i3]*(l==k) + 1.e-8)/
	    ((lambda_sqr[i1] - lambda_sqr[i2] + 1.e-8)*(lambda_sqr[i1] - lambda_sqr[i3] + 1.e-8));
	}
      }
    }
  }
  */
  /*
  if ((std::abs(lambda_sqr[0] - lambda_sqr[1]) > 1e-12) && (std::abs(lambda_sqr[0] - lambda_sqr[2]) > 1e-12)){
    // 3 distinct eigenvalues
    for (int i1=0;i1<3; ++i1){
      int i2 = (i1+1)%3;
      int i3 = (i1+2)%3;
      for (int j=0; j<3; ++j){
	for (int k=0; k<3; ++k){
	  nn[i1][j][k] = 0.;
	  for (int l=0; l<3; ++l){
	    nn[i1][j][k] += (be[j][l] - lambda_sqr[i2]*(j==l))*(be[l][k] - lambda_sqr[i3]*(l==k))/
	      ((lambda_sqr[i1] - lambda_sqr[i2])*(lambda_sqr[i1] - lambda_sqr[i3]));
	  }
	}
      }
    }
  }
  else if ((std::abs(lambda_sqr[0] - lambda_sqr[1]) <= 1e-12) && (std::abs(lambda_sqr[0] - lambda_sqr[2]) <= 1e-12)){
    // All equal
    for (int k=0; k<3; ++k){
      for (int i=0; i<3; ++i){
	for (int j=0; j<3; ++j){
	  nn[k][i][j] = (k==i)*(k==j);
	}
      }
    }
  }
  else{ //one unique (make it the first in lambda_sqr)
    int i1; // Index of unique eigval
    if (std::abs(lambda_sqr[0] - lambda_sqr[1]) <= 1e-12){
      i1 = 2;
    }
    else if (std::abs(lambda_sqr[0] - lambda_sqr[2]) <= 1e-12){
      i1 = 1;
    }
    else{
      i1 = 0;
    }
    int i2 = (i1+1)%3;
    int i3 = (i1+2)%3;
    for (int j=0; j<3; ++j){
      for (int k=0; k<3; ++k){
	nn[i1][j][k] = 0.;
	for (int l=0; l<3; ++l){
	  nn[i1][j][k] += (be[j][l] - lambda_sqr[i2]*(j==l))*(be[l][k] - lambda_sqr[i2]*(l==k))/
	    std::pow(lambda_sqr[i1] - lambda_sqr[i2],2);
	}
	nn[i2][j][k] = 0.5*((j==k) - nn[i1][j][k]); // Note: need to update this to something more accurate, but this works 
	nn[i3][j][k] = nn[i2][j][k];
      }
    }
    
  }
  */
  
}

template class Residual<Sacado::Fad::DFad<double>, 1>;
template class Residual<Sacado::Fad::DFad<double>, 2>;
template class Residual<Sacado::Fad::DFad<double>, 3>;
