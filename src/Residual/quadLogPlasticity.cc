#include "../../include/Residual.h"
#include "../../include/supplementary/supplementaryFunctions.h"

template <class T, int dim>
void Residual<T,dim>::Dev(const dealii::Table<1, T >& beta, dealii::Table<1, T >& Dev_beta)
{
  //evaluate the norm of the deviatoric part of the vector of principal components
  // Dev(beta_i) = beta_i - avg(beta)
  Dev_beta.reinit(dealii::TableIndices<1>(dim));
  T avg_beta = 0.;
  for (unsigned int i=0; i<dim; ++i){
    avg_beta += 1./dim*beta[i];
  }
  for (unsigned int i=0; i<dim; ++i){
    Dev_beta[i] = beta[i] - avg_beta;
  }
}

template <class T, int dim>
T Residual<T,dim>::Trace(const dealii::Table<1, T >& tensor_1)
{
  T out = 0.;
  for (int i=0; i<dim; ++i){
    out += tensor_1[i];
  }
  return out;
}

template <class T, int dim>
void Residual<T,dim>::evaluateQuadLogStress(dealii::Table<3, T >& P_FpinvT, const dealii::Table<3, T > &F, const dealii::Table<3, T >& Cpinv_conv, const dealii::Table<1,T> &alpha_conv, dealii::Table<3, T >& Cpinv, dealii::Table<1,T> &alpha)
{
  //determine first Piola-Kirchhoff stress tensor P, using return mapping algorithm
  double kappa = lambda + 2.*mu/3.; //Bulk modulus
  dealii::Table<2,T> be(dim,dim), tau(dim,dim);
  dealii::Table<1,T> lambda_sqr(dim);
  dealii::Table<3,T> nn(dim,dim,dim);  
  dealii::Table<1,T> eps(dim), beta(dim), Dev_eps(dim);
  T Tr_eps, q_hard = 0.;
  unsigned int n_q_points= P_FpinvT.size(0);
  double K = 0.; // linear hardening coefficient; perfect plasticity for now
  double tau_y = 1.e15; // super high for testing

  for (unsigned int q=0; q<n_q_points; ++q){
    // Get trial values
    //std::cout << "be:\n";
    for (int i=0; i<dim; ++i){
      for (int j=0; j<dim; ++j){
	be[i][j] = 0.;
	for (int k=0; k<dim; ++k){
	  be[i][j] += F[q][i][k]*F[q][j][k];
	    //for (int l=0; l<dim; ++l){
	    //be[i][j] += F[q][i][k]*Cpinv_conv[q][k][l]*F[q][j][l];
	    //}
	}
	//std::cout << be[i][j].val() << " ";
      }
      //std::cout << "\n";
    }
    //std::cout << "\n";

    if (false){//be[0][1]*be[0][1] + be[0][2]*be[0][2] + be[1][2]*be[1][2] == 0){
    //if (true){
      //std::cout << "diagonal\n";
      T det_be = be[0][0]*(be[1][1]*be[2][2] - be[2][1]*be[1][2]) -
	be[0][1]*(be[1][0]*be[2][2] - be[2][0]*be[1][2]) +
	be[0][2]*(be[1][0]*be[2][1] - be[2][0]*be[1][1]);
      for (int i=0; i<dim; ++i){
	for (int j=0; j<dim; ++j){
	  tau[i][j] = mu*(be[i][j] - (i==j)) + 0.5*lambda*std::log(det_be)*(i==j);
	//std::cout << tau[i][j] << " ";
	}
      }
      //std::cout << "\n";
    }
    else{
      spec_decomp_sym3x3(be,lambda_sqr,nn); // get the spectral decomposition of be (be careful with nn)
      //std::cout << lambda_sqr[0].dx(0) << " " << lambda_sqr[1].dx(0) << " " << lambda_sqr[2].dx(0) << "\n";
      //std::cout << "Eigenvalues: " << lambda_sqr[0] << " " << lambda_sqr[1] << " " << lambda_sqr[2] << "\n";
      dealii::Table<2,T> be_check(dim,dim);
      for (int i=0; i<dim; ++i){
	for (int j=0; j<dim; ++j){
	  be_check[i][j] = 0.;
	  for (int k=0; k<dim; ++k){
	    be_check[i][j] += lambda_sqr[k]*nn[k][i][j];
	  }
	  if (false){//std::abs(be_check[i][j] - be[i][j]) > 1e-12){
	    std::cout << "Spectral decomposition: " << be_check[i][j] << "\n----------------------: " << be[i][j] << std::endl;
	  }
	}
      }      
      for (int i=0; i<dim; ++i){
	eps[i] = 0.5*std::log(lambda_sqr[i]);
      }
      Tr_eps = Trace(eps);
      Dev(eps,Dev_eps);
      for (int i=0; i<dim; ++i){
	beta[i] = kappa*Tr_eps + 2.*mu*Dev_eps[i];
      }
      q_hard = -K*alpha_conv[q];
    
      // Evaluated the yield function
      T f_tr = evaluateYield(beta,tau_y,q);

      if (f_tr < 0){
	// No additional plastic flow, trial state holds
	alpha[q] = alpha_conv[q];
      }
      else{
	// Plastic flow
	/*
	  T gamma_dt = f_tr/(2*mu + 2./3.*K); // mu:Lame, K:linear hardening coefficient
	  dealii::Table<1,T> Dev_beta(dim);
	  Dev(beta,Dev_beta); 
	  T Dev_beta_norm = norm(Dev_beta);
	  for (int i=0; i<dim; ++i){
	  eps[i] -= gamma_dt*Dev_beta[i]/Dev_beta_norm;
	  lambda_sqr[i] = std::exp(2.*eps[i]);
	  }
	  // With updated stretches and log strain, update beta (stress) and be
	  Tr_eps = Trace(eps);
	  Dev(eps,Dev_eps);
	  for (int i=0; i<dim; ++i){
	  beta[i] = kappa*Tr_eps + 2.*mu*Dev_eps[i];
	  }
	  for (int i=0; i<dim; ++i){
	  for (int j=0; j<dim; ++j){
	  be[i][j] = 0.;
	  for (int k=0; k<dim; ++k){
	  be[i][j] += lambda_sqr[k]*nn[k][i][j];
	  }
	  }
	  }
	  alpha[q] = alpha_conv[q] + gamma_dt*std::sqrt(2./3.);
	*/
      } // end plastic flow specific update

      // Form Kirchhoff stress tau
      for (int i=0; i<dim; ++i){
	for (int j=0; j<dim; ++j){
	  tau[i][j] = 0.;
	  for (int k=0; k<dim; ++k){
	    tau[i][j] += beta[k]*nn[k][i][j];
	  }
	}
      }
      /*
      T trb = be_check[0][0] + be_check[1][1] + be_check[2][2];
      for (unsigned int i=0; i<dim; ++i){
	for (unsigned int j=0; j<dim; ++j){
	  P_FpinvT[q][i][j] = 0.5*lambda*(trb - 3.)*F[q][i][j];
	  for (unsigned int k=0; k<dim; ++k){
	    P_FpinvT[q][i][j] += mu*(be_check[i][k] - (i==k))*F[q][k][j];
	  }
	}
      }
      // */
      //*
      T trb = be[0][0] + be[1][1] + be[2][2];
      for (unsigned int i=0; i<dim; ++i){
	for (unsigned int j=0; j<dim; ++j){
	  P_FpinvT[q][i][j] = 0.5*lambda*(trb - 3.)*F[q][i][j];
	  for (unsigned int k=0; k<dim; ++k){
	    P_FpinvT[q][i][j] += mu*(be[i][k] - (i==k))*F[q][k][j];
	  }
	}
      }
      // */

    }

    dealii::Table<2, T > Finv(dim,dim), Ftmp(dim,dim);
    for (unsigned int i=0; i<dim; ++i){
      for (unsigned int j=0; j<dim; ++j){
	Ftmp[i][j] = F[q][i][j];
      }
    }
    /*
    getInverse<T,dim>(Ftmp,Finv);
    // Get P_FpinvT: P*F_p^{-T} = tau*F^{-T}
    for (unsigned int i=0; i<dim; ++i){
      for (unsigned int j=0; j<dim; ++j){
	P_FpinvT[q][i][j]=0.0;
	for (unsigned int k=0; k<dim; ++k){
	  P_FpinvT[q][i][j] += tau[i][k]*Finv[j][k];
	}
      }
    }
    */
    /*
    // Update Cp_inv
    for (int i=0; i<dim; ++i){
      for (int j=0; j<dim; ++j){
	Cpinv[q][i][j] = 0.;
	for (int k=0; k<dim; ++k){
	  for (int l=0; l<dim; ++l){
	    Cpinv[q][i][j] += Finv[i][k]*be[k][l]*Finv[j][l];
	  }
	}
      }
    }
    */
  } // End loop over quad points

}

template class Residual<Sacado::Fad::DFad<double>, 1>;
template class Residual<Sacado::Fad::DFad<double>, 2>;
template class Residual<Sacado::Fad::DFad<double>, 3>;
