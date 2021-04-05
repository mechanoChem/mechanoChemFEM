
#include <deal.II/base/config.h>

#include <deal.II/base/exceptions.h>
#include <deal.II/base/symmetric_tensor.h>

#include <deal.II/differentiation/ad.h>

#include <deal.II/physics/transformations.h>
#include "../../include/Residual.h"
#include "../../include/supplementary/supplementaryFunctions.h"

template <typename Number>
Tensor<2, 1, Number>
dediagonalize_tensor(const dealii::SymmetricTensor<2, 1, Number> &T,
		     const double /*rotation_angle*/,
		     const unsigned int /*axis*/ = 0)
{
  AssertThrow(false, ExcNotImplemented());
  return Tensor<2, 1, Number>({{T[0][0]}});
}


template <typename Number>
Tensor<2, 2, Number>
dediagonalize_tensor(const dealii::SymmetricTensor<2, 2, Number> &T,
		     const double rotation_angle,
		     const unsigned int /*axis*/ = 0)
{
  const Tensor<2, 2> R =
    dealii::Physics::Transformations::Rotations::rotation_matrix_2d(
								    rotation_angle);
  return R * T;
}


template <typename Number>
Tensor<2, 3, Number>
dediagonalize_tensor(const dealii::SymmetricTensor<2, 3, Number> &T,
		     const double       rotation_angle,
		     const unsigned int axis = 0)
{
  AssertIndexRange(axis, 3);

  Tensor<2, 3> R;
  switch (axis)
    {
    case (0):
      R = dealii::Physics::Transformations::Rotations::rotation_matrix_3d(
									  {1, 0, 0}, rotation_angle);
      break;
    case (1):
      R = dealii::Physics::Transformations::Rotations::rotation_matrix_3d(
									  {0, 1, 0}, rotation_angle);
      break;
    case (2):
      R = dealii::Physics::Transformations::Rotations::rotation_matrix_3d(
									  {0, 0, 1}, rotation_angle);
      break;
    default:
      AssertThrow(false, ExcNotImplemented());
      break;
    }
  return R * T;
}


template <int dim, typename Number>
std::array<std::pair<Number, Tensor<1, dim, Number>>, dim>
jacobi2(dealii::SymmetricTensor<2, dim, Number> A)
{
  static_assert(numbers::NumberTraits<Number>::is_complex == false,
		"This implementation of the Jacobi algorithm does "
		"not support complex numbers");

  // Sums of diagonal resp. off-diagonal elements
  Number sd, so;
  // sin(phi), cos(phi), tan(phi) and temporary storage
  Number s, c, t;
  // More temporary storage
  Number g, h, z, theta;
  // Threshold value
  Number thresh;

  // Initialize the transformation matrix as the
  // identity tensor
  dealii::Tensor<2, dim, Number> Q(
				   dealii::unit_symmetric_tensor<dim, Number>());

  // The diagonal elements of the tridiagonal matrix;
  // this will ultimately store the eigenvalues
  std::array<Number, dim> w;
  for (int i = 0; i < dim; i++)
    w[i] = A[i][i];

  // Calculate (tr(A))^{2}
  sd = trace(A);
  sd *= sd;

  // Number of iterations
  const unsigned int max_n_it = 150;
  for (unsigned int it = 0; it <= max_n_it; it++)
    {
      // Test for convergence
      so = 0.0;
      for (int p = 0; p < dim; p++)
	for (int q = p + 1; q < dim; q++)
	  so += std::abs(A[p][q]);
      if (so == 0.0)
	break;

      // Throw if no convergence is achieved within a
      // stipulated number of iterations
      if (it == max_n_it)
	{
	  AssertThrow(
		      false,
		      ExcMessage(
				 "No convergence in iterative Jacobi eigenvector algorithm.")) return std::
	    array<std::pair<Number, Tensor<1, dim, Number>>, dim>();
	}

      // Compute threshold value which dictates whether or
      // not a Jacobi rotation is performed
      const unsigned int n_it_skip = 4;
      if (it < n_it_skip)
	thresh = 0.2 * so / (dim * dim);
      else
	thresh = 0.0;

      // Perform sweep
      for (int p = 0; p < dim; p++)
	for (int q = p + 1; q < dim; q++)
	  {
	    g = 100.0 * std::abs(A[p][q]);

	    // After a given number of iterations the
	    // rotation is skipped if the off-diagonal
	    // element is small
	    if (it > n_it_skip && std::abs(w[p]) + g == std::abs(w[p]) &&
		std::abs(w[q]) + g == std::abs(w[q]))
	      {
		A[p][q] = 0.0;
	      }
	    else if (std::abs(A[p][q]) > thresh)
	      {
		// Calculate Jacobi transformation
		h = w[q] - w[p];

		// Compute surrogate for angle theta resulting from
		// angle transformation and subsequent smallest solution
		// of quadratic equation
		if (std::abs(h) + g == std::abs(h))
		  {
		    // Prevent overflow for large theta^2. This computation
		    // is the algebraic equivalent of t = 1/(2*theta).
		    t = A[p][q] / h;
		  }
		else
		  {
		    theta = 0.5 * h / A[p][q];
		    if (theta < 0.0)
		      t = -1.0 / (std::sqrt(1.0 + theta * theta) - theta);
		    else
		      t = 1.0 / (std::sqrt(1.0 + theta * theta) + theta);
		  }

		// Compute trigonometric functions for rotation
		// in such a way as to prevent overflow for
		// large theta.
		c = 1.0 / std::sqrt(1.0 + t * t);
		s = t * c;
		z = t * A[p][q];

		// Apply Jacobi transformation...
		A[p][q] = 0.0;
		w[p] -= z;
		w[q] += z;
		// ... by executing the various rotations in sequence
		for (int r = 0; r < p; r++)
		  {
		    t       = A[r][p];
		    A[r][p] = c * t - s * A[r][q];
		    A[r][q] = s * t + c * A[r][q];
		  }
		for (int r = p + 1; r < q; r++)
		  {
		    t       = A[p][r];
		    A[p][r] = c * t - s * A[r][q];
		    A[r][q] = s * t + c * A[r][q];
		  }
		for (int r = q + 1; r < dim; r++)
		  {
		    t       = A[p][r];
		    A[p][r] = c * t - s * A[q][r];
		    A[q][r] = s * t + c * A[q][r];
		  }

		// Update the eigenvectors
		for (int r = 0; r < dim; r++)
		  {
		    t       = Q[r][p];
		    Q[r][p] = c * t - s * Q[r][q];
		    Q[r][q] = s * t + c * Q[r][q];
		  }
	      }
	  }
    }

  // Structure the data to be outputted
  std::array<std::pair<Number, Tensor<1, dim, Number>>, dim> eig_vals_vecs;
  for (unsigned int e = 0; e < dim; ++e)
    {
      eig_vals_vecs[e].first = w[e];

      // The column "e" of Q contains the non-normalized
      // eigenvector associated with the eigenvalue "e"
      for (unsigned int a = 0; a < dim; ++a)
	{
	  eig_vals_vecs[e].second[a] = Q[a][e];
	}

      // Normalize
      Assert(eig_vals_vecs[e].second.norm() != 0.0, ExcDivideByZero());
      eig_vals_vecs[e].second /= eig_vals_vecs[e].second.norm();
    }
  return eig_vals_vecs;
}

template <int dim, typename Number>
std::array<std::pair<Number, Tensor<1, dim, Number>>,
           std::integral_constant<int, dim>::value>
eigenvectors2(const SymmetricTensor<2, dim, Number> &T,
             const SymmetricTensorEigenvectorMethod method)
{
  // Not much to do when there's only a single entry
  if (dim == 1)
    return jacobi2(T);

  std::array<std::pair<Number, Tensor<1, dim, Number>>, dim> eig_vals_vecs;

  if (dealii::Differentiation::AD::is_ad_number<Number>::value && dim > 1)
    {
      // If the tensor is diagonal, then we have a bit on an issue when using
      // auto-differentiable numbers. The reason for this is that all of the
      // algorithms have shortcuts by which to return result in this case.
      // This artificially decouples the eigenvalues/vectors, which upon
      // differentiation leads to the wrong result (each is, incorrectly,
      // insensitive with respect to the other). To work around this manipulate
      // tensor @p T in an objective manner: through an infinitesimal rotation we
      // make it non-diagonal (although we introduce some numerical error).
      bool is_diagonal = true;
      for (unsigned int i = 0; i < dim; ++i)
        for (unsigned int j = i + 1; j < dim; ++j)
          if (T[i][j] != 0.0)
            {
              is_diagonal = false;
              break;
            }

      // If our tensor is not diagonal, then just carry on as per usual.
      if (!is_diagonal)
        eig_vals_vecs = jacobi2(T);
      else
        {
          Assert(
            method != dealii::SymmetricTensorEigenvectorMethod::hybrid,
            ExcMessage(
              "The hybrid method cannot be used with auto-differentiable numbers "
              "when the tensor upon which an eigen-decomposition is being performed "
              "is diagonal. This is because the hybrid method immediately assumes "
              "the values of the eigenvectors (since the characteristic polynomial) "
              "is not solved, and therefore the sensitivity of the eigenvalues with "
              "respect to one another is not resolved."));

          // These parameters are heuristicaly chosen through "rigorous"
          // eye-ball analysis of the errors of tests based on
          // ad-common-tests/symmetric_tensor_functions_03.h. This checks the
          // first and second derivatives of the representation of a Neo-Hookean
          // material described by an Ogden-type model (i.e. using an
          // eigen-decomposition of the right Cauchy-Green tensor). Using this
          // comparison between the well-understood result expected from the
          // Neo-Hookean model and its Ogden equivalent, these parameters are a
          // first approximation to those required to collectively minimize
          // error in the energy values, first and second derivatives. What's
          // apparent is that all AD numbers and eigen-decomposition algorithms
          // are not made equal!
          double sf = 1.0;
          if (dealii::Differentiation::AD::is_taped_ad_number<Number>::value)
            {
              // ADOL-C taped
              if (method ==
                  dealii::SymmetricTensorEigenvectorMethod::ql_implicit_shifts)
                sf = 2e11;
              else if (method == SymmetricTensorEigenvectorMethod::jacobi)
                sf = (dim == 2 ? 1e6 : 1e9);
              else
                AssertThrow(false, ExcNotImplemented());
            }
          else if (dealii::Differentiation::AD::is_sacado_rad_number<Number>::value)
            {
              // Sacado::Rad
              if (method ==
                  dealii::SymmetricTensorEigenvectorMethod::ql_implicit_shifts)
                sf = (dim == 2 ? 1e8 : 1e9);
              else if (method == SymmetricTensorEigenvectorMethod::jacobi)
                sf = (dim == 2 ? 1e8 : 1e9);
              else
                AssertThrow(false, ExcNotImplemented());
            }
          else
            {
              // Everything else
              Assert(dealii::Differentiation::AD::is_tapeless_ad_number<Number>::value,
                     ExcInternalError());
              Assert(
                dealii::Differentiation::AD::is_sacado_dfad_number<Number>::value ||
                  dealii::Differentiation::AD::is_adolc_tapeless_number<Number>::value,
                ExcInternalError());

              if (method ==
                  dealii::SymmetricTensorEigenvectorMethod::ql_implicit_shifts)
                sf = (dim == 2 ? 1e7 : 2.5e7);
              else if (method == dealii::SymmetricTensorEigenvectorMethod::jacobi)
                sf = (dim == 2 ? 1e2 : 1e7);
              else
                AssertThrow(false, ExcNotImplemented());
            }

          using scalar_type =
            typename dealii::Differentiation::AD::ADNumberTraits<Number>::scalar_type;
          const double delta = sf * std::numeric_limits<scalar_type>::epsilon();
          //const double rotation_angle = delta * dealii::numbers::PI / 180.0;
          const double rotation_angle = 0.1 * dealii::numbers::PI / 180.0;

          if (dim == 2)
            {
              const Tensor<2, dim, Number> T_prime_ns =
                dediagonalize_tensor(
                  T, rotation_angle);

              // We can't symmetrize the tensor, otherwise the sensitivities
              // cancel out. So we take the upper triangle as an approximation
              // instead.
              // TODO[JPP]: Perform the eigen-decomposition on the non-symmetric
              //            T_prime_ns. This is, however, nontrivial to
              //            implement in this context. See:
              //            http://www.alglib.net/eigen/nonsymmetric/nonsymmetricevd.php
              //            https://groups.google.com/forum/#!topic/stan-users/QJe1TNioiyg
              SymmetricTensor<2, dim, Number> T_prime;
              for (unsigned int i = 0; i < dim; ++i)
                for (unsigned int j = i; j < dim; ++j)
                  T_prime[i][j] = T_prime_ns[i][j];

              eig_vals_vecs = jacobi2(T_prime);
            }
          else
            {
              Assert(dim == 3, ExcDimensionMismatch(dim, 3));

              SymmetricTensor<2, dim, Number> T_prime;
              Tensor<2, dim, Number>          T_prime_ns;
              for (unsigned int i = 0; i < dim; ++i)
                {
                  // This is a little bit hacky, so here's a brief explanation
                  // as to what the principal of this operation is: What we're
                  // trying to do here is perturb our tensor such that the
                  // sensitivity of the eigenvectors with respect to each other
                  // can be established. So, one at a time, we compute the
                  // perturbation of the input tensor such that the maximal
                  // number of off-diagonal entries are non-zero for any given
                  // "i". This means that we rotation not about the "ith" axis,
                  // but rather some offset of it. Note: This does NOT lead to
                  // an exact value or derivative of the eigendata being
                  // computed, so one should be aware that for this case (where
                  // the eigenvalues are equal), the linearization of any
                  // resulting quantities is only approximate.
                  const unsigned int axis = (i + 2) % 3;
                  T_prime_ns = dediagonalize_tensor(T, rotation_angle, axis);

                  // We can't symmetrize the tensor, otherwise the sensitivities
                  // cancel out. So we take the upper triangle as an
                  // approximation instead.
                  // TODO[JPP]: Keep the full row and perform the
                  //            eigen-decomposition on the
                  //            non-symmetric T_prime_ns. See the related
                  //            comment above in the 2d case.
                  for (unsigned int j = i; j < dim; ++j)
                    T_prime[i][j] = T_prime_ns[i][j];
                }
              eig_vals_vecs = jacobi2(T_prime);
            }
        }
    }
  else
    eig_vals_vecs = jacobi2(T);

  // Sort in descending order before output.
  std::sort(
    eig_vals_vecs.begin(),
    eig_vals_vecs.end(),
    dealii::internal::SymmetricTensorImplementation::SortEigenValuesVectors<dim,
                                                                    Number>());
  return eig_vals_vecs;
}

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
void Residual<T,dim>::evaluateQuadLogStress(dealii::Table<3, T >& P_FpinvT, const dealii::Table<3, T > &F, const dealii::Table<3, double>& Cpinv_conv, const dealii::Table<1,double> &alpha_conv, dealii::Table<3, T >& Cpinv, dealii::Table<1,T> &alpha, double tau_y, int currentIteration)
{
  //determine first Piola-Kirchhoff stress tensor P, using return mapping algorithm
  double kappa = lambda + 2.*mu/3.; //Bulk modulus
  dealii::SymmetricTensor<2,dim,T> be_tmp;
  dealii::SymmetricTensor<2,dim,Sacado::Fad::DFad<T> > be_tmp2;	    
  dealii::Table<2,T> tau(dim,dim), be(dim,dim), be_check(dim,dim), be_check2_tmp(dim,dim), be_check2(dim,dim), P_FpinvT_check(dim,dim);
  dealii::Table<1,T> lambda_sqr(dim);
  dealii::Table<3,T> nn(dim,dim,dim);  
  dealii::Table<1,T> eps(dim), beta(dim), Dev_eps(dim);
  T Tr_eps, q_hard = 0.;
  unsigned int n_q_points= P_FpinvT.size(0);
  double K = 1.e6; // linear hardening coefficient; perfect plasticity for now
  //double tau_y = 300.e6; // super high for testing

  for (unsigned int q=0; q<n_q_points; ++q){
    // Get trial values
    //std::cout << "be:" << std::endl;
    for (int i=0; i<dim; ++i){
      for (int j=0; j<dim; ++j){
	be[i][j] = 0.;
	for (int k=0; k<dim; ++k){
	  for (int l=0; l<dim; ++l){
	    be[i][j] += F[q][i][k]*Cpinv_conv[q][k][l]*F[q][j][l];
	  }
	}
	//std::cout << be[i][j] << std::endl;
      }
    }
    //std::cout << std::endl;
    //*
    for (int i=0; i<dim; ++i){
      for (int j=0; j<dim; ++j){
	be_tmp[i][j] = be[i][j];
	//be_tmp2[i][j] = be[i][j];
	//be_tmp2[i][j].diff(dim*i+j,9);
      }
    }
    // */
    std::array<std::pair<T,Tensor<1,dim,T> >,std::integral_constant<int,dim>::value> eig;//, eig2;
    //std::array<std::pair<Sacado::Fad::DFad<T>,Tensor<1,dim,Sacado::Fad::DFad<T> > >,std::integral_constant<int,dim>::value> eig2;
    //eig2 = eigenvectors2(be_tmp2,dealii::SymmetricTensorEigenvectorMethod::jacobi);
    //std::array<Sacado::Fad::DFad<T>,dim> eig2;
    //eig2 = eigenvalues(be_tmp2);
    //eig = dealii::eigenvectors(be_tmp);
    eig = dealii::eigenvectors(be_tmp,dealii::SymmetricTensorEigenvectorMethod::jacobi);
    /*
      if (currentIteration < 2){
      eig = dealii::eigenvectors(be_tmp,dealii::SymmetricTensorEigenvectorMethod::jacobi);
      }
      else{
      eig = dealii::eigenvectors(be_tmp);
      }
    */
    /*
    // be from eigenvalue decomposition derivative
    for (int i=0; i<dim; ++i){
      for (int j=0; j<dim; ++j){
	be_check2_tmp[i][j] = 0.;
	for (int k=0; k<dim; ++k){
	  be_check2_tmp[i][j] += eig2[k].first.val()*eig2[k].first.dx(dim*i+j);
	  //be_check2_tmp[i][j] += eig2[k].val()*eig2[k].dx(dim*i+j);
	}
      }
    }
    for (int i=0; i<dim; ++i){
      for (int j=0; j<dim; ++j){
	be_check2[i][j] = 0.5*(be_check2_tmp[i][j] + be_check2_tmp[j][i]);
      }
    }
    // */
    /*
      try{
      eig = dealii::eigenvectors(be,dealii::SymmetricTensorEigenvectorMethod::jacobi);
      //eig = dealii::eigenvectors(be);
      }
      catch(...){
      //eig = dealii::eigenvectors(be);
      //eig = dealii::eigenvectors(be,dealii::SymmetricTensorEigenvectorMethod::jacobi);
      for (int i=0; i<dim; ++i){
      for (int j=0; j<dim; ++j){
      //std::cout << be[i][j] << " ";
      std::cout << F[q][i][j] << " ";
      }
      std::cout << std::endl;
      }
      std::cout << std::endl;
      eig = dealii::eigenvectors(be,dealii::SymmetricTensorEigenvectorMethod::jacobi);
      }
    */
    //auto eig = dealii::eigenvectors(be);
    for (int i=0; i<dim; ++i){
      eps[i] = 0.5*std::log(eig[i].first);
    }
    Tr_eps = Trace(eps);
    Dev(eps,Dev_eps);
    for (int i=0; i<dim; ++i){
      beta[i] = kappa*Tr_eps + 2.*mu*Dev_eps[i];
    }
    q_hard = -K*alpha_conv[q];
    //*
    // Evaluated the yield function
    T f_tr = evaluateYield(beta,tau_y,q);

    if (f_tr < 0){
      // No additional plastic flow, trial state holds
      alpha[q] = alpha_conv[q];
    }
    else{
      // Plastic flow
      T gamma_dt = f_tr/(2.*mu + 2./3.*K); // mu:Lame, K:linear hardening coefficient
      //T gamma_dt = (1.-1.e-8)*f_tr/(2.*mu + 2./3.*K); // mu:Lame, K:linear hardening coefficient, with a little extra pullback for stability
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
	    be[i][j] += lambda_sqr[k]*eig[k].second[i]*eig[k].second[j];
	  }
	}
      }
      alpha[q] = alpha_conv[q] + gamma_dt*std::sqrt(2./3.);
    } // end plastic flow specific update
    // */
    // Form Kirchhoff stress tau
    //*
    for (int i=0; i<dim; ++i){
      for (int j=0; j<dim; ++j){
	tau[i][j] = 0.;
	for (int k=0; k<dim; ++k){
	  tau[i][j] += beta[k]*eig[k].second[i]*eig[k].second[j];
	}
      }
    }
    // */
    /*
    // be from eigenvalue decomposition
    for (int i=0; i<dim; ++i){
      for (int j=0; j<dim; ++j){
	be_check[i][j] = 0.;
	for (int k=0; k<dim; ++k){
	  be_check[i][j] += eig[k].first*eig[k].second[i]*eig[k].second[j];
	}
      }
    }
    // */
    /*
      std::cout << "be:" << std::endl;
      for (int i=0; i<dim; ++i){
      for (int j=0; j<dim; ++j){
      std::cout << be[i][j] << std::endl;
      }
      } 
      std::cout << std::endl;
      // */
    /*
      std::cout << "be_check:" << std::endl;
      for (int i=0; i<dim; ++i){
      for (int j=0; j<dim; ++j){
      std::cout << be_check[i][j] << std::endl;
      }
      } 
      std::cout << std::endl;
      // */
    /*
      std::cout << "be_check2:" << std::endl;
      for (int i=0; i<dim; ++i){
      for (int j=0; j<dim; ++j){
      std::cout << be_check2[i][j] << std::endl;
      }
      } 
      std::cout << std::endl;
      // */
    /* St. Venant-Kirchhoff, for debugging
    T trb_check = be_check[0][0] + be_check[1][1] + be_check[2][2];
    for (unsigned int i=0; i<dim; ++i){
      for (unsigned int j=0; j<dim; ++j){
	P_FpinvT[q][i][j] = 0.5*lambda*(trb_check - 3.)*F[q][i][j];
	for (unsigned int k=0; k<dim; ++k){
	  P_FpinvT[q][i][j] += mu*(be_check[i][k] - (i==k))*F[q][k][j];
	}
      }
    }
    // */
    /* St. Venant-Kirchhoff, for debugging
    T trb_check2 = be_check2[0][0] + be_check2[1][1] + be_check2[2][2];
    for (unsigned int i=0; i<dim; ++i){
      for (unsigned int j=0; j<dim; ++j){
	P_FpinvT[q][i][j] = 0.5*lambda*(trb_check2 - 3.)*F[q][i][j];
	for (unsigned int k=0; k<dim; ++k){
	  P_FpinvT[q][i][j] += mu*(be_check2[i][k] - (i==k))*F[q][k][j];
	}
	
      }
    }
    // */
    /* St. Venant-Kirchhoff, for debugging
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
    /*
      std::cout << "P:" << std::endl;
      for (int i=0; i<dim; ++i){
      for (int j=0; j<dim; ++j){
      std::cout << P_FpinvT[q][i][j] << std::endl;
      }
      } 
      std::cout << std::endl;
      // */
    /*
      std::cout << "P_check:" << std::endl;
      for (int i=0; i<dim; ++i){
      for (int j=0; j<dim; ++j){
      std::cout << P_FpinvT_check[i][j] << std::endl;
      }
      } 
      std::cout << std::endl;
      // */
    dealii::Table<2, T > Finv(dim,dim), Ftmp(dim,dim);
    for (unsigned int i=0; i<dim; ++i){
      for (unsigned int j=0; j<dim; ++j){
	Ftmp[i][j] = F[q][i][j];
      }
    }
    getInverse<T,dim>(Ftmp,Finv);
    //*
    // Get P_FpinvT: P*F_p^{-T} = tau*F^{-T}
    for (unsigned int i=0; i<dim; ++i){
      for (unsigned int j=0; j<dim; ++j){
	P_FpinvT[q][i][j]=0.0;
	for (unsigned int k=0; k<dim; ++k){
	  P_FpinvT[q][i][j] += tau[i][k]*Finv[j][k];
	}
	if (std::isinf(P_FpinvT[q][i][j].val()) || std::isnan(P_FpinvT[q][i][j].val())){
	  std::cout << "nan or inf in stress\n";
	}
      }
    }
    // */
    //*
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
    // */
  } // End loop over quad points

}

template class Residual<Sacado::Fad::DFad<double>, 1>;
template class Residual<Sacado::Fad::DFad<double>, 2>;
template class Residual<Sacado::Fad::DFad<double>, 3>;
