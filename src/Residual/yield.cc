#include "../../include/Residual.h"

template <class T, int dim>
T Residual<T,dim>::norm(const dealii::Table<1, T >& tensor_1)
{
  T tensor_norm = 0.;
  for (int i=0; i<dim; ++i){
    tensor_norm += tensor_1[i]*tensor_1[i];
  }
  return std::sqrt(tensor_norm);
  
}

template <class T, int dim>
T Residual<T,dim>::evaluateYield(const dealii::Table<1, T >& beta,double tau_y,T q)
{
  //evaluate the yield function
  dealii::Table<1,T> Dev_beta(dim);
  Dev(beta,Dev_beta); // get deviatoric part of the beta vector (principal stresses)
  T norm_Dev_beta = norm(Dev_beta); // take the norm

  return norm_Dev_beta - std::sqrt(2./3.)*(tau_y - q);
}

template class Residual<Sacado::Fad::DFad<double>, 1>;
template class Residual<Sacado::Fad::DFad<double>, 2>;
template class Residual<Sacado::Fad::DFad<double>, 3>;
