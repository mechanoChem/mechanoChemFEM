#include"../../include/Residual.h"

template <class T, int dim>
void Residual<T,dim>::evaluateNeoHookeanStress(dealii::Table<3, T >& P, dealii::Table<3, T > &Fe)
{
	//determine second Piola-Kirchhoff stress tensor, S, first Piola-Kirchhoff stress tensor P
	unsigned int n_q_points= P.size(0);
	for (unsigned int q=0; q<n_q_points; ++q){
		
		Table<2, T > C (dim, dim);
		Table<2, T > S(dim, dim); 
		for (unsigned int i=0; i<dim; ++i){
			 for (unsigned int j=0; j<dim; ++j){
				 C[i][j] = 0.0;
				 for (unsigned int k=0; k<dim; ++k){
					 C[i][j] += Fe[q][k][i]*Fe[q][k][j];
				 }
			 }
		}
		 
  	T detC = 0.0, logdetC = 0.0;
  	Table<2, T> invC(dim, dim); 
		
  	for (unsigned int i = 0; i < dim; ++i){
    	for (unsigned int j = 0; j < dim; ++j)
      {
				invC[i][j] = 0.0;
			}
		}
  	getInverse<T, dim>(C, invC, detC);
  	for (unsigned int i = 0; i < dim; ++i){
			for (unsigned int j = 0; j < dim; ++j){
				S[i][j] += 0.5*lambda*detC*invC[i][j] - (0.5*lambda + mu)*invC[i][j] + mu*(i==j);
    	}
 		}
		 
		//P
  	for (unsigned int i=0; i<dim; ++i){
			for (unsigned int j=0; j<dim; ++j){
				P[q][i][j]=0.0;
				for (unsigned int k=0; k<dim; ++k){
					P[q][i][j]+=Fe[q][i][k]*S[k][j];
				}
			}
		}
	}
}

template class Residual<Sacado::Fad::DFad<double>, 1>;
template class Residual<Sacado::Fad::DFad<double>, 2>;
template class Residual<Sacado::Fad::DFad<double>, 3>;
template class Residual<double, 1>;
template class Residual<double, 2>;
template class Residual<double, 3>;