#include"../../include/Residual.h"

template <class T, int dim>
double Residual<T,dim>::SVK3D(unsigned int i, unsigned int j, unsigned int k, unsigned int l)
{
  return lambda*(i==j)*(k==l) + mu*((i==k)*(j==l) + (i==l)*(j==k));
}

template <class T, int dim>
double Residual<T,dim>::SVK2D(unsigned int i, unsigned int j)
{
  if (i==j && i<2) return lambda + 2*mu;
  else if (i==2 && j==2) return mu;
  else if ((i+j)==1) return lambda;
  else return 0.0;
}

template <class T, int dim>
void Residual<T,dim>::evaluateSaint_Venant_KirchhoffStress(dealii::Table<3, T >& P, dealii::Table<3, T > &Fe, dealii::Table<3, T >& E)
{
	unsigned int n_q_points= P.size(0);
	for (unsigned int q=0; q<n_q_points; ++q){
		Table<2, T > S (dim, dim);
		if(dim==3){
			for (unsigned int i=0; i<dim; ++i){
	  		for (unsigned int j=0; j<dim; ++j){
	    		S[i][j]=0;
	    		for (unsigned int k=0; k<dim; ++k){
	      		for (unsigned int l=0; l<dim; ++l){
	        		S[i][j] += SVK3D(i, j, k, l)*E[q][k][l];
	      		}
	    		}
	  		}
			}
		}
		else if(dim==2){
			S[0][0]=SVK2D(0,0)*E[q][0][0]+SVK2D(0,1)*E[q][1][1]+SVK2D(0,2)*(E[q][0][1]+E[q][1][0]);
  		S[1][1]=SVK2D(1,0)*E[q][0][0]+SVK2D(1,1)*E[q][1][1]+SVK2D(1,2)*(E[q][0][1]+E[q][1][0]);
  		S[0][1]=S[1][0]=SVK2D(2,0)*E[q][0][0]+SVK2D(2,1)*E[q][1][1]+SVK2D(2,2)*(E[q][0][1]+E[q][1][0]);
		}
		else if(dim==1){
 	 	 S[0][0]=lambda*(1+mu)*(1-2*mu)/mu*E[q][0][0];
		}
	
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