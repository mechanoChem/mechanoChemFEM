#include"../../include/Residual.h"

//default evaluateStress No thermal chemical coupling
template <class T, int dim>
void Residual<T, dim>::setLameParametersByYoungsModulusPoissonRatio(double youngsModulus, double poissonRatio)
{
	lambda=(youngsModulus*poissonRatio)/((1+poissonRatio)*(1-2*poissonRatio));
	mu=youngsModulus/(2*(1+poissonRatio));
}

template <class T, int dim>
void Residual<T,dim>::evaluateStrain(dealii::Table<3, T >& Fe, dealii::Table<3, T >& E, deformationMap<T, dim>& defMap, bool infinitesimal_strain_indicator)
{
	unsigned int n_q_points= Fe.size(0);

	for (unsigned int q=0; q<n_q_points; ++q){   		
		for (unsigned int i=0; i<dim; ++i){
			for (unsigned int j=0; j<dim; ++j){
			  if(infinitesimal_strain_indicator) {
			    E[q][i][j]=(defMap.F[q][i][j] - (i==j)+defMap.F[q][j][i] - (i==j))/2;
			  }
			  else{
			    E[q][i][j] = -0.5*(i==j);
					for (unsigned int k=0; k<dim; ++k){
				  	E[q][i][j] += 0.5*Fe[q][k][i]*Fe[q][k][j];
					}
			  }
			}	
		}
	}
}

template <class T, int dim>
void Residual<T,dim>::residualForMechanics(const FEValues<dim>& fe_values, unsigned int DOF, Table<1, T >& R,dealii::Table<3, T > P)
{
  unsigned int dofs_per_cell= fe_values.dofs_per_cell;
  unsigned int n_q_points= fe_values.n_quadrature_points;

  //evaluate ResidualULocal
  for (unsigned int i=0; i<dofs_per_cell; ++i) {
    const unsigned int ck = fe_values.get_fe().system_to_component_index(i).first - DOF;
    if (ck>=0 && ck<dim){
      for (unsigned int q=0; q<n_q_points; ++q){
				for (unsigned int d = 0; d < dim; d++){
	  			R[i] +=  fe_values.shape_grad(i, q)[d]*P[q][ck][d]*fe_values.JxW(q);
				}
      }
  	}
  }
}

template <class T, int dim>
void Residual<T,dim>::residualForNeummanBC(const FEValues<dim>& fe_values, const FEFaceValues<dim>& fe_face_values, unsigned int DOF, dealii::Table<1, T >& R, dealii::Table<2,T >& gradn){
  unsigned int dofs_per_cell= fe_values.dofs_per_cell;
	unsigned int n_face_q_points=fe_face_values.n_quadrature_points;
	unsigned int ck;
  for (unsigned int i=0; i<dofs_per_cell; ++i) {
    ck = fe_values.get_fe().system_to_component_index(i).first-DOF;
		if(ck>=0 && ck<dim){
			for (unsigned int q=0; q<n_face_q_points; ++q){
				R[i] += fe_face_values.shape_value(i, q)*gradn[q][ck]*fe_face_values.JxW(q);
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
