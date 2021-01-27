#include"../../include/Residual.h"
/*
*residualForDiffusionEq
*/
template <class T, int dim>
void Residual<T,dim>::residualForGeneralWeakForm(const FEValues<dim>& fe_values, unsigned int DOF, dealii::Table<1,T >& R, dealii::Table<1, T >& scalar_term, dealii::Table<2, T >& vector_term)
{
  unsigned int dofs_per_cell= fe_values.dofs_per_cell;
  unsigned int n_q_points= fe_values.n_quadrature_points;  
  //evaluate Residual
  for (unsigned int i=0; i<dofs_per_cell; ++i) {
    const unsigned int ck = fe_values.get_fe().system_to_component_index(i).first-DOF;
		if (ck==0){
    	for (unsigned int q=0; q<n_q_points; ++q){
				R[i] +=  fe_values.shape_value(i, q)* scalar_term[q]*fe_values.JxW(q);
				for (unsigned int j = 0; j < dim; j++){
	 			 	R[i] += -fe_values.shape_grad(i, q)[j]*vector_term[q][j]*fe_values.JxW(q);
				}
      }
    }
  }		
}

template <class T, int dim>
void Residual<T,dim>::residualForGeneralWeakForm(const FEValues<dim>& fe_values, unsigned int DOF, dealii::Table<1,T >& R, deformationMap<T, dim>& defMap, dealii::Table<1, T >& scalar_term, dealii::Table<2, T >& vector_term)
{
  unsigned int dofs_per_cell= fe_values.dofs_per_cell;
  unsigned int n_q_points= fe_values.n_quadrature_points;  
  //evaluate Residual
  for (unsigned int i=0; i<dofs_per_cell; ++i) {
    const unsigned int ck = fe_values.get_fe().system_to_component_index(i).first-DOF;
		if (ck==0){
    	for (unsigned int q=0; q<n_q_points; ++q){
				dealii::Table<1,T > shapeGradSpat(dim); 
				for (unsigned int k=0; k<dim; ++k){
					shapeGradSpat[k]=0;
				}
		  	for (unsigned int j=0; j<dim; ++j){
					for(unsigned int k=0; k<dim; ++k){
		      	shapeGradSpat[j] += fe_values.shape_grad(i, q)[k]*defMap.invF[q][k][j];
					}
				}
				
				R[i] +=  fe_values.shape_value(i, q)* scalar_term[q]*fe_values.JxW(q)*defMap.detF[q];
				for (unsigned int j = 0; j < dim; j++){
	 			 	R[i] += -shapeGradSpat[j]*vector_term[q][j]*fe_values.JxW(q)*defMap.detF[q];
				}
      }
    }
  }		
}

template <class T, int dim>
void Residual<T,dim>::residualForEqualityEq(const FEValues<dim>& fe_values, unsigned int DOF, dealii::Table<1, T >& R, dealii::Table<1,T >& value)
{
  unsigned int dofs_per_cell= fe_values.dofs_per_cell;
  unsigned int n_q_points= fe_values.n_quadrature_points;  
  //evaluate Residual
  for (unsigned int i=0; i<dofs_per_cell; ++i) {
    const unsigned int ck = fe_values.get_fe().system_to_component_index(i).first-DOF;
		if (ck==0){
    	for (unsigned int q=0; q<n_q_points; ++q){
				R[i] +=  fe_values.shape_value(i, q)* value[q]*fe_values.JxW(q);
      }
    }
  }		
}

template <class T, int dim>
void Residual<T,dim>::residualForEqualityEq(const FEValues<dim>& fe_values, unsigned int DOF, dealii::Table<1, T >& R, dealii::Table<2,T >& value)
{
  unsigned int dofs_per_cell= fe_values.dofs_per_cell;
  unsigned int n_q_points= fe_values.n_quadrature_points;  
  //evaluate Residual
  for (unsigned int i=0; i<dofs_per_cell; ++i) {
    const unsigned int ck = fe_values.get_fe().system_to_component_index(i).first-DOF;
		if (ck>=0 && ck<dim){
    	for (unsigned int q=0; q<n_q_points; ++q){
				R[i] +=  fe_values.shape_value(i, q)* value[q][ck]*fe_values.JxW(q);
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