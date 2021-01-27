#include"../../include/Residual.h"


template <class T, int dim>
void Residual<T,dim>::residualForStokesEq(const FEValues<dim>& fe_values, unsigned int DOF, Table<1, T >& R,  Table<3, T>& gradV, dealii::Table<1, T >& pressure)
{
  unsigned int dofs_per_cell= fe_values.dofs_per_cell;
  unsigned int n_q_points= fe_values.n_quadrature_points;
	
	//dealii::Table<3, T > P_stoke(n_q_points,dim,dim);
	//evaluateStokesStress(P_stoke,fe_values,gradV,pressure);

  for (unsigned int i=0; i<dofs_per_cell; ++i) {
    const unsigned int ck = fe_values.get_fe().system_to_component_index(i).first - DOF;
    if (ck>=0 && ck<dim){
      for (unsigned int q=0; q<n_q_points; ++q){
				for (unsigned int j = 0; j < dim; j++){
	  			R[i] += fe_values.shape_grad(i, q)[j]*(viscosity*gradV[q][ck][j])*fe_values.JxW(q);				
				}
					R[i] +=fe_values.shape_grad(i, q)[ck]*-pressure[q]*fe_values.JxW(q);
      }
  	}	
	}
}


template <class T, int dim>
void Residual<T,dim>::residualForStokesEq(const FEValues<dim>& fe_values, unsigned int DOF, Table<1, T >& R, deformationMap<T, dim>& defMap, Table<3, T>& gradV, dealii::Table<1, T >& pressure)
	{
	  unsigned int dofs_per_cell= fe_values.dofs_per_cell;
	  unsigned int n_q_points= fe_values.n_quadrature_points;
	
		//dealii::Table<3, T > P_stoke(n_q_points,dim,dim);
		//evaluateStokesStress(P_stoke,fe_values,gradV,pressure);

	  for (unsigned int i=0; i<dofs_per_cell; ++i) {
	    const unsigned int ck = fe_values.get_fe().system_to_component_index(i).first - DOF;
	    if (ck>=0 && ck<dim){
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
					
					for (unsigned int j = 0; j < dim; j++){
		  			R[i] += fe_values.shape_grad(i, q)[j]*(viscosity*gradV[q][ck][j])*fe_values.JxW(q)*defMap.detF[q];				
					}
						R[i] +=fe_values.shape_grad(i, q)[ck]*-pressure[q]*fe_values.JxW(q)*defMap.detF[q];
	      }
	  	}	
		}
		
	}



template <class T, int dim>
void Residual<T,dim>:: residualForNavier_StokesEq(const FEValues<dim>& fe_values, unsigned int DOF, Table<1, T >& R, deformationMap<T, dim>& defMap, Table<2, T>& V, dealii::Table<2,double>& V_conv,Table<3, T>& gradV, dealii::Table<1, T >& pressure)
{
  unsigned int dofs_per_cell= fe_values.dofs_per_cell;
  unsigned int n_q_points= fe_values.n_quadrature_points;

	//dealii::Table<3, T > P_stoke(n_q_points,dim,dim);
	//evaluateStokesStress(P_stoke,fe_values,gradV,pressure);

  for (unsigned int i=0; i<dofs_per_cell; ++i) {
    const unsigned int ck = fe_values.get_fe().system_to_component_index(i).first - DOF;
    if (ck>=0 && ck<dim){
      for (unsigned int q=0; q<n_q_points; ++q){
				R[i] +=fe_values.shape_value(i, q)*density*(V[q][ck]-V_conv[q][ck])/dt*fe_values.JxW(q)*defMap.detF[q];
				for (unsigned int j = 0; j < dim; j++){
					R[i] +=fe_values.shape_value(i, q)*density*V[q][j]*gradV[q][ck][j]*fe_values.JxW(q)*defMap.detF[q];
					
	  			R[i] += fe_values.shape_grad(i, q)[j]*(viscosity*gradV[q][ck][j])*fe_values.JxW(q)*defMap.detF[q];				
				}
					R[i] +=fe_values.shape_grad(i, q)[ck]*-pressure[q]*fe_values.JxW(q)*defMap.detF[q];
      }
  	}	
	}
	
}




template <class T, int dim>
void Residual<T,dim>::evaluateStokesStress(dealii::Table<3, T > P_stoke, const FEValues<dim>& fe_values, Table<3, T>& gradV, dealii::Table<1, T >& pressure)
{
	unsigned int n_q_points= fe_values.n_quadrature_points;
	//Stokes's stress constitutive equation
	for (unsigned int q=0; q<n_q_points; ++q){
		for (unsigned int i=0; i<dim; ++i){
	  	for (unsigned int j=0; j<dim; ++j){
				P_stoke[q][i][j]=viscosity*(gradV[q][i][j]+gradV[q][j][i]);
			}
			//P_stoke[q][i][i]+=-pressure[q];
		}
	}
}

template <class T, int dim>
void Residual<T,dim>::residualForContinuityEq(const FEValues<dim>& fe_values, unsigned int DOF, Table<1, T >& R, Table<1, T >&divVelocity)
{
  unsigned int dofs_per_cell= fe_values.dofs_per_cell;
  unsigned int n_q_points= fe_values.n_quadrature_points;
	//continuity equation without time term
  for (unsigned int i=0; i<dofs_per_cell; ++i) {
    const unsigned int ck = fe_values.get_fe().system_to_component_index(i).first - DOF;
    for (unsigned int q=0; q<n_q_points; ++q){
	  	R[i] += -fe_values.shape_value(i, q)*density*divVelocity[q]*fe_values.JxW(q);
    }
	}
}

template <class T, int dim>
void Residual<T,dim>::residualForContinuityEq(const FEValues<dim>& fe_values, unsigned int DOF, Table<1, T >& R, deformationMap<T, dim>& defMap,  Table<1, T >&divVelocity)
{
  unsigned int dofs_per_cell= fe_values.dofs_per_cell;
  unsigned int n_q_points= fe_values.n_quadrature_points;
	//continuity equation without time term
  for (unsigned int i=0; i<dofs_per_cell; ++i) {
    const unsigned int ck = fe_values.get_fe().system_to_component_index(i).first - DOF;
    for (unsigned int q=0; q<n_q_points; ++q){
	  	R[i] += -fe_values.shape_value(i, q)*density*divVelocity[q]*fe_values.JxW(q)*defMap.detF[q];
    }
	}
}


template class Residual<Sacado::Fad::DFad<double>, 1>;
template class Residual<Sacado::Fad::DFad<double>, 2>;
template class Residual<Sacado::Fad::DFad<double>, 3>;
template class Residual<double, 1>;
template class Residual<double, 2>;
template class Residual<double, 3>;