#include"../../include/Residual.h"
/*
*residualForDiffusionEq
*/
template <class T, int dim>
void Residual<T,dim>::residualForDiffusionEq(const FEValues<dim>& fe_values, unsigned int DOF, dealii::Table<1, T >& R, dealii::Table<1, T >& c, dealii::Table<1,double>& c_conv, dealii::Table<2,T >& flux)
{
  unsigned int dofs_per_cell= fe_values.dofs_per_cell;
  unsigned int n_q_points= fe_values.n_quadrature_points;
  
  //evaluate Residual
  for (unsigned int i=0; i<dofs_per_cell; ++i) {
    const unsigned int ck = fe_values.get_fe().system_to_component_index(i).first-DOF;
    for (unsigned int q=0; q<n_q_points; ++q){
      if (ck==0){
				R[i] +=  fe_values.shape_value(i, q)*((c[q]-c_conv[q])/dt)*fe_values.JxW(q);
				for (unsigned int j = 0; j < dim; j++){
	 			 	R[i] += -fe_values.shape_grad(i, q)[j]*flux[q][j]*fe_values.JxW(q);
				}
      }
    }
  }	
}

template <class T, int dim>
void Residual<T,dim>::residualForDiffusionEq(const FEValues<dim>& fe_values, unsigned int DOF, dealii::Table<1, T >& R,deformationMap<T, dim>& defMap, dealii::Table<1, T >& c, dealii::Table<1,double>& c_conv, dealii::Table<2,T >& flux)
{
  unsigned int dofs_per_cell= fe_values.dofs_per_cell;
  unsigned int n_q_points= fe_values.n_quadrature_points;
  
  //evaluate Residual
  for (unsigned int i=0; i<dofs_per_cell; ++i) {
    const unsigned int ck = fe_values.get_fe().system_to_component_index(i).first-DOF;
    for (unsigned int q=0; q<n_q_points; ++q){
      if (ck==0){
				dealii::Table<1,T > shapeGradSpat(dim); 
				for (unsigned int k=0; k<dim; ++k){
					shapeGradSpat[k]=0;
				}
		  	for (unsigned int j=0; j<dim; ++j){
					for(unsigned int k=0; k<dim; ++k){
		      	shapeGradSpat[j] += fe_values.shape_grad(i, q)[k]*defMap.invF[q][k][j];
					}
				}
				
				R[i] +=  fe_values.shape_value(i, q)*((c[q]-c_conv[q])/dt)*fe_values.JxW(q)*defMap.detF[q];
				for (unsigned int j = 0; j < dim; j++){
	 			 	R[i] += -shapeGradSpat[j]*flux[q][j]*fe_values.JxW(q)*defMap.detF[q];
				}
      }
    }
  }	
}

/*
*residualForDiff_ReacEq
*/
template <class T, int dim>
void Residual<T,dim>::residualForDiff_ReacEq(const FEValues<dim>& fe_values, unsigned int DOF, dealii::Table<1, T >& R, dealii::Table<1, T >& c, dealii::Table<1,double>& c_conv, dealii::Table<2,T >& flux, dealii::Table<1, T >& reaction)
{
  unsigned int dofs_per_cell= fe_values.dofs_per_cell;
  unsigned int n_q_points= fe_values.n_quadrature_points;
  
  //evaluate Residual
  for (unsigned int i=0; i<dofs_per_cell; ++i) {
    const unsigned int ck = fe_values.get_fe().system_to_component_index(i).first-DOF;
    for (unsigned int q=0; q<n_q_points; ++q){
      if (ck==0){
				R[i] +=  fe_values.shape_value(i, q)*((c[q]-c_conv[q])/dt)*fe_values.JxW(q);
				for (unsigned int j = 0; j < dim; j++){
	 			 	R[i] += -fe_values.shape_grad(i, q)[j]*flux[q][j]*fe_values.JxW(q);
				}
				R[i] +=  -fe_values.shape_value(i, q)*reaction[q]*fe_values.JxW(q);
      }
    }
  }	
}
template <class T, int dim>
void Residual<T,dim>::residualForDiff_ReacEq(const FEValues<dim>& fe_values, unsigned int DOF, dealii::Table<1, T >& R,deformationMap<T, dim>& defMap, dealii::Table<1, T >& c, dealii::Table<1,double>& c_conv, dealii::Table<2,T >& flux, dealii::Table<1, T >& reaction)
{
  unsigned int dofs_per_cell= fe_values.dofs_per_cell;
  unsigned int n_q_points= fe_values.n_quadrature_points;
  
  //evaluate Residual
  for (unsigned int i=0; i<dofs_per_cell; ++i) {
    const unsigned int ck = fe_values.get_fe().system_to_component_index(i).first-DOF;
    for (unsigned int q=0; q<n_q_points; ++q){
      if (ck==0){
				
				dealii::Table<1,T > shapeGradSpat(dim); 
				for (unsigned int k=0; k<dim; ++k){
					shapeGradSpat[k]=0;
				}
		  	for (unsigned int j=0; j<dim; ++j){
					for(unsigned int k=0; k<dim; ++k){
		      	shapeGradSpat[j] += fe_values.shape_grad(i, q)[k]*defMap.invF[q][k][j];
					}
				}
				
				R[i] +=  fe_values.shape_value(i, q)*((c[q]-c_conv[q])/dt)*fe_values.JxW(q)*defMap.detF[q];
				for (unsigned int j = 0; j < dim; j++){
	 			 	R[i] += -shapeGradSpat[j]*flux[q][j]*fe_values.JxW(q)*defMap.detF[q];
				}
				R[i] +=  -fe_values.shape_value(i, q)*reaction[q]*fe_values.JxW(q)*defMap.detF[q];
      }
    }
  }	
}

/*
*residualForPoissonEq
*/
template <class T, int dim>
void Residual<T,dim>::residualForPoissonEq(const FEValues<dim>& fe_values, unsigned int DOF, dealii::Table<1,T >& R, dealii::Table<2, T >& phi_grad, dealii::Table<1, T >& rhs)
{
  unsigned int dofs_per_cell= fe_values.dofs_per_cell;
  unsigned int n_q_points= fe_values.n_quadrature_points;
  
  //evaluate Residual
  for (unsigned int i=0; i<dofs_per_cell; ++i) {
    const unsigned int ck = fe_values.get_fe().system_to_component_index(i).first-DOF;
    for (unsigned int q=0; q<n_q_points; ++q){
      if (ck==0){
				R[i] += -fe_values.shape_value(i, q)*rhs[q]*fe_values.JxW(q);							
				for (unsigned int j = 0; j < dim; j++){
	  			R[i] += -fe_values.shape_grad(i, q)[j]*phi_grad[q][j]*fe_values.JxW(q);
				}
      }
    }
  }	
}
template <class T, int dim>
void Residual<T,dim>::residualForPoissonEq(const FEValues<dim>& fe_values, unsigned int DOF, dealii::Table<1,T >& R,  deformationMap<T, dim>& defMap, dealii::Table<2, T >& phi_grad, dealii::Table<1, T >& rhs)
{
  unsigned int dofs_per_cell= fe_values.dofs_per_cell;
  unsigned int n_q_points= fe_values.n_quadrature_points;
  
  //evaluate Residual
  for (unsigned int i=0; i<dofs_per_cell; ++i) {
    const unsigned int ck = fe_values.get_fe().system_to_component_index(i).first-DOF;
    for (unsigned int q=0; q<n_q_points; ++q){
      if (ck==0){
				dealii::Table<1,T > shapeGradSpat(dim); 
				for (unsigned int k=0; k<dim; ++k){
					shapeGradSpat[k]=0;
				}
		  	for (unsigned int j=0; j<dim; ++j){
					for(unsigned int k=0; k<dim; ++k){
		      	shapeGradSpat[j] += fe_values.shape_grad(i, q)[k]*defMap.invF[q][k][j];
					}
				}
				
				R[i] += -fe_values.shape_value(i, q)*rhs[q]*fe_values.JxW(q)*defMap.detF[q];							
				for (unsigned int j = 0; j < dim; j++){
	  			R[i] += -shapeGradSpat[j]*phi_grad[q][j]*fe_values.JxW(q)*defMap.detF[q];
				}
      }
    }
  }	
}

/*
*residualForPoissonEq for vector variabls
*/
template <class T, int dim>
void Residual<T,dim>::residualForPoissonEq(const FEValues<dim>& fe_values, unsigned int DOF, dealii::Table<1,T >& R, dealii::Table<3, T >& phi_grad, dealii::Table<2, T >& rhs)
{
  unsigned int dofs_per_cell= fe_values.dofs_per_cell;
  unsigned int n_q_points= fe_values.n_quadrature_points;
  
  //evaluate Residual
  for (unsigned int i=0; i<dofs_per_cell; ++i) {
    const unsigned int ck = fe_values.get_fe().system_to_component_index(i).first-DOF;
		if (ck>=0 && ck<dim){
    	for (unsigned int q=0; q<n_q_points; ++q){
				R[i] += -fe_values.shape_value(i, q)*rhs[q][ck]*fe_values.JxW(q);							
				for (unsigned int j = 0; j < dim; j++){
	  			R[i] += -fe_values.shape_grad(i, q)[j]*phi_grad[q][ck][j]*fe_values.JxW(q);
				}
      }
    }
  }	
}
template <class T, int dim>
void Residual<T,dim>::residualForPoissonEq(const FEValues<dim>& fe_values, unsigned int DOF, dealii::Table<1,T >& R,  deformationMap<T, dim>& defMap, dealii::Table<3, T >& phi_grad, dealii::Table<2, T >& rhs)
{
  unsigned int dofs_per_cell= fe_values.dofs_per_cell;
  unsigned int n_q_points= fe_values.n_quadrature_points;
  
  //evaluate Residual
  for (unsigned int i=0; i<dofs_per_cell; ++i) {
    const unsigned int ck = fe_values.get_fe().system_to_component_index(i).first-DOF;
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
				
				R[i] += -fe_values.shape_value(i, q)*rhs[q][ck]*fe_values.JxW(q)*defMap.detF[q];							
				for (unsigned int j = 0; j < dim; j++){
	  			R[i] += -shapeGradSpat[j]*phi_grad[q][ck][j]*fe_values.JxW(q)*defMap.detF[q];
				}
      }
    }
  }	
}






template <class T, int dim>
void Residual<T,dim>::residualForNeummanBC(const FEValues<dim>& fe_values, const FEFaceValues<dim>& fe_face_values, unsigned int DOF, dealii::Table<1, T >& R, dealii::Table<1,T >& gradn)
{
  unsigned int dofs_per_cell= fe_values.dofs_per_cell;
	unsigned int n_face_q_points=fe_face_values.n_quadrature_points;
	unsigned int ck;
  for (unsigned int i=0; i<dofs_per_cell; ++i) {
    ck = fe_values.get_fe().system_to_component_index(i).first-DOF;
		if(ck==0){
			for (unsigned int q=0; q<n_face_q_points; ++q){
				R[i] += fe_face_values.shape_value(i, q)*gradn[q]*fe_face_values.JxW(q);
			}
		}
	}	
}



template <class T, int dim>
void Residual<T,dim>::residualForNeummanBC(const FEValues<dim>& fe_values, const FEFaceValues<dim>& fe_face_values, unsigned int DOF, dealii::Table<1, T >& R, deformationMap<T, dim>& defMap_face, dealii::Table<1,T >& gradn)
{
  unsigned int dofs_per_cell= fe_values.dofs_per_cell;
	unsigned int n_face_q_points=fe_face_values.n_quadrature_points;
	unsigned int ck;

	
  for (unsigned int i=0; i<dofs_per_cell; ++i) {
    ck = fe_values.get_fe().system_to_component_index(i).first-DOF;
		if(ck==0){
			for (unsigned int q=0; q<n_face_q_points; ++q){
				dealii::Table<1, T > surfaceNormal(dim);
				T l2_norm=0;
				const Tensor<1,dim> normal=fe_face_values.normal_vector(q);
				for(unsigned int i=0;i<dim;i++){
					surfaceNormal[i]=0;
					for(unsigned int j=0;j<dim;j++){
						//(F^-T)*normal;
						surfaceNormal[i]=surfaceNormal[i]+defMap_face.invF[q][j][i]*normal[j];
					}
					l2_norm=l2_norm+surfaceNormal[i]*surfaceNormal[i];
				}
				l2_norm=std::sqrt(l2_norm);
				
				R[i] += fe_face_values.shape_value(i, q)*gradn[q]*fe_face_values.JxW(q)*l2_norm*defMap_face.detF[q];
			}
		}
	}	
}




template <class T, int dim>
void Residual<T,dim>::residualForNeummanBC(const FEValues<dim>& fe_values, const FEFaceValues<dim>& fe_face_values, unsigned int DOF, dealii::Table<1, T >& R, double gradn)
{
  unsigned int dofs_per_cell= fe_values.dofs_per_cell;
	unsigned int n_face_q_points=fe_face_values.n_quadrature_points;
	unsigned int ck;
  for (unsigned int i=0; i<dofs_per_cell; ++i) {
    ck = fe_values.get_fe().system_to_component_index(i).first-DOF;
		if(ck==0){
			for (unsigned int q=0; q<n_face_q_points; ++q){
				R[i] += fe_face_values.shape_value(i, q)*gradn*fe_face_values.JxW(q);
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