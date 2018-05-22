#include"../../include/Residual.h"

template <class T, int dim>
void Residual<T,dim>::scalling(const FEValues<dim>& fe_values, unsigned int DOF, Table<1, Sacado::Fad::DFad<double> >& R, double scallingFactor)
{
	unsigned int dofs_per_cell= fe_values.dofs_per_cell;
  for (unsigned int i=0; i<dofs_per_cell; ++i) {
    const unsigned int ck = fe_values.get_fe().system_to_component_index(i).first - DOF ;
    if (ck==0){
			R[i] = R[i]*scallingFactor;
    }
  }
}


/**
*volume integration
*/
template <class T, int dim>
double Residual<T,dim>::volumeIntegration(const FEValues<dim>& fe_values, Table<1, Sacado::Fad::DFad<double> >& value_quad)
{
	double value=0;
	unsigned int n_q_points= fe_values.n_quadrature_points;
  for (unsigned int q=0; q<n_q_points; ++q){
		value += value_quad[q].val()*fe_values.JxW(q);							
  }
	return value;
}

template <class T, int dim>
double Residual<T,dim>::volumeIntegration(const FEValues<dim>& fe_values, Table<1, Sacado::Fad::DFad<double> >& value_quad, deformationMap<T, dim>& defMap)
{
	double value=0;
	unsigned int n_q_points= fe_values.n_quadrature_points;
  for (unsigned int q=0; q<n_q_points; ++q){
		value += value_quad[q].val()*fe_values.JxW(q)*defMap.detF[q].val();							
  }
	return value;
}

//============================================================================================
template <class T, int dim>
double Residual<T,dim>::volumeIntegration(const FEValues<dim>& fe_values, double value_quad)
{
	double value=0;
	unsigned int n_q_points= fe_values.n_quadrature_points;
  for (unsigned int q=0; q<n_q_points; ++q){
		value += value_quad*fe_values.JxW(q);							
  }
	return value;
}
template <class T, int dim>
double Residual<T,dim>::volumeIntegration(const FEValues<dim>& fe_values, double value_quad, deformationMap<T, dim>& defMap)
{
	double value=0;
	unsigned int n_q_points= fe_values.n_quadrature_points;
  for (unsigned int q=0; q<n_q_points; ++q){
		value += value_quad*fe_values.JxW(q)*defMap.detF[q].val();							
  }
	return value;
}
//============================================================================================
template <class T, int dim>
double Residual<T,dim>::volumeIntegration(const FEValues<dim>& fe_values, Table<1, double >& value_quad)
{
	double value=0;
	unsigned int n_q_points= fe_values.n_quadrature_points;
  for (unsigned int q=0; q<n_q_points; ++q){
		value += value_quad[q]*fe_values.JxW(q);							
  }
	return value;
}

template <class T, int dim>
double Residual<T,dim>::volumeIntegration(const FEValues<dim>& fe_values, Table<1, double >& value_quad, deformationMap<T, dim>& defMap)
{
	double value=0;
	unsigned int n_q_points= fe_values.n_quadrature_points;
  for (unsigned int q=0; q<n_q_points; ++q){
		value += value_quad[q]*fe_values.JxW(q)*defMap.detF[q].val();							
  }
	return value;
}
//============================================================================================
/**
*Surface integration
*/
template <class T, int dim>
double Residual<T,dim>::surfaceIntegration(const FEFaceValues<dim>& fe_face_values, double value_quad)
{
	double value=0;
	unsigned int n_face_q_points=fe_face_values.n_quadrature_points;
  for (unsigned int q=0; q<n_face_q_points; ++q){
		value += value_quad*fe_face_values.JxW(q);							
  }
	return value;
}


template <class T, int dim>
double Residual<T,dim>::surfaceIntegration(const FEFaceValues<dim>& fe_face_values, double value_quad, deformationMap<T, dim>& defMap_face)
{
	double value=0;
	unsigned int n_face_q_points=fe_face_values.n_quadrature_points;
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
		value += value_quad*fe_face_values.JxW(q)*l2_norm.val()*defMap_face.detF[q].val();							
  }
	return value;
}

//============================================================================================

template <class T, int dim>
double Residual<T,dim>::surfaceIntegration(const FEFaceValues<dim>& fe_face_values, dealii::Table<1, Sacado::Fad::DFad<double> >& value_quad)
{
	double value=0;
	unsigned int n_face_q_points=fe_face_values.n_quadrature_points;
  for (unsigned int q=0; q<n_face_q_points; ++q){
		value += value_quad[q].val()*fe_face_values.JxW(q);							
  }
	return value;
}

template <class T, int dim>
double Residual<T,dim>::surfaceIntegration(const FEFaceValues<dim>& fe_face_values, dealii::Table<1, Sacado::Fad::DFad<double> >& value_quad, deformationMap<T, dim>& defMap_face)
{
	double value=0;
	unsigned int n_face_q_points=fe_face_values.n_quadrature_points;
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
		value += value_quad[q].val()*fe_face_values.JxW(q)*l2_norm.val()*defMap_face.detF[q].val();							
  }
	return value;
}


//============================================================================================
template <class T, int dim>
double Residual<T,dim>::surfaceIntegration(const FEFaceValues<dim>& fe_face_values, dealii::Table<1, double >& value_quad)
{
	double value=0;
	unsigned int n_face_q_points=fe_face_values.n_quadrature_points;
  for (unsigned int q=0; q<n_face_q_points; ++q){
		value += value_quad[q]*fe_face_values.JxW(q);							
  }
	return value;
}

template <class T, int dim>
double Residual<T,dim>::surfaceIntegration(const FEFaceValues<dim>& fe_face_values, dealii::Table<1, double >& value_quad, deformationMap<T, dim>& defMap_face)
{
	double value=0;
	unsigned int n_face_q_points=fe_face_values.n_quadrature_points;
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
		value += value_quad[q]*fe_face_values.JxW(q)*l2_norm.val()*defMap_face.detF[q].val();							
  }
	return value;
}


template class Residual<Sacado::Fad::DFad<double>, 1>;
template class Residual<Sacado::Fad::DFad<double>, 2>;
template class Residual<Sacado::Fad::DFad<double>, 3>;