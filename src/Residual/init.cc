#include"../../include/Residual.h"


template <class T, int dim>
Residual<T,dim>::Residual(){}

template <class T, int dim>
Residual<T, dim>::~Residual (){}


template <class T, int dim>
void Residual<T, dim>::reinit(dealii::ParameterHandler& _params)
{
	params=&_params;
}

template <class T, int dim>
void Residual<T, dim>::setLameParametersByYoungsModulusPoissonRatio(double youngsModulus, double poissonRatio)
{
	lambda=(youngsModulus*poissonRatio)/((1+poissonRatio)*(1-2*poissonRatio));
	mu=youngsModulus/(2*(1+poissonRatio));
}

template class Residual<Sacado::Fad::DFad<double>, 1>;
template class Residual<Sacado::Fad::DFad<double>, 2>;
template class Residual<Sacado::Fad::DFad<double>, 3>;