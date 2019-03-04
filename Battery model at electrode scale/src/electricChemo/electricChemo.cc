#include"electricChemo.h"

template <class T,int dim>
ElectricChemo<T,dim>::ElectricChemo(){}

template <class T,int dim>
ElectricChemo<T,dim>::ElectricChemo(dealii::ParameterHandler& _params):params(&_params)
{declare_parameters();}

template <class T,int dim>
ElectricChemo<T,dim>::~ElectricChemo(){}

template class ElectricChemo<Sacado::Fad::DFad<double>,1>;
template class ElectricChemo<Sacado::Fad::DFad<double>,2>;
template class ElectricChemo<Sacado::Fad::DFad<double>,3>;

template class ElectricChemo<double,1>;
template class ElectricChemo<double,2>;
template class ElectricChemo<double,3>;