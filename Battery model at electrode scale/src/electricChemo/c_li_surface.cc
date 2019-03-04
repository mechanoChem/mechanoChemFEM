#include"electricChemo.h"

template <class T,int dim>
T ElectricChemo<T,dim>::c_li_surface_parabolic(T c_li, T jn, T D_s, T R_s_0)
{
	return c_li-R_s_0*jn/5.0/D_s;
}


template class ElectricChemo<Sacado::Fad::DFad<double>,1>;
template class ElectricChemo<Sacado::Fad::DFad<double>,2>;
template class ElectricChemo<Sacado::Fad::DFad<double>,3>;

template class ElectricChemo<double,1>;
template class ElectricChemo<double,2>;
template class ElectricChemo<double,3>;