#include"electricChemo.h"

template <class T,int dim>
T ElectricChemo<T,dim>::formula_dUdt(T UnitC)
{
	T dUdt;
	if(UnitC<=0.2)  dUdt=0.01442*UnitC*UnitC-0.00291*UnitC-0.000138;
	else if(UnitC<=0.4) dUdt=0.00634*UnitC*UnitC*UnitC-0.006625*UnitC*UnitC+0.002635*UnitC-0.0004554;
	else if(UnitC<=0.5) dUdt=0.001059*UnitC-0.0004793;
	else if(UnitC<=0.7) dUdt=0.00025*UnitC-7.5e-5;
	else if(UnitC<=0.8) dUdt=-0.001*UnitC+0.0008;
	else if(UnitC<=0.85) dUdt=0.0333*UnitC*UnitC-0.057*UnitC+0.02427;
	else if(UnitC<=0.95) dUdt=0.002*UnitC*UnitC-0.0039*UnitC+0.00177;
  else if (UnitC<=1) dUdt=-0.0014*UnitC+0.0012;
	return dUdt;
}


template class ElectricChemo<Sacado::Fad::DFad<double>,1>;
template class ElectricChemo<Sacado::Fad::DFad<double>,2>;
template class ElectricChemo<Sacado::Fad::DFad<double>,3>;

template class ElectricChemo<double,1>;
template class ElectricChemo<double,2>;
template class ElectricChemo<double,3>;
