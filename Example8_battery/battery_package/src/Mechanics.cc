/*
zhenlin wang 2019
*module transportation
*/
#include "../include/Mechanics.h"

template <int dim>
Mechanics<dim>::Mechanics(){}

template <int dim>
Mechanics<dim>::Mechanics(Battery_fields<dim>& _fields, Residual<Sacado::Fad::DFad<double>,dim>& _ResidualEq)
{
	Solution=&_fields;
	ResidualEq=&_ResidualEq;
}

template <int dim>
Mechanics<dim>::Mechanics(Battery_fields<dim>& _fields,Residual<Sacado::Fad::DFad<double>,dim>& _ResidualEq, int _primiary_dof)
{
	Solution=&_fields;
	ResidualEq=&_ResidualEq;
	primiary_dof=_primiary_dof;
}

template <int dim>
void Mechanics<dim>::set_up_fields(Battery_fields<dim>& _battery_fields,Residual<Sacado::Fad::DFad<double>,dim>& _ResidualEq, int _primiary_dof)
{
	Solution=&_battery_fields;
	ResidualEq=&_ResidualEq;
	primiary_dof=_primiary_dof;
}

template <int dim>
void Mechanics<dim>::set_primiary_dof(int _primiary_dof)
{
	primiary_dof=_primiary_dof;
}

template <int dim>
void Mechanics<dim>::r_get_residual(const FEValues<dim>& fe_values, Table<1, Sacado::Fad::DFad<double> >& R, Table<1, Sacado::Fad::DFad<double>>& ULocal, Table<1, double >& ULocalConv)
{
//evaluate primary fields	
	ResidualEq->setLameParametersByYoungsModulusPoissonRatio(youngsModulus, poissonRatio);
	unsigned int n_q_points= fe_values.n_quadrature_points;
	dealii::Table<3, Sacado::Fad::DFad<double> >Fe(n_q_points,dim,dim);
	deformationMap<Sacado::Fad::DFad<double>, dim> defMap(n_q_points); 
	getDeformationMap<Sacado::Fad::DFad<double>, dim>(fe_values, primiary_dof, ULocal, defMap);
	dealii::Table<3, Sacado::Fad::DFad<double> > P(n_q_points,dim,dim);
	set_stress(defMap.F, P);
	  
	//chemo
	ResidualEq->residualForMechanics(fe_values, primiary_dof, R, P);
	
}

template <int dim>
void Mechanics<dim>::set_stress(dealii::Table<3,Sacado::Fad::DFad<double> >& F, dealii::Table<3,Sacado::Fad::DFad<double> >& P)
{
	ResidualEq->evaluateNeoHookeanStress(P, F);
}


template class Mechanics<1>;
template class Mechanics<2>;
template class Mechanics<3>;