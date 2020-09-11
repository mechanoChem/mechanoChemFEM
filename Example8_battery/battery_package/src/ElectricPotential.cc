/*
zhenlin wang 2019
*module transportation
*/
#include "../include/ElectricPotential.h"
#include "../include/Battery_fields.h"

template <int dim>
ElectricPotential<dim>::ElectricPotential(){}

template <int dim>
ElectricPotential<dim>::ElectricPotential(Battery_fields& _fields, Residual<Sacado::Fad::DFad<double>,dim>& _ResidualEq, int _primiary_dof)
{
	Solution=&_fields;
	ResidualEq=&_ResidualEq;
	primiary_dof=_primiary_dof;
}

template <int dim>
void ElectricPotential<dim>::set_up_fields(Battery_fields& _fields, Residual<Sacado::Fad::DFad<double>,dim>& _ResidualEq, int _primiary_dof)
{
	Solution=&_fields;
	ResidualEq=&_ResidualEq;
	primiary_dof=_primiary_dof;
}

template <int dim>
void ElectricPotential<dim>::r_get_residual(const typename hp::DoFHandler<dim>::active_cell_iterator &cell, const FEValues<dim>& fe_values, Table<1, Sacado::Fad::DFad<double> >& R, Table<1, Sacado::Fad::DFad<double>>& ULocal, Table<1, double >& ULocalConv)
{
//evaluate primary fields
	unsigned int n_q_points= fe_values.n_quadrature_points;
	  
	//chemo
	TableIndices<1> tableIndex1(n_q_points);
	TableIndices<2> tableIndex2(n_q_points,n_q_points);
	Solution->quad_fields[primiary_dof].value_conv.reinit(tableIndex1);
	Solution->quad_fields[primiary_dof].value.reinit(tableIndex1);
	Solution->quad_fields[primiary_dof].value_grad.reinit(tableIndex2);
	
	dealii::Table<2,Sacado::Fad::DFad<double> > current(n_q_points,dim);
	dealii::Table<1,Sacado::Fad::DFad<double> > zeros(n_q_points);
	evaluateScalarFunction<double,dim>(fe_values, primiary_dof, ULocalConv, Solution->quad_fields[primiary_dof].value_conv);
	evaluateScalarFunction<Sacado::Fad::DFad<double>,dim>(fe_values, primiary_dof, ULocal, Solution->quad_fields[primiary_dof].value);
	evaluateScalarFunctionGradient<Sacado::Fad::DFad<double>,dim>(fe_values, primiary_dof, ULocal, Solution->quad_fields[primiary_dof].value_grad);
	set_current(current);
	
	//call residual functions
	ResidualEq->residualForPoissonEq(fe_values, primiary_dof, R, current, zeros);
	
}

template <int dim>
void ElectricPotential<dim>::set_current(dealii::Table<2,Sacado::Fad::DFad<double> >& current)
{
	unsigned int n_q_points= current.size(0);
	for (unsigned int q=0; q<n_q_points; ++q){
		for (unsigned int i=0; i<dim; ++i) current[q][i]=-Solution->quad_fields[primiary_dof].value_grad[q][i];
	}
}

template class ElectricPotential<1>;
template class ElectricPotential<2>;
template class ElectricPotential<3>;