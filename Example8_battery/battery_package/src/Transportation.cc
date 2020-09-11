/*
zhenlin wang 2020
*module transportation
*/
#include "../include/Transportation.h"

template <int dim>
Transportation<dim>::Transportation(){}

template <int dim>
Transportation<dim>::Transportation(Battery_fields& _fields, Residual<Sacado::Fad::DFad<double>,dim>& _ResidualEq)
{
	Solution=&_fields;
	ResidualEq=&_ResidualEq;
}

template <int dim>
Transportation<dim>::Transportation(Battery_fields& _fields, Residual<Sacado::Fad::DFad<double>,dim>& _ResidualEq, int _primiary_dof)
{
	Solution=&_fields;
	ResidualEq=&_ResidualEq;
	primiary_dof=_primiary_dof;
}

template <int dim>
void Transportation<dim>::set_up_fields(Battery_fields& _battery_fields, Residual<Sacado::Fad::DFad<double>,dim>& _ResidualEq, int _primiary_dof)
{
	Solution=&_battery_fields;
	ResidualEq=&_ResidualEq;
	primiary_dof=_primiary_dof;
}

template <int dim>
void Transportation<dim>::set_primiary_dof(int _primiary_dof)
{
	primiary_dof=_primiary_dof;
}

template <int dim>
void Transportation<dim>::r_get_residual(const typename hp::DoFHandler<dim>::active_cell_iterator &cell, const FEValues<dim>& fe_values, Table<1, Sacado::Fad::DFad<double> >& R, Table<1, Sacado::Fad::DFad<double>>& ULocal, Table<1, double >& ULocalConv)
{
//evaluate primary fields
	unsigned int n_q_points= fe_values.n_quadrature_points;
	  
	//chemo
	TableIndices<1> tableIndex1(n_q_points);
	TableIndices<2> tableIndex2(n_q_points,n_q_points);
	Solution->quad_fields[primiary_dof].value_conv.reinit(tableIndex1);
	Solution->quad_fields[primiary_dof].value.reinit(tableIndex1);
	Solution->quad_fields[primiary_dof].value_grad.reinit(tableIndex2);
	
	dealii::Table<1,Sacado::Fad::DFad<double> > react(n_q_points);
	dealii::Table<2,Sacado::Fad::DFad<double> > diffu(n_q_points, dim);
	evaluateScalarFunction<double,dim>(fe_values, primiary_dof, ULocalConv, Solution->quad_fields[primiary_dof].value_conv);
	evaluateScalarFunction<Sacado::Fad::DFad<double>,dim>(fe_values, primiary_dof, ULocal, Solution->quad_fields[primiary_dof].value);
	evaluateScalarFunctionGradient<Sacado::Fad::DFad<double>,dim>(fe_values, primiary_dof, ULocal, Solution->quad_fields[primiary_dof].value_grad);
	set_reaction_term(react);
	set_diffusion_term(diffu);
	
	//call residual functions
	ResidualEq->residualForDiff_ReacEq(fe_values,primiary_dof, R,Solution->quad_fields[primiary_dof].value, Solution->quad_fields[primiary_dof].value_conv, diffu,react);
	
}

template <int dim>
void Transportation<dim>::set_reaction_term(dealii::Table<1,Sacado::Fad::DFad<double> >& react)
{
	unsigned int n_q_points= react.size(0);
	for (unsigned int q=0; q<n_q_points; ++q) react[q]=0;
}

template <int dim>
void Transportation<dim>::set_diffusion_term(dealii::Table<2,Sacado::Fad::DFad<double> >& diffu)
{
	unsigned int n_q_points= diffu.size(0);
	for (unsigned int q=0; q<n_q_points; ++q){
		for (unsigned int i=0; i<dim; ++i) diffu[q][i]=-Solution->quad_fields[primiary_dof].value_grad[q][i];
	}
}

template class Transportation<1>;
template class Transportation<2>;
template class Transportation<3>;