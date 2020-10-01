/*
zhenlin wang 2019
*module transportation
*/
#include "../include/PoissonEquation.h"
#include "../include/Battery_fields.h"

template <int dim>
PoissonEquation<dim>::PoissonEquation(){}

template <int dim>
PoissonEquation<dim>::PoissonEquation(Battery_fields& _fields, Residual<Sacado::Fad::DFad<double>,dim>& _ResidualEq, int _primiary_dof)
{
	Solution=&_fields;
	ResidualEq=&_ResidualEq;
	primiary_dof=_primiary_dof;
}

template <int dim>
void PoissonEquation<dim>::set_up_fields(Battery_fields& _fields, Residual<Sacado::Fad::DFad<double>,dim>& _ResidualEq, int _primiary_dof)
{
	Solution=&_fields;
	ResidualEq=&_ResidualEq;
	primiary_dof=_primiary_dof;
}

template <int dim>
void PoissonEquation<dim>::r_get_residual(const typename hp::DoFHandler<dim>::active_cell_iterator &cell, const FEValues<dim>& fe_values, Table<1, Sacado::Fad::DFad<double> >& R, Table<1, Sacado::Fad::DFad<double>>& ULocal, Table<1, double >& ULocalConv)
{
//evaluate primary fields
	unsigned int n_q_points= fe_values.n_quadrature_points;
	
	dealii::Table<2,Sacado::Fad::DFad<double> > field(n_q_points,dim);
	dealii::Table<1,Sacado::Fad::DFad<double> > source(n_q_points);
	set_field_and_source_term(field,source);
	
	//call residual functions
	ResidualEq->residualForPoissonEq(fe_values, primiary_dof, R, field, source);
	
}

template <int dim>
void PoissonEquation<dim>::set_field_and_source_term(dealii::Table<2,Sacado::Fad::DFad<double> >& field, dealii::Table<1,Sacado::Fad::DFad<double> >& source)
{
	field=Solution->quad_fields[primiary_dof].value_grad;
	source=table_scaling<1,Sacado::Fad::DFad<double> > (source,0);
}

template class PoissonEquation<1>;
template class PoissonEquation<2>;
template class PoissonEquation<3>;