/*
zhenlin wang 2020
*module transportation
*/
#ifndef PoissonEquation_h
#define PoissonEquation_h

#include "mechanoChemFEM.h"
#include <deal.II/base/table_indices.h>
#include "Battery_fields.h"
#include "ElectricChemo.h"

template <int dim>
class PoissonEquation
{
	public:
		PoissonEquation();
		PoissonEquation(Battery_fields<dim>& _fields, Residual<Sacado::Fad::DFad<double>,dim>& _ResidualEq);
		PoissonEquation(Battery_fields<dim>& _fields, Residual<Sacado::Fad::DFad<double>,dim>& _ResidualEq, int _primiary_dof);
		void set_up_fields(Battery_fields<dim>& _battery_fields, Residual<Sacado::Fad::DFad<double>,dim>& ResidualEq, int _primiary_dof);
		void set_up_fields(Battery_fields<dim>& _battery_fields, ElectricChemo<dim,Sacado::Fad::DFad<double>>& _electricChemoFormula, Residual<Sacado::Fad::DFad<double>,dim>& ResidualEq, int _primiary_dof);
		void set_primiary_dof(int _primiary_dof);
		virtual void set_field_and_source_term(dealii::Table<2,Sacado::Fad::DFad<double> >& field, dealii::Table<1,Sacado::Fad::DFad<double> >& source);
		void r_get_residual(const FEValues<dim>& fe_values, Table<1, Sacado::Fad::DFad<double> >& R, Table<1, Sacado::Fad::DFad<double>>& ULocal, Table<1, double >& ULocalConv);
		int primiary_dof;	
		
		Battery_fields<dim>* battery_fields; 
		ElectricChemo<dim, Sacado::Fad::DFad<double>>* electricChemoFormula;
		Residual<Sacado::Fad::DFad<double>,dim>* ResidualEq;
};
#endif