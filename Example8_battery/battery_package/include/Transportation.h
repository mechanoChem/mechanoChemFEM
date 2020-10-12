/*
zhenlin wang 2020
*module transportation
*/
#ifndef Transportation_h
#define Transportation_h

#include "mechanoChemFEM.h"
#include <deal.II/base/table_indices.h>
#include "Battery_fields.h"
template <int dim>
class Transportation
{
	public:
		Transportation();
		Transportation(Battery_fields<dim>& _fields, Residual<Sacado::Fad::DFad<double>,dim>& ResidualEq);
		Transportation(Battery_fields<dim>& _fields, Residual<Sacado::Fad::DFad<double>,dim>& ResidualEq, int _primiary_dof);
		//this is a overloaded function 
		void set_up_fields(Battery_fields<dim>& _battery_fields, Residual<Sacado::Fad::DFad<double>,dim>& ResidualEq, int _primiary_dof);
		void set_primiary_dof(int _primiary_dof);
		virtual void set_reaction_term(dealii::Table<1,Sacado::Fad::DFad<double> >& react);
		virtual void set_diffusion_term(dealii::Table<2,Sacado::Fad::DFad<double> >& diffu);
		void r_get_residual(const FEValues<dim>& fe_values, Table<1, Sacado::Fad::DFad<double> >& R, Table<1, Sacado::Fad::DFad<double>>& ULocal, Table<1, double >& ULocalConv);
		int primiary_dof;	
		
		Battery_fields<dim>* Solution; 
		Residual<Sacado::Fad::DFad<double>,dim>* ResidualEq;
		
};
#endif