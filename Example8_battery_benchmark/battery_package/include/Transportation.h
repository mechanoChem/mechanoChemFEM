/*
zhenlin wang 2020
*module transportation
*/
#ifndef Transportation_h
#define Transportation_h

#include "mechanoChemFEM.h"
#include <deal.II/base/table_indices.h>
#include "Battery_fields.h"
#include "ElectricChemo.h"
#include "SDdata.h"
template <int dim>
class Transportation
{
	public:
		Transportation();
		Transportation(Battery_fields<dim>& _fields, Residual<Sacado::Fad::DFad<double>,dim>& ResidualEq);
		Transportation(Battery_fields<dim>& _fields, Residual<Sacado::Fad::DFad<double>,dim>& ResidualEq, int _primiary_dof);
		//this is a overloaded function 
		void set_up_fields(Battery_fields<dim>& _battery_fields, Residual<Sacado::Fad::DFad<double>,dim>& ResidualEq, int _primiary_dof);
		void set_up_fields(Battery_fields<dim>& _battery_fields, ElectricChemo<dim,Sacado::Fad::DFad<double>>& _electricChemoFormula, Residual<Sacado::Fad::DFad<double>,dim>& ResidualEq, int _primiary_dof);
		void set_primiary_dof(int _primiary_dof);
		virtual void set_diffusion_reaction_term(dealii::Table<2,Sacado::Fad::DFad<double> >& diffu, dealii::Table<1,Sacado::Fad::DFad<double> >& react);
		void r_get_residual(const FEValues<dim>& fe_values, Table<1, Sacado::Fad::DFad<double> >& R, Table<1, Sacado::Fad::DFad<double>>& ULocal, Table<1, double >& ULocalConv);
		void r_get_residual_with_interface(const typename hp::DoFHandler<dim>::active_cell_iterator &cell, const FEValues<dim>& fe_values, Table<1, Sacado::Fad::DFad<double> >& R, Table<1, Sacado::Fad::DFad<double>>& ULocal, Table<1, double >& ULocalConv, std::vector<SDdata<dim>> &cell_SDdata);
		int primiary_dof;	
		
		Battery_fields<dim>* battery_fields; 
		ElectricChemo<dim, Sacado::Fad::DFad<double>>* electricChemoFormula;
		Residual<Sacado::Fad::DFad<double>,dim>* ResidualEq;
		
};
#endif
