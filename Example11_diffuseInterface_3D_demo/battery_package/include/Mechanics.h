/*
zhenlin wang 2019
*module transportation
*/
#ifndef Mechanics_h
#define Mechanics_h

#include "mechanoChemFEM.h"
#include <deal.II/base/table_indices.h>
#include "Battery_fields.h"
#include "ElectricChemo.h"

template <int dim>
class Mechanics
{
	public:
		Mechanics();
		Mechanics(Battery_fields<dim>& _fields, Residual<Sacado::Fad::DFad<double>,dim>& ResidualEq);
		Mechanics(Battery_fields<dim>& _fields, Residual<Sacado::Fad::DFad<double>,dim>& ResidualEq, int _primiary_dof);
		//this is a overloaded function 
		void set_up_fields(Battery_fields<dim>& _battery_fields, Residual<Sacado::Fad::DFad<double>,dim>& ResidualEq, int _primiary_dof);
		void set_up_fields(Battery_fields<dim>& _battery_fields, ElectricChemo<dim,Sacado::Fad::DFad<double>>& _electricChemoFormula, Residual<Sacado::Fad::DFad<double>,dim>& ResidualEq, int _primiary_dof);
		void set_primiary_dof(int _primiary_dof);
		virtual void set_stress(dealii::Table<3,Sacado::Fad::DFad<double> >& F, dealii::Table<3,Sacado::Fad::DFad<double> >& P);
		void r_get_residual(const FEValues<dim>& fe_values, Table<1, Sacado::Fad::DFad<double> >& R, Table<1, Sacado::Fad::DFad<double>>& ULocal, Table<1, double >& ULocalConv);
		int primiary_dof;	
		
		Battery_fields<dim>* battery_fields; 
		ElectricChemo<dim, Sacado::Fad::DFad<double>>* electricChemoFormula;
		Residual<Sacado::Fad::DFad<double>,dim>* ResidualEq;
		double youngsModulus=1,	poissonRatio=0.3;
		
};
#endif