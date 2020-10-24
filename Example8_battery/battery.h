/*
zhenlin wang 2019
*coupled diffusion reaction
*/
#ifndef battery_h
#include "mechanoChemFEM.h"
#include "battery_package/include/battery_components.h"

template <int dim>
class battery: public mechanoChemFEM<dim>
{
	public:
		battery(std::string parameter_file_Dir);
		//this is a overloaded function 
		void apply_boundary_condition();
		void apply_Neumann_boundary_condition();
		void declare_parameters();
		void get_residual(const typename hp::DoFHandler<dim>::active_cell_iterator &cell, const FEValues<dim>& fe_values, Table<1, Sacado::Fad::DFad<double> >& R, Table<1, Sacado::Fad::DFad<double>>& ULocal, Table<1, double >& ULocalConv);
		void setMultDomain();
		void output_w_domain();
		void make_grid();
		void define_battery_fields();
		
		
		ParameterHandler* params;		
		nlohmann::json* params_json;
		ConstraintMatrix* constraints;
		Battery_fields<dim> battery_fields;
		ElectricChemo<dim, Sacado::Fad::DFad<double>> electricChemoFormula;
		Lithium<dim> lithium;
		Lithium_phaseField<dim> lithium_mu;
		Lithium_cation<dim> lithium_cation;
		Electrode_potential<dim> phi_s;
		Electrolyte_potential<dim> phi_e;
		Displacement<dim> displacement;
		
};

#endif