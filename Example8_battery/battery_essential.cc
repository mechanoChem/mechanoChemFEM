/*
zhenlin wang 2019
*coupled diffusion reaction
*/
#include "battery.h"

template <int dim>
battery<dim>::battery(std::string parameter_file_Dir)
{
	//pass the pointer to "constraints" in that defined in mechanoChemFEM
	constraints=this->constraints_mechanoChemFEM;
	//This let you use one params to get all parameters pre-defined in the mechanoChemFEM
	params=this->params_mechanoChemFEM;
	lithium.declare_parameters(*params);
	displacement.declare_parameters(*params);
	lithium_mu.declare_parameters(*params);
	
	params_json=this->params_mechanoChemFEM_json;
	lithium.declare_parameters(*params_json);
	displacement.declare_parameters(*params_json);
	lithium_mu.declare_parameters(*params_json);
	
	//Declear the parameters before load it
	this->load_parameters(parameter_file_Dir);		

	this->define_primary_fields();

	
	//set active battery fields
	battery_fields.set_up_active_fields(this->primary_variables,dim);	
	if(battery_fields.active_fields_index["Lithium"]>-1) lithium.set_up_fields(battery_fields, this->ResidualEq, battery_fields.active_fields_index["Lithium"]);
	if(battery_fields.active_fields_index["Lithium_phaseField"]>-1) lithium_mu.set_up_fields(battery_fields, this->ResidualEq, battery_fields.active_fields_index["Lithium_phaseField"]);
	if(battery_fields.active_fields_index["Displacement"]>-1) displacement.set_up_fields(battery_fields,this->ResidualEq, battery_fields.active_fields_index["Displacement"]);

	this->init_ibvp();
}


template <int dim>
void battery<dim>::get_residual(const typename hp::DoFHandler<dim>::active_cell_iterator &cell, const FEValues<dim>& fe_values, Table<1, Sacado::Fad::DFad<double> >& R, Table<1, Sacado::Fad::DFad<double>>& ULocal, Table<1, double >& ULocalConv)
{	
	battery_fields.update_fields(fe_values,ULocal, ULocalConv);
	if(battery_fields.active_fields_index["Lithium"]>-1) lithium.r_get_residual(cell,fe_values,R,ULocal, ULocalConv);
	if(battery_fields.active_fields_index["Lithium_phaseField"]>-1) lithium_mu.r_get_residual(cell,fe_values,R,ULocal, ULocalConv);
	if(battery_fields.active_fields_index["Displacement"]>-1) displacement.r_get_residual(cell,fe_values,R,ULocal, ULocalConv);
	apply_Neumann_boundary_condition();
}

template class battery<1>;
template class battery<2>;
template class battery<3>;