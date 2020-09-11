/*
zhenlin wang 2019
*module transportation
*/
#include "../include/Battery_fields.h"
Battery_fields::Battery_fields():pcout (std::cout, (Utilities::MPI::this_mpi_process(MPI_COMM_WORLD) == 0))
{
	active_fields_index["Lithium"]=-1;
	active_fields_index["Lithium_phaseField"]=-1;
	active_fields_index["Lithium_cation"]=-1;
	active_fields_index["Electrode_potential"]=-1;
	active_fields_index["Electrolyte_potential"]=-1;
	active_fields_index["Temperature"]=-1;
	active_fields_index["Displacement"]=-1;
	//declare_parameters_batteryFields();
}

void Battery_fields::set_up_active_fields(std::vector<std::vector<std::string> > primary_variables, int dim)
{
	int active_dof=0;
	std::map<std::string,int>::iterator it;
	for(unsigned int i=0;i<primary_variables.size();i++){
		it = active_fields_index.find(primary_variables[i][0]);
		if (it == active_fields_index.end())this->pcout<<primary_variables[i][0]<<" is not a battery field"<<std::endl;
		else {it->second=active_dof; this->pcout<<primary_variables[i][0]<<" is a battery field, its index is assigned to be "<<active_dof<<std::endl;}
		if(std::strcmp(primary_variables[i][1].c_str(),"component_is_scalar")==0) active_dof++;
		else if(std::strcmp(primary_variables[i][1].c_str(),"component_is_vector")==0) active_dof+=dim;
	}


	quad_fields.resize(active_dof);
	this->pcout<<"=========== Battery Fields Defined ==========="<<std::endl;
}

void Battery_fields::declare_parameters_batteryFields()
{
	params_battery_fileds->enter_subsection("battery_fields");
	params_battery_fileds->declare_entry("lithium_index","-1",Patterns::Integer() );
	params_battery_fileds->declare_entry("Lithium_phaseField","-1",Patterns::Integer() );
	params_battery_fileds->declare_entry("Lithium_cation","-1",Patterns::Integer() );
	params_battery_fileds->declare_entry("Electrode_potential","-1",Patterns::Integer() );
	params_battery_fileds->declare_entry("Electrolyte_potential","-1",Patterns::Integer() );
	params_battery_fileds->declare_entry("Temperature","-1",Patterns::Integer() );
	params_battery_fileds->declare_entry("Displacement","-1",Patterns::Integer() );
	params_battery_fileds->leave_subsection();	
}