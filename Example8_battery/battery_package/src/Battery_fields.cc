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

void Battery_fields::set_up_active_fields(std::vector<std::vector<std::string> > primary_variables, int _dim)
{
	int active_dof=0;
	dim=_dim;
	std::map<std::string,int>::iterator it;
	for(unsigned int i=0;i<primary_variables.size();i++){
		it = active_fields_index.find(primary_variables[i][0]);
		if (it == active_fields_index.end())this->pcout<<primary_variables[i][0]<<" is not a battery field"<<std::endl;
		else {it->second=active_dof; this->pcout<<primary_variables[i][0]<<" is a battery field, its index is assigned to be "<<active_dof<<std::endl;}
		if(std::strcmp(primary_variables[i][1].c_str(),"component_is_scalar")==0) active_dof++;
		else if(std::strcmp(primary_variables[i][1].c_str(),"component_is_vector")==0) active_dof+=dim;
	}


	quad_fields.resize(active_dof);
	this->pcout<<std::endl;
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

void Battery_fields::update_fields(const FEValues<1>& fe_values, Table<1, Sacado::Fad::DFad<double>>& ULocal, Table<1, double >& ULocalConv)
{
	unsigned int n_q_points= fe_values.n_quadrature_points;
	TableIndices<1> tableIndex1(n_q_points);
	TableIndices<2> tableIndex2(n_q_points,n_q_points);
	
	for (std::map<std::string,int>::iterator it=active_fields_index.begin(); it!=active_fields_index.end(); ++it){
		if(it->second>-1 and std::strcmp(it->first.c_str(),"Displacement")!=0){
			int primiary_dof=it->second;
			quad_fields[primiary_dof].value_conv.reinit(tableIndex1);
			quad_fields[primiary_dof].value.reinit(tableIndex1);
			quad_fields[primiary_dof].value_grad.reinit(tableIndex2);
			evaluateScalarFunction<double,1>(fe_values, primiary_dof, ULocalConv, quad_fields[primiary_dof].value_conv);
			evaluateScalarFunction<Sacado::Fad::DFad<double>,1>(fe_values, primiary_dof, ULocal, quad_fields[primiary_dof].value);
			evaluateScalarFunctionGradient<Sacado::Fad::DFad<double>,1>(fe_values, primiary_dof, ULocal, quad_fields[primiary_dof].value_grad);
		}
	}
}

void Battery_fields::update_fields(const FEValues<2>& fe_values, Table<1, Sacado::Fad::DFad<double>>& ULocal, Table<1, double >& ULocalConv)
{
	unsigned int n_q_points= fe_values.n_quadrature_points;
	TableIndices<1> tableIndex1(n_q_points);
	TableIndices<2> tableIndex2(n_q_points,n_q_points);
	
	for (std::map<std::string,int>::iterator it=active_fields_index.begin(); it!=active_fields_index.end(); ++it){
		if(it->second>-1 and std::strcmp(it->first.c_str(),"Displacement")!=0){
			int primiary_dof=it->second;
			quad_fields[primiary_dof].value_conv.reinit(tableIndex1);
			quad_fields[primiary_dof].value.reinit(tableIndex1);
			quad_fields[primiary_dof].value_grad.reinit(tableIndex2);
			evaluateScalarFunction<double,2>(fe_values, primiary_dof, ULocalConv, quad_fields[primiary_dof].value_conv);
			evaluateScalarFunction<Sacado::Fad::DFad<double>,2>(fe_values, primiary_dof, ULocal, quad_fields[primiary_dof].value);
			evaluateScalarFunctionGradient<Sacado::Fad::DFad<double>,2>(fe_values, primiary_dof, ULocal, quad_fields[primiary_dof].value_grad);
		}
	}
}

void Battery_fields::update_fields(const FEValues<3>& fe_values, Table<1, Sacado::Fad::DFad<double>>& ULocal, Table<1, double >& ULocalConv)
{
	unsigned int n_q_points= fe_values.n_quadrature_points;
	TableIndices<1> tableIndex1(n_q_points);
	TableIndices<2> tableIndex2(n_q_points,n_q_points);
	
	for (std::map<std::string,int>::iterator it=active_fields_index.begin(); it!=active_fields_index.end(); ++it){
		if(it->second>-1 and std::strcmp(it->first.c_str(),"Displacement")!=0){
			int primiary_dof=it->second;
			quad_fields[primiary_dof].value_conv.reinit(tableIndex1);
			quad_fields[primiary_dof].value.reinit(tableIndex1);
			quad_fields[primiary_dof].value_grad.reinit(tableIndex2);
			evaluateScalarFunction<double,3>(fe_values, primiary_dof, ULocalConv, quad_fields[primiary_dof].value_conv);
			evaluateScalarFunction<Sacado::Fad::DFad<double>,3>(fe_values, primiary_dof, ULocal, quad_fields[primiary_dof].value);
			evaluateScalarFunctionGradient<Sacado::Fad::DFad<double>,3>(fe_values, primiary_dof, ULocal, quad_fields[primiary_dof].value_grad);
		}
	}
}






