/*
zhenlin wang 2020
*battery package_battery_fields
*/
#include "../include/Battery_fields.h"
template <int dim>
Battery_fields<dim>::Battery_fields():pcout (std::cout, (Utilities::MPI::this_mpi_process(MPI_COMM_WORLD) == 0))
{
	active_fields_index["Lithium"]=-1;
	active_fields_index["Lithium_phaseField"]=-1;
	active_fields_index["Lithium_cation"]=-1;
	active_fields_index["Electrode_potential"]=-1;
	active_fields_index["Electrolyte_potential"]=-1;
	active_fields_index["Temperature"]=-1;
	active_fields_index["Displacement"]=-1;
}

template <int dim>
void Battery_fields<dim>::set_up_active_fields(std::vector<std::vector<std::string> > primary_variables)
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
	this->pcout<<std::endl;
	this->pcout<<"=========== Battery Fields Defined ==========="<<std::endl;
}

template <int dim>
void Battery_fields<dim>::declare_parameters(nlohmann::json& _params)
{
	params_json=&_params;
	(*params_json)["Geometry"]["pos_electrode_len"]=0.0;
	(*params_json)["Geometry"]["num_elem_pos_electrode"]=1;
	(*params_json)["Geometry"]["separator_len"]=0;
	(*params_json)["Geometry"]["num_elem_speparator"]=1;
	(*params_json)["Geometry"]["neg_electrode_len"]=0.0;
	(*params_json)["Geometry"]["num_elem_neg_electrode"]=1;
	(*params_json)["Problem"]["battery_variables"]={"None"};
}

template <int dim>
void Battery_fields<dim>::update_fields(const typename hp::DoFHandler<dim>::active_cell_iterator &cell, const FEValues<dim>& fe_values, Table<1, Sacado::Fad::DFad<double>>& ULocal, Table<1, double >& ULocalConv)
{
	current_cell=&cell;
	current_domain_id=cell->material_id();
	unsigned int n_q_points= fe_values.n_quadrature_points;
	TableIndices<1> tableIndex1(n_q_points);
	TableIndices<2> tableIndex2(n_q_points,n_q_points);
	
	for (std::map<std::string,int>::iterator it=active_fields_index.begin(); it!=active_fields_index.end(); ++it){
		if(it->second>-1 and std::strcmp(it->first.c_str(),"Displacement")!=0){
			int primiary_dof=it->second;
			quad_fields[primiary_dof].value_conv.reinit(tableIndex1);
			quad_fields[primiary_dof].value.reinit(tableIndex1);
			quad_fields[primiary_dof].value_grad.reinit(tableIndex2);
			evaluateScalarFunction<double,dim>(fe_values, primiary_dof, ULocalConv, quad_fields[primiary_dof].value_conv);
			evaluateScalarFunction<Sacado::Fad::DFad<double>,dim>(fe_values, primiary_dof, ULocal, quad_fields[primiary_dof].value);
			evaluateScalarFunctionGradient<Sacado::Fad::DFad<double>,dim>(fe_values, primiary_dof, ULocal, quad_fields[primiary_dof].value_grad);
		}
	}
}

template class Battery_fields<1>;
template class Battery_fields<2>;
template class Battery_fields<3>;