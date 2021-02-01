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
	params_json=this->params_mechanoChemFEM_json;
	
	electricChemoFormula.declare_parameters(*params_json);
	battery_fields.declare_parameters(*params_json);
	lithium.declare_parameters(*params_json);
	lithium_mu.declare_parameters(*params_json);
	lithium_cation.declare_parameters(*params_json);
	phi_s.declare_parameters(*params_json);
	phi_e.declare_parameters(*params_json);
	displacement.declare_parameters(*params_json);
	
	//Declear the parameters before load it
	this->load_parameters(parameter_file_Dir);		
	electricChemoFormula.init(battery_fields);
	define_battery_fields();
	this->define_primary_fields();
	//set active battery fields
	battery_fields.set_up_active_fields(this->primary_variables);	
	if(battery_fields.active_fields_index["Lithium"]>-1) lithium.set_up_fields(battery_fields, electricChemoFormula, this->ResidualEq, battery_fields.active_fields_index["Lithium"]);
	if(battery_fields.active_fields_index["Lithium_phaseField"]>-1) lithium_mu.set_up_fields(battery_fields, electricChemoFormula, this->ResidualEq, battery_fields.active_fields_index["Lithium_phaseField"]);
	if(battery_fields.active_fields_index["Lithium_cation"]>-1) lithium_cation.set_up_fields(battery_fields, electricChemoFormula, this->ResidualEq, battery_fields.active_fields_index["Lithium_cation"]);
	if(battery_fields.active_fields_index["Electrode_potential"]>-1) phi_s.set_up_fields(battery_fields, electricChemoFormula, this->ResidualEq, battery_fields.active_fields_index["Electrode_potential"]);
	if(battery_fields.active_fields_index["Electrolyte_potential"]>-1) phi_e.set_up_fields(battery_fields, electricChemoFormula, this->ResidualEq, battery_fields.active_fields_index["Electrolyte_potential"]);
	if(battery_fields.active_fields_index["Displacement"]>-1) displacement.set_up_fields(battery_fields, electricChemoFormula, this->ResidualEq, battery_fields.active_fields_index["Displacement"]);
	if(battery_fields.active_fields_index["Diffuse_interface"]>-1) diffuse_interface.set_up_fields(battery_fields, electricChemoFormula, this->ResidualEq, battery_fields.active_fields_index["Diffuse_interface"]);
	this->init_ibvp();
	
	//set up projection fields
	std::vector<std::vector<std::string> > computed_primary_variables={ {"VonMises", "component_is_scalar"}};
	computedNodalField.setup(battery_fields,*params_json);
	computedNodalField.setupComputedField(computed_primary_variables);
	computedNodalField.primary_variables_dof=this->primary_variables_dof;
}

template <int dim>
battery<dim>::~battery(){}
template <int dim>
void battery<dim>::define_battery_fields()
{
	std::vector<std::string> battery_fileds=(*params_json)["Problem"]["battery_variables"];
	
	std::vector<std::string> battery_fileds_s;
	std::vector<int> FE_support_list_v;
	if (battery_fileds[0]!="None"){
		for(unsigned int i=0;i<battery_fileds.size();i++){
			battery_fileds_s.push_back(battery_fileds[i]);
			if(battery_fileds[i]=="Displacement") battery_fileds_s.push_back("component_is_vector");
			else battery_fileds_s.push_back("component_is_scalar");
			FE_support_list_v.push_back(1);
		}
		(*params_json)["Problem"]["primary_variables_list"]=battery_fileds_s;
		(*params_json)["Problem"]["FE_support_list"]=FE_support_list_v;
	}
}

template <int dim>
void battery<dim>::declare_parameters(){}

template <int dim>
void battery<dim>::get_residual(const typename hp::DoFHandler<dim>::active_cell_iterator &cell, const FEValues<dim>& fe_values, Table<1, Sacado::Fad::DFad<double> >& R, Table<1, Sacado::Fad::DFad<double>>& ULocal, Table<1, double >& ULocalConv)
{	
  int cell_id = cell->active_cell_index();
	battery_fields.update_fields(cell, fe_values, ULocal, ULocalConv);

  // update reaction rate at the interface 
	double tem=(*params_json)["ElectroChemo"]["jn_react"];
	double fliptime=(*params_json)["ElectroChemo"]["flip_time"];
	
  cell_SDdata[cell_id].reaction_rate = tem;
	if(this->current_time>fliptime) cell_SDdata[cell_id].reaction_rate=-tem;

  if (cell_SDdata[cell_id].is_interface_element){
	  if(battery_fields.active_fields_index["Lithium"]>-1) lithium.r_get_residual_with_interface(cell, fe_values, R, ULocal, ULocalConv, cell_SDdata);
	  if(battery_fields.active_fields_index["Lithium_cation"]>-1) lithium_cation.r_get_residual_with_interface(cell, fe_values, R, ULocal, ULocalConv, cell_SDdata);
  }
  else{
	  if(battery_fields.active_fields_index["Lithium"]>-1) lithium.r_get_residual(fe_values, R, ULocal, ULocalConv);
	  if(battery_fields.active_fields_index["Lithium_cation"]>-1) lithium_cation.r_get_residual(fe_values, R, ULocal, ULocalConv);
  }
	if(battery_fields.active_fields_index["Lithium_phaseField"]>-1) lithium_mu.r_get_residual(fe_values, R, ULocal, ULocalConv);
	if(battery_fields.active_fields_index["Electrode_potential"]>-1) phi_s.r_get_residual(fe_values, R, ULocal, ULocalConv);
	if(battery_fields.active_fields_index["Displacement"]>-1) displacement.r_get_residual(fe_values, R, ULocal, ULocalConv);
	if(battery_fields.active_fields_index["Diffuse_interface"]>-1) diffuse_interface.r_get_residual(fe_values, R, ULocal, ULocalConv);
	
	apply_Neumann_boundary_condition();
}

template <int dim>
void battery<dim>::run()
{
  identify_diffuse_interface();

	this->pcout<<std::endl<<std::endl;
	this->pcout<<"======== RUNNING... ========"<<std::endl;	
	clock_t t_solve;	
	t_solve = clock();
  for (; this->current_time<=this->total_time; this->current_time+=this->current_dt){
    this->current_increment++;
		PetscPrintf(this->mpi_communicator,"************");
		PetscPrintf(this->mpi_communicator,"current increment=%d, current time= %f",this->current_increment, this->current_time);
		PetscPrintf(this->mpi_communicator,"************\n");
		this->solve_ibvp();
		
	  t_solve = clock() - t_solve;
		this->pcout<<"It took me"<< ((float)t_solve)/CLOCKS_PER_SEC<<" seconds for this solve"<<std::endl<<std::endl;
		
		this->FEMdata_out.clear_data_vectors();
		Vector<double> localized_U(this->solution_prev);
		this->FEMdata_out.data_out.add_data_vector (localized_U, computedNodalField);	
		std::string output_path = this->output_directory+"output-"+std::to_string(this->current_increment+this->off_output_index)+".vtk";
		this->FEMdata_out.write_vtk(this->solution_prev, output_path);	
    //this->output_results();
	}
	this->pcout<<"Finish running!!"<<std::endl;
}

template class battery<1>;
template class battery<2>;
template class battery<3>;
