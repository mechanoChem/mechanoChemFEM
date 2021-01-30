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
	
	this->init_ibvp();
	output_diffuse_interface();
	//setMultDomain();
	//output_w_domain();

	//exit(-1);
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
			for(unsigned int j=0;j<3;j++) FE_support_list_v.push_back(1);
		}
		(*params_json)["Problem"]["primary_variables_list"]=battery_fileds_s;
		(*params_json)["Problem"]["FE_support_list"]=FE_support_list_v;
	}
}

template <int dim>
void battery<dim>::setup_diffuse_interface_FEM()
{
	this->pcout<<"setup_diffuse_interface"<<std::endl;
	dof_handler_interface=new hp::DoFHandler<dim>(this->triangulation);
	std::vector<unsigned int > variables_dof_tem;
	variables_interface.resize(1);
  variables_interface[0].push_back("diffuse_interface");
	variables_interface[0].push_back("component_is_scalar");
	std::vector<std::vector<int> > FE_support_interface(1);
	FE_support_interface[0].push_back(1);	
	this->setup_FeSystem(fe_system_interface,fe_collection_interface, q_collection_interface, variables_dof_tem,variables_interface,FE_support_interface,*(this->volume_quadrature) );
	dof_handler_interface->distribute_dofs (fe_collection_interface);
	DoFRenumbering::component_wise (*dof_handler_interface);
	const types::global_dof_index n_local_dofs = DoFTools::count_dofs_with_subdomain_association(*dof_handler_interface, this->this_mpi_process);
	const types::global_dof_index n_total_dofs=dof_handler_interface->n_dofs();
										
	diffuse_interface.reinit (this->mpi_communicator,n_total_dofs,n_local_dofs); 
	setup_diffuse_interface();
  identify_diffuse_interface();
}

template <int dim>
void battery<dim>::output_diffuse_interface(){
	this->pcout<<"output_diffuse_interface"<<std::endl;
	dealii::Vector<double> diffuse_interface_(diffuse_interface);
	std::string output_path = this->output_directory+"output_w_domain-"+std::to_string(this->current_increment+this->off_output_index)+".vtk";
	this->FEMdata_out.clear_data_vectors();
	this->FEMdata_out.data_out.add_data_vector(diffuse_interface_, "diffuse_interface");
	this->FEMdata_out.write_vtk(this->solution_prev, output_path);	
}

template <int dim>
void battery<dim>::output_w_domain(){
	Vector<float> material_id(this->triangulation.n_active_cells()); 
  typename hp::DoFHandler<dim>::active_cell_iterator elem = this->dof_handler.begin_active(), endc = this->dof_handler.end();             
  unsigned int j = 0;                                                                                                      
  for (;elem!=endc; ++elem){                                                                                                
		material_id(j++) = elem->material_id();
	}
	std::string output_path = this->output_directory+"output_w_domain-"+std::to_string(this->current_increment+this->off_output_index)+".vtk";
	this->FEMdata_out.clear_data_vectors();
	this->FEMdata_out.data_out.add_data_vector(material_id, "material");
	this->FEMdata_out.write_vtk(this->solution_prev, output_path);	
}

template <int dim>
void battery<dim>::declare_parameters(){}

template <int dim>
void battery<dim>::get_residual(const typename hp::DoFHandler<dim>::active_cell_iterator &cell, const FEValues<dim>& fe_values, Table<1, Sacado::Fad::DFad<double> >& R, Table<1, Sacado::Fad::DFad<double>>& ULocal, Table<1, double >& ULocalConv)
{	
  int cell_id = cell->active_cell_index();
	battery_fields.update_fields(cell, fe_values, ULocal, ULocalConv);

  // update reaction rate at the interface 
  cell_SDdata[cell_id].reaction_rate = 0.01;

  if (cell_SDdata[cell_id].is_interface_element){
	  if(battery_fields.active_fields_index["Lithium"]>-1) lithium.r_get_residual_with_interface(cell, fe_values, R, ULocal, ULocalConv, cell_SDdata);
  }
  else{
	  if(battery_fields.active_fields_index["Lithium"]>-1) lithium.r_get_residual(fe_values, R, ULocal, ULocalConv);
  }

	if(battery_fields.active_fields_index["Lithium_phaseField"]>-1) lithium_mu.r_get_residual(fe_values, R, ULocal, ULocalConv);
	if(battery_fields.active_fields_index["Electrode_potential"]>-1) phi_s.r_get_residual(fe_values, R, ULocal, ULocalConv);
	if(battery_fields.active_fields_index["Displacement"]>-1) displacement.r_get_residual(fe_values, R, ULocal, ULocalConv);
	
	
	apply_Neumann_boundary_condition();

}

template class battery<1>;
template class battery<2>;
template class battery<3>;
