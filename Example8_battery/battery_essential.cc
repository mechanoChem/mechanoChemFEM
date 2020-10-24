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
	output_w_domain();
}

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
void battery<dim>::make_grid()
{
	this->pcout<<"make grid..."<<std::endl;
	double X_0,Y_0,Z_0,X_end,Y_end,Z_end;
	int element_div_x, element_div_y,element_div_z;
	
	X_0=(*params_json)["Geometry"]["x_min"];
	Y_0=(*params_json)["Geometry"]["y_min"];
	Z_0=(*params_json)["Geometry"]["z_min"];

	X_end=(*params_json)["Geometry"]["x_max"];
	Y_end=(*params_json)["Geometry"]["y_max"];
	Z_end=(*params_json)["Geometry"]["z_max"];

	element_div_x=(*params_json)["Geometry"]["num_elem_x"].get<int>();
	element_div_y=(*params_json)["Geometry"]["num_elem_y"].get<int>();
	element_div_z=(*params_json)["Geometry"]["num_elem_z"].get<int>();
	
	double pos_electrode_len=(*params_json)["Geometry"]["pos_electrode_len"];
	double separator_len=(*params_json)["Geometry"]["separator_len"];
	double neg_electrode_len=(*params_json)["Geometry"]["neg_electrode_len"];
	int num_elem_pos_electrode=(*params_json)["Geometry"]["num_elem_pos_electrode"];
	int num_elem_separator=(*params_json)["Geometry"]["num_elem_separator"];
	int num_elem_neg_electrode=(*params_json)["Geometry"]["num_elem_neg_electrode"];
	bool colorize = false;
  std::vector< std::vector< double > > step_sizes;
  step_sizes.resize(dim);
	
	if(std::abs(pos_electrode_len+neg_electrode_len+separator_len-X_end+X_0)>1.0e-3){
    throw std::invalid_argument( "pos_electrode_len+neg_electrode_len+separator_len!=(X_end-X_0) ");
	}
	
  for (unsigned int j = 0; j < num_elem_pos_electrode; ++j) step_sizes[0].push_back(pos_electrode_len/num_elem_pos_electrode); 
	for (unsigned int j = 0; j < num_elem_separator; ++j) step_sizes[0].push_back(separator_len/num_elem_separator); 
	for (unsigned int j = 0; j < num_elem_neg_electrode; ++j) step_sizes[0].push_back(neg_electrode_len/num_elem_neg_electrode); 
	
  for (unsigned int j = 0; j < element_div_y; ++j) step_sizes[1].push_back((Y_end-Y_0)/element_div_y);
	if(dim==3)	for (unsigned int j = 0; j < element_div_z; ++j) step_sizes[2].push_back((Z_end-Z_0)/element_div_z);
	
  if(dim==2) GridGenerator::subdivided_hyper_rectangle (this->triangulation, step_sizes, Point<dim>(X_0,Y_0), Point<dim>(X_end,Y_end), colorize);
	else GridGenerator::subdivided_hyper_rectangle (this->triangulation, step_sizes, Point<dim>(X_0,Y_0,Z_0), Point<dim>(X_end,Y_end,Z_end), colorize);
}

template <int dim>
void battery<dim>::setMultDomain()
{
	this->pcout<<"setMultDomain... neg_electrode | separator | pos_electrode"<<std::endl;
	double X_0,Y_0,Z_0,X_end,Y_end,Z_end;
	int element_div_x, element_div_y,element_div_z;
	
	X_0=(*params_json)["Geometry"]["x_min"];
	Y_0=(*params_json)["Geometry"]["y_min"];
	Z_0=(*params_json)["Geometry"]["z_min"];

	X_end=(*params_json)["Geometry"]["x_max"];
	Y_end=(*params_json)["Geometry"]["y_max"];
	Z_end=(*params_json)["Geometry"]["z_max"];
	
	double pos_electrode_len=(*params_json)["Geometry"]["pos_electrode_len"];
	double separator_len=(*params_json)["Geometry"]["separator_len"];
	double neg_electrode_len=(*params_json)["Geometry"]["neg_electrode_len"];
  for (typename Triangulation<dim>::active_cell_iterator cell = this->dof_handler.begin_active(); cell != this->dof_handler.end(); ++cell){
    Point<dim> cell_center = cell->center();
		if(cell_center[0]<X_0+neg_electrode_len) cell->set_material_id(battery_fields.neg_electrode_domain_id);
		else if(cell_center[0]<X_0+neg_electrode_len+separator_len) cell->set_material_id(battery_fields.separator_domain_id);
		else cell->set_material_id(battery_fields.pos_electrode_domain_id);
	}
	//assign Fe_system to corresponding cells
	this->set_active_fe_indices (this->FE_support, this->dof_handler);
}

template <int dim>
void battery<dim>::get_residual(const typename hp::DoFHandler<dim>::active_cell_iterator &cell, const FEValues<dim>& fe_values, Table<1, Sacado::Fad::DFad<double> >& R, Table<1, Sacado::Fad::DFad<double>>& ULocal, Table<1, double >& ULocalConv)
{	
	battery_fields.update_fields(cell, fe_values, ULocal, ULocalConv);
	if(cell->material_id()==battery_fields.separator_domain_id) {
		if(battery_fields.active_fields_index["Lithium_cation"]>-1) lithium_cation.r_get_residual(fe_values, R, ULocal, ULocalConv);
		if(battery_fields.active_fields_index["Electrolyte_potential"]>-1) phi_e.r_get_residual(fe_values, R, ULocal, ULocalConv);
	}
	else{
		if(battery_fields.active_fields_index["Lithium"]>-1) lithium.r_get_residual(fe_values, R, ULocal, ULocalConv);
		if(battery_fields.active_fields_index["Lithium_phaseField"]>-1) lithium_mu.r_get_residual(fe_values, R, ULocal, ULocalConv);
		if(battery_fields.active_fields_index["Electrode_potential"]>-1) phi_s.r_get_residual(fe_values, R, ULocal, ULocalConv);
		if(battery_fields.active_fields_index["Displacement"]>-1) displacement.r_get_residual(fe_values, R, ULocal, ULocalConv);
	}
	apply_Neumann_boundary_condition();
}

template class battery<1>;
template class battery<2>;
template class battery<3>;