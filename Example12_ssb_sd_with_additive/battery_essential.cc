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

  //std::cout << "before init_ibvp..." << std::endl;
	this->init_ibvp();
  //std::cout << "after init_ibvp..." << std::endl;
	
	//set up projection fields
	std::vector<std::vector<std::string> > computed_primary_variables={ {"VonMises", "component_is_scalar"}};
	computedNodalField.setup(battery_fields,*params_json);
	computedNodalField.setupComputedField(computed_primary_variables);
	computedNodalField.primary_variables_dof=this->primary_variables_dof;
}

template <int dim>
battery<dim>::~battery(){}

template <int dim>
void battery<dim>::make_grid()
{
  try
  {
    std::string mesh_directory=(*params_json)["Problem"]["mesh"];
    this->pcout << "reading external  mesh:"<<mesh_directory<<std::endl;
    GridIn<dim> gridin;
    gridin.attach_triangulation(this->triangulation);
    std::ifstream f(mesh_directory);
    gridin.read_abaqus(f);
    this->pcout << "reading external mesh (done):"<<mesh_directory<<std::endl;
  }
  catch (...)
  {
    double X_0,Y_0,Z_0,X_end,Y_end,Z_end;
    int element_div_x, element_div_y,element_div_z;
    
    X_0=(*params_json)["Geometry"]["x_min"];
    Y_0=(*params_json)["Geometry"]["y_min"];
    Z_0=(*params_json)["Geometry"]["z_min"];
    
    X_end=(*params_json)["Geometry"]["x_max"];
    Y_end=(*params_json)["Geometry"]["y_max"];
    Z_end=(*params_json)["Geometry"]["z_max"];;
    
    element_div_x=(*params_json)["Geometry"]["num_elem_x"].get<int>();
    element_div_y=(*params_json)["Geometry"]["num_elem_y"].get<int>();
    element_div_z=(*params_json)["Geometry"]["num_elem_z"].get<int>();
    
    bool colorize = false;
    std::vector< std::vector< double > > step_sizes;
    step_sizes.resize(dim);
    
    for (unsigned int j = 0; j < element_div_x; ++j) step_sizes[0].push_back((X_end-X_0)/element_div_x); 
    for (unsigned int j = 0; j < element_div_y; ++j) step_sizes[1].push_back((Y_end-Y_0)/element_div_y);
    if(dim==3)	for (unsigned int j = 0; j < element_div_z; ++j) step_sizes[2].push_back((Z_end-Z_0)/element_div_z);
    if(dim==2) GridGenerator::subdivided_hyper_rectangle (this->triangulation, step_sizes, Point<dim>(X_0,Y_0), Point<dim>(X_end,Y_end), colorize);
    else GridGenerator::subdivided_hyper_rectangle (this->triangulation, step_sizes, Point<dim>(X_0,Y_0,Z_0), Point<dim>(X_end,Y_end,Z_end), colorize);
  }
}
template <int dim>
void battery<dim>::define_battery_fields()
{
	std::vector<std::string> battery_fields=(*params_json)["Problem"]["battery_variables"];
	
	std::vector<std::string> battery_fields_s;
	std::vector<int> FE_support_list_v(battery_fields.size()*7); // 7 different material blocks
	if (battery_fields[0]!="None"){
		for(unsigned int i=0;i<battery_fields.size();i++){
			battery_fields_s.push_back(battery_fields[i]);
			if(battery_fields[i]=="Displacement") battery_fields_s.push_back("component_is_vector");
			else battery_fields_s.push_back("component_is_scalar");

		  //int  electrolyte_id=0, active_particle_id=1, interface_id=2, li_metal_id=3, additive_id=4, li_metal_interface_id=5, additive_interface_id=6;
     //"battery_variables":["Lithium","Lithium_cation","Electrode_potential","Electrolyte_potential","Diffuse_interface","Displacement"],

		  // define all at interface 
			//FE_support_list_v[i+2*battery_fields.size()]=1;

			if(battery_fields[i]=="Lithium" or battery_fields[i]=="Lithium_phaseField"){
				//FE_support_list_v[i]=1;
				//FE_support_list_v[i+battery_fields.size()]=0;

				FE_support_list_v[i+battery_fields.size()*electrolyte_id]=0;
				FE_support_list_v[i+battery_fields.size()*active_particle_id]=1;
				FE_support_list_v[i+battery_fields.size()*interface_id]=1;
				FE_support_list_v[i+battery_fields.size()*li_metal_id]=1;
				FE_support_list_v[i+battery_fields.size()*additive_id]=0;
				FE_support_list_v[i+battery_fields.size()*li_metal_interface_id]=1;
				FE_support_list_v[i+battery_fields.size()*additive_interface_id]=1;
			}
      else if(battery_fields[i]=="Electrode_potential"){
				//FE_support_list_v[i]=1;
				//FE_support_list_v[i+battery_fields.size()]=0;

				FE_support_list_v[i+battery_fields.size()*electrolyte_id]=0;
				FE_support_list_v[i+battery_fields.size()*active_particle_id]=1;
				FE_support_list_v[i+battery_fields.size()*interface_id]=1;
				FE_support_list_v[i+battery_fields.size()*li_metal_id]=1;
				FE_support_list_v[i+battery_fields.size()*additive_id]=1;
				FE_support_list_v[i+battery_fields.size()*li_metal_interface_id]=1; // change it to 0
				FE_support_list_v[i+battery_fields.size()*additive_interface_id]=1;
			}
			else if(battery_fields[i]=="Lithium_cation" or battery_fields[i]=="Electrolyte_potential"){
				//FE_support_list_v[i]=0;
				//FE_support_list_v[i+battery_fields.size()]=1;

				FE_support_list_v[i+battery_fields.size()*electrolyte_id]=1;
				FE_support_list_v[i+battery_fields.size()*active_particle_id]=0;
				FE_support_list_v[i+battery_fields.size()*interface_id]=1;
				FE_support_list_v[i+battery_fields.size()*li_metal_id]=0;
				FE_support_list_v[i+battery_fields.size()*additive_id]=0;
				FE_support_list_v[i+battery_fields.size()*li_metal_interface_id]=1;
				FE_support_list_v[i+battery_fields.size()*additive_interface_id]=0; // change it to 0
			}
			else if(battery_fields[i]=="Displacement"){
				FE_support_list_v[i+battery_fields.size()*electrolyte_id]=1;
				FE_support_list_v[i+battery_fields.size()*active_particle_id]=1;
				FE_support_list_v[i+battery_fields.size()*interface_id]=1;
				FE_support_list_v[i+battery_fields.size()*li_metal_id]=1;
				FE_support_list_v[i+battery_fields.size()*additive_id]=1;
				FE_support_list_v[i+battery_fields.size()*li_metal_interface_id]=1;
				FE_support_list_v[i+battery_fields.size()*additive_interface_id]=1;
			}
      else if(battery_fields[i]=="Diffuse_interface"){
				//FE_support_list_v[i]=1;
				//FE_support_list_v[i+battery_fields.size()]=0;
				FE_support_list_v[i+battery_fields.size()*electrolyte_id]=0;
				FE_support_list_v[i+battery_fields.size()*active_particle_id]=0;
				FE_support_list_v[i+battery_fields.size()*interface_id]=1;
				FE_support_list_v[i+battery_fields.size()*li_metal_id]=0;
				FE_support_list_v[i+battery_fields.size()*additive_id]=0;
				FE_support_list_v[i+battery_fields.size()*li_metal_interface_id]=1;
				FE_support_list_v[i+battery_fields.size()*additive_interface_id]=1;
			}
			//turn on all fields everywhere
			//FE_support_list_v[i]=1;
			//FE_support_list_v[i+battery_fields.size()]=1;
			//FE_support_list_v[i+2*battery_fields.size()]=1;
		}
		(*params_json)["Problem"]["primary_variables_list"]=battery_fields_s;
		(*params_json)["Problem"]["FE_support_list"]=FE_support_list_v;
	}
}

template <int dim>
void battery<dim>::declare_parameters(){}


template <int dim>
void battery<dim>::output_w_domain(){
	Vector<float> material_id(this->triangulation.n_active_cells()); 
  typename hp::DoFHandler<dim>::active_cell_iterator elem = this->dof_handler.begin_active(), endc = this->dof_handler.end();             
  unsigned int j = 0;                                                                                                      
  for (;elem!=endc; ++elem){                                                                                                
		material_id(j++) = elem->material_id();
	}
	std::string output_path = this->output_directory+"output_w_domain-"+std::to_string(this->current_increment)+".vtk";
	this->FEMdata_out.clear_data_vectors();
	this->FEMdata_out.data_out.add_data_vector(material_id, "material");
	this->FEMdata_out.write_vtk(this->solution_prev, output_path);	
}


template <int dim>
void battery<dim>::get_residual(const typename hp::DoFHandler<dim>::active_cell_iterator &cell, const FEValues<dim>& fe_values, Table<1, Sacado::Fad::DFad<double> >& R, Table<1, Sacado::Fad::DFad<double>>& ULocal, Table<1, double >& ULocalConv)
{	
  int cell_id = cell->active_cell_index();
	Point<dim> center=cell->center();
	double separator_line=(*params_json)["ElectroChemo"]["separator_line"];
	int orientation=(*params_json)["ElectroChemo"]["orientation"];
  //std::cout << " cell_id " << cell_id << " mat_id " << cell->material_id() << std::endl;

  if (cell->material_id()==interface_id){		 // 0
    //std::cout << " residual(interface): cell_id " << cell_id << std::endl;
    get_residual_at_diffuse_interface(cell, fe_values, R, ULocal, ULocalConv);
  }
  else if (cell->material_id()==li_metal_interface_id){		 // 5
    get_residual_at_li_metal_interface(cell, fe_values, R, ULocal, ULocalConv);
  }
  else if (cell->material_id()==additive_interface_id){		// 6
    get_residual_at_additive_interface(cell, fe_values, R, ULocal, ULocalConv);
  }
	else if (cell->material_id()==li_metal_id){ // 3
	  battery_fields.update_fields(cell, fe_values, ULocal, ULocalConv);
		if(battery_fields.active_fields_index["Lithium"]>-1) 
    { 
      lithium.set_cell_id(cell_id);
      lithium.r_get_residual(fe_values, R, ULocal, ULocalConv, pressure_old);
    }
		if(battery_fields.active_fields_index["Electrode_potential"]>-1) phi_s.r_get_residual(fe_values, R, ULocal, ULocalConv);
	  if(battery_fields.active_fields_index["Displacement"]>-1) 
    {
      displacement.set_cell_id(cell_id);
      displacement.r_get_residual(fe_values, R, ULocal, ULocalConv, pressure);
      //std::cout << "-particle- center " << center << " pressure "  << pressure[cell_id][0] << std::endl;
    }
	}
	else if (cell->material_id()==active_particle_id){ // 1
	  battery_fields.update_fields(cell, fe_values, ULocal, ULocalConv);
		if(battery_fields.active_fields_index["Lithium"]>-1) 
    { 
      lithium.set_cell_id(cell_id);
      lithium.r_get_residual(fe_values, R, ULocal, ULocalConv, pressure_old);
    }
		if(battery_fields.active_fields_index["Lithium_phaseField"]>-1) lithium_mu.r_get_residual(fe_values, R, ULocal, ULocalConv);
		if(battery_fields.active_fields_index["Electrode_potential"]>-1) phi_s.r_get_residual(fe_values, R, ULocal, ULocalConv);
	  if(battery_fields.active_fields_index["Displacement"]>-1) 
    {
      displacement.set_cell_id(cell_id);
      displacement.r_get_residual(fe_values, R, ULocal, ULocalConv, pressure);
      //std::cout << "-particle- center " << center << " pressure "  << pressure[cell_id][0] << std::endl;
    }
	}
	else if (cell->material_id()==electrolyte_id){ //0
	  battery_fields.update_fields(cell, fe_values, ULocal, ULocalConv);

	  if(battery_fields.active_fields_index["Lithium_cation"]>-1)
    {
      lithium_cation.set_cell_id(cell_id);
      lithium_cation.r_get_residual(fe_values, R, ULocal, ULocalConv, pressure_old);
    }
		if(battery_fields.active_fields_index["Electrolyte_potential"]>-1) phi_e.r_get_residual(fe_values, R, ULocal, ULocalConv);
	  if(battery_fields.active_fields_index["Displacement"]>-1)
    {
      displacement.set_cell_id(cell_id);
      displacement.r_get_residual(fe_values, R, ULocal, ULocalConv, pressure);
      //std::cout << "-electrolyte- center " << center << " pressure "  << pressure[cell_id][0] << std::endl;
    }
	}
	else if (cell->material_id()==additive_id){ // 4
	  battery_fields.update_fields(cell, fe_values, ULocal, ULocalConv);
		if(battery_fields.active_fields_index["Electrode_potential"]>-1) phi_s.r_get_residual(fe_values, R, ULocal, ULocalConv);
	  if(battery_fields.active_fields_index["Displacement"]>-1) 
    {
      displacement.set_cell_id(cell_id);
      displacement.r_get_residual(fe_values, R, ULocal, ULocalConv, pressure);
    }
	}

  pressure_gp0[cell_id] = pressure[cell_id][0];
  pressure_gp1[cell_id] = pressure[cell_id][1];
  pressure_gp2[cell_id] = pressure[cell_id][2];
  pressure_gp3[cell_id] = pressure[cell_id][3];
  jn[cell_id] = cell_SDdata[cell_id].reaction_rate_li.val();
  crack_id[cell_id] = cell_SDdata[cell_id].crack_id;
  T_n[cell_id] = cell_SDdata[cell_id].T_n;

  // the flux sign will be updated based on flip time
  double flux_sign = (*params_json)["ElectroChemo"]["flux_sign"];
  double fliptime=(*params_json)["ElectroChemo"]["flip_time"];
  //if (this->current_time >= fliptime)
  //{
    //flux_sign = flux_sign * -1.0;
    //(*params_json)["ElectroChemo"]["flux_sign"] = flux_sign;
  //}
  if (cell_id == 0) std::cout << " flux_sign " << flux_sign << " fliptime " << fliptime << std::endl;

	//apply_Neumann_boundary_condition();
	//BC
	for (unsigned int faceID=0; faceID<2*dim; faceID++){
		if(cell->face(faceID)->boundary_id()==1+orientation or cell->face(faceID)->boundary_id()==dim+1+orientation ){
			double current_IpA=(*params_json)["ElectroChemo"]["applied_current"];
      if(this->current_increment<=0){current_IpA=0; }
      else if(this->current_increment<=1){current_IpA=1*current_IpA/10; }
      else if(this->current_increment<=2){current_IpA=2*current_IpA/10; }
      else if(this->current_increment<=3){current_IpA=4*current_IpA/10; }
      else if(this->current_increment<=4){current_IpA=8*current_IpA/10; }

      if ( this->current_time >= fliptime + 0.0*this->current_dt and this->current_time < fliptime + 1.0*this->current_dt) {current_IpA=1*current_IpA/10; }
      if ( this->current_time >= fliptime + 1.0*this->current_dt and this->current_time < fliptime + 2.0*this->current_dt) {current_IpA=2*current_IpA/10; }
      if ( this->current_time >= fliptime + 2.0*this->current_dt and this->current_time < fliptime + 3.0*this->current_dt) {current_IpA=4*current_IpA/10; }
      if ( this->current_time >= fliptime + 3.0*this->current_dt and this->current_time < fliptime + 4.0*this->current_dt) {current_IpA=8*current_IpA/10; }

			if (cell->face(faceID)->boundary_id()==1+orientation) current_IpA=-current_IpA;
      //std::cout << " faceID= " << faceID << " BC_id= " << cell->face(faceID)->boundary_id() << " I= " << current_IpA << " dt= " << this->current_dt << std::endl;
		  FEFaceValues<dim> fe_face_values(fe_values.get_fe(), *(this->common_face_quadrature), update_values | update_quadrature_points | update_JxW_values);
			fe_face_values.reinit(cell,faceID);
			this->ResidualEq.residualForNeummanBC(fe_values, fe_face_values, battery_fields.active_fields_index["Electrode_potential"], R, current_IpA * flux_sign);
		}
	}
}

template <int dim>
void battery<dim>::output_results()
{
  //std::cout << " # of cells: " << this->triangulation.n_active_cells() << std::endl;
	Vector<float> material_id(this->triangulation.n_active_cells()); 
	Vector<float> subdomain_id(this->triangulation.n_active_cells()); 
  typename hp::DoFHandler<dim>::active_cell_iterator elem = this->dof_handler.begin_active(), endc = this->dof_handler.end();             
  unsigned int j = 0;                                                                                                      
  unsigned int _j = 0;                                                                                                      
  for (;elem!=endc; ++elem){                                                                                                
      //if (elem->subdomain_id() == this->this_mpi_process) {
		    material_id(j++) = elem->material_id();
		    subdomain_id(_j++) = elem->subdomain_id();
      //}
	}
  //std::cout << " after elem id " << std::endl;
	Vector<double> _jump_n(this->triangulation.n_active_cells()); 
	Vector<double> _jump_m(this->triangulation.n_active_cells()); 
	Vector<double> _jump_w(this->triangulation.n_active_cells()); 
  Utilities::MPI::sum(jump_n, MPI_COMM_WORLD, _jump_n);
  Utilities::MPI::sum(jump_m, MPI_COMM_WORLD, _jump_m);
  Utilities::MPI::sum(jump_w, MPI_COMM_WORLD, _jump_w);

  Utilities::MPI::sum(jn, MPI_COMM_WORLD, jn);
  Utilities::MPI::sum(pressure_gp0, MPI_COMM_WORLD, pressure_gp0);
  Utilities::MPI::sum(pressure_gp1, MPI_COMM_WORLD, pressure_gp1);
  Utilities::MPI::sum(pressure_gp2, MPI_COMM_WORLD, pressure_gp2);
  Utilities::MPI::sum(pressure_gp3, MPI_COMM_WORLD, pressure_gp3);
  Utilities::MPI::sum(crack_id, MPI_COMM_WORLD, crack_id);
  Utilities::MPI::sum(T_n, MPI_COMM_WORLD, T_n);

	//write vtk and snapshot for solution
	if(this->save_output){ 
		std::string output_path = this->output_directory+"output-new-"+std::to_string(this->current_increment)+".vtk";
		this->FEMdata_out.clear_data_vectors();
	  this->FEMdata_out.data_out.add_data_vector(material_id, "mat_id");
    this->FEMdata_out.data_out.add_data_vector(subdomain_id, "sub_id");

	  this->FEMdata_out.data_out.add_data_vector(_jump_n, "jump_n");
	  this->FEMdata_out.data_out.add_data_vector(_jump_m, "jump_m");
	  this->FEMdata_out.data_out.add_data_vector(_jump_w, "jump_w");

	  this->FEMdata_out.data_out.add_data_vector(jn, "jn");
	  this->FEMdata_out.data_out.add_data_vector(pressure_gp0, "p_gp0");
	  this->FEMdata_out.data_out.add_data_vector(pressure_gp1, "p_gp1");
	  this->FEMdata_out.data_out.add_data_vector(pressure_gp2, "p_gp2");
	  this->FEMdata_out.data_out.add_data_vector(pressure_gp3, "p_gp3");
	  this->FEMdata_out.data_out.add_data_vector(crack_id, "crack_id");
	  this->FEMdata_out.data_out.add_data_vector(T_n, "T_n");
		if(this->current_increment%this->skip_output==0) this->FEMdata_out.write_vtk(this->solution_prev, output_path);	

	}
	if(this->save_snapshot){
		std::string snapshot_path = this->snapshot_directory+"snapshot-"+std::to_string(this->current_increment+this->off_output_index)+".dat";
		this->FEMdata_out.create_vector_snapshot(this->solution, snapshot_path);
    save_sd_data();
	}
}

template <int dim>
void battery<dim>::solve_ibvp()
{
	bool to_flip=(*params_json)["ElectroChemo"]["to_flip"];
	bool is_converged = this->nonlinearSolve(this->solution);
  std::cout << " is converged? = " << is_converged << std::endl;
  this->solution_k = this->solution_prev;
  if (not is_converged)
  {
    to_flip = true;
    (*params_json)["ElectroChemo"]["to_flip"] = true;

    double flux_sign = (*params_json)["ElectroChemo"]["flux_sign"];
    flux_sign = flux_sign * -1.0;
    (*params_json)["ElectroChemo"]["flux_sign"] = flux_sign;
  }
  else
  {
    (*params_json)["ElectroChemo"]["to_flip"] = false;
  }

  std::cout << " to flip = " << (*params_json)["ElectroChemo"]["to_flip"] << std::endl;

	//update
  if (not to_flip)
  {
	  this->solution_prev=this->solution;
    std::cout << " solution updated " << std::endl;
  }
  else
  {
    this->solution_prev=this->solution_k;
    this->solution=this->solution_prev;
    std::cout << " solution NOT updated " << std::endl;
  }

}


template <int dim>
void battery<dim>::run()
{
  crack_id.reinit(this->triangulation.n_active_cells());
  jump_n.reinit(this->triangulation.n_active_cells());
  jump_m.reinit(this->triangulation.n_active_cells());
  jump_w.reinit(this->triangulation.n_active_cells());
  jn.reinit(this->triangulation.n_active_cells());
  pressure_gp0.reinit(this->triangulation.n_active_cells());
  pressure_gp1.reinit(this->triangulation.n_active_cells());
  pressure_gp2.reinit(this->triangulation.n_active_cells());
  pressure_gp3.reinit(this->triangulation.n_active_cells());
  T_n.reinit(this->triangulation.n_active_cells());
  pressure.resize(this->triangulation.n_active_cells()); // pressure needs to be loaded from restart files
  pressure_old.resize(this->triangulation.n_active_cells()); // pressure_old needs to be loaded from restart files
  is_new_step.resize(this->triangulation.n_active_cells());

  for (unsigned int count = 0; count < pressure.size(); count++)
  {
    pressure[count].resize(N_GPs);
    pressure_old[count].resize(N_GPs);

    for (unsigned int q = 0; q < N_GPs; ++q){
      pressure[count][q] = 0.0;
      pressure_old[count][q] = 0.0;
    }
  }

  //std::cout << "before run..." << std::endl;

	bool resuming_from_snapshot=(*params_json)["Problem"]["resuming_from_snapshot"];
	if(resuming_from_snapshot) {
		apply_initial_condition();		
    std::string snapfile=(*params_json)["Problem"]["snapshot_file"];
	  this->pcout<<"resuming from snapshot"<<std::endl;
	  this->FEMdata_out.resume_vector_from_snapshot(this->solution,snapfile);
	  std::string output_path = this->output_directory+"output-resume.vtk";
		this->FEMdata_out.clear_data_vectors();
    this->FEMdata_out.write_vtk(this->solution, output_path);
	  this->solution_prev=this->solution;
    load_sd_data();

    double resume_at_time=(*params_json)["Problem"]["resume_at_time"];
    this->current_time += resume_at_time;
    this->current_increment += this->off_output_index;

    for (unsigned int count = 0; count < cell_SDdata.size(); count++)
    {
      jn[count] = cell_SDdata[count].reaction_rate_li.val();
      crack_id[count] = cell_SDdata[count].crack_id;
      T_n[count] = cell_SDdata[count].T_n;
    }
	}
  else
  {
    this->off_output_index = 0;
  }


  std::cout << "before output w..." << std::endl;
	output_w_domain();
  std::cout << "after output w..." << std::endl;
  this->output_results();

	this->pcout<<std::endl<<std::endl;
	this->pcout<<"======== RUNNING... ========"<<std::endl;	

	clock_t t_solve;	
  double input_dt_from_param = this->current_dt;
  this->current_dt = std::min(1.0, input_dt_from_param);
  //(*params_json)["ElectroChemo"]["flip_time"] = 10.0;
	//double new_fliptime=(*params_json)["ElectroChemo"]["flip_time"];
  //std::cout << "fliptime " << fliptime << " new " << new_fliptime << " and " << (*params_json)["ElectroChemo"]["flip_time"] << std::endl;

  for (; this->current_time<=this->total_time; ){

	  double fliptime=(*params_json)["ElectroChemo"]["flip_time"];

	  t_solve = clock();

    // To use a small time step for the first 10s or the fliptime+10s.
    if (this->current_time>=10.0)
    {
      this->current_dt = input_dt_from_param;
    }
    // also make the ramp up with a small dt
    if (this->current_time + this->current_dt >= fliptime and this->current_time + this->current_dt < fliptime+10.0 )
    {
      this->current_dt = std::min(1.0, input_dt_from_param);
    }
    else if (this->current_time + this->current_dt >= fliptime+10.0 )
    {
      this->current_dt = input_dt_from_param;
    }
    this->current_time+=this->current_dt;
    this->current_increment++;

    // for parallel computing purpose
    jn.reinit(this->triangulation.n_active_cells());
    pressure_gp0.reinit(this->triangulation.n_active_cells());
    pressure_gp1.reinit(this->triangulation.n_active_cells());
    pressure_gp2.reinit(this->triangulation.n_active_cells());
    pressure_gp3.reinit(this->triangulation.n_active_cells());
    crack_id.reinit(this->triangulation.n_active_cells());
    T_n.reinit(this->triangulation.n_active_cells());

		PetscPrintf(this->mpi_communicator,"************");
		PetscPrintf(this->mpi_communicator,"current increment=%d, current time= %f",this->current_increment, this->current_time);
		PetscPrintf(this->mpi_communicator,"************\n");

    {
      typename hp::DoFHandler<dim>::active_cell_iterator elem = this->dof_handler.begin_active(), endc = this->dof_handler.end();             
      unsigned int j = 0;                                                                                                      
      for (;elem!=endc; ++elem){                                                                                                
        int cell_id = elem->active_cell_index();
        is_new_step[cell_id] = true;
	    }
    }

		this->solve_ibvp();
	  bool to_flip=(*params_json)["ElectroChemo"]["to_flip"]; // has to be here. to_flip will be updated in solve_ibvp()
    if (to_flip)
    {
      (*params_json)["ElectroChemo"]["flip_time"] = this->current_time;
      std::cout << " flip discharge/charge sign at: " << this->current_time << " (s) " << std::endl;
    }
    else
    {
      // update history variables
      for (unsigned i = 0; i < cell_SDdata.size(); ++i) {
        cell_SDdata[i].xi_conv = cell_SDdata[i].xi_old;
        cell_SDdata[i].xi_conv_c_e = cell_SDdata[i].xi_old_c_e;
        cell_SDdata[i].xi_conv_phi_s = cell_SDdata[i].xi_old_phi_s;
        cell_SDdata[i].xi_conv_phi_e = cell_SDdata[i].xi_old_phi_e;
        cell_SDdata[i].C_Li_plus_old = cell_SDdata[i].C_Li_plus_new;

        cell_SDdata[i].xi_conv_u_sd = cell_SDdata[i].xi_old_u_sd;

        cell_SDdata[i].Tn_old = cell_SDdata[i].Tn_new;
        // WARNING: should not use the following, as many cell_SDdata are not initialized.
        //std::cout << "C_Li_plus_old[0]" << cell_SDdata[i].C_Li_plus_old[0] << std::endl;
        //for (int q=0; q<4; q++) cell_SDdata[i].C_Li_plus_old[q] = cell_SDdata[i].C_Li_plus_new[q];
        //std::cout << "C_Li_plus_old[3]" << cell_SDdata[i].C_Li_plus_old[3] << std::endl;
      }
      std::cout << " SD updated! " << std::endl;
    }
		
	  t_solve = clock() - t_solve;
		this->pcout<<"It took me "<< ((float)t_solve)/CLOCKS_PER_SEC<<" seconds for this solve"<<std::endl<<std::endl;

		// std::string snapfile="snapshot_phase_2/snapshot-"+std::to_string(this->current_increment)+".dat";
		// 	  this->FEMdata_out.resume_vector_from_snapshot(this->solution,snapfile);
		// 	  this->solution_prev=this->solution;
		
     //this->FEMdata_out.clear_data_vectors();
     //Vector<double> localized_U(this->solution_prev);
     //this->FEMdata_out.data_out.add_data_vector (localized_U, computedNodalField);
     //std::string output_path = this->output_directory+"output-"+std::to_string(this->current_increment)+".vtk";
     //this->FEMdata_out.write_vtk(this->solution_prev, output_path);

    if (not to_flip)
    {
      this->output_results();
      pressure_old = pressure;
      std::cout << " pressure updated " << std::endl;
    }
    else
    {
      int resume_non_conv_back_at = 7;
      (*params_json)["Problem"]["sd_data_file"] = this->snapshot_directory+"restart-"+std::to_string(this->current_increment +this->off_output_index - resume_non_conv_back_at)+".dat";

      std::cout << "restart from " << this->current_increment << " " <<  this->off_output_index << " " <<  resume_non_conv_back_at << " " << this->current_increment +this->off_output_index - resume_non_conv_back_at << std::endl;

      std::cout << " load sd data from " <<  (*params_json)["Problem"]["sd_data_file"] << std::endl;
      (*params_json)["Problem"]["snapshot_file"] = this->snapshot_directory+"snapshot-"+std::to_string(this->current_increment+this->off_output_index - resume_non_conv_back_at)+".dat";

      this->current_increment = this->current_increment - resume_non_conv_back_at + 2;

      std::string snapfile=(*params_json)["Problem"]["snapshot_file"];
	  this->pcout<<"resuming from snapshot " << snapfile <<std::endl;
	  this->FEMdata_out.resume_vector_from_snapshot(this->solution,snapfile);
	  this->solution_prev=this->solution;
      load_sd_data();
      for (unsigned int count = 0; count < cell_SDdata.size(); count++)
      {
        jn[count] = cell_SDdata[count].reaction_rate_li.val();
        crack_id[count] = cell_SDdata[count].crack_id;
        T_n[count] = cell_SDdata[count].T_n;
      }
    } // not to_flip
  } // time loop
  this->pcout<<"Finish running!!"<<std::endl;
}

template class battery<1>;
template class battery<2>;
template class battery<3>;
