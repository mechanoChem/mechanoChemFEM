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
void battery<dim>::make_grid()
{
	std::string mesh_directory=(*params_json)["Problem"]["mesh"];
	this->pcout << "reading external  mesh:"<<mesh_directory<<std::endl;
  GridIn<dim> gridin;
  gridin.attach_triangulation(this->triangulation);
  std::ifstream f(mesh_directory);
  gridin.read_abaqus(f);
}
template <int dim>
void battery<dim>::define_battery_fields()
{
	std::vector<std::string> battery_fileds=(*params_json)["Problem"]["battery_variables"];
	
	std::vector<std::string> battery_fileds_s;
	std::vector<int> FE_support_list_v(battery_fileds.size()*3);
	if (battery_fileds[0]!="None"){
		for(unsigned int i=0;i<battery_fileds.size();i++){
			battery_fileds_s.push_back(battery_fileds[i]);
			if(battery_fileds[i]=="Displacement") battery_fileds_s.push_back("component_is_vector");
			else battery_fileds_s.push_back("component_is_scalar");
		  // define all at interface 
			FE_support_list_v[i+2*battery_fileds.size()]=1;
			if(battery_fileds[i]=="Lithium" or battery_fileds[i]=="Lithium_phaseField" or battery_fileds[i]=="Electrode_potential"){
				FE_support_list_v[i]=1;
				FE_support_list_v[i+battery_fileds.size()]=0;
			}
			else if(battery_fileds[i]=="Lithium_cation" or battery_fileds[i]=="Electrolyte_potential"){
				FE_support_list_v[i]=0;
				FE_support_list_v[i+battery_fileds.size()]=1;
			}
			else if(battery_fileds[i]=="Displacement"){
				FE_support_list_v[i]=1;
				FE_support_list_v[i+battery_fileds.size()]=1;
			}
			//turn on all fileds everywhere
			//FE_support_list_v[i]=1;
			//FE_support_list_v[i+battery_fileds.size()]=1;
			//FE_support_list_v[i+2*battery_fileds.size()]=1;
		}
		(*params_json)["Problem"]["primary_variables_list"]=battery_fileds_s;
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
	std::string output_path = this->output_directory+"output_w_domain-"+std::to_string(this->current_increment+this->off_output_index)+".vtk";
	this->FEMdata_out.clear_data_vectors();
	this->FEMdata_out.data_out.add_data_vector(material_id, "material");
	this->FEMdata_out.write_vtk(this->solution_prev, output_path);	
}


template <int dim>
void battery<dim>::get_residual(const typename hp::DoFHandler<dim>::active_cell_iterator &cell, const FEValues<dim>& fe_values, Table<1, Sacado::Fad::DFad<double> >& R, Table<1, Sacado::Fad::DFad<double>>& ULocal, Table<1, double >& ULocalConv)
{	
  int cell_id = cell->active_cell_index();
	double separator_line=(*params_json)["ElectroChemo"]["separator_line"];
	int orientation=(*params_json)["ElectroChemo"]["orientation"];

  if (cell->material_id()==interface_id){		
    get_residual_at_diffuse_interface(cell, fe_values, R, ULocal, ULocalConv);
  }
	else if (cell->material_id()==active_particle_id){
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
    }
	}
	else if (cell->material_id()==electrolyte_id){
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
    }
	}
	

  double flux_sign = 1;
	double fliptime=(*params_json)["ElectroChemo"]["flip_time"];
  if (this->current_time >= fliptime)
  {
    flux_sign = -1;
  }

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
			if (cell->face(faceID)->boundary_id()==1+orientation) current_IpA=-current_IpA;
		  FEFaceValues<dim> fe_face_values(fe_values.get_fe(), *(this->common_face_quadrature), update_values | update_quadrature_points | update_JxW_values);
			fe_face_values.reinit(cell,faceID);
			this->ResidualEq.residualForNeummanBC(fe_values, fe_face_values, battery_fields.active_fields_index["Electrode_potential"], R, current_IpA * flux_sign);
		}
	}
}

template <int dim>
void battery<dim>::output_results()
{
	Vector<float> material_id(this->triangulation.n_active_cells()); 
  typename hp::DoFHandler<dim>::active_cell_iterator elem = this->dof_handler.begin_active(), endc = this->dof_handler.end();             
  unsigned int j = 0;                                                                                                      
  for (;elem!=endc; ++elem){                                                                                                
		material_id(j++) = elem->material_id();
	}

	//write vtk and snapshot for solution
	if(this->save_output){ 
		std::string output_path = this->output_directory+"output-new-"+std::to_string(this->current_increment+this->off_output_index)+".vtk";
		this->FEMdata_out.clear_data_vectors();
	  this->FEMdata_out.data_out.add_data_vector(material_id, "mat_id");
	  this->FEMdata_out.data_out.add_data_vector(crack_id, "crack_id");
	  this->FEMdata_out.data_out.add_data_vector(jump_n, "jump_n");
	  this->FEMdata_out.data_out.add_data_vector(jump_m, "jump_m");
	  this->FEMdata_out.data_out.add_data_vector(jump_w, "jump_w");
	  this->FEMdata_out.data_out.add_data_vector(T_n, "T_n");
		if(this->current_increment%this->skip_output==0) this->FEMdata_out.write_vtk(this->solution_prev, output_path);	

	}
	if(this->save_snapshot){
		std::string snapshot_path = this->snapshot_directory+"snapshot-"+std::to_string(this->current_increment+this->off_output_index)+".dat";
		this->FEMdata_out.create_vector_snapshot(this->solution, snapshot_path);
	}
}



template <int dim>
void battery<dim>::run()
{
  crack_id.reinit(this->triangulation.n_active_cells());
  jump_n.reinit(this->triangulation.n_active_cells());
  jump_m.reinit(this->triangulation.n_active_cells());
  jump_w.reinit(this->triangulation.n_active_cells());
  T_n.reinit(this->triangulation.n_active_cells());
  pressure.resize(this->triangulation.n_active_cells());
  pressure_old.resize(this->triangulation.n_active_cells());
  is_new_step.resize(this->triangulation.n_active_cells());

	output_w_domain();

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

    {
      typename hp::DoFHandler<dim>::active_cell_iterator elem = this->dof_handler.begin_active(), endc = this->dof_handler.end();             
      unsigned int j = 0;                                                                                                      
      for (;elem!=endc; ++elem){                                                                                                
        int cell_id = elem->active_cell_index();
        is_new_step[cell_id] = true;
	    }
    }

    // update history variables
    for (unsigned i = 0; i < cell_SDdata.size(); ++i) {
      cell_SDdata[i].xi_conv = cell_SDdata[i].xi_old;
      cell_SDdata[i].xi_conv_c_e = cell_SDdata[i].xi_old_c_e;
      cell_SDdata[i].xi_conv_phi_s = cell_SDdata[i].xi_old_phi_s;
      cell_SDdata[i].xi_conv_phi_e = cell_SDdata[i].xi_old_phi_e;
      cell_SDdata[i].C_Li_plus_old = cell_SDdata[i].C_Li_plus_new;
      // WARNING: should not use the following, as many cell_SDdata are not initialized.
      //std::cout << "C_Li_plus_old[0]" << cell_SDdata[i].C_Li_plus_old[0] << std::endl;
      //for (int q=0; q<4; q++) cell_SDdata[i].C_Li_plus_old[q] = cell_SDdata[i].C_Li_plus_new[q];
      //std::cout << "C_Li_plus_old[3]" << cell_SDdata[i].C_Li_plus_old[3] << std::endl;
    }
		
	  t_solve = clock() - t_solve;
		this->pcout<<"It took me "<< ((float)t_solve)/CLOCKS_PER_SEC<<" seconds for this solve"<<std::endl<<std::endl;

		// std::string snapfile="snapshot_phase_2/snapshot-"+std::to_string(this->current_increment+this->off_output_index)+".dat";
		// 	  this->FEMdata_out.resume_vector_from_snapshot(this->solution,snapfile);
		// 	  this->solution_prev=this->solution;
		
     //this->FEMdata_out.clear_data_vectors();
     //Vector<double> localized_U(this->solution_prev);
     //this->FEMdata_out.data_out.add_data_vector (localized_U, computedNodalField);
     //std::string output_path = this->output_directory+"output-"+std::to_string(this->current_increment+this->off_output_index)+".vtk";
     //this->FEMdata_out.write_vtk(this->solution_prev, output_path);
    this->output_results();

    pressure_old = pressure;
	}
	this->pcout<<"Finish running!!"<<std::endl;
}

template <int dim>
void battery<dim>::identify_diffuse_interface()
{
  int primary_dof = -1;
  int opposite_flux_dof_li = -1;
  int opposite_flux_dof_potential = -1;
	if(battery_fields.active_fields_index["Diffuse_interface"]>-1) primary_dof=battery_fields.active_fields_index["Diffuse_interface"];
	if(battery_fields.active_fields_index["Lithium_cation"]>-1) opposite_flux_dof_li=battery_fields.active_fields_index["Lithium_cation"];
	if(battery_fields.active_fields_index["Electrolyte_potential"]>-1) opposite_flux_dof_potential=battery_fields.active_fields_index["Electrolyte_potential"];

	
  std::cout << "---------- primary dof for diffusive interface ------ " << primary_dof  << " opposite dof " << opposite_flux_dof_li << " "<< opposite_flux_dof_potential << std::endl;

  hp::FEValues<dim> hp_fe_values (this->fe_collection, this->q_collection, update_values | update_quadrature_points  | update_JxW_values | update_gradients);	

  Vector<double> localized_U(this->solution_prev);
  int total_cell_num = this->triangulation.n_active_cells();
  cell_SDdata.resize(total_cell_num);
  typename hp::DoFHandler<dim>::active_cell_iterator cell = this->dof_handler.begin_active(), endc=this->dof_handler.end();
  for (;cell!=endc; ++cell){
		if (cell->subdomain_id() == this->this_mpi_process){
      if (cell->material_id()==interface_id)
      {	
				hp_fe_values.reinit (cell);
	    	const FEValues<dim> &fe_values = hp_fe_values.get_present_fe_values();

	      int cell_id = cell->active_cell_index();
        //std::cout <<  " cell_id " << cell_id << " cell_SDdata.size() " << cell_SDdata.size() << " interface_id " << interface_id << std::endl;
	      cell_SDdata[cell_id].cell_id = cell_id;
	      cell_SDdata[cell_id].opposite_flux_dof_li = opposite_flux_dof_li;
	      cell_SDdata[cell_id].opposite_flux_dof_potential = opposite_flux_dof_potential;

	      const unsigned int dofs_per_cell = cell->get_fe().dofs_per_cell;

	      std::vector<unsigned int> local_dof_indices (dofs_per_cell);
	      cell->get_dof_indices (local_dof_indices);
	      std::vector<double> local_diffuse_interface;
				
	      for (unsigned int i=0; i<dofs_per_cell; ++i) {
	        unsigned int ck = fe_values.get_fe().system_to_component_index(i).first - primary_dof;
	        if (ck == 0) {
            //std::cout 
              //<< "---------- primary dof for diffusive interface ------ i = " << i 
              //<< " val = " << localized_U(local_dof_indices[i])
              //<< std::endl;
	          local_diffuse_interface.push_back(localized_U(local_dof_indices[i]));
	        }
	      }
				
        cell_SDdata[cell_id].is_interface_element = true;

        cell_SDdata[cell_id].rlocal.reinit(1);
        cell_SDdata[cell_id].rlocal(0) = 0.0;
        cell_SDdata[cell_id].xi_old.reinit(1);
        cell_SDdata[cell_id].xi_old(0) = 0.0;
        cell_SDdata[cell_id].xi_conv.reinit(1);
        cell_SDdata[cell_id].xi_conv(0) = 0.0;
        cell_SDdata[cell_id].Kcc.reinit(4,4);
        cell_SDdata[cell_id].Kcxi.reinit(4,1);
        cell_SDdata[cell_id].Kxic.reinit(1,4);
        cell_SDdata[cell_id].Kxixi_inv.reinit(1,1);

        cell_SDdata[cell_id].rlocal_c_e.reinit(1);
        cell_SDdata[cell_id].rlocal_c_e(0) = 0.0;
        cell_SDdata[cell_id].xi_old_c_e.reinit(1);
        cell_SDdata[cell_id].xi_old_c_e(0) = 0.0;
        cell_SDdata[cell_id].xi_conv_c_e.reinit(1);
        cell_SDdata[cell_id].xi_conv_c_e(0) = 0.0;
        cell_SDdata[cell_id].Kcc_c_e.reinit(4,4);
        cell_SDdata[cell_id].Kcxi_c_e.reinit(4,1);
        cell_SDdata[cell_id].Kxic_c_e.reinit(1,4);
        cell_SDdata[cell_id].Kxixi_inv_c_e.reinit(1,1);

        cell_SDdata[cell_id].C_Li_plus_old.reinit(4); // size of gps
        cell_SDdata[cell_id].C_Li_plus_new.reinit(4);

        cell_SDdata[cell_id].rlocal_phi_s.reinit(1);
        cell_SDdata[cell_id].rlocal_phi_s(0) = 0.0;
        cell_SDdata[cell_id].xi_old_phi_s.reinit(1);
        cell_SDdata[cell_id].xi_old_phi_s(0) = 0.0;
        cell_SDdata[cell_id].xi_conv_phi_s.reinit(1);
        cell_SDdata[cell_id].xi_conv_phi_s(0) = 0.0;
        cell_SDdata[cell_id].Kcc_phi_s.reinit(4,4);
        cell_SDdata[cell_id].Kcxi_phi_s.reinit(4,1);
        cell_SDdata[cell_id].Kxic_phi_s.reinit(1,4);
        cell_SDdata[cell_id].Kxixi_inv_phi_s.reinit(1,1);

        cell_SDdata[cell_id].rlocal_phi_e.reinit(1);
        cell_SDdata[cell_id].rlocal_phi_e(0) = 0.0;
        cell_SDdata[cell_id].xi_old_phi_e.reinit(1);
        cell_SDdata[cell_id].xi_old_phi_e(0) = 0.0;
        cell_SDdata[cell_id].xi_conv_phi_e.reinit(1);
        cell_SDdata[cell_id].xi_conv_phi_e(0) = 0.0;
        cell_SDdata[cell_id].Kcc_phi_e.reinit(4,4);
        cell_SDdata[cell_id].Kcxi_phi_e.reinit(4,1);
        cell_SDdata[cell_id].Kxic_phi_e.reinit(1,4);
        cell_SDdata[cell_id].Kxixi_inv_phi_e.reinit(1,1);

        cell_SDdata[cell_id].rlocal_u_sd.reinit(3);
        cell_SDdata[cell_id].rlocal_u_sd = 0.0;
        cell_SDdata[cell_id].xi_old_u_sd.reinit(3);
        cell_SDdata[cell_id].xi_old_u_sd = 0.0;
        cell_SDdata[cell_id].xi_conv_u_sd.reinit(3);
        cell_SDdata[cell_id].xi_conv_u_sd = 0.0;
        cell_SDdata[cell_id].Kuu_sd.reinit(4*dim,4*dim);
        cell_SDdata[cell_id].Kuxi_sd.reinit(4*dim,3);
        cell_SDdata[cell_id].Kxiu_sd.reinit(3,4*dim);
        cell_SDdata[cell_id].Kxixi_inv_u_sd.reinit(3,3);

        cell_SDdata[cell_id].ULocal_k.reinit(40);

        unsigned int n_q_points = fe_values.n_quadrature_points;
        for (unsigned int q = 0; q < n_q_points; ++q) {
          cell_SDdata[cell_id].area_elem += fe_values.JxW(q);
        }

        unsigned int count_larger_c = 0, count_smaller_c = 0, count_equal_c = 0;

        // get the side of the local and global node number
        for (unsigned int i=0; i<local_diffuse_interface.size(); ++i) {
          if (local_diffuse_interface[i] >= iso_value){
            cell_SDdata[cell_id].lnode_plus.push_back(i);
            cell_SDdata[cell_id].one_plus_node = cell->vertex(i);
            if (std::abs(local_diffuse_interface[i] - iso_value) < 1e-12)
            {
              count_equal_c += 1;
            }
            else
            {
              count_larger_c += 1;
            }
          };
          if (local_diffuse_interface[i] < iso_value){
            cell_SDdata[cell_id].lnode_minus.push_back(i);
            count_smaller_c += 1;
          };
          //std::cout << " --- ** -- " << i << std::endl;
        }

        std::vector<types::global_dof_index> local_face_dof_indices(this->fe_system[interface_id]->dofs_per_face);
        int count = 0;
        for (unsigned int f = 0; f < GeometryInfo<dim>::faces_per_cell; ++f) {
          cell->face(f)->get_dof_indices(local_face_dof_indices, interface_id);
          double c_1 = 0.0;
          double c_2 = 0.0;
          std::vector<double> local_local_diffuse_interface_face;
          for (unsigned int i = 0; i < local_face_dof_indices.size(); ++i) {
            const unsigned int ck = this->fe_system[interface_id]->face_system_to_component_index(i).first - primary_dof;
            if (ck == 0) local_local_diffuse_interface_face.push_back(localized_U(local_face_dof_indices[i]));
          }

          c_1 = local_local_diffuse_interface_face[0];
          c_2 = local_local_diffuse_interface_face[1];
          //std::cout << " count " << count << " c_1 " << c_1 << " c_2 " << c_2 
            //<< " larger " << count_larger_c 
            //<< " smaller " << count_smaller_c 
            //<< " equal " << count_equal_c 
            //<< std::endl;
          // slightly perturb the iso_value to avoid node cut
          if (count_larger_c == 0 and count_equal_c > 0)
          {
            if (std::abs(c_1 - iso_value) < 1e-12) c_1 = iso_value + iso_value * 0.001;
            if (std::abs(c_2 - iso_value) < 1e-12) c_2 = iso_value + iso_value * 0.001;
          }
          else if (count_smaller_c == 0 and count_equal_c > 0)
          {
            if (std::abs(c_1 - iso_value) < 1e-12) c_1 = iso_value - iso_value * 0.001;
            if (std::abs(c_2 - iso_value) < 1e-12) c_2 = iso_value - iso_value * 0.001;
          }
          else
          {
            if (count_equal_c > 0)
            {
              // all set to larger value
              if (std::abs(c_1 - iso_value) < 1e-12) c_1 = iso_value + iso_value * 0.001;
              if (std::abs(c_2 - iso_value) < 1e-12) c_2 = iso_value + iso_value * 0.001;
            }
          }

          if (c_1 == iso_value and c_2 == iso_value){
            // for the case where the edge of element is aligned with the contour. 
            cell_SDdata[cell_id].edge1_node = cell->face(f)->vertex(0) ;
            cell_SDdata[cell_id].edge2_node = cell->face(f)->vertex(1) ;
          }
          else if ((c_1 >= iso_value and c_2 < iso_value) || (c_1 <= iso_value and c_2 > iso_value)){
            // Without equal sign between c_2 vs iso_value can prevent assigning the same node with c=iso_value to both edge1_node and edge2_node
            if (count == 0){
              cell_SDdata[cell_id].edge1_node1 = cell->face(f)->vertex(0) ;
              cell_SDdata[cell_id].edge1_node2 = cell->face(f)->vertex(1) ;
              cell_SDdata[cell_id].edge1_local_s = (c_1 - iso_value)/(c_1 - c_2);
              cell_SDdata[cell_id].edge1_node =  cell_SDdata[cell_id].edge1_node1 - cell_SDdata[cell_id].edge1_local_s * (cell_SDdata[cell_id].edge1_node1 - cell_SDdata[cell_id].edge1_node2);
            }
            else if (count == 1){
              cell_SDdata[cell_id].edge2_node1 = cell->face(f)->vertex(0) ;
              cell_SDdata[cell_id].edge2_node2 = cell->face(f)->vertex(1) ;
              cell_SDdata[cell_id].edge2_local_s = (c_1 - iso_value)/(c_1 - c_2);
              cell_SDdata[cell_id].edge2_node =  cell_SDdata[cell_id].edge2_node1 - cell_SDdata[cell_id].edge2_local_s * (cell_SDdata[cell_id].edge2_node1 - cell_SDdata[cell_id].edge2_node2);
            }
            count++;
          }
        }
        double elem_length = cell_SDdata[cell_id].edge1_node.distance(cell_SDdata[cell_id].edge2_node);
        //std::cout << "elem_length " << elem_length << " " << cell_SDdata[cell_id].edge1_node << " "<< cell_SDdata[cell_id].edge2_node << std::endl;
        cell_SDdata[cell_id].interface_length = elem_length;

        double dx = cell_SDdata[cell_id].edge1_node[0] - cell_SDdata[cell_id].edge2_node[0];
        double dy = cell_SDdata[cell_id].edge1_node[1] - cell_SDdata[cell_id].edge2_node[1];
        double mid_x = 0.5*(cell_SDdata[cell_id].edge1_node[0] + cell_SDdata[cell_id].edge2_node[0]);
        double mid_y = 0.5*(cell_SDdata[cell_id].edge1_node[1] + cell_SDdata[cell_id].edge2_node[1]);
        //std::cout << "----- p1 ---- " << cell_SDdata[cell_id].edge1_node <<  " p2 " << cell_SDdata[cell_id].edge2_node << " dx " << dx <<  " dy " << dy  << std::endl;
        // two possible normal directions
        //std::cout << "----- normal ---- " << -dy <<  " " << dx << " or " << dy <<  " " << -dx  <<  " plus_node " << cell_SDdata[cell_id].one_plus_node<< std::endl;

        // correct outward normal for the plus region
        if ( (-dy * (cell_SDdata[cell_id].one_plus_node[0]-mid_x) + dx * (cell_SDdata[cell_id].one_plus_node[1] -mid_y)) < 0)
        {
            cell_SDdata[cell_id].crk_n[0] = -dy / sqrt(dy*dy+dx*dx);
            cell_SDdata[cell_id].crk_n[1] = dx / sqrt(dy*dy+dx*dx);
        }
        else
        {
            cell_SDdata[cell_id].crk_n[0] = dy / sqrt(dy*dy+dx*dx);
            cell_SDdata[cell_id].crk_n[1] = -dx / sqrt(dy*dy+dx*dx);
        }
        //std::cout << "----- final normal ---- " << cell_SDdata[cell_id].crk_n[0] <<  " " << cell_SDdata[cell_id].crk_n[1]  << " length :" << elem_length << std::endl;

        /// update the area_elem for the actual sizes
        /// should not do the following. As the crack length is smaller if the cutting region is changed. This is reflected in the local residual function.
        if ( cell_SDdata[cell_id].lnode_plus.size() == 1)
        {
          // https://www.mathopenref.com/coordtrianglearea.html
          Point<dim, double> A = cell->vertex(cell_SDdata[cell_id].lnode_plus[0]);
          Point<dim, double> B = cell_SDdata[cell_id].edge1_node;
          Point<dim, double> C = cell_SDdata[cell_id].edge2_node;
          double area = 0.5 * std::abs(A[0]*(B[1]-C[1]) + B[0]*(C[1]-A[1]) +  C[0]*(A[1]-B[1]));
          //std::cout << "----- plus node size :" << cell_SDdata[cell_id].lnode_plus.size() << " new area " << area << " old area: " << cell_SDdata[cell_id].area_elem << std::endl;
          //std::cout << " A " << A << std::endl;
          //std::cout << " B " << B << std::endl;
          //std::cout << " C " << C << std::endl;
          if (abs(area) < 1.e-12) area = 1.0e-10;
          cell_SDdata[cell_id].computed_area = area;
        }

        if ( cell_SDdata[cell_id].lnode_plus.size() == 3)
        {
          // https://www.mathopenref.com/coordtrianglearea.html
          Point<dim, double> A = cell->vertex(cell_SDdata[cell_id].lnode_minus[0]);
          Point<dim, double> B = cell_SDdata[cell_id].edge1_node;
          Point<dim, double> C = cell_SDdata[cell_id].edge2_node;
          double area = 0.5 * std::abs(A[0]*(B[1]-C[1]) + B[0]*(C[1]-A[1]) +  C[0]*(A[1]-B[1]));
          //std::cout << "----- plus node size :" << cell_SDdata[cell_id].lnode_plus.size() << " new area " << cell_SDdata[cell_id].area_elem - area << " old area: " << cell_SDdata[cell_id].area_elem << std::endl;
          if (abs(area) < 1.e-12) area = 1.0e-10;
          cell_SDdata[cell_id].computed_area = cell_SDdata[cell_id].area_elem - area;
        }

        if ( cell_SDdata[cell_id].lnode_plus.size() == 2)
        {
          // https://www.mathopenref.com/coordtrianglearea.html
          double area_1 =0, area_2 =0, area_3 =0, area_4 = 0;
          {
            Point<dim, double> A = cell->vertex(cell_SDdata[cell_id].lnode_plus[0]);
            Point<dim, double> B = cell_SDdata[cell_id].edge1_node;
            Point<dim, double> C = cell_SDdata[cell_id].edge2_node;
            area_1 = 0.5 * std::abs(A[0]*(B[1]-C[1]) + B[0]*(C[1]-A[1]) +  C[0]*(A[1]-B[1]));
            if (abs(area_1) < 1.e-12) area_1 = 1.0e-10;
          }

          {
            Point<dim, double> A = cell->vertex(cell_SDdata[cell_id].lnode_plus[1]);
            Point<dim, double> B = cell_SDdata[cell_id].edge1_node;
            Point<dim, double> C = cell_SDdata[cell_id].edge2_node;
            area_2 = 0.5 * std::abs(A[0]*(B[1]-C[1]) + B[0]*(C[1]-A[1]) +  C[0]*(A[1]-B[1]));
            if (abs(area_2) < 1.e-12) area_2 = 1.0e-10;
          }

          {
            Point<dim, double> A = cell->vertex(cell_SDdata[cell_id].lnode_plus[1]);
            Point<dim, double> B = cell->vertex(cell_SDdata[cell_id].lnode_plus[0]);
            Point<dim, double> C = cell_SDdata[cell_id].edge2_node;
            area_3 = 0.5 * std::abs(A[0]*(B[1]-C[1]) + B[0]*(C[1]-A[1]) +  C[0]*(A[1]-B[1]));
            if (abs(area_3) < 1.e-12) area_3 = 1.0e-10;
          }
          {
            Point<dim, double> A = cell->vertex(cell_SDdata[cell_id].lnode_plus[1]);
            Point<dim, double> B = cell->vertex(cell_SDdata[cell_id].lnode_plus[0]);
            Point<dim, double> C = cell_SDdata[cell_id].edge1_node;
            area_4 = 0.5 * std::abs(A[0]*(B[1]-C[1]) + B[0]*(C[1]-A[1]) +  C[0]*(A[1]-B[1]));
            if (abs(area_4) < 1.e-12) area_4 = 1.0e-10;
          }
          double area = 0.5 * (area_1 + area_2 + area_3 + area_4 );
          //std::cout << "----- plus node size :" << cell_SDdata[cell_id].lnode_plus.size() << " new area " << area << " old area: " << cell_SDdata[cell_id].area_elem << std::endl;
          cell_SDdata[cell_id].computed_area = area;
        }

        //
        //

        Triangulation<1> triangulation_1d;
        GridGenerator::hyper_cube	(	triangulation_1d, 0.,  elem_length);
        DoFHandler<1>      dof_handler(triangulation_1d);
        int problem_dof = 1;
        int poly_order = 1;
        int quad_order = 2;
        FESystem<1> fe (FE_Q<1>(poly_order), problem_dof);
        dof_handler.distribute_dofs (fe);

        QGauss<1>  quadrature_formula(quad_order);
        FEValues<1> fe_values_1d (fe, quadrature_formula,
                                 update_values   | update_gradients |
                                 update_quadrature_points | update_JxW_values);

        typename DoFHandler<1>::active_cell_iterator cell_1d = dof_handler.begin_active(),
                                                       endc_1d = dof_handler.end();

        double vol = 0.0;
        cell_SDdata[cell_id].shape_value_1d.reinit(2,quadrature_formula.size());
        cell_SDdata[cell_id].jxw_1d.reinit(quadrature_formula.size());

        for (; cell_1d!=endc_1d; ++cell_1d)
        {
            fe_values_1d.reinit (cell_1d);
            for (unsigned q=0; q< quadrature_formula.size(); ++q)
            {
              //std::cout << " q " << q << std::endl;
              vol += fe_values_1d.JxW(q);
              for (unsigned int i = 0; i < 2; ++i) {
                cell_SDdata[cell_id].shape_value_1d(i, q) = fe_values_1d.shape_value(i, q);
                cell_SDdata[cell_id].jxw_1d(q) = fe_values_1d.JxW(q);
                //r_local[i] += fe_values_1d.shape_value(i, q) * dRc * fe_values_1d.JxW(q);
              } // q_point
            }
        } // cell_1d
        //std::cout <<  " total elem # = " << triangulation_1d.n_active_cells() << " length: " << elem_length << " vol " << vol<< std::endl;
        //for (auto p0: cell_SDdata[cell_id].lnode_plus) std::cout << " plus node: " << p0 << std::endl;
        //for (auto p0: cell_SDdata[cell_id].lnode_minus) std::cout << " minus node: " << p0 << std::endl;
        //std::cout << " crk_n: " << cell_SDdata[cell_id].crk_n[0] << "\t" << cell_SDdata[cell_id].crk_n[1] << std::endl;
      } // interface id
		}		// this process
	} // cell

  {
    typename hp::DoFHandler<dim>::active_cell_iterator cell = this->dof_handler.begin_active(), endc=this->dof_handler.end();
    for (;cell!=endc; ++cell){
	  	if (cell->subdomain_id() == this->this_mpi_process){
        int cell_id = cell->active_cell_index();
        if  (cell_SDdata[cell_id].is_interface_element)
        {
          Point<dim> cell_center=cell->center();
          //std::cout 
            //<< " cell_id = " << cell_id << " center = " 
            //<< cell_center << " has interface!  R = " 
            //<< std::sqrt((cell_center(0)-2)*(cell_center(0)-2) + (cell_center(1)-2)*(cell_center(1)-2))
            //<< std::endl;
        }
      }
    }
  }
}

template class battery<1>;
template class battery<2>;
template class battery<3>;
