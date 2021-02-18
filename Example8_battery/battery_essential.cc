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
	battery_fields.update_fields(cell, fe_values, ULocal, ULocalConv);

  // update reaction rate at the interface 
	double tem=(*params_json)["ElectroChemo"]["jn_react"];
	double fliptime=(*params_json)["ElectroChemo"]["flip_time"];
	int Li_index=battery_fields.active_fields_index["Lithium"];
	int Li_plus_index=battery_fields.active_fields_index["Lithium_cation"];
  if (cell->material_id()==interface_id){		
		if(battery_fields.active_fields_index["Diffuse_interface"]>-1) diffuse_interface.r_get_residual(fe_values, R, ULocal, ULocalConv);
		if(battery_fields.active_fields_index["Lithium"]>-1) lithium.r_get_residual(fe_values, R, ULocal, ULocalConv);
		if(battery_fields.active_fields_index["Lithium_phaseField"]>-1) lithium_mu.r_get_residual(fe_values, R, ULocal, ULocalConv);
		if(battery_fields.active_fields_index["Electrode_potential"]>-1) phi_s.r_get_residual(fe_values, R, ULocal, ULocalConv);
	  if(battery_fields.active_fields_index["Lithium_cation"]>-1) lithium_cation.r_get_residual(fe_values, R, ULocal, ULocalConv);
		if(battery_fields.active_fields_index["Electrolyte_potential"]>-1) phi_e.r_get_residual(fe_values, R, ULocal, ULocalConv);
  }
	else if (cell->material_id()==active_particle_id){
		if(battery_fields.active_fields_index["Lithium"]>-1) lithium.r_get_residual(fe_values, R, ULocal, ULocalConv);
		if(battery_fields.active_fields_index["Lithium_phaseField"]>-1) lithium_mu.r_get_residual(fe_values, R, ULocal, ULocalConv);
		if(battery_fields.active_fields_index["Electrode_potential"]>-1) phi_s.r_get_residual(fe_values, R, ULocal, ULocalConv);
	}
	else if (cell->material_id()==electrolyte_id){
	  if(battery_fields.active_fields_index["Lithium_cation"]>-1) lithium_cation.r_get_residual(fe_values, R, ULocal, ULocalConv);
		if(battery_fields.active_fields_index["Electrolyte_potential"]>-1) phi_e.r_get_residual(fe_values, R, ULocal, ULocalConv);
	}
	if(battery_fields.active_fields_index["Displacement"]>-1) displacement.r_get_residual(fe_values, R, ULocal, ULocalConv);
	
	apply_Neumann_boundary_condition();
}

template <int dim>
void battery<dim>::run()
{
  //identify_diffuse_interface();
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
		
	  t_solve = clock() - t_solve;
		this->pcout<<"It took me"<< ((float)t_solve)/CLOCKS_PER_SEC<<" seconds for this solve"<<std::endl<<std::endl;
		
		// this->FEMdata_out.clear_data_vectors();
		// Vector<double> localized_U(this->solution_prev);
		// this->FEMdata_out.data_out.add_data_vector (localized_U, computedNodalField);
		// std::string output_path = this->output_directory+"output-"+std::to_string(this->current_increment+this->off_output_index)+".vtk";
		// this->FEMdata_out.write_vtk(this->solution_prev, output_path);
    this->output_results();
	}
	this->pcout<<"Finish running!!"<<std::endl;
}

template <int dim>
void battery<dim>::identify_diffuse_interface()
{
  int primary_dof = -1;
  int opposite_flux_dof = -1;
	if(battery_fields.active_fields_index["Diffuse_interface"]>-1) primary_dof=battery_fields.active_fields_index["Diffuse_interface"];
	if(battery_fields.active_fields_index["Lithium_cation"]>-1) opposite_flux_dof=battery_fields.active_fields_index["Lithium_cation"];
	
  std::cout << "---------- primary dof for diffusive interface ------ " << primary_dof  << " opposite dof " << opposite_flux_dof << std::endl;

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
	      cell_SDdata[cell_id].cell_id = cell_id;
	      cell_SDdata[cell_id].opposite_flux_dof = opposite_flux_dof;

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

        // get the side of the local and global node number
        for (unsigned int i=0; i<local_diffuse_interface.size(); ++i) {
          if (local_diffuse_interface[i] >= iso_value){
            cell_SDdata[cell_id].lnode_plus.push_back(i);
          };
          if (local_diffuse_interface[i] < iso_value){
            cell_SDdata[cell_id].lnode_minus.push_back(i);
          };
        }

        std::vector<types::global_dof_index> local_face_dof_indices(this->fe_system[interface_id]->dofs_per_face);
        int count = 0;
        for (unsigned int f = 0; f < GeometryInfo<dim>::faces_per_cell; ++f) {
          cell->face(f)->get_dof_indices(local_face_dof_indices, 0);
          double c_1 = 0.0;
          double c_2 = 0.0;
          std::vector<double> local_local_diffuse_interface_face;
          for (unsigned int i = 0; i < local_face_dof_indices.size(); ++i) {
            const unsigned int ck = this->fe_system[interface_id]->face_system_to_component_index(i).first - primary_dof;
            if (ck == 0)
              local_local_diffuse_interface_face.push_back(localized_U(local_face_dof_indices[i]));
          }

          c_1 = local_local_diffuse_interface_face[0];
          c_2 = local_local_diffuse_interface_face[1];

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
            for (unsigned int i = 0; i < 2; ++i) {
                for (unsigned q=0; q< quadrature_formula.size(); ++q)
                {
                    //std::cout << " q " << q << std::endl;
                    vol += fe_values_1d.JxW(q);

                    cell_SDdata[cell_id].shape_value_1d(i, q) = fe_values_1d.shape_value(i, q);
                    cell_SDdata[cell_id].jxw_1d(q) = fe_values_1d.JxW(q);
                    //r_local[i] += fe_values_1d.shape_value(i, q) * dRc * fe_values_1d.JxW(q);
                } // q_point
            }
        }
        //std::cout <<  " total elem # = " << triangulation_1d.n_active_cells() << " length: " << elem_length << " vol " << vol<< std::endl;
      }
			//}
		}		
	}

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
