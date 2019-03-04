#include"initBoundValProbs.h"
template <int dim>
void initBoundValProbs<dim>::run()
{
	Residual<Sacado::Fad::DFad<double>,dim> _ResidualEq;
	ResidualEq=&_ResidualEq;
	
	ElectricChemo<Sacado::Fad::DFad<double>,dim> _electricChemoFormula(*params);
	electricChemoFormula=&_electricChemoFormula;
	declare_parameters();
	
	params->read_input ("../parameters.prm");
	electricChemoFormula->setParametersFromHandler();
	
	params->enter_subsection("Problem");
	bool printParameter=params->get_bool("print_parameter");
	std::string output_directory=params->get("output_directory");
	std::string snapshot_directory=params->get("snapshot_directory");
	const QGauss<dim> _volume_quadrature(params->get_integer("volume_quadrature"));
	volume_quadrature=&_volume_quadrature;
	const QGauss<dim-1> _face_quadrature(params->get_integer("face_quadrature"));
	common_face_quadrature=&_face_quadrature;
	
	active_material_id=params->get_integer("active_material_id");
	electrolyte_id=params->get_integer("electrolyte_id");
	current_collector_id=params->get_integer("current_collector_id");
	binder_id=params->get_integer("binder_id");
	solid_id=params->get_integer("solid_id");
	
	active_material_fe=params->get_integer("active_material_fe");
	electrolyte_fe=params->get_integer("electrolyte_fe");
	current_collector_fe=params->get_integer("current_collector_fe");
	binder_fe=2params->get_integer("binder_fe");
	solid_fe=3params->get_integer("solid_fe");
	
	
	
	double dt_0=params->get_double("dt");
	double total_time=params->get_double("totalTime");
	bool step_load=params->get_bool("step_load");
	double IpA_origin=params->get_double("IpA");
	params->leave_subsection();	
	
	params->enter_subsection("Geometry");
	int particle_num=params->get_integer("particle_number");
	electrode_Y1=params->get_double("electrode_Y1");
	electrode_Y2=params->get_double("electrode_Y2");
  params->leave_subsection();
	
	RVEdata.resize	(particle_num*5,0);
	if(printParameter) {
		if(this_mpi_process == 0) params->print_parameters (std::cout, ParameterHandler::Text);
	}
	
	FEMdata<dim,PETScWrappers::MPI::Vector> _FEMdata_out(this->dof_handler);
	FEMdata_out=&_FEMdata_out;
	FEMdata_out->set_output_name(primary_variables);

	this->setup_FeSystem(fe_system, fe_collection, q_collection, primary_variables_dof,primary_variables,FE_support,*volume_quadrature);
	
	//set dof_id
	u_dof=primary_variables_dof[0];
  v_dof=primary_variables_dof[1];
  p_dof=primary_variables_dof[2];
	c_li_plus_dof=primary_variables_dof[3];
	phi_e_dof=primary_variables_dof[4];
	c_li_dof=primary_variables_dof[5];
	phi_s_dof=primary_variables_dof[6];
	T_dof=primary_variables_dof[7];
	u_m_dof=primary_variables_dof[8];
	
	
	current_dt=dt_0;
	make_grid();
	//setMultDomain();
	//refine_grid();
  mark_boundary();
  this->set_active_fe_indices (FE_support, hpFEM<dim>::dof_handler,1);
  setup_system();	
	setup_constraints();
	apply_initial_condition();
	//output initial
	Vector<float> material_id(this->triangulation.n_active_cells()); 
  typename hp::DoFHandler<dim>::active_cell_iterator elem = this->dof_handler.begin_active(), endc = this->dof_handler.end();             
  unsigned int j = 0;                                                                                                      
  for (;elem!=endc; ++elem){                                                                                                
		material_id(j++) = elem->material_id();
	}

	std::string output_path = output_directory+"output-0.vtk";
	
  FEMdata_out->data_out.add_data_vector(material_id, "material");
  FEMdata_out->write_vtk(solution_prev, output_path);		
	//parameters
	double small_dt=0.001;
	current_increment=0;
	current_time=0;
	double initial_time=0;	
       	setup_constraints();
	if(step_load){
		current_dt=small_dt;	
		setup_constraints();
	}
	bool constrain=true;
  for (current_time=initial_time; current_time<=total_time; current_time+=current_dt){
    current_increment++;
		//step loading
    if(step_load){
      if(current_increment<=0){current_IpA=IpA_origin/20; }
      else if(current_increment<=1){current_IpA=1*IpA_origin/10; }
      else if(current_increment<=2){current_IpA=2*IpA_origin/10; }
      else if(current_increment<=3){current_IpA=4*IpA_origin/10; }
      else if(current_increment<=4){current_IpA=8*IpA_origin/10; }
    	else  {	
	  		current_dt=dt_0;
	  		if (constrain){ setup_constraints(); constrain=false;}
				current_IpA=1*IpA_origin;
    	} 
    } 
		
		clock_t t_solve;	
  	t_solve = clock();
		
		dU=0;
		this->nonlinearSolve(solution,dU);
		
	  t_solve = clock() - t_solve;
		PetscPrintf(this->mpi_communicator,"************");
		PetscPrintf(this->mpi_communicator,"current increment=%d, current time= %f,  It took me %f seconds for one solve.\n ",current_increment, current_time,((float)t_solve)/CLOCKS_PER_SEC);
		
		
		//update
		solution_prev=solution;
		
		
		//write vtk and snapshot for solution
		
		//putput v
		if (this_mpi_process == 0 ){
		  std::ofstream myfile;
		  myfile.open ( (output_directory+"RVEdata.txt").c_str(),std::ios_base::app);
			if(!myfile.is_open()) {std::cout<<"file failed to open!"; exit(1);}
			for(unsigned int i=0;i<particle_num*5;i++) {
			  myfile <<RVEdata[i]<<" ";
			}
			myfile <<"\n";
		}
		
		std::string output_path = output_directory+"output-"+std::to_string(current_increment)+".vtk";
		FEMdata_out->clear_data_vectors();
		FEMdata_out->data_out.add_data_vector(material_id, "material");
	  FEMdata_out->write_vtk(solution_prev, output_path);		
		
		//std::string snapshot_path = snapshot_directory+"snapshot-"+std::to_string(current_increment)+".txt";
		//FEMdata_out->create_vector_snapshot(solution, snapshot_path);
		
		//reset velocity and pressure to zero. 
	  typename hp::DoFHandler<dim>::active_cell_iterator cell = hpFEM<dim>::dof_handler.begin_active(), endc=hpFEM<dim>::dof_handler.end();
	  for (;cell!=endc; ++cell){
			if (cell->subdomain_id() == this_mpi_process){
	    	hp::FEValues<dim> hp_fe_values (fe_collection, q_collection, update_values | update_quadrature_points  | update_JxW_values | update_gradients);
	    	hp_fe_values.reinit (cell);
	    	const FEValues<dim> &fe_values = hp_fe_values.get_present_fe_values();
	    	const unsigned int dofs_per_cell = cell->get_fe().dofs_per_cell;
				std::vector<unsigned int> local_dof_indices (dofs_per_cell);
	    	cell->get_dof_indices (local_dof_indices);
				
				if(cell->material_id()==electrolyte_id ){
					for (unsigned int i=0; i<dofs_per_cell; ++i) {
						const unsigned int ck = fe_values.get_fe().system_to_component_index(i).first-v_dof;
						if(ck>=0 and ck<dim){
							solution_prev(local_dof_indices[i])=0;
							solution_prev(local_dof_indices[i])=0;
						}
					}
				}
			}
	  }

		solution_prev.compress(VectorOperation::insert);
		solution_prev.compress(VectorOperation::insert);
		
	}
}


template class initBoundValProbs<1>;
template class initBoundValProbs<2>;
template class initBoundValProbs<3>;
