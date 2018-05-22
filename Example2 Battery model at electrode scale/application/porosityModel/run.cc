#include"initBoundValProbs.h"
#include"nodalField.h"

template <int dim>
void initBoundValProbs<dim>::run(std::vector<std::vector<std::string> > _primary_variables, std::vector<std::vector<int> > _FE_support)
{
	primary_variables=_primary_variables;
	FE_support=_FE_support;
	
	Residual<Sacado::Fad::DFad<double>,dim> _ResidualEq;
	ResidualEq=&_ResidualEq;
	
	electricChemoFormula->setParametersFromHandler();
	
	params->enter_subsection("Problem");
	bool printParameter=params->get_bool("print_parameter");
	std::string output_directory=params->get("output_directory");
	std::string snapshot_directory=params->get("snapshot_directory");
	bool resuming_from_snapshot=params->get_bool("resuming_from_snapshot");
	std::string snapfile=params->get("snapshot_file");
	
	const QGauss<dim> _volume_quadrature(params->get_integer("volume_quadrature"));
	volume_quadrature=&_volume_quadrature;
	const QGauss<dim-1> _face_quadrature(params->get_integer("face_quadrature"));
	common_face_quadrature=&_face_quadrature;
	
	int first_domain_id=params->get_integer("first_domain_id");
	separator_id=params->get_integer("separator_id");
	electrode_id=params->get_integer("electrode_id");
	separator_fe=params->get_integer("separator_fe");
	electrode_fe=params->get_integer("electrode_fe");
	
	double dt_0=params->get_double("dt");
	double total_time=params->get_double("totalTime");
	bool step_load=params->get_bool("step_load");
	double IpA_origin=params->get_double("IpA");
	params->leave_subsection();	
	
	params->enter_subsection("Geometry");
	electrode_Y1=params->get_double("electrode_Y1");
	electrode_Y2=params->get_double("electrode_Y2");
  params->leave_subsection();
	
	if(printParameter) {
		if(this_mpi_process == 0) params->print_parameters (std::cout, ParameterHandler::Text);
	}
	
	
	FEMdata<dim,PETScWrappers::MPI::Vector> _FEMdata_out(this->dof_handler);
	FEMdata_out=&_FEMdata_out;
	FEMdata_out->set_output_name(primary_variables);

	this->setup_FeSystem(fe_system, fe_collection, q_collection, primary_variables_dof,primary_variables,FE_support,*volume_quadrature);
	
	//set dof_id
	u_dof=primary_variables_dof[0];
	c_li_plus_dof=primary_variables_dof[1];
	phi_e_dof=primary_variables_dof[2];
	c_li_dof=primary_variables_dof[3];
	phi_s_dof=primary_variables_dof[4];
	T_dof=primary_variables_dof[5];
		
	current_dt=dt_0;
	make_grid();
	setMultDomain();
	//refine_grid();
  mark_boundary();
  this->set_active_fe_indices (FE_support, hpFEM<dim>::dof_handler,first_domain_id);
  setup_system();	
	setup_constraints();
	//output initial

	if(!resuming_from_snapshot) {
	  apply_initial_condition();
	  std::string output_path = output_directory+"output-0.vtk";
	  FEMdata_out->write_vtk(solution_prev, output_path);
	  solution=solution_prev;
	  std::string snapshot_path = snapshot_directory+"snapshot-"+std::to_string(0)+".txt";
	  FEMdata_out->create_vector_snapshot(solution_prev, snapshot_path);
	}
	else {
	  pcout<<"resuming from snapshot"<<std::endl;
	  FEMdata_out->resume_vector_from_snapshot(solution,snapfile);
	  std::string output_path = output_directory+"output-resume.vtk";
          FEMdata_out->write_vtk(solution, output_path);
	  solution_prev=solution;
	}
	pcout<<"initial condition applied"<<std::endl;
	
	//parameters
	double small_dt=0.001;
	current_increment=0;
	current_time=0;
	double initial_time=0;	
	
	if(step_load){
		current_dt=small_dt;	
	}
	
	nodalField<dim> computedNodalField(*params);
	std::vector<std::vector<std::string> > computed_primary_variables={{"jn", "component_is_scalar"}};
	//computed_primary_variables[0].push_back("jn"); computed_primary_variables[0].push_back("component_is_scalar");
	computedNodalField.setupComputedField(computed_primary_variables);
	computedNodalField.primary_variables_dof=primary_variables_dof;
	
	
	bool constrain=true;
	current_IpA=IpA_origin;
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
				current_IpA=1*IpA_origin;
    	} 
    } 
		clock_t t_solve;	
  	t_solve = clock();
		this->nonlinearSolve(solution);
	  t_solve = clock() - t_solve;
		PetscPrintf(this->mpi_communicator,"************");
		PetscPrintf(this->mpi_communicator,"current increment=%d, current time= %f,  It took me %f seconds for one solve.\n ",current_increment, current_time,((float)t_solve)/CLOCKS_PER_SEC);
	
		//update
		solution_prev=solution;
		//write vtk and snapshot for solution
		FEMdata_out->clear_data_vectors();
		Vector<double> localized_U(solution);
		
		FEMdata_out->data_out.add_data_vector (localized_U, computedNodalField);
		std::string output_path = output_directory+"output-"+std::to_string(current_increment)+".vtk";
	  FEMdata_out->write_vtk(solution_prev, output_path);	
		
		//std::string snapshot_path = snapshot_directory+"snapshot-"+std::to_string(current_increment)+".txt";
		//FEMdata_out->create_vector_snapshot(solution, snapshot_path);
		
	}
}


template class initBoundValProbs<1>;
template class initBoundValProbs<2>;
template class initBoundValProbs<3>;
