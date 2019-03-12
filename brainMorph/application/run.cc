#include"initBoundValProbs.h"
#include <deal.II/dofs/dof_tools.h>
template <int dim>
void initBoundValProbs<dim>::run()
{	
	Residual<Sacado::Fad::DFad<double>,dim> _ResidualEq;
	ResidualEq=&_ResidualEq;
	
	FEMdata<dim,PETScWrappers::MPI::Vector> _FEMdata_out(this->dof_handler);
	FEMdata_out=&_FEMdata_out;
	FEMdata_out->set_output_name(primary_variables);
	
	//
	declare_parameters();
	params->read_input ("parameters.prm");	
	params->enter_subsection("Problem");
	bool printParameter=params->get_bool("print_parameter");
	bool with_additional_data=params->get_bool("with_additional_data");
	bool resuming_additional_data_from_snapshot=params->get_bool("resuming_additional_data_from_snapshot");
	
	std::string output_directory=params->get("output_directory");
	std::string snapshot_directory=params->get("snapshot_directory");
	std::string snapfile=params->get("snapshot_file");
	
	std::string additional_data_snapshot=params->get("additional_data_snapshot");
	
	
	current_dt=params->get_double("dt");
	current_increment=params->get_integer("current_increment");
	current_time=params->get_double("current_time");
	bool resuming_from_snapshot=params->get_bool("resuming_from_snapshot");
	double total_time=params->get_double("totalTime");
	
	const QGauss<dim> _volume_quadrature(params->get_integer("volume_quadrature"));
	volume_quadrature=&_volume_quadrature;
	const QGauss<dim-1> _face_quadrature(params->get_integer("face_quadrature"));
	common_face_quadrature=&_face_quadrature;
	
	int first_domain_id=params->get_integer("first_domain_id");
	Cortex_id=params->get_integer("Cortex_id");
  Subcortex_id=params->get_integer("Subcortex_id");
	Ventricle_id=params->get_integer("Ventricle_id");
	
	Cortex_fe=params->get_integer("Cortex_fe");
  Subcortex_fe=params->get_integer("Subcortex_fe");
	Ventricle_fe=params->get_integer("Ventricle_fe");
	params->leave_subsection();	
	
	params->enter_subsection("Concentration");
	std::string advection_type=params->get("advection_type");
	params->leave_subsection();	
	if(printParameter) {
		if(this_mpi_process == 0) params->print_parameters (std::cout, ParameterHandler::Text);
	}
	this->setup_FeSystem(fe_system, fe_collection, q_collection, primary_variables_dof,primary_variables,FE_support,*volume_quadrature);
	u_dof=primary_variables_dof[0];
	c1_dof=primary_variables_dof[1];
	c2_dof=primary_variables_dof[2];
		
	make_grid();
	//setMultDomain();
	//refine_grid();
  mark_boundary();
  this->set_active_fe_indices (FE_support, hpFEM<dim>::dof_handler,first_domain_id);
  setup_system();
	setup_constraints();
	apply_initial_condition();
	if(!resuming_from_snapshot) {
	  // apply_initial_condition();
	  std::string output_path = output_directory+"output-0.vtk";
	  FEMdata_out->write_vtk(solution_0, output_path);
	  solution_prev=solution_0;
	  solution=solution_0;
	}
	else {
	  pcout<<"resuming from snapshot"<<std::endl;
	  FEMdata_out->resume_vector_from_snapshot(solution,snapfile);
	  std::string output_path = output_directory+"output-resume.vtk";
          FEMdata_out->write_vtk(solution, output_path);
	  solution_prev=solution;
	}
	pcout<<"initial condition applied"<<std::endl;
	

	if(with_additional_data){
		dof_handler_add=new hp::DoFHandler<dim>(this->triangulation);

		FEMdata<dim,PETScWrappers::MPI::Vector> FEMdata_out_add(*dof_handler_add);
		FEMdata_out_add.set_output_name(variables_add); 
		
		setup_addtionalSystem();
		
		if(std::strcmp(advection_type.c_str(),"fromFile")==0){
		  if(resuming_additional_data_from_snapshot){FEMdata_out_add.resume_vector_from_snapshot(additional_data,additional_data_snapshot);}
			else {
				apply_advection_direction();
				std::string snapshot_path = "advection.dat";
				FEMdata_out_add.create_vector_snapshot(additional_data, snapshot_path);
			}
		}
 
		Vector<float> material_id(this->triangulation.n_active_cells()); 
  	typename hp::DoFHandler<dim>::active_cell_iterator elem = dof_handler_add->begin_active(), endc = dof_handler_add->end();             
 	 	unsigned int j = 0;                                                                                                      
  	for (;elem!=endc; ++elem){                                                                                                
			material_id(j++) = elem->material_id();
		}
	  FEMdata_out_add.data_out.add_data_vector(material_id, "material");

		std::string output_path = output_directory+"additional_info.vtk";
		FEMdata_out_add.write_vtk(additional_data, output_path);
	}

	//parameters
	dataStack.resize	(3,0);
	double initial_time=0;
	//initialize solution, solution_prev
	//solution_prev=solution_0;
	//solution=solution_0;
	double timeCheck1 = 14999, timeCheck2 = 15000;
	for (; current_time<=total_time; current_time+=current_dt){

    current_increment++;
		if(current_time>=timeCheck1 and current_time<timeCheck2 ) current_dt=1;
		else if(current_increment>=timeCheck2) current_dt=0.1;
		//else if(current_increment>=30 and current_increment<40 ) current_dt=2;
		
		clock_t t_solve;	
  	t_solve = clock();
		this->nonlinearSolve(solution);
		
	  t_solve = clock() - t_solve;
		PetscPrintf(this->mpi_communicator,"************");
		PetscPrintf(this->mpi_communicator,"current increment=%d, current time= %f,  It took me %f seconds for one solve.\n ",current_increment, current_time,((float)t_solve)/CLOCKS_PER_SEC);
		
		//update
		solution_prev=solution;
		
		//write vtk and snapshot for solution
		std::string output_path = output_directory+"output-"+std::to_string(current_increment)+".vtk";
		
		FEMdata_out->clear_data_vectors();
	  FEMdata_out->write_vtk(solution_prev, output_path);		
		
		std::string snapshot_path = snapshot_directory+"snapshot-"+std::to_string(current_increment)+".dat";
		FEMdata_out->create_vector_snapshot(solution, snapshot_path);
		//write surface area
		if (this_mpi_process == 0 ){
		  std::ofstream myfile;
		  myfile.open ( (output_directory+"suface_area.dat").c_str(),std::ios_base::app);
			if(!myfile.is_open()) {std::cout<<"file failed to open!"; exit(1);}
			for(unsigned int i=0;i<dataStack.size();i++) {
			  myfile <<dataStack[i]<<" ";
			}
			myfile <<"\n";
		}
	}
}


template class initBoundValProbs<1>;
template class initBoundValProbs<2>;
template class initBoundValProbs<3>;
