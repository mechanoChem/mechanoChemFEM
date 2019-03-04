/*
zhenlin wang 2019
*/

#include"../../include/initBoundValProbs.h"
template <int dim>
void initBoundValProbs<dim>::run()
{		
	pre_run();
		
	make_grid();
	refine_grid();
	setMultDomain();
  mark_boundary();
  setup_system();
	setup_constraints();
	if(!resuming_from_snapshot) {
	  apply_initial_condition();
	  std::string output_path = output_directory+"output-0.vtk";
	  FEMdata_out.write_vtk(solution, output_path);
	  solution_prev=solution;
	  if(save_snapshot){
			std::string snapshot_path = snapshot_directory+"snapshot-"+std::to_string(0)+".dat";
	  	FEMdata_out.create_vector_snapshot(solution, snapshot_path);
		}
	}
	else {
	  pcout<<"resuming from snapshot"<<std::endl;
	  FEMdata_out.resume_vector_from_snapshot(solution,snapfile);
	  std::string output_path = output_directory+"output-resume.vtk";
    FEMdata_out.write_vtk(solution, output_path);
	  solution_prev=solution;
	}
	pcout<<"initial condition applied"<<std::endl;
	PetscPrintf(this->mpi_communicator,"running....\n\n");
	clock_t t_solve;	
	t_solve = clock();
  for (; current_time<=total_time; current_time+=current_dt){
    current_increment++;
		PetscPrintf(this->mpi_communicator,"************");
		PetscPrintf(this->mpi_communicator,"current increment=%d, current time= %f",current_increment, current_time);
		PetscPrintf(this->mpi_communicator,"************\n");
		this->nonlinearSolve(solution);
		
	  t_solve = clock() - t_solve;
		PetscPrintf(this->mpi_communicator,"It took me %f seconds for this solve.\n ",((float)t_solve)/CLOCKS_PER_SEC);
		PetscPrintf(this->mpi_communicator,"\n\n");
		//update
		solution_prev=solution;
		
		//write vtk and snapshot for solution
		std::string output_path = output_directory+"output-"+std::to_string(current_increment+off_output_index)+".vtk";
		
		FEMdata_out.clear_data_vectors();
		if(current_increment%skip_output==0) FEMdata_out.write_vtk(solution_prev, output_path);	
		if(save_snapshot){
			std::string snapshot_path = snapshot_directory+"snapshot-"+std::to_string(current_increment+off_output_index)+".dat";
			FEMdata_out.create_vector_snapshot(solution, snapshot_path);
		}
	}
	PetscPrintf(this->mpi_communicator,"Finish running!!\n");
}


template class initBoundValProbs<1>;
template class initBoundValProbs<2>;
template class initBoundValProbs<3>;
