/*
zhenlin wang 2019
*/

#include"../../include/mechanoChemFEM.h"
template <int dim>
void mechanoChemFEM<dim>::run()
{		
	apply_boundary_condition();
  pcout << "   Number of active cells:       " << hpFEM<dim>::triangulation.n_active_cells() << std::endl;
  pcout << "   Number of degrees of freedom: " << hpFEM<dim>::dof_handler.n_dofs() << std::endl; 
	
	if(!resuming_from_snapshot) {
		apply_initial_condition();		
	  std::string output_path = output_directory+"output-0.vtk";
	  FEMdata_out.write_vtk(solution, output_path);
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
		solve_ibvp();
		
	  t_solve = clock() - t_solve;
		PetscPrintf(this->mpi_communicator,"It took me %f seconds for this solve.\n ",((float)t_solve)/CLOCKS_PER_SEC);
		PetscPrintf(this->mpi_communicator,"\n\n");

		//call post TS update function (e.g. for adaptive mesh refinement)
		update_post_TS ();
		
		output_results();
	}
	PetscPrintf(this->mpi_communicator,"Finish running!!\n");
}


template class mechanoChemFEM<1>;
template class mechanoChemFEM<2>;
template class mechanoChemFEM<3>;
