/*
zhenlin wang 2019
*/

#include"../../include/initBoundValProbs.h"

template <int dim>
void initBoundValProbs<dim>::apply_initial_condition()
{ 
	pcout << "applying initial condition\n";
	if(!resuming_from_snapshot) {
		int totalDOF=this->totalDOF(primary_variables);	
	  VectorTools::interpolate(this->dof_handler, InitialConditions<dim>(totalDOF, *params), solution_prev); 
	
		solution_prev.compress(VectorOperation::insert);
		solution=solution_prev;
		
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
}


template class initBoundValProbs<1>;
template class initBoundValProbs<2>;
template class initBoundValProbs<3>;
