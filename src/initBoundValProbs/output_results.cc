/*
zhenlin wang 2019
*/
#include"../../include/mechanoChemFEM.h"

template <int dim>
void mechanoChemFEM<dim>::output_results()
{
	
	//write vtk and snapshot for solution
	if(save_output){ 
		std::string output_path = output_directory+"output-"+std::to_string(current_increment+off_output_index)+".vtk";
		FEMdata_out.clear_data_vectors();
		if(current_increment%skip_output==0) FEMdata_out.write_vtk(solution_prev, output_path);	
	}
	if(save_snapshot){
		std::string snapshot_path = snapshot_directory+"snapshot-"+std::to_string(current_increment+off_output_index)+".dat";
		FEMdata_out.create_vector_snapshot(solution, snapshot_path);
	}
}


template class mechanoChemFEM<1>;
template class mechanoChemFEM<2>;
template class mechanoChemFEM<3>;
