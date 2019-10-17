/*
zhenlin wang 2019
*/
#include"../../include/initBoundValProbs.h"

template <int dim>
void initBoundValProbs<dim>::output()
{
	
	//write vtk and snapshot for solution
	std::string output_path = output_directory+"output-"+std::to_string(current_increment+off_output_index)+".vtk";
	
	FEMdata_out.clear_data_vectors();
	if(current_increment%skip_output==0) FEMdata_out.write_vtk(solution_prev, output_path);	
	if(save_snapshot){
		std::string snapshot_path = snapshot_directory+"snapshot-"+std::to_string(current_increment+off_output_index)+".dat";
		FEMdata_out.create_vector_snapshot(solution, snapshot_path);
	}
}


template class initBoundValProbs<1>;
template class initBoundValProbs<2>;
template class initBoundValProbs<3>;