#include "../../include/FEMdata.h"

template  <int dim, class vectorType>
void FEMdata<dim, vectorType>::write_vtk(vectorType& solution, std::string path)
{
	dealii::Vector<double> Un(solution);
	int this_mpi_process;
	MPI_Comm_rank(mpi_communicator, &this_mpi_process);
	if(this_mpi_process == 0){
		std::ofstream file(path.c_str());
		//Add nodal DOF data
		data_out.add_data_vector (Un, nodal_solution_names, dealii::DataOut<dim,dealii::hp::DoFHandler<dim> >::type_dof_data, nodal_data_component_interpretation);
	
		data_out.build_patches (); data_out.write_vtk (file); file.close();
	}
}

template class FEMdata<1, dealii::Vector<double> >;
template class FEMdata<2, dealii::Vector<double> >;
template class FEMdata<3, dealii::Vector<double> >;

template class FEMdata<1, dealii::PETScWrappers::MPI::Vector>;
template class FEMdata<2, dealii::PETScWrappers::MPI::Vector>;
template class FEMdata<3, dealii::PETScWrappers::MPI::Vector>;
