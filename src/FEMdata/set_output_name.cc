#include "../../include/FEMdata.h"

template  <int dim, class vectorType>
void FEMdata<dim, vectorType>::set_output_name(std::vector<std::vector<std::string> > primary_variables)
{
	
	for(unsigned int i=0;i<primary_variables.size();i++){
		if(std::strcmp(primary_variables[i][1].c_str(),"component_is_scalar")==0) {
			nodal_solution_names.push_back(primary_variables[i][0]);
			nodal_data_component_interpretation.push_back(dealii::DataComponentInterpretation::component_is_scalar);
		}
		else if(std::strcmp(primary_variables[i][1].c_str(),"component_is_vector")==0){
		  for (unsigned int j=0; j<dim; ++j){
		    nodal_solution_names.push_back(primary_variables[i][0]);
				nodal_data_component_interpretation.push_back(dealii::DataComponentInterpretation::component_is_part_of_vector);
		  }
		}	
		else{
			PetscPrintf (mpi_communicator,"primary_variables component type does not support \n"); exit(1);
		}
	}
}

template class FEMdata<1, dealii::Vector<double> >;
template class FEMdata<2, dealii::Vector<double> >;
template class FEMdata<3, dealii::Vector<double> >;

template class FEMdata<1, dealii::PETScWrappers::MPI::Vector>;
template class FEMdata<2, dealii::PETScWrappers::MPI::Vector>;
template class FEMdata<3, dealii::PETScWrappers::MPI::Vector>;