/*
zhenlin wang 2019
*/

#include"../../include/mechanoChemFEM.h"
template <int dim>
void mechanoChemFEM<dim>::run()
{			
	pcout<<std::endl<<std::endl;
	pcout<<"======== RUNNING... ========"<<std::endl;	
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
		
		output_results();
	}
	PetscPrintf(this->mpi_communicator,"Finish running!!\n");
}


template class mechanoChemFEM<1>;
template class mechanoChemFEM<2>;
template class mechanoChemFEM<3>;
