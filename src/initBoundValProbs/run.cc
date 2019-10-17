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
	apply_initial_condition();


	PetscPrintf(this->mpi_communicator,"running....\n\n");
	clock_t t_solve;	
	t_solve = clock();
  for (; current_time<=total_time; current_time+=current_dt){
    current_increment++;
		PetscPrintf(this->mpi_communicator,"************");
		PetscPrintf(this->mpi_communicator,"current increment=%d, current time= %f",current_increment, current_time);
		PetscPrintf(this->mpi_communicator,"************\n");
		solve();
		
	  t_solve = clock() - t_solve;
		PetscPrintf(this->mpi_communicator,"It took me %f seconds for this solve.\n ",((float)t_solve)/CLOCKS_PER_SEC);
		PetscPrintf(this->mpi_communicator,"\n\n");
		
		output();
	}
	PetscPrintf(this->mpi_communicator,"Finish running!!\n");
}


template class initBoundValProbs<1>;
template class initBoundValProbs<2>;
template class initBoundValProbs<3>;
