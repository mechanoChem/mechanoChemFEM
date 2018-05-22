#include "solve_test_PETSc.h"


template <int dim>
void solve_test_PETSc<dim>::nonlinearSolvePetsc_test()
{
	
	mpi_communicator=MPI_COMM_WORLD;
	MPI_Comm_rank(mpi_communicator, &this_mpi_process);
	MPI_Comm_size(mpi_communicator, &n_mpi_processes);
	
	parametersClass<dim> params;
	
	params.setString("nonLinearScheme", "classicNewton");
	params.setString("solver","PETScMUMPS");
	//params.setString("system_matrix_symmetricFlag", false); //defalut is false
	params.setDouble("relative_norm_tolerance",1.0e-15);
	params.setDouble("absolute_norm_tolerance",1.0e-12);
	params.setInt("maxIteration",50);
	
	this->reinitParams(params);
	
	n_total_dofs=2;
	n_local_dofs=n_total_dofs/n_mpi_processes;
	max_couplings_between_dofs=2;
	if(this_mpi_process==n_mpi_processes-1) n_local_dofs=n_total_dofs-(n_mpi_processes-1)*n_local_dofs;
	
	std::cout<<"this_mpi_process="<<this_mpi_process<<"; n_local_dofs="<<n_local_dofs<<std::endl;
	
	
	solution.reinit (mpi_communicator,n_total_dofs,n_local_dofs); 
	solution=1;
	this->setupLinearSystem(n_total_dofs, n_local_dofs, max_couplings_between_dofs);
	this->nonlinearSolve(solution);
	double solution_norm=solution.l2_norm();
	
	assert(solution_norm<1.0e-5);
	PetscPrintf (mpi_communicator,"solution_norm=%f, nonlinearSolvePetsc_test: PASS\n", solution_norm);
}

template class solve_test_PETSc<1>;
template class solve_test_PETSc<2>;
template class solve_test_PETSc<3>;