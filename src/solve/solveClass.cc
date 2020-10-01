#include "../../include/solveClass.h"


template <int dim, class matrixType, class vectorType>
solveClass<dim, matrixType, vectorType>::solveClass():mpi_communicator(MPI_COMM_WORLD)
{
	params_solve=&this->params_ref;
	params_solve_json=&this->params_ref_json;
	declare_parameters_solveClass();
}

template <int dim, class matrixType, class vectorType>
solveClass<dim, matrixType, vectorType>::~solveClass()
{
	system_matrix.clear();
}


template <int dim, class matrixType, class vectorType>
void solveClass<dim, matrixType, vectorType>::reinitLinearSystem()
{
	system_matrix=0;
	system_rhs=0;
}

template <int dim, class matrixType, class vectorType>
void solveClass<dim, matrixType, vectorType>::distribute_local_to_global(dealii::FullMatrix<double>& local_matrix, dealii::Vector<double>& local_rhs, std::vector<types::global_dof_index> local_dof_indices)
{
	constraints_solver.distribute_local_to_global (local_matrix, local_rhs, local_dof_indices, system_matrix, system_rhs);
}

template <int dim, class matrixType, class vectorType>
void solveClass<dim, matrixType, vectorType>::apply_dU_constrain(vectorType& dU)
{}

template <int dim, class matrixType, class vectorType>
void solveClass<dim, matrixType, vectorType>::solveLinearSystem(PETScWrappers::MPI::Vector& dU){}

template <int dim, class matrixType, class vectorType>
void solveClass<dim, matrixType, vectorType>::updateLinearSystem(){}

template class solveClass<1, dealii::SparseMatrix<double>, dealii::Vector<double> >;
template class solveClass<2, dealii::SparseMatrix<double>, dealii::Vector<double> >;
template class solveClass<3, dealii::SparseMatrix<double>, dealii::Vector<double> >;
template class solveClass<1, PETScWrappers::MPI::SparseMatrix, PETScWrappers::MPI::Vector>;
template class solveClass<2, PETScWrappers::MPI::SparseMatrix, PETScWrappers::MPI::Vector>;
template class solveClass<3, PETScWrappers::MPI::SparseMatrix, PETScWrappers::MPI::Vector>;