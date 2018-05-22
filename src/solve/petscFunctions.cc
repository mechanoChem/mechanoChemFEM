#include "../../include/solveClass.h"

template <int dim, class matrixType, class vectorType>
void solveClass<dim, matrixType, vectorType>::LinearSystemCompressAdd()
{
  system_matrix.compress(VectorOperation::add);
  system_rhs.compress(VectorOperation::add);
}


template <int dim, class matrixType, class vectorType>
void solveClass<dim, matrixType, vectorType>::setupLinearSystem(int n_total_dofs, int n_local_dofs, int max_couplings_between_dofs)
{
  system_matrix.reinit (mpi_communicator, n_total_dofs, n_total_dofs, n_local_dofs, n_local_dofs, max_couplings_between_dofs);
	system_rhs.reinit (mpi_communicator, n_total_dofs, n_local_dofs); 
}

template class solveClass<1, PETScWrappers::MPI::SparseMatrix, PETScWrappers::MPI::Vector>;
template class solveClass<2, PETScWrappers::MPI::SparseMatrix, PETScWrappers::MPI::Vector>;
template class solveClass<3, PETScWrappers::MPI::SparseMatrix, PETScWrappers::MPI::Vector>;