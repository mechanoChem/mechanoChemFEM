#include "../../include/solveClass.h"

template <int dim, class matrixType, class vectorType>
solveClass<dim, matrixType, vectorType>::~solveClass()
{
	system_matrix.clear();
	system_rhs.clear();
}


template class solveClass<1, PETScWrappers::MPI::SparseMatrix, PETScWrappers::MPI::Vector>;
template class solveClass<2, PETScWrappers::MPI::SparseMatrix, PETScWrappers::MPI::Vector>;
template class solveClass<3, PETScWrappers::MPI::SparseMatrix, PETScWrappers::MPI::Vector>;