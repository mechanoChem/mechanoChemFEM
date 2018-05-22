#include"../../include/supplementary/CustomerPreconditioner.h"

template<class matrixType, class vectorType>
CustomerPreconditioner<matrixType, vectorType>::CustomerPreconditioner(const matrixType &A) : system_matrix (&A)
{}

template<class matrixType, class vectorType>	
void CustomerPreconditioner<matrixType, vectorType>::vmult (vectorType &dst, vectorType &src) const
{
	system_matrix->vmult(dst,src);
}

template class CustomerPreconditioner<PETScWrappers::MPI::SparseMatrix, PETScWrappers::MPI::Vector>;
template class CustomerPreconditioner<SparseMatrix<double> , Vector<double> >;