#include "../../include/solveClass.h"

template <int dim, class matrixType, class vectorType>
void solveClass<dim, matrixType, vectorType>::solveLinearSystem_default_direct(dealii::Vector<double>& dU)
{	
}

template class solveClass<1, dealii::SparseMatrix<double>, dealii::Vector<double> >;
template class solveClass<2, dealii::SparseMatrix<double>, dealii::Vector<double> >;
template class solveClass<3, dealii::SparseMatrix<double>, dealii::Vector<double> >;