#include "../../include/solveClass.h"

template <int dim, class matrixType, class vectorType>
void solveClass<dim, matrixType, vectorType>::solveLinearSystem_default_direct(dealii::Vector<double>& dU)
{
	params_solve->enter_subsection("Linear_solver");
	std::string solver=params_solve->get("solver_method");
	params_solve->leave_subsection();	
	
}

template class solveClass<1, dealii::SparseMatrix<double>, dealii::Vector<double> >;
template class solveClass<2, dealii::SparseMatrix<double>, dealii::Vector<double> >;
template class solveClass<3, dealii::SparseMatrix<double>, dealii::Vector<double> >;