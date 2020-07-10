#include "../../include/solveClass.h"

template <int dim, class matrixType, class vectorType>
void solveClass<dim, matrixType, vectorType>::declare_parameters_solveClass()
{
	params_solve->enter_subsection("Nonlinear_solver");
	params_solve->declare_entry("nonLinear_method","classicNewton", Patterns::Selection("classicNewton|NewtonLS"),"support method: classicNewton, NewtonLS");
	params_solve->declare_entry("relative_norm_tolerance", "1.0e-12", Patterns::Double() );
	params_solve->declare_entry("absolute_norm_tolerance", "1.0e-12", Patterns::Double() );
	params_solve->declare_entry("max_iterations", "5", Patterns::Integer() );
	params_solve->declare_entry("Line_search_scheme", "NONE", Patterns::Selection("NONE|Backtracking"),"support shceme: Backtracking");
	params_solve->leave_subsection();	
	
	params_solve->enter_subsection("Backtracking");
	params_solve->declare_entry("Backtracking_gamma","0",Patterns::Double());
	params_solve->declare_entry("Backtracking_beta","0",Patterns::Double());
	params_solve->declare_entry("Backtracking_max_iterations","0",Patterns::Integer() );
	params_solve->leave_subsection();	
											 
	params_solve->enter_subsection("Linear_solver");
	params_solve->declare_entry("solver_method","PETScsuperLU",
		 										Patterns::Selection("PETScMUMPS|PETScsuperLU|own_solver"),"linear solver");
	params_solve->declare_entry("system_matrix_symmetricFlag", "false", Patterns::Bool() );
	params_solve->leave_subsection();	
}

template class solveClass<1, dealii::SparseMatrix<double>, dealii::Vector<double> >;
template class solveClass<2, dealii::SparseMatrix<double>, dealii::Vector<double> >;
template class solveClass<3, dealii::SparseMatrix<double>, dealii::Vector<double> >;
template class solveClass<1, PETScWrappers::MPI::SparseMatrix, PETScWrappers::MPI::Vector>;
template class solveClass<2, PETScWrappers::MPI::SparseMatrix, PETScWrappers::MPI::Vector>;
template class solveClass<3, PETScWrappers::MPI::SparseMatrix, PETScWrappers::MPI::Vector>;