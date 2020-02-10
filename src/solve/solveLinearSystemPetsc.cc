#include "../../include/solveClass.h"

template <int dim, class matrixType, class vectorType>
void solveClass<dim, matrixType, vectorType>::solveLinearSystem_default_direct(PETScWrappers::MPI::Vector& dU)
{
	params_solve->enter_subsection("Linear_solver");
	std::string solver_method=params_solve->get("solver_method");
	
	
	if(std::strcmp(solver_method.c_str(),"PETScsuperLU")==0){
    SolverControl solver_control;
    PETScWrappers::SolverPreOnly solver(solver_control, mpi_communicator);
    PETScWrappers::PreconditionLU preconditioner(system_matrix);
    solver.solve(system_matrix, dU, system_rhs, preconditioner);
	}
	else if(std::strcmp(solver_method.c_str(),"PETScMUMPS")==0){
		bool symmetricFlag=params_solve->get_bool("system_matrix_symmetricFlag");
    SolverControl cn;
    PETScWrappers::SparseDirectMUMPS solver(cn, mpi_communicator);
    if(symmetricFlag) solver.set_symmetric_mode(true);
    solver.solve (system_matrix, dU, system_rhs);
	}
	else if(std::strcmp(solver_method.c_str(),"own_solver")==0){
		solveLinearSystem(dU);		
	}
	
	PETScWrappers::Vector localized_dU(dU);
	constraints_solver.distribute (localized_dU);
	dU=localized_dU;
	params_solve->leave_subsection();	
}

template class solveClass<1, PETScWrappers::MPI::SparseMatrix, PETScWrappers::MPI::Vector>;
template class solveClass<2, PETScWrappers::MPI::SparseMatrix, PETScWrappers::MPI::Vector>;
template class solveClass<3, PETScWrappers::MPI::SparseMatrix, PETScWrappers::MPI::Vector>;
