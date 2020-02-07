#include "../../include/solveClass.h"

template <int dim, class matrixType, class vectorType>
bool solveClass<dim, matrixType, vectorType>::nonlinearSolve(vectorType& U, bool converge_flag)
{	
	params_solve->enter_subsection("Nonlinear_solver");
	std::string nonLinearScheme=params_solve->get("nonLinear_method");	
	
	vectorType U_initial_0;
	if(!converge_flag) U_initial_0=U;
	if(std::strcmp(nonLinearScheme.c_str(),"classicNewton")==0){
		

		vectorType dU;
		dU.reinit(U);
		
	  double tol=params_solve->get_double("relative_norm_tolerance");
		double abs_tol=params_solve->get_double("absolute_norm_tolerance");
		int maxIteration=params_solve->get_integer("max_iterations");
		params_solve->leave_subsection();
		double initial_norm=0;
		double current_norm=0;
	  double res=1.0;
	  int currentIteration=0;
		
	  while (true){
	    if (currentIteration>=maxIteration) {
				PetscPrintf (mpi_communicator,"Maximum number of iterations reached without convergence. \n");
				if(!converge_flag) U=U_initial_0;
				return false;
				break;
			}
	    if (current_norm>1/std::pow(tol,2)){
				PetscPrintf (mpi_communicator,"\nNorm is too high. \n\n");
				if(!converge_flag) U=U_initial_0;
				return false;
				break;
			}
	    
			updateLinearSystem();
			
			current_norm=system_rhs.l2_norm();
			initial_norm=std::max(initial_norm, current_norm);
	    res=current_norm/initial_norm;
	    PetscPrintf (mpi_communicator,"Iter:%2u. Residual norm: %10.2e. Relative norm: %10.2e \n", currentIteration, current_norm, res); 
	    if (res<tol || current_norm< abs_tol){
				PetscPrintf (mpi_communicator,"Residual converged in %u iterations.\n", currentIteration);
			  return true;
			 break;
		 }
			
			solveLinearSystem_default_direct(dU);
			
			U+=dU;
			currentIteration++;
		}
	}
	
	if(std::strcmp(nonLinearScheme.c_str(),"NewtonLS")==0){
		
		vectorType dU;
		dU.reinit(U);
		
	  double tol=params_solve->get_double("relative_norm_tolerance");
		double abs_tol=params_solve->get_double("absolute_norm_tolerance");
		int maxIteration=params_solve->get_integer("max_iterations");
		std::string Line_search_scheme=params_solve->get("Line_search_scheme");
		params_solve->leave_subsection();
		
		double initial_norm=0;
		double current_norm=0;
	  double res=1.0;
	  int currentIteration=0;
		
	  while (true){
	    if (currentIteration>=maxIteration) {
				PetscPrintf (mpi_communicator,"Maximum number of iterations reached without convergence. \n");
				if(!converge_flag) U=U_initial_0;
				return false;
				 break;
			 }
	    if (current_norm>1/std::pow(tol,2)){
				PetscPrintf (mpi_communicator,"\nNorm is too high. \n\n");
				if(!converge_flag) U=U_initial_0;
				return false;
				 break;
			}
	    
			updateLinearSystem();
			
			current_norm=system_rhs.l2_norm();
			initial_norm=std::max(initial_norm, current_norm);
	    res=current_norm/initial_norm;
	    PetscPrintf (mpi_communicator,"Iter:%2u. Residual norm: %10.2e. Relative norm: %10.2e \n", currentIteration, current_norm, res); 
	    if (res<tol || current_norm< abs_tol){
				PetscPrintf (mpi_communicator,"Residual converged in %u iterations.\n", currentIteration);
				return true;
				 break;
			 }
			solveLinearSystem_default_direct(dU);
			
			double alpha=1;
			if(std::strcmp(Line_search_scheme.c_str(),"Backtracking")==0){
				vectorType U_initial;
				U_initial.reinit(U);
				params_solve->enter_subsection("Backtracking");
				double gamma=params_solve->get_double("Backtracking_gamma");
				double beta=params_solve->get_double("Backtracking_beta");
				int max_backtracking_iterations=params_solve->get_double("Backtracking_max_iterations");
				params_solve->leave_subsection();
				U+=dU;
				updateLinearSystem();
				double norm_0=system_rhs.l2_norm();
				double alpha_1=alpha*gamma;
				
				U=U_initial;
				U.add(alpha_1,dU);			
				updateLinearSystem();
				double norm_1=system_rhs.l2_norm(); 
				for(unsigned int line_i=0;line_i<max_backtracking_iterations;line_i++){
					if(norm_1>norm_0) {
						break;
					}
					alpha=alpha_1;
					norm_0=norm_1;
					alpha_1=alpha*beta;
					U=U_initial;
					U.add(alpha_1,dU);			
					updateLinearSystem();
					norm_1=system_rhs.l2_norm(); 	
				}
				U=U_initial;
			}

			U.add(alpha,dU);
			currentIteration++;
		}
	}
	
}

template <int dim, class matrixType, class vectorType>
bool solveClass<dim, matrixType, vectorType>::nonlinearSolve(vectorType& U,vectorType& dU,bool converge_flag)
{	
	vectorType U_initial_0;
	if(!converge_flag) U_initial_0=U;
	
	params_solve->enter_subsection("Nonlinear_solver");
	
	std::string nonLinearScheme=params_solve->get("nonLinear_method");
  double tol=params_solve->get_double("relative_norm_tolerance");
	double abs_tol=params_solve->get_double("absolute_norm_tolerance");
	int maxIteration=params_solve->get_integer("max_iterations");
	
	if(std::strcmp(nonLinearScheme.c_str(),"classicNewton")==0){	

		params_solve->leave_subsection();	
		double initial_norm=0;
		double current_norm=0;
	  double res=1.0;
	  int currentIteration=0;
		
	  while (true){
	    if (currentIteration>=maxIteration) {
				PetscPrintf (mpi_communicator,"Maximum number of iterations reached without convergence. \n");
				if(!converge_flag) U=U_initial_0;
				return false;
				 break;
			 }
	    if (current_norm>1/std::pow(tol,2)){
				PetscPrintf (mpi_communicator,"\nNorm is too high. \n\n");
				if(!converge_flag) U=U_initial_0;
				return false;
				break; 
			}
	    
			updateLinearSystem();
			
			current_norm=system_rhs.l2_norm();
			initial_norm=std::max(initial_norm, current_norm);
	    res=current_norm/initial_norm;
	    PetscPrintf (mpi_communicator,"Iter:%2u. Residual norm: %10.2e. Relative norm: %10.2e \n", currentIteration, current_norm, res); 
	    if (res<tol || current_norm< abs_tol){
				PetscPrintf (mpi_communicator,"Residual converged in %u iterations.\n", currentIteration);
				return true;
				break;
			 }
			
			solveLinearSystem_default_direct(dU);
			apply_dU_constrain(dU);
			
			
			U+=dU;
			currentIteration++;
		}
	}
	
	
	if(std::strcmp(nonLinearScheme.c_str(),"NewtonLS")==0){
		std::string Line_search_scheme=params_solve->get("Line_search_scheme");
		params_solve->leave_subsection();
		
		double initial_norm=0;
		double current_norm=0;
	  double res=1.0;
	  int currentIteration=0;
		
	  while (true){
	    if (currentIteration>=maxIteration) {
				PetscPrintf (mpi_communicator,"Maximum number of iterations reached without convergence. \n");
				if(!converge_flag) U=U_initial_0;
				return false;
				 break;
			 }
	    if (current_norm>1/std::pow(tol,2)){
				PetscPrintf (mpi_communicator,"\nNorm is too high. \n\n");
				if(!converge_flag) U=U_initial_0;
				return false;
				break;
			}
	    
			updateLinearSystem();
			
			current_norm=system_rhs.l2_norm();
			initial_norm=std::max(initial_norm, current_norm);
	    res=current_norm/initial_norm;
	    PetscPrintf (mpi_communicator,"Iter:%2u. Residual norm: %10.2e. Relative norm: %10.2e \n", currentIteration, current_norm, res); 
	    if (res<tol || current_norm< abs_tol){
				PetscPrintf (mpi_communicator,"Residual converged in %u iterations.\n", currentIteration);
				return true;
		    break;
			}
			solveLinearSystem_default_direct(dU);
			apply_dU_constrain(dU);
			
			double alpha=1;
			if(std::strcmp(Line_search_scheme.c_str(),"Backtracking")==0){
				vectorType U_initial;
				U_initial.reinit(U);
				params_solve->enter_subsection("Backtracking");
				double gamma=params_solve->get_double("Backtracking_gamma");
				double beta=params_solve->get_double("Backtracking_beta");
				int max_backtracking_iterations=params_solve->get_double("Backtracking_max_iterations");
				params_solve->leave_subsection();
				U+=dU;
				updateLinearSystem();
				double norm_0=system_rhs.l2_norm();
				double alpha_1=alpha*gamma;
				
				U=U_initial;
				U.add(alpha_1,dU);			
				updateLinearSystem();
				double norm_1=system_rhs.l2_norm(); 
				for(unsigned int line_i=0;line_i<max_backtracking_iterations;line_i++){
					if(norm_1>norm_0) {
						break;
					}
					alpha=alpha_1;
					norm_0=norm_1;
					alpha_1=alpha*beta;
					U=U_initial;
					U.add(alpha_1,dU);			
					updateLinearSystem();
					norm_1=system_rhs.l2_norm(); 	
				}
				U=U_initial;
			}
			U.add(alpha,dU);
			currentIteration++;
		}
	}	
}

template class solveClass<1, dealii::SparseMatrix<double>, dealii::Vector<double> >;
template class solveClass<2, dealii::SparseMatrix<double>, dealii::Vector<double> >;
template class solveClass<3, dealii::SparseMatrix<double>, dealii::Vector<double> >;
template class solveClass<1, PETScWrappers::MPI::SparseMatrix, PETScWrappers::MPI::Vector>;
template class solveClass<2, PETScWrappers::MPI::SparseMatrix, PETScWrappers::MPI::Vector>;
template class solveClass<3, PETScWrappers::MPI::SparseMatrix, PETScWrappers::MPI::Vector>;
