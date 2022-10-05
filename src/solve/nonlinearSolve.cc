#include "../../include/solveClass.h"

template <int dim, class matrixType, class vectorType>
int solveClass<dim, matrixType, vectorType>::nonlinearSolve(vectorType& U, bool converge_flag)
{	
	vectorType U_initial_0;
	if(!converge_flag) U_initial_0=U;
	
	vectorType dU;
	dU.reinit(U);
	
	std::string nonLinearScheme;
  double tol;
	double abs_tol;
	int maxIteration;
	if(this->use_ParameterHandler){
		params_solve->enter_subsection("Nonlinear_solver");
		nonLinearScheme=params_solve->get("nonLinear_method");	
	  tol=params_solve->get_double("relative_norm_tolerance");
		abs_tol=params_solve->get_double("absolute_norm_tolerance");
		maxIteration=params_solve->get_integer("max_iterations");
		params_solve->leave_subsection();
	}
	if(this->use_ParameterJson){
		nonLinearScheme=(*params_solve_json)["Nonlinear_solver"]["nonLinear_method"];	
	  tol=(*params_solve_json)["Nonlinear_solver"]["relative_norm_tolerance"];
		abs_tol=(*params_solve_json)["Nonlinear_solver"]["absolute_norm_tolerance"];
		maxIteration=(*params_solve_json)["Nonlinear_solver"]["max_iterations"];
	}
	
	if(std::strcmp(nonLinearScheme.c_str(),"classicNewton")==0){
		double initial_norm=0;
		double current_norm=0;
		double min_norm=1e10;
	  double res=1.0;
	  int currentIteration=0;
		
	  while (true){
	    if (currentIteration>=maxIteration) {
				PetscPrintf (mpi_communicator,"Maximum number of iterations reached without convergence. \n");
				if(!converge_flag) U=U_initial_0;
					//return false;
					return -1;
					break;
			}
	    if (current_norm>1/std::pow(tol,2)){
				PetscPrintf (mpi_communicator,"\nNorm is too high. \n\n");
				if(!converge_flag) U=U_initial_0;
					//return false;
					return -1;
					break;
			}
		if (std::log10(current_norm/min_norm) > 3){
				PetscPrintf (mpi_communicator,"\nNorm is increasing. \n\n");
				if(!converge_flag) U=U_initial_0;
					//return false;
					return -1;
					break;
      		}
      	if (isnan(current_norm)){
				PetscPrintf (mpi_communicator,"\nNorm is nan. \n\n");
				if(!converge_flag) U=U_initial_0;
					//return false;
					return -1;
					break;
      		}
	    
			updateLinearSystem();
			
			current_norm=system_rhs.l2_norm();
			initial_norm=std::max(initial_norm, current_norm);
			min_norm=std::min(min_norm, current_norm);
	    res=current_norm/initial_norm;
	    PetscPrintf (mpi_communicator,"Iter:%2u. Residual norm: %10.2e. Relative norm: %10.2e \n", currentIteration, current_norm, res); 
	    if (res<tol || current_norm< abs_tol){
				PetscPrintf (mpi_communicator,"Residual converged in %u iterations.\n", currentIteration);
			//return true;
			return currentIteration;
			break;
		 }
			
			solveLinearSystem_default_direct(dU);
			
			U+=dU;
			currentIteration++;
		}
	}
	
	if(std::strcmp(nonLinearScheme.c_str(),"NewtonLS")==0){
		std::string Line_search_scheme;
		if(this->use_ParameterHandler){
			params_solve->enter_subsection("Nonlinear_solver");
			Line_search_scheme=params_solve->get("Line_search_scheme");
			params_solve->leave_subsection();
		}
		if(this->use_ParameterJson){
			Line_search_scheme=(*params_solve_json)["Nonlinear_solver"]["Line_search_scheme"]["method"];
		}
		
		double initial_norm=0;
		double min_norm=1e10;
		double current_norm=0;
	  double res=1.0;
	  int currentIteration=0;
		
	  while (true){
	    if (currentIteration>=maxIteration) {
				PetscPrintf (mpi_communicator,"Maximum number of iterations reached without convergence. \n");
				if(!converge_flag) U=U_initial_0;
					//return false;
					return -1;
					break;
			 }
	    if (current_norm>1/std::pow(tol,2)){
				PetscPrintf (mpi_communicator,"\nNorm is too high. \n\n");
				if(!converge_flag) U=U_initial_0;
					//return false; 
					return -1;
					break;
			}
		if (std::log10(current_norm/min_norm) > 3){
				PetscPrintf (mpi_communicator,"\nNorm is increasing. \n\n");
				if(!converge_flag) U=U_initial_0;
					//return false;
					return -1;
					break;
			}
		if (isnan(current_norm)){
				PetscPrintf (mpi_communicator,"\nNorm is nan. \n\n");
				if(!converge_flag) U=U_initial_0;
					//return false;
					return -1;
					break;
			}
	    
			updateLinearSystem();
			
			current_norm=system_rhs.l2_norm();
			initial_norm=std::max(initial_norm, current_norm);
			min_norm=std::min(min_norm, current_norm);
	    res=current_norm/initial_norm;
	    PetscPrintf (mpi_communicator,"Iter:%2u. Residual norm: %10.2e. Relative norm: %10.2e \n", currentIteration, current_norm, res); 
	    if (res<tol || current_norm< abs_tol){
				PetscPrintf (mpi_communicator,"Residual converged in %u iterations.\n", currentIteration);
				//return true;
				return currentIteration;
				break;
			}
			solveLinearSystem_default_direct(dU);
			
			double alpha=1;
			if(std::strcmp(Line_search_scheme.c_str(),"Backtracking")==0){
				vectorType U_initial;
				U_initial.reinit(U);
				
				double gamma;
				double beta;
				int max_backtracking_iterations;
				
				if(this->use_ParameterHandler){
					params_solve->enter_subsection("Backtracking");
					gamma=params_solve->get_double("Backtracking_gamma");
					beta=params_solve->get_double("Backtracking_beta");
					max_backtracking_iterations=params_solve->get_double("Backtracking_max_iterations");
					params_solve->leave_subsection();
				}
				
				if(this->use_ParameterJson){
					gamma=(*params_solve_json)["Nonlinear_solver"]["Line_search_scheme"]["Backtracking"]["Backtracking_gamma"];
					beta=(*params_solve_json)["Nonlinear_solver"]["Line_search_scheme"]["Backtracking"]["Backtracking_beta"];
					max_backtracking_iterations=(*params_solve_json)["Nonlinear_solver"]["Line_search_scheme"]["Backtracking"]["Backtracking_max_iterations"];
				}
				
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
int solveClass<dim, matrixType, vectorType>::nonlinearSolve(vectorType& U,vectorType& dU,bool converge_flag)
{	
	vectorType U_initial_0;
	if(!converge_flag) U_initial_0=U;
	
	std::string nonLinearScheme;
  double tol;
	double abs_tol;
	int maxIteration;
	if(this->use_ParameterHandler){
		params_solve->enter_subsection("Nonlinear_solver");
		nonLinearScheme=params_solve->get("nonLinear_method");	
	  tol=params_solve->get_double("relative_norm_tolerance");
		abs_tol=params_solve->get_double("absolute_norm_tolerance");
		maxIteration=params_solve->get_integer("max_iterations");
		params_solve->leave_subsection();
	}
	if(this->use_ParameterJson){
		nonLinearScheme=(*params_solve_json)["Nonlinear_solver"]["nonLinear_method"];	
	  tol=(*params_solve_json)["Nonlinear_solver"]["relative_norm_tolerance"];
		abs_tol=(*params_solve_json)["Nonlinear_solver"]["absolute_norm_tolerance"];
		maxIteration=(*params_solve_json)["Nonlinear_solver"]["max_iterations"];
	}
	
	if(std::strcmp(nonLinearScheme.c_str(),"classicNewton")==0){
		double initial_norm=0;
		double current_norm=0;
		double min_norm=1e10;
	  double res=1.0;
	  int currentIteration=0;
		
	  while (true){
	    if (currentIteration>=maxIteration) {
				PetscPrintf (mpi_communicator,"Maximum number of iterations reached without convergence. \n");
				if(!converge_flag) U=U_initial_0;
					//return false;
					return -1;
					break;
			}
	    if (current_norm>1/std::pow(tol,2)){
				PetscPrintf (mpi_communicator,"\nNorm is too high. \n\n");
				if(!converge_flag) U=U_initial_0;
					//return false;
					return -1;
					break;
			}
		if (std::log10(current_norm/min_norm) > 3){
				PetscPrintf (mpi_communicator,"\nNorm is increasing. \n\n");
					if(!converge_flag) U=U_initial_0;
					//return false;
					return -1;
					break;
      		}
		if (isnan(current_norm)){
				PetscPrintf (mpi_communicator,"\nNorm is nan. \n\n");
					if(!converge_flag) U=U_initial_0;
					//return false;
					return -1;
					break;
      		}
	    
			updateLinearSystem();
			
			current_norm=system_rhs.l2_norm();
			initial_norm=std::max(initial_norm, current_norm);
			min_norm=std::min(min_norm, current_norm);
	    res=current_norm/initial_norm;
	    PetscPrintf (mpi_communicator,"Iter:%2u. Residual norm: %10.2e. Relative norm: %10.2e \n", currentIteration, current_norm, res); 
	    if (res<tol || current_norm< abs_tol){
				PetscPrintf (mpi_communicator,"Residual converged in %u iterations.\n", currentIteration);
				//return true;
				return currentIteration;
				break;
			}
			
			solveLinearSystem_default_direct(dU);
			
			U+=dU;
			currentIteration++;
		}
	}
	
	if(std::strcmp(nonLinearScheme.c_str(),"NewtonLS")==0){
		std::string Line_search_scheme;
		if(this->use_ParameterHandler){
			params_solve->enter_subsection("Nonlinear_solver");
			Line_search_scheme=params_solve->get("Line_search_scheme");
			params_solve->leave_subsection();
		}
		if(this->use_ParameterJson){
			Line_search_scheme=(*params_solve_json)["Nonlinear_solver"]["Line_search_scheme"]["method"];
		}
		
		double initial_norm=0;
		double current_norm=0;
		double min_norm=1e10;
	  double res=1.0;
	  int currentIteration=0;
		
	  while (true){
	    if (currentIteration>=maxIteration) {
				PetscPrintf (mpi_communicator,"Maximum number of iterations reached without convergence. \n");
				if(!converge_flag) U=U_initial_0;
					//return false; 
					return -1;
					break;
			 }
	    if (current_norm>1/std::pow(tol,2)){
				PetscPrintf (mpi_communicator,"\nNorm is too high. \n\n");
				if(!converge_flag) U=U_initial_0;
					//return false;
					return -1;
					break;
			}
		if (std::log10(current_norm/min_norm) > 3){
				PetscPrintf (mpi_communicator,"\nNorm is increasing. \n\n");
					if(!converge_flag) U=U_initial_0;
					//return false;
					return -1;
					break;
      		}
		if (isnan(current_norm)){
				PetscPrintf (mpi_communicator,"\nNorm is nan. \n\n");
					if(!converge_flag) U=U_initial_0;
					//return false;
					return -1;
					break;
      		}
	    
			updateLinearSystem();
			
			current_norm=system_rhs.l2_norm();
			initial_norm=std::max(initial_norm, current_norm);
			min_norm=std::min(min_norm, current_norm);
	    res=current_norm/initial_norm;
	    PetscPrintf (mpi_communicator,"Iter:%2u. Residual norm: %10.2e. Relative norm: %10.2e \n", currentIteration, current_norm, res); 
	    if (res<tol || current_norm< abs_tol){
				PetscPrintf (mpi_communicator,"Residual converged in %u iterations.\n", currentIteration);
				//return true;
				return currentIteration;
				break;
			}
			solveLinearSystem_default_direct(dU);
			
			double alpha=1;
			if(std::strcmp(Line_search_scheme.c_str(),"Backtracking")==0){
				vectorType U_initial;
				U_initial.reinit(U);
				
				double gamma;
				double beta;
				int max_backtracking_iterations;
				
				if(this->use_ParameterHandler){
					params_solve->enter_subsection("Backtracking");
					gamma=params_solve->get_double("Backtracking_gamma");
					beta=params_solve->get_double("Backtracking_beta");
					max_backtracking_iterations=params_solve->get_double("Backtracking_max_iterations");
					params_solve->leave_subsection();
				}
				
				if(this->use_ParameterJson){
					gamma=(*params_solve_json)["Nonlinear_solver"]["Line_search_scheme"]["Backtracking"]["Backtracking_gamma"];
					beta=(*params_solve_json)["Nonlinear_solver"]["Line_search_scheme"]["Backtracking"]["Backtracking_beta"];
					max_backtracking_iterations=(*params_solve_json)["Nonlinear_solver"]["Line_search_scheme"]["Backtracking"]["Backtracking_max_iterations"];
				}
				
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
