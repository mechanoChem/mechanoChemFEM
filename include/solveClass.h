
#ifndef solve_h
#define solve_h

#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <string.h>
#include <ctime>
#include <cmath>
#include <vector>
#include <map>
#include <mpi.h>

#include <deal.II/lac/block_vector.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/block_sparse_matrix.h>
#include <deal.II/lac/sparse_direct.h>
#include <deal.II/lac/precondition.h>
#include <deal.II/lac/constraint_matrix.h>
#include <deal.II/lac/solver_gmres.h>
#include <deal.II/lac/petsc_vector.h>
#include <deal.II/lac/petsc_parallel_vector.h>
#include <deal.II/lac/petsc_parallel_sparse_matrix.h>
#include <deal.II/lac/petsc_solver.h>
#include <deal.II/lac/petsc_precondition.h>
#include <deal.II/base/parameter_handler.h>

#include "hpFEM.h"


template <int dim, class matrixType, class vectorType>
class solveClass
{
public:
	/**
	*constructor without arguments
	*/
	solveClass();
	/**
	*constructor with abstract classhpFEM<dim>
	*/
	solveClass(hpFEM<dim>& _FEMbase);
	/**
	*constructor with abstract classhpFEM<dim> and parametersClass<dim>
	*/
  solveClass(hpFEM<dim>& _FEMbase, dealii::ParameterHandler& _params);
  ~solveClass();
	
	/**
	*set parametersClass
	*/
	void reinitParams(dealii::ParameterHandler& _params);
	
	/**
	*re-set abstract class hpFEM<dim>. Though it should be passed by constructor  
	*/
	void reinitMethod(hpFEM<dim>& _FEMbase);
	
	void declare_parameters();
	
	/**
	*set up size and partition of system_matrix and system_rhs
	*/
	void setupLinearSystem(int n_total_dofs, int n_local_dofs, int max_couplings_between_dofs);
	
	/**
	*set system_matrix=0 and system_rhs=0
	*/
	void reinitLinearSystem();

	/**
	*linear Solve
	*/
	void linearSolve(vectorType& U);
	
	/**
	*nonlinear Solve
	*/
	void nonlinearSolve(vectorType& U);
	
	/**
	*nonlinear Solve with access of incremental solution dU
	*/
	void nonlinearSolve(vectorType& U,vectorType& dU);
	
	/**
	*solve Ax=b with dealii matrix and vector type by default_direct solver
	*/
	void solveLinearSystem_default_direct(dealii::Vector<double>& dU);
	/**
	*solve Ax=b with PETSc matrix and vector type by default_direct solver
	*/
  void solveLinearSystem_default_direct(PETScWrappers::MPI::Vector& dU);
	/**
	*solve Ax=b with PETSc matrix and vector type by customer solver
	*/
	virtual void solveLinearSystem(PETScWrappers::MPI::Vector& dU);
	
	/**
	*wrapper for dealii::constraints.distribute_local_to_global
	*/
	void distribute_local_to_global(dealii::FullMatrix<double>& local_matrix, dealii::Vector<double>& local_rhs, std::vector<types::global_dof_index> local_dof_indices);
	/**
	*do PETScWrappers::VectorBase::compress(VectorOperation::add) for system_matrix and system_rhs;
	*/
	void LinearSystemCompressAdd();
	
	virtual void updateLinearSystem();
	
	/**
	*do PETScWrappers::VectorBase::compress(VectorOperation::add) for system_matrix and system_rhs;
	*/
	virtual void apply_dU_constrain(vectorType& dU);
	
	hpFEM<dim>* FEMbase;
	dealii::ParameterHandler* params;
	
  matrixType system_matrix;
  vectorType system_rhs;
	
	MPI_Comm mpi_communicator;
	

};

#endif