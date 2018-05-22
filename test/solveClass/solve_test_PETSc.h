#ifndef solve_test_h
#define solve_test_h

#include <Sacado.hpp>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <math.h>
#include <mpi.h>

#include <cstdlib>
#include <ctime>

#include "solveClass.h"

template <int dim>
class solve_test_PETSc: public solveClass<dim, PETScWrappers::MPI::SparseMatrix, PETScWrappers::MPI::Vector> 
{
  public:
     solve_test_PETSc ();
    ~ solve_test_PETSc();

	
		void updateLinearSystem();
		  
		/*
		* test nonlinearSolve
		*test exmaple:drive residual vector to zero
		*test redisual vector: [(x_1^2,..., x_i^2,..., x_50^2), i=2:49 ]
		*/
		void nonlinearSolvePetsc_test();
		
		PETScWrappers::MPI::Vector solution;
		
		int n_total_dofs;
		int n_local_dofs;
		int max_couplings_between_dofs;
		
	  MPI_Comm mpi_communicator;
		
		int start_support_point;
	  int n_mpi_processes;
	  int this_mpi_process;
		
};

#endif