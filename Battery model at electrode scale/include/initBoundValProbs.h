#ifndef initBoundValProbs_h
#define initBoundValProbs_h

//more headers in deal.ii

#include <deal.II/base/tensor_function.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/logstream.h>
#include <deal.II/base/function.h>
#include <deal.II/base/conditional_ostream.h>
#include <deal.II/base/utilities.h>
#include <deal.II/base/timer.h>
#include <deal.II/base/parameter_handler.h>

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

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/grid_in.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/grid/tria_boundary_lib.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/grid_refinement.h>
#include <deal.II/grid/manifold_lib.h>

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_renumbering.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_dgq.h>
#include <deal.II/fe/fe_nothing.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/fe_values.h>

#include <deal.II/hp/dof_handler.h>
#include <deal.II/hp/q_collection.h>
#include <deal.II/hp/fe_values.h>
#include <deal.II/hp/fe_collection.h>

#include <deal.II/numerics/vector_tools.h>
#include <deal.II/numerics/matrix_tools.h>
#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/solution_transfer.h>

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


#include "hpFEM.h"
#include "Residual.h"
#include "solveClass.h"
#include "FEMdata.h"
#include "electricChemo.h"
#include "supplementary/functionEvaluations.h"


using namespace dealii;

template <int dim>
class initBoundValProbs: public solveClass<dim, PETScWrappers::MPI::SparseMatrix, PETScWrappers::MPI::Vector>, public hpFEM<dim>
{
  public:
    initBoundValProbs (ParameterHandler& _params);
    ~ initBoundValProbs();
		
		/**
		*domain_id;
		*/
		unsigned int 
			separator_id,
			electrode_id;
		
		unsigned int 
			separator_fe,
			electrode_fe;
		/**
		*dof_id;
		*/
   unsigned int
			u_dof,
			c_li_plus_dof,
			phi_e_dof,
			c_li_dof,
			phi_s_dof,
			T_dof;
		
		ElectricChemo<Sacado::Fad::DFad<double>,dim>* electricChemoFormula;
		
		/**
		*global variables
		*/
		double current_dt;
		double current_IpA;
		double electrode_Y1, electrode_Y2;
		
		void step_load();
			
		/*
		--------------------------------------------------------------------------------------------
		* functions and variables defined below should be generic to any initBoundValProbs
		--------------------------------------------------------------------------------------------
		*/
		
		ParameterHandler* params;
		void declare_parameters();
		
		Residual<Sacado::Fad::DFad<double>,dim>* ResidualEq;
		FEMdata<dim, PETScWrappers::MPI::Vector>* FEMdata_out;
		/**
		*strickly using this name without any argument!
		*/
		void updateLinearSystem();		
		/**
		*function need to be modified for speific problem.
		*/
    void run (std::vector<std::vector<std::string> > _primary_variables, std::vector<std::vector<int> > _FE_support);	
	  void setup_system();	
		void make_grid();
		void mark_boundary();
		void apply_initial_condition();
		void refine_grid();
		void setMultDomain();
		void setup_constraints();

		std::vector<std::vector<std::string> > primary_variables;
		std::vector<std::vector<int> > FE_support;
		
		std::vector<unsigned int > primary_variables_dof;
		std::vector<std::shared_ptr<FESystem<dim>> > fe_system;
		
		PETScWrappers::MPI::Vector solution;
		PETScWrappers::MPI::Vector solution_prev;
		
    const QGauss<dim>* volume_quadrature;
		const QGauss<dim-1>* common_face_quadrature;
    hp::FECollection<dim> fe_collection;
    hp::QCollection<dim>  q_collection;
    ConstraintMatrix      constraints;
		
		int current_increment;
		double current_time;
		
	  MPI_Comm mpi_communicator;
	  const unsigned int n_mpi_processes;
	  const unsigned int this_mpi_process;
		TimerOutput computing_timer;
		ConditionalOStream pcout;		
};

#endif