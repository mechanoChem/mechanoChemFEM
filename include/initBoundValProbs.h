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
#include "supplementary/functionEvaluations.h"
#include "supplementary/parameters.h"
#include "supplementary/InitialConditions.h"

template <int dim>
class initBoundValProbs: public solveClass<dim, PETScWrappers::MPI::SparseMatrix, PETScWrappers::MPI::Vector>, public hpFEM<dim>
{
  public:
    initBoundValProbs (std::vector<std::vector<std::string> > _primary_variables, std::vector<std::vector<int> > _FE_support, ParameterHandler& _params);
    ~ initBoundValProbs();
		
		double current_dt;
				
		ParameterHandler* params;
		
		Residual<Sacado::Fad::DFad<double>,dim> ResidualEq;
		FEMdata<dim, PETScWrappers::MPI::Vector> FEMdata_out;
		
		/**
		*declare generic paramters
		*/
		void declare_parameters_initBoundValProbs();
		
		/**
		*this function should not be modified. 
		*/
		void setup_system();
		
		/**
		*overload function from solveClass
		*/
		virtual void updateLinearSystem();
		
		/**
		*abstract function
		*/
		virtual void get_residual(const typename hp::DoFHandler<dim>::active_cell_iterator &cell, const FEValues<dim>& fe_values, Table<1, Sacado::Fad::DFad<double> >& R, Table<1, Sacado::Fad::DFad<double>>& ULocal, Table<1, double >& ULocalConv);
		
		
		/**
		*generic function calls before running the simulations
		*/
		virtual void pre_run();
		/**
		*generic function calls to runthe simulations
		*/
    virtual void run ();	
		/**
		*default:make hyper_rectangle mesh
		*/
		virtual void make_grid();
		/**
		*mark boundary
		*/
		virtual void mark_boundary();
		/**
		*default: set all solution vector to be zero 
		*/
		virtual void apply_initial_condition();
		/**
		*default: pure virtual function
		*/
		virtual void refine_grid();
		/**
		*default: pure virtual function
		*/
		virtual void setMultDomain();
		/**
		*default: set up constrains
		*/
		virtual void setup_constraints();

		std::vector<std::vector<std::string> > primary_variables;
		std::vector<unsigned int > primary_variables_dof;
		
		std::vector<std::vector<int> > FE_support;
		std::vector<std::shared_ptr<FESystem<dim>> > fe_system;
		
		PETScWrappers::MPI::Vector solution;
		PETScWrappers::MPI::Vector solution_prev;
		
    const QGauss<dim>* volume_quadrature;
		const QGauss<dim-1>* common_face_quadrature;
    hp::FECollection<dim> fe_collection;
    hp::QCollection<dim>  q_collection;
    ConstraintMatrix      constraints;
		
		double domain_X, domain_Y, domain_Z;
		
		int current_increment;
		double current_time;
		double total_time;
		
	  MPI_Comm mpi_communicator;
	  const unsigned int n_mpi_processes;
	  const unsigned int this_mpi_process;
		ConditionalOStream pcout;
		
		std::string output_directory;
		std::string snapshot_directory;
		int skip_output;
		std::string snapfile;
		bool resuming_from_snapshot;
		bool save_snapshot;
		int off_output_index;
		
};

#endif