/*
zhenlin wang 2019
*/

#include"../../include/mechanoChemFEM.h"

template <int dim>
mechanoChemFEM<dim>::mechanoChemFEM(std::vector<std::vector<std::string> > _primary_variables, std::vector<std::vector<int> > _FE_support, ParameterHandler& _params)
	:solveClass<dim, PETScWrappers::MPI::SparseMatrix, PETScWrappers::MPI::Vector>(*this, _params),FEMdata_out(this->dof_handler),params(&_params),
	primary_variables(_primary_variables),FE_support(_FE_support),
	mpi_communicator(MPI_COMM_WORLD), n_mpi_processes (Utilities::MPI::n_mpi_processes(mpi_communicator)), 
  this_mpi_process (Utilities::MPI::this_mpi_process(mpi_communicator)), pcout (std::cout, (Utilities::MPI::this_mpi_process(MPI_COMM_WORLD) == 0))
{
	pcout<<"mechanoChemFEM initiated"<<std::endl;
	
	pre_run();
		
	make_grid();
	refine_grid();
	setMultDomain();
  mark_boundary();
	setup_linear_system();
		
}


template <int dim>
mechanoChemFEM<dim>::~mechanoChemFEM (){
	this->dof_handler.clear ();
	this->system_matrix.clear();
	this->system_rhs.clear();
	solution.clear();
	solution_prev.clear();
	delete volume_quadrature;
	delete common_face_quadrature;
}

template <int dim>
void mechanoChemFEM<dim>::refine_grid(){}

template <int dim>
void mechanoChemFEM<dim>::setMultDomain(){}

template <int dim>
void mechanoChemFEM<dim>::ini_updateLinearSystem(){}

template <int dim>
void mechanoChemFEM<dim>:: get_residual(const typename hp::DoFHandler<dim>::active_cell_iterator &cell, const FEValues<dim>& fe_values, Table<1, Sacado::Fad::DFad<double> >& R, Table<1, Sacado::Fad::DFad<double>>& ULocal, Table<1, double >& ULocalConv){}


template class mechanoChemFEM<1>;
template class mechanoChemFEM<2>;
template class mechanoChemFEM<3>;