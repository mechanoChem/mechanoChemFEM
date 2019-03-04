/*
zhenlin wang 2019
*/

#include"../../include/initBoundValProbs.h"

template <int dim>
initBoundValProbs<dim>::initBoundValProbs(std::vector<std::vector<std::string> > _primary_variables, std::vector<std::vector<int> > _FE_support, ParameterHandler& _params)
	:solveClass<dim, PETScWrappers::MPI::SparseMatrix, PETScWrappers::MPI::Vector>(*this, _params),FEMdata_out(this->dof_handler),params(&_params),
	primary_variables(_primary_variables),FE_support(_FE_support),
	mpi_communicator(MPI_COMM_WORLD), n_mpi_processes (Utilities::MPI::n_mpi_processes(mpi_communicator)), 
  this_mpi_process (Utilities::MPI::this_mpi_process(mpi_communicator)), pcout (std::cout, (Utilities::MPI::this_mpi_process(MPI_COMM_WORLD) == 0))
{
	pcout<<"initBoundValProbs initiated"<<std::endl;
}


template <int dim>
initBoundValProbs<dim>::~initBoundValProbs (){
	this->dof_handler.clear ();
	this->system_matrix.clear();
	this->system_rhs.clear();
	solution.clear();
	solution_prev.clear();
}

template <int dim>
void initBoundValProbs<dim>::refine_grid(){}

template <int dim>
void initBoundValProbs<dim>::setMultDomain(){}

template <int dim>
void initBoundValProbs<dim>:: get_residual(const typename hp::DoFHandler<dim>::active_cell_iterator &cell, const FEValues<dim>& fe_values, Table<1, Sacado::Fad::DFad<double> >& R, Table<1, Sacado::Fad::DFad<double>>& ULocal, Table<1, double >& ULocalConv){}


template class initBoundValProbs<1>;
template class initBoundValProbs<2>;
template class initBoundValProbs<3>;