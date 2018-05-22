#include "initBoundValProbs.h"

template <int dim>
initBoundValProbs<dim>::initBoundValProbs(std::vector<std::vector<std::string> > _primary_variables, std::vector<std::vector<int> > _FE_support, ParameterHandler& _params)
	:solveClass<dim, PETScWrappers::MPI::SparseMatrix, PETScWrappers::MPI::Vector>(*this, _params),params(&_params),
	primary_variables(_primary_variables),FE_support(_FE_support),
	mpi_communicator(MPI_COMM_WORLD), n_mpi_processes (Utilities::MPI::n_mpi_processes(mpi_communicator)), 
  this_mpi_process (Utilities::MPI::this_mpi_process(mpi_communicator)), pcout (std::cout, (Utilities::MPI::this_mpi_process(MPI_COMM_WORLD) == 0)),
  computing_timer (pcout, TimerOutput::summary, TimerOutput::wall_times)
{
}

template <int dim>
initBoundValProbs<dim>::
initBoundValProbs (std::vector<std::vector<std::string> > _primary_variables, std::vector<std::vector<int> > _FE_support,std::vector<std::vector<std::string> > _variables_add, std::vector<std::vector<int> > _FE_support_add, ParameterHandler& _params)
	:solveClass<dim, PETScWrappers::MPI::SparseMatrix, PETScWrappers::MPI::Vector>(*this, _params),params(&_params),
	primary_variables(_primary_variables),FE_support(_FE_support), variables_add(_variables_add), FE_support_add(_FE_support_add),
	mpi_communicator(MPI_COMM_WORLD), n_mpi_processes (Utilities::MPI::n_mpi_processes(mpi_communicator)), 
  this_mpi_process (Utilities::MPI::this_mpi_process(mpi_communicator)), pcout (std::cout, (Utilities::MPI::this_mpi_process(MPI_COMM_WORLD) == 0)),
  computing_timer (pcout, TimerOutput::summary, TimerOutput::wall_times)
{
}

template <int dim>
initBoundValProbs<dim>::~initBoundValProbs (){}

template class initBoundValProbs<1>;
template class initBoundValProbs<2>;
template class initBoundValProbs<3>;