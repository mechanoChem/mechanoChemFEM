#include "../../include/FEMdata.h"

template  <int dim, class vectorType>
FEMdata<dim, vectorType>::FEMdata():mpi_communicator(MPI_COMM_WORLD)
{}
	
template  <int dim, class vectorType>
FEMdata<dim, vectorType>::FEMdata(dealii::hp::DoFHandler<dim>& _dof_handler):dof_handler(&_dof_handler),mpi_communicator(MPI_COMM_WORLD)
{
	data_out.attach_dof_handler (*dof_handler);
}

template  <int dim, class vectorType>
void FEMdata<dim, vectorType>::clear_data_vectors()
{
	data_out.clear_data_vectors();
}

template  <int dim, class vectorType>
FEMdata<dim, vectorType>::~FEMdata(){}

template class FEMdata<1, dealii::Vector<double> >;
template class FEMdata<2, dealii::Vector<double> >;
template class FEMdata<3, dealii::Vector<double> >;

template class FEMdata<1, dealii::PETScWrappers::MPI::Vector>;
template class FEMdata<2, dealii::PETScWrappers::MPI::Vector>;
template class FEMdata<3, dealii::PETScWrappers::MPI::Vector>;