#include"initBoundValProbs.h"

template <int dim>
void initBoundValProbs<dim>::setup_system()
{
	GridTools::partition_triangulation (n_mpi_processes, hpFEM<dim>::triangulation);
	hpFEM<dim>::dof_handler.distribute_dofs (fe_collection);
	//DoFRenumbering::subdomain_wise(dof_handler);
	DoFRenumbering::component_wise (hpFEM<dim>::dof_handler);
	const types::global_dof_index n_local_dofs = DoFTools::count_dofs_with_subdomain_association(hpFEM<dim>::dof_handler, this_mpi_process);
	const types::global_dof_index n_total_dofs=hpFEM<dim>::dof_handler.n_dofs();
	const types::global_dof_index max_couplings_between_dofs=hpFEM<dim>::dof_handler.max_couplings_between_dofs();
	
	this->setupLinearSystem(n_total_dofs, n_local_dofs, max_couplings_between_dofs);
  									
	solution.reinit (mpi_communicator,n_total_dofs,n_local_dofs); 
	solution_prev.reinit (mpi_communicator,n_total_dofs,n_local_dofs);
}

template class initBoundValProbs<1>;
template class initBoundValProbs<2>;
template class initBoundValProbs<3>;