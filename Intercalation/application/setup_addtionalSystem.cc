#include "initBoundValProbs.h"

template <int dim>
void initBoundValProbs<dim>::setup_addtionalSystem()
{
	//adjoint Fe_system with same volume_quadrature	
	params->enter_subsection("Problem");
		int first_domain_id=params->get_integer("first_domain_id");
	params->leave_subsection();	
	
	std::vector<unsigned int > variables_dof_tem;
	this->setup_FeSystem(fe_system_add,fe_collection_add, q_collection_add, variables_dof_tem,variables_add,FE_support_add,*volume_quadrature);

	this->set_active_fe_indices (FE_support_add, *dof_handler_add, first_domain_id);
	dof_handler_add->distribute_dofs (fe_collection_add);
	
	DoFRenumbering::component_wise (*dof_handler_add);
	const types::global_dof_index n_local_dofs = DoFTools::count_dofs_with_subdomain_association(*dof_handler_add, this_mpi_process);
	const types::global_dof_index n_total_dofs=dof_handler_add->n_dofs();
										
	additional_data.reinit (mpi_communicator,n_total_dofs,n_local_dofs);
 
	additional_data=0;
}
template class initBoundValProbs<1>;
template class initBoundValProbs<2>;
template class initBoundValProbs<3>;
