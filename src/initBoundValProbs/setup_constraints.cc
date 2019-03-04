/*
zhenlin wang 2019
*/
#include"../../include/initBoundValProbs.h"

template <int dim>
void initBoundValProbs<dim>::setup_constraints()
{	
	hpFEM<dim>::constraints.clear ();
  DoFTools::make_hanging_node_constraints (hpFEM<dim>::dof_handler, hpFEM<dim>::constraints);
  hpFEM<dim>::constraints.close ();
  pcout << "   Number of active cells:       " << hpFEM<dim>::triangulation.n_active_cells() << std::endl;
  pcout << "   Number of degrees of freedom: " << hpFEM<dim>::dof_handler.n_dofs() << std::endl; 
  
}


template class initBoundValProbs<1>;
template class initBoundValProbs<2>;
template class initBoundValProbs<3>;
