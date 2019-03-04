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
}


template class initBoundValProbs<1>;
template class initBoundValProbs<2>;
template class initBoundValProbs<3>;
