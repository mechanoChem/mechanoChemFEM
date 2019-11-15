/*
zhenlin wang 2019
*/
#include"../../include/mechanoChemFEM.h"

template <int dim>
void mechanoChemFEM<dim>::apply_boundary_condition()
{	
	hpFEM<dim>::constraints.clear ();
  DoFTools::make_hanging_node_constraints (hpFEM<dim>::dof_handler, hpFEM<dim>::constraints);
  hpFEM<dim>::constraints.close ();  
}


template class mechanoChemFEM<1>;
template class mechanoChemFEM<2>;
template class mechanoChemFEM<3>;
