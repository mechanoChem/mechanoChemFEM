/*
zhenlin wang 2019
*/
#include"../../include/mechanoChemFEM.h"

template <int dim>
void mechanoChemFEM<dim>::apply_boundary_condition()
{	
	constraints_mechanoChemFEM->clear ();
  DoFTools::make_hanging_node_constraints (this->dof_handler, *constraints_mechanoChemFEM);
  constraints_mechanoChemFEM->close ();  
}


template class mechanoChemFEM<1>;
template class mechanoChemFEM<2>;
template class mechanoChemFEM<3>;
