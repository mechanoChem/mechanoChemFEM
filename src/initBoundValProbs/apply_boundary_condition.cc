/*
zhenlin wang 2019
*/
#include"../../include/mechanoChemFEM.h"

template <int dim>
void mechanoChemFEM<dim>::apply_boundary_condition()
{	
	solveClass<dim, PETScWrappers::MPI::SparseMatrix, PETScWrappers::MPI::Vector>::constraints.clear ();
  DoFTools::make_hanging_node_constraints (this->dof_handler, this->constraints);
  solveClass<dim, PETScWrappers::MPI::SparseMatrix, PETScWrappers::MPI::Vector>::constraints.close ();  
}


template class mechanoChemFEM<1>;
template class mechanoChemFEM<2>;
template class mechanoChemFEM<3>;
