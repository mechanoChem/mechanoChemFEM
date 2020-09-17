/*
zhenlin wang 2019
*/

#include"../../include/mechanoChemFEM.h"

template <int dim>
void mechanoChemFEM<dim>::apply_initial_condition()
{ 
	pcout << "applying initial condition\n";
	int totalDOF=this->totalDOF(primary_variables);	
  VectorTools::interpolate(this->dof_handler, InitialConditions<dim>(totalDOF, primary_variables, primary_variables_dof, *params_mechanoChemFEM), solution_prev); 
	
	solution_prev.compress(VectorOperation::insert);
	solution=solution_prev;
}


template class mechanoChemFEM<1>;
template class mechanoChemFEM<2>;
template class mechanoChemFEM<3>;
