/*
zhenlin wang 2019
*/

#include"../../include/initBoundValProbs.h"

template <int dim>
void initBoundValProbs<dim>::apply_initial_condition()
{ 
	pcout << "applying initial condition\n";
	int totalDOF=this->totalDOF(primary_variables);	
  VectorTools::interpolate(this->dof_handler, InitialConditions<dim>(totalDOF, *params), solution_prev); 
	
	solution_prev.compress(VectorOperation::insert);
	solution=solution_prev;
}


template class initBoundValProbs<1>;
template class initBoundValProbs<2>;
template class initBoundValProbs<3>;
