/*
zhenlin wang 2019
*/

#include"../../include/mechanoChemFEM.h"
template <int dim>
void mechanoChemFEM<dim>::solve_ibvp()
{		
	this->nonlinearSolve(solution);
	//update
	solution_prev=solution;
}
template class mechanoChemFEM<1>;
template class mechanoChemFEM<2>;
template class mechanoChemFEM<3>;