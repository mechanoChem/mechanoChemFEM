/*
zhenlin wang 2019
*/

#include"../../include/initBoundValProbs.h"
template <int dim>
void initBoundValProbs<dim>::solve()
{		
	this->nonlinearSolve(solution);
}
template class initBoundValProbs<1>;
template class initBoundValProbs<2>;
template class initBoundValProbs<3>;