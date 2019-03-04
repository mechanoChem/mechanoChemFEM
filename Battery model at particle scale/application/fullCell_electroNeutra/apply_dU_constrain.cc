#include"initBoundValProbs.h"

template <int dim>
void initBoundValProbs<dim>::apply_dU_constrain(PETScWrappers::MPI::Vector& dU)
{
	PETScWrappers::Vector localized_dU(this->dU);	
  for(unsigned int i=0;i<constrain_index[0].size();i++){
		localized_dU(constrain_index[0][i])=localized_dU(constrain_index[1][i]);
	}
	dU = localized_dU;
}

template class initBoundValProbs<1>;
template class initBoundValProbs<2>;
template class initBoundValProbs<3>;