#include "initBoundValProbs.h"

template <int dim>
void initBoundValProbs<dim>::refine_grid (){
	//nothing need here currently
}

template class initBoundValProbs<1>;
template class initBoundValProbs<2>;
template class initBoundValProbs<3>;
