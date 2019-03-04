#include"initBoundValProbs.h"

template <int dim>
void initBoundValProbs<dim>::setMultDomain()
{
	typename Triangulation<dim>::active_cell_iterator cell = hpFEM<dim>::dof_handler.begin_active(), endc = hpFEM<dim>::dof_handler.end();
	for (;cell!=endc; ++cell){
		const Point<dim> cell_center = cell->center();
		if(cell_center[1]<=electrode_Y1 or cell_center[1]>=electrode_Y2){
			cell->set_material_id (electrode_id);
	  }
	  else cell->set_material_id (separator_id);
	}
}

template class initBoundValProbs<1>;
template class initBoundValProbs<2>;
template class initBoundValProbs<3>;
