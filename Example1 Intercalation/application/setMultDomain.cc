#include "initBoundValProbs.h"

template <int dim>
void initBoundValProbs<dim>::setMultDomain()
{
	params->enter_subsection("Problem");
	double x=params->get_double("Shell_center_X");
	double y=params->get_double("Shell_center_Y");
	double z=params->get_double("Shell_center_Z");
	double inner_radius=params->get_double("inner_radius");
	double outer_radius=params->get_double("outer_radius");
	params->leave_subsection();	
	const Point<dim> center(x,y,z);
   for (typename Triangulation<dim>::active_cell_iterator cell = this->dof_handler.begin_active(); cell != this->dof_handler.end(); ++cell){
     Point<dim> cellcenter = cell->center();
      if (cellcenter.distance(center) >= 0.1*inner_radius+0.9*outer_radius)
       cell->set_material_id(Cortex_id);
      else if(cellcenter.distance(center)>=1.05*inner_radius)
       cell->set_material_id(Subcortex_id);
      //else cell->set_material_id(Ventricle_id);
      else cell->set_material_id(Subcortex_id);   
   }
}

template class initBoundValProbs<1>;
template class initBoundValProbs<2>;
template class initBoundValProbs<3>;
