#include "initBoundValProbs.h"

template <int dim>
void initBoundValProbs<dim>::mark_boundary()
{
	params->enter_subsection("Problem");
	double x=params->get_double("Shell_center_X");
	double y=params->get_double("Shell_center_Y");
	double z=params->get_double("Shell_center_Z");
	double inner_radius=params->get_double("inner_radius");
	double outer_radius=params->get_double("outer_radius");
	params->leave_subsection();	
	const Point<dim> center(x,y,z);
	
  typename  Triangulation<dim>::active_cell_iterator cell = hpFEM<dim>::triangulation.begin_active(), endc = hpFEM<dim>::triangulation.end();
  for (;cell!=endc; ++cell){
		for (unsigned int f=0; f<GeometryInfo<dim>::faces_per_cell; ++f){
			if (cell->face(f)->at_boundary()){
			  const Point<dim> face_center = cell->face(f)->center();
      
			  if (std::abs(face_center[0])<1.0e-6){
			    cell->face(f)->set_boundary_id (1); //X-                                                                                            
			  }
			  else if(face_center.distance(center)<inner_radius+2) cell->face(f)->set_boundary_id (2);
			  else cell->face(f)->set_boundary_id (0);
		        }
		}
	}
}

template class initBoundValProbs<1>;
template class initBoundValProbs<2>;
template class initBoundValProbs<3>;
