#include "initBoundValProbs.h"

template <int dim>
void initBoundValProbs<dim>::mark_boundary()
{
	params->enter_subsection("Geometry");
	double X_0=params->get_double("X_0");
	double Y_0=params->get_double("Y_0");
	double Z_0=params->get_double("Z_0");
	
	double X_end=params->get_double("X_end");
	double Y_end=params->get_double("Y_end");
	double Z_end=params->get_double("Z_end");
	params->leave_subsection();	
	
  typename  Triangulation<dim>::active_cell_iterator cell = hpFEM<dim>::triangulation.begin_active(), endc = hpFEM<dim>::triangulation.end();
  for (;cell!=endc; ++cell){
    if(dim==2){
    	for (unsigned int f=0; f<GeometryInfo<dim>::faces_per_cell; ++f){
      	if (cell->face(f)->at_boundary()){
					cell->face(f)->set_boundary_id (0);
					const Point<dim> face_center = cell->face(f)->center();
					if (std::abs(face_center[0]-X_0)<1.0e-6){
	  				cell->face(f)->set_boundary_id (1); //X-
					}
					if (std::abs(face_center[1]-Y_0)<1.0e-6){
						cell->face(f)->set_boundary_id (2); //Y-
			 		}
			 		if (std::abs(face_center[0] - X_end)<1.0e-6){
						cell->face(f)->set_boundary_id (3); //X+
					}
			 		if (std::abs(face_center[1] - Y_end)<1.0e-6){
				 	 cell->face(f)->set_boundary_id (4); //Y+
			 		}
  			}
			}
		}
		else if(dim==3){
  		for (unsigned int f=0; f<GeometryInfo<dim>::faces_per_cell; ++f){
				if (cell->face(f)->at_boundary()){
					cell->face(f)->set_boundary_id (0);
					const Point<dim> face_center = cell->face(f)->center();
					if (std::abs(face_center[0]-X_0)<1.0e-8){
	  				cell->face(f)->set_boundary_id (1); //X-
					}
					if (std::abs(face_center[1]-Y_0)<1.0e-8){
	  				cell->face(f)->set_boundary_id (2); //Y-
					}
					if (std::abs(face_center[2]-Z_0)<1.0e-8){
	  				cell->face(f)->set_boundary_id (3); //Z-
					}
					if (std::abs(face_center[0] - X_end)<1.0e-8){
	  				cell->face(f)->set_boundary_id (4); //X+
					}
					if (std::abs(face_center[1] - Y_end)<1.0e-8){
	  				cell->face(f)->set_boundary_id (5); //Y+
					}
					if (std::abs(face_center[2] - Z_end)<1.0e-8){
	  				cell->face(f)->set_boundary_id (6); //Z+
					}
      	}
    	}
		}
	}
}

template class initBoundValProbs<1>;
template class initBoundValProbs<2>;
template class initBoundValProbs<3>;
