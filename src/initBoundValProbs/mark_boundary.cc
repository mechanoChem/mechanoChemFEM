/*
zhenlin wang 2019
*/

#include"../../include/mechanoChemFEM.h"

template <int dim>
void mechanoChemFEM<dim>::mark_boundary()
{
	double X_0,Y_0,Z_0,X_end,Y_end,Z_end;
	
	if(this->use_ParameterHandler){
		params_mechanoChemFEM->enter_subsection("Geometry");
		X_0=params_mechanoChemFEM->get_double("x_min");
		Y_0=params_mechanoChemFEM->get_double("y_min");
		Z_0=params_mechanoChemFEM->get_double("z_min");
	
		X_end=params_mechanoChemFEM->get_double("x_max");
		Y_end=params_mechanoChemFEM->get_double("y_max");
		Z_end=params_mechanoChemFEM->get_double("z_max");

		params_mechanoChemFEM->leave_subsection();	
	}
	if(this->use_ParameterJson){
		X_0=(*params_mechanoChemFEM_json)["Geometry"]["x_min"];
		Y_0=(*params_mechanoChemFEM_json)["Geometry"]["y_min"];
		Z_0=(*params_mechanoChemFEM_json)["Geometry"]["z_min"];
	
		X_end=(*params_mechanoChemFEM_json)["Geometry"]["x_max"];
		Y_end=(*params_mechanoChemFEM_json)["Geometry"]["y_max"];
		Z_end=(*params_mechanoChemFEM_json)["Geometry"]["z_max"];
	}	
	
  typename  Triangulation<dim>::active_cell_iterator cell = hpFEM<dim>::triangulation.begin_active(), endc = hpFEM<dim>::triangulation.end();
  for (;cell!=endc; ++cell){
    if(dim==2){
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
			 		if (std::abs(face_center[0] - X_end)<1.0e-8){
						cell->face(f)->set_boundary_id (3); //X+
					}
			 		if (std::abs(face_center[1] - Y_end)<1.0e-8){
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

template class mechanoChemFEM<1>;
template class mechanoChemFEM<2>;
template class mechanoChemFEM<3>;
