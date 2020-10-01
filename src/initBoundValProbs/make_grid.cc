/*
zhenlin wang 2019
*/

#include"../../include/mechanoChemFEM.h"

template <int dim>
void mechanoChemFEM<dim>::make_grid()
{	
	double X_0,Y_0,Z_0,X_end,Y_end,Z_end;
	int element_div_x, element_div_y,element_div_z;
	
	if(this->use_ParameterHandler){
		params_mechanoChemFEM->enter_subsection("Geometry");
		X_0=params_mechanoChemFEM->get_double("x_min");
		Y_0=params_mechanoChemFEM->get_double("y_min");
		Z_0=params_mechanoChemFEM->get_double("z_min");
	
		X_end=params_mechanoChemFEM->get_double("x_max");
		Y_end=params_mechanoChemFEM->get_double("y_max");
		Z_end=params_mechanoChemFEM->get_double("z_max");
	
		element_div_x=params_mechanoChemFEM->get_integer("num_elem_x");
		element_div_y=params_mechanoChemFEM->get_integer("num_elem_y");
		element_div_z=params_mechanoChemFEM->get_integer("num_elem_z");
		params_mechanoChemFEM->leave_subsection();	
	}
	if(this->use_ParameterJson){
		X_0=(*params_mechanoChemFEM_json)["Geometry"]["x_min"];
		Y_0=(*params_mechanoChemFEM_json)["Geometry"]["y_min"];
		Z_0=(*params_mechanoChemFEM_json)["Geometry"]["z_min"];
	
		X_end=(*params_mechanoChemFEM_json)["Geometry"]["x_max"];
		Y_end=(*params_mechanoChemFEM_json)["Geometry"]["y_max"];
		Z_end=(*params_mechanoChemFEM_json)["Geometry"]["z_max"];;
	
		element_div_x=(*params_mechanoChemFEM_json)["Geometry"]["num_elem_x"].get<int>();
		element_div_y=(*params_mechanoChemFEM_json)["Geometry"]["num_elem_y"].get<int>();
		element_div_z=(*params_mechanoChemFEM_json)["Geometry"]["num_elem_z"].get<int>();
	}
	
	bool colorize = false;
  std::vector< std::vector< double > > step_sizes;
  step_sizes.resize(dim);
	
  for (unsigned int j = 0; j < element_div_x; ++j) step_sizes[0].push_back((X_end-X_0)/element_div_x); 
  for (unsigned int j = 0; j < element_div_y; ++j) step_sizes[1].push_back((Y_end-Y_0)/element_div_y);
	if(dim==3)	for (unsigned int j = 0; j < element_div_z; ++j) step_sizes[2].push_back((Z_end-Z_0)/element_div_z);
  if(dim==2) GridGenerator::subdivided_hyper_rectangle (this->triangulation, step_sizes, Point<dim>(X_0,Y_0), Point<dim>(X_end,Y_end), colorize);
	else GridGenerator::subdivided_hyper_rectangle (this->triangulation, step_sizes, Point<dim>(X_0,Y_0,Z_0), Point<dim>(X_end,Y_end,Z_end), colorize);
	
}


template class mechanoChemFEM<1>;
template class mechanoChemFEM<2>;
template class mechanoChemFEM<3>;
