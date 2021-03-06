/*
zhenlin wang 2019
*/

#include"../../include/mechanoChemFEM.h"

template <int dim>
void mechanoChemFEM<dim>::make_grid()
{	
	
	params_mechanoChemFEM->enter_subsection("Geometry");
	double X_0=params_mechanoChemFEM->get_double("X_0");
	double Y_0=params_mechanoChemFEM->get_double("Y_0");
	double Z_0=params_mechanoChemFEM->get_double("Z_0");
	
	double X_end=params_mechanoChemFEM->get_double("X_end");
	double Y_end=params_mechanoChemFEM->get_double("Y_end");
	double Z_end=params_mechanoChemFEM->get_double("Z_end");
	
	int element_div_x=params_mechanoChemFEM->get_double("element_div_x");
	int element_div_y=params_mechanoChemFEM->get_double("element_div_y");
	int element_div_z=params_mechanoChemFEM->get_double("element_div_z");
	params_mechanoChemFEM->leave_subsection();	
	
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
