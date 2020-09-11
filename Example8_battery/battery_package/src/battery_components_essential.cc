/*
zhenlin wang 2019
*coupled diffusion reaction
*/
#include "../include/battery_components.h"
template <int dim>
void Lithium<dim>::declare_parameters(ParameterHandler& _params){
	params=&_params;
	params->enter_subsection("parameters_Lithium");
	params->declare_entry("lithium_diffusivity","0",Patterns::Double() );
	params->declare_entry("lithium_ini","1",Patterns::Double() );
	params->leave_subsection();	
}

template class Lithium<1>;
template class Lithium<2>;
template class Lithium<3>;

template <int dim>
void Displacement<dim>::declare_parameters(ParameterHandler& _params){
	params=&_params;
	params->enter_subsection("parameters_Mechanics");
	params->declare_entry("youngsModulus","0",Patterns::Double() );
	params->declare_entry("poissonRatio","0",Patterns::Double() );
	params->declare_entry("lithium_0","1",Patterns::Double() );
	params->leave_subsection();	
}

template class Displacement<1>;
template class Displacement<2>;
template class Displacement<3>;