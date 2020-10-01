/*
zhenlin wang 2019
*coupled diffusion reaction
*/
#include "../include/battery_components.h"

/*
*Lithium
*/
template <int dim>
void Lithium<dim>::declare_parameters(ParameterHandler& _params){
	params=&_params;
	params->enter_subsection("parameters_Lithium");
	params->declare_entry("lithium_diffusivity","0",Patterns::Double() );
	params->declare_entry("lithium_ini","1",Patterns::Double() );
	params->leave_subsection();	
}

template <int dim>
void Lithium<dim>::declare_parameters(nlohmann::json& _params){
	params_json=&_params;
	(*params_json)["Lithium"]["lithium_diffusivity"]=0;
	(*params_json)["Lithium"]["lithium_ini"]=0;
}

template class Lithium<1>;
template class Lithium<2>;
template class Lithium<3>;


/*
*Lithium phaseField
*/
template <int dim>
void Lithium_phaseField<dim>::declare_parameters(ParameterHandler& _params){
	params=&_params;
	params->enter_subsection("parameters_Lithium_phaseField");
	params->declare_entry("kappa","0",Patterns::Double() );
	params->declare_entry("c_alpha","0",Patterns::Double() );
	params->declare_entry("c_beta","0",Patterns::Double() );
	params->declare_entry("omega","0",Patterns::Double() );
	params->leave_subsection();
}

template <int dim>
void Lithium_phaseField<dim>::declare_parameters(nlohmann::json& _params){
	params_json=&_params;
	(*params_json)["Lithium_phaseField"]["kappa"]=0;
	(*params_json)["Lithium_phaseField"]["c_alpha"]=0;
	(*params_json)["Lithium_phaseField"]["c_beta"]=0;
	(*params_json)["Lithium_phaseField"]["omega"]=0;
}

template class  Lithium_phaseField<1>;
template class  Lithium_phaseField<2>;
template class  Lithium_phaseField<3>;


template <int dim>
void Displacement<dim>::declare_parameters(ParameterHandler& _params){
	params=&_params;
	params->enter_subsection("parameters_Mechanics");
	params->declare_entry("youngs_modulus","0",Patterns::Double() );
	params->declare_entry("poisson_ratio","0",Patterns::Double() );
	params->declare_entry("lithium_a","0",Patterns::Double() );
	params->declare_entry("lithium_b","1",Patterns::Double() );
	
  params->declare_entry("Feiga_11","1",Patterns::Double() );
  params->declare_entry("Feiga_22","1",Patterns::Double() );
  params->declare_entry("Feiga_33","1",Patterns::Double() );

  params->declare_entry("Feigb_11","1",Patterns::Double() );
  params->declare_entry("Feigb_22","1",Patterns::Double() );
  params->declare_entry("Feigb_33","1",Patterns::Double() );
	params->leave_subsection();	
}

template <int dim>
void Displacement<dim>::declare_parameters(nlohmann::json& _params){
	params_json=&_params;
	(*params_json)["Mechanics"]["youngs_modulus"]=0;
	(*params_json)["Mechanics"]["poisson_ratio"]=0;
	(*params_json)["Mechanics"]["lithium_a"]=0;
	(*params_json)["Mechanics"]["lithium_b"]=1;
	(*params_json)["Mechanics"]["Feiga_11"]=0;
	(*params_json)["Mechanics"]["Feiga_22"]=0;
	(*params_json)["Mechanics"]["Feiga_33"]=0;
	(*params_json)["Mechanics"]["Feigb_11"]=0;
	(*params_json)["Mechanics"]["Feigb_22"]=0;
	(*params_json)["Mechanics"]["Feigb_33"]=0;
}
template class Displacement<1>;
template class Displacement<2>;
template class Displacement<3>;