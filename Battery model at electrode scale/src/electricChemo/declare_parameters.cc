#include"electricChemo.h"

template <class T,int dim>
void ElectricChemo<T,dim>::declare_parameters()
{
	//declare paramter for Butler-Volmer equations
	params->enter_subsection("ElectroChemo");
	params->declare_entry("F","0",dealii::Patterns::Double() );
	params->declare_entry("Rr","0",dealii::Patterns::Double() );
	params->declare_entry("alpha_a","0",dealii::Patterns::Double() );
	params->declare_entry("alpha_c","0",dealii::Patterns::Double() );
  params->declare_entry("k_neg","0",dealii::Patterns::Double() );
	params->declare_entry("k_pos","0",dealii::Patterns::Double() );
  params->declare_entry("c_max_neg","0",dealii::Patterns::Double() );
	params->declare_entry("c_max_pos","0",dealii::Patterns::Double() );
	params->leave_subsection();
}

template <class T,int dim>
void ElectricChemo<T,dim>::setParametersFromHandler()
{
	params->enter_subsection("ElectroChemo");
	F=params->get_double("F");
	Rr=params->get_double("Rr");
	k_neg=params->get_double("k_neg");
	k_pos=params->get_double("k_pos");
	alpha_pos=params->get_double("alpha_a");	
	alpha_neg=params->get_double("alpha_c");
	c_max_neg=params->get_double("c_max_neg");
	c_max_pos=params->get_double("c_max_pos");	
	params->leave_subsection();
}
	
	
	
	
	
template class ElectricChemo<Sacado::Fad::DFad<double>,1>;
template class ElectricChemo<Sacado::Fad::DFad<double>,2>;
template class ElectricChemo<Sacado::Fad::DFad<double>,3>;

template class ElectricChemo<double,1>;
template class ElectricChemo<double,2>;
template class ElectricChemo<double,3>;