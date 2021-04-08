/*
zhenlin wang 2020
*Electrode_potential
*/
#include "../../include/battery_components.h"


template <int dim>
void Electrode_potential<dim>::declare_parameters(nlohmann::json& _params){
	params_json=&_params;
	(*params_json)["ElectroChemo"]["sigma_s_neg"]=0;
	(*params_json)["ElectroChemo"]["sigma_s_pos"]=0;
}

template <int dim>
void Electrode_potential<dim>::set_field_and_source_term(dealii::Table<2,Sacado::Fad::DFad<double> >& field, dealii::Table<1,Sacado::Fad::DFad<double> >& source)
{
	double sigma=(*params_json)["ElectroChemo"]["sigma_s_neg"];
	
	field=table_scaling<2,Sacado::Fad::DFad<double> > (this->battery_fields->quad_fields[this->primiary_dof].value_grad,-sigma);
    //std::cout << " electrode field bulk " << field[0][0] << std::endl;
  //std::cout << " ---------- in electrode set field and source " << std::endl;
}


template <int dim>
void Electrode_potential<dim>::set_field_and_source_term_interface(dealii::Table<2,Sacado::Fad::DFad<double> >& field, dealii::Table<1,Sacado::Fad::DFad<double> >& source, dealii::Table<2, Sacado::Fad::DFad<double>> &grad)
{
	double sigma=(*params_json)["ElectroChemo"]["sigma_s_neg"];
	
	field=table_scaling<2,Sacado::Fad::DFad<double> > (grad,-sigma);
    //std::cout << " electrode field interface " << field[0][0] << std::endl;
}

template class  Electrode_potential<1>;
template class  Electrode_potential<2>;
template class  Electrode_potential<3>;

