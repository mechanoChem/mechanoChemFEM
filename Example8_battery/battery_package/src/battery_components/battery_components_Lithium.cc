/*
zhenlin wang 2019
*coupled diffusion reaction
*/
#include "../../include/battery_components.h"

/*
*Lithium
*/
template <int dim>
void Lithium<dim>::declare_parameters(nlohmann::json& _params){
	params_json=&_params;
}

template <int dim>
void Lithium<dim>::set_diffusion_term(dealii::Table<2,Sacado::Fad::DFad<double> >& diffu){
	
	double M=(*params_json)["ElectroChemo"]["D_li_neg"];
	//M=(*params_json);
	if(this->battery_fields->current_domain_id==this->battery_fields->pos_electrode_domain_id) M=(*params_json)["ElectroChemo"]["D_li_pos"];

	int phaesField_index=this->battery_fields->active_fields_index["Lithium_phaseField"];
	if(phaesField_index==-1) diffu=table_scaling<2,Sacado::Fad::DFad<double>,double >(this->battery_fields->quad_fields[this->primiary_dof].value_grad,-M);
	else{
		diffu=table_scaling<2,Sacado::Fad::DFad<double>,double >(this->battery_fields->quad_fields[phaesField_index].value_grad,-M);
	}
}

template class Lithium<1>;
template class Lithium<2>;
template class Lithium<3>;


/*
*Lithium phaseField
*/
template <int dim>
void Lithium_phaseField<dim>::declare_parameters(nlohmann::json& _params){
	params_json=&_params;
}

/*
*Lithium phaseField
*/
template <int dim>
void Lithium_phaseField<dim>::set_field_and_source_term(dealii::Table<2,Sacado::Fad::DFad<double> >& field, dealii::Table<1,Sacado::Fad::DFad<double> >& source)
{
	int Lithium_index=this->battery_fields->active_fields_index["Lithium_phaseField"];
	if (Lithium_index==-1){std::cout<<"Lithium_phaseField is defined, but Lithium is not defined"<<std::endl; exit(-1);}
	double kappa=(*params_json)["Lithium_phaseField"]["kappa"];
	double C_alpha=(*params_json)["Lithium_phaseField"]["c_alpha"];
	double C_beta=(*params_json)["Lithium_phaseField"]["c_beta"];
	double omega=(*params_json)["Lithium_phaseField"]["omega"];
	
	field=table_scaling<2,Sacado::Fad::DFad<double>,double >(this->battery_fields->quad_fields[Lithium_index].value_grad,-kappa);
	unsigned int n_q_points= source.size(0);
	dealii::Table<1,Sacado::Fad::DFad<double> > C_li(this->battery_fields->quad_fields[Lithium_index].value);
	for (unsigned int q=0;q<n_q_points;q++){
		source[q]=2*omega*(C_li[q]-C_alpha)*(C_li[q]-C_beta)*(2*C_li[q]-C_alpha-C_beta)-this->battery_fields->quad_fields[this->primiary_dof].value[q];
	}
}

template class  Lithium_phaseField<1>;
template class  Lithium_phaseField<2>;
template class  Lithium_phaseField<3>;

