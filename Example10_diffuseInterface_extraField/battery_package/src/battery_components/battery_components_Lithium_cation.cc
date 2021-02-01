/*
zhenlin wang 2019
*coupled diffusion reaction
*/
#include "../../include/battery_components.h"

/*
*Lithium
*/
template <int dim>
void Lithium_cation<dim>::declare_parameters(nlohmann::json& _params){
	params_json=&_params;
}

template <int dim>
void Lithium_cation<dim>::set_diffusion_reaction_term(dealii::Table<2,Sacado::Fad::DFad<double> >& diffu, dealii::Table<1,Sacado::Fad::DFad<double> >& react)
{	
	// double Rr=(*params_json)["ElectroChemo"]["Rr"];
	// double t_0=(*params_json)["ElectroChemo"]["t_0"];
	// double F=(*params_json)["ElectroChemo"]["F"];
	// double Temp=(*params_json)["ElectroChemo"]["T_0"];
	
	// unsigned int n_q_points= diffu.size(0);
	// int phi_e_index=this->battery_fields->active_fields_index["Electrolyte_potential"];
	// int c_li_plus_index=this->battery_fields->active_fields_index["Lithium_cation"];
	// dealii::Table<2,Sacado::Fad::DFad<double>> phi_e_grad=this->battery_fields->quad_fields[phi_e_index].value_grad;
	// dealii::Table<1,Sacado::Fad::DFad<double>> c_li_plus=this->battery_fields->quad_fields[c_li_plus_index].value;
	// dealii::Table<2,Sacado::Fad::DFad<double>> c_li_plus_grad=this->battery_fields->quad_fields[c_li_plus_index].value_grad;
	//
	// dealii::Table<1,Sacado::Fad::DFad<double> > D_li_plus=this->electricChemoFormula->D_li_plus();
	// dealii::Table<1,Sacado::Fad::DFad<double>> sigma_e = this->electricChemoFormula->sigma_e();
	//
	// dealii::Table<2,Sacado::Fad::DFad<double> > i_phi_e(n_q_points, dim);
	// for(unsigned int q=0; q<n_q_points;q++){
	// 	for(unsigned int i=0; i<dim;i++){
	// 		i_phi_e[q][i]-sigma_e[q]*phi_e_grad[q][i]+2*Rr*Temp/F*sigma_e[q]*(1-t_0)/c_li_plus[q]*c_li_plus_grad[q][i];
	// 		diffu[q][i]=-D_li_plus[q]*c_li_plus_grad[q][i]+t_0/F*i_phi_e[q][i];
	// 	}
	// }
	double M=(*params_json)["ElectroChemo"]["D_li_plus"];
	double jn_react=(*params_json)["ElectroChemo"]["jn_react"];
	double eps_0=1.0e-5;
	int interface_index=this->battery_fields->active_fields_index["Diffuse_interface"];
	if(this->battery_fields->quad_fields[interface_index].value[0]>=1-eps_0 ) M=0;	 
	else{	 
		unsigned int n_q_points= react.size(0);
		//for (unsigned int q=0; q<n_q_points; ++q) react[q]=-jn_react;
	}
	
	int phaesField_index=this->battery_fields->active_fields_index["Lithium_phaseField"];
	diffu=table_scaling<2,Sacado::Fad::DFad<double>,double >(this->battery_fields->quad_fields[this->primiary_dof].value_grad,-M);
}

template class Lithium_cation<1>;
template class Lithium_cation<2>;
template class Lithium_cation<3>;



