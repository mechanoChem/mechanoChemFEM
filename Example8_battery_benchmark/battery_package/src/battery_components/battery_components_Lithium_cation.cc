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
	double Rr=(*params_json)["ElectroChemo"]["Rr"];
	double t_0=(*params_json)["ElectroChemo"]["t_0"];
	double F=(*params_json)["ElectroChemo"]["F"];
	double Temp=(*params_json)["ElectroChemo"]["T_0"];
	Sacado::Fad::DFad<double> Temp_s=Temp;

	unsigned int n_q_points= diffu.size(0);
	int phi_e_index=this->battery_fields->active_fields_index["Electrolyte_potential"];
	int c_li_plus_index=this->battery_fields->active_fields_index["Lithium_cation"];
	dealii::Table<2,Sacado::Fad::DFad<double>> phi_e_grad=this->battery_fields->quad_fields[phi_e_index].value_grad;
	dealii::Table<1,Sacado::Fad::DFad<double>> c_li_plus=this->battery_fields->quad_fields[c_li_plus_index].value;
	dealii::Table<1,double> c_li_plus_prev=this->battery_fields->quad_fields[c_li_plus_index].value_conv;
	
	dealii::Table<2,Sacado::Fad::DFad<double>> c_li_plus_grad=this->battery_fields->quad_fields[c_li_plus_index].value_grad;

	dealii::Table<1,double > D_li_plus(n_q_points);
	//dealii::Table<1,Sacado::Fad::DFad<double> > D_li_plus(n_q_points);
	dealii::Table<1,Sacado::Fad::DFad<double>>sigma_e(n_q_points);

	dealii::Table<2,Sacado::Fad::DFad<double> > i_phi_e(n_q_points, dim);
	for(unsigned int q=0; q<n_q_points;q++){
		Sacado::Fad::DFad<double> c_li_plus_prev_val=c_li_plus_prev[q];
		D_li_plus[q]=this->electricChemoFormula->D_li_plus(c_li_plus_prev_val, Temp_s).val()/1000;
		sigma_e[q]=this->electricChemoFormula->sigma_e(c_li_plus_prev_val, Temp_s).val();
		for(unsigned int i=0; i<dim;i++){
			i_phi_e[q][i]=-sigma_e[q]*phi_e_grad[q][i]-2*Rr*Temp/F*sigma_e[q]*(1-t_0)/c_li_plus[q]*c_li_plus_grad[q][i];
			diffu[q][i]=-D_li_plus[q]*c_li_plus_grad[q][i]+t_0/F*i_phi_e[q][i];
		}
	}	
}

template class Lithium_cation<1>;
template class Lithium_cation<2>;
template class Lithium_cation<3>;



