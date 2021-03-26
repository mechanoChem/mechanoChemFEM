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
	// double M=(*params_json)["ElectroChemo"]["D_li_plus"];
	// double jn_react=(*params_json)["ElectroChemo"]["jn_react"];
	// double eps_0=1.0e-5;
	// int interface_index=this->battery_fields->active_fields_index["Diffuse_interface"];
	// diffu=table_scaling<2,Sacado::Fad::DFad<double>,double >(this->battery_fields->quad_fields[this->primiary_dof].value_grad,-M);
	
  //std::cout << " ---------- in lithium cation set diffusion reaction " << std::endl;
}

template <int dim>
void Lithium_cation<dim>::set_diffusion_reaction_term_interface(dealii::Table<2,Sacado::Fad::DFad<double> >& diffu, dealii::Table<1,Sacado::Fad::DFad<double> >& react, dealii::Table<2, Sacado::Fad::DFad<double>> &phi_e_grad, dealii::Table<2, Sacado::Fad::DFad<double>> &c_li_plus_grad, dealii::Table<1, Sacado::Fad::DFad<double>> &c_li_plus, dealii::Table<1, double> &c_li_plus_old)
{	
	double Rr=(*params_json)["ElectroChemo"]["Rr"];
	double t_0=(*params_json)["ElectroChemo"]["t_0"];
	double F=(*params_json)["ElectroChemo"]["F"];
	double Temp=(*params_json)["ElectroChemo"]["T_0"];

	unsigned int n_q_points= diffu.size(0);
	int phi_e_index=this->battery_fields->active_fields_index["Electrolyte_potential"];
	int c_li_plus_index=this->battery_fields->active_fields_index["Lithium_cation"];

	//dealii::Table<2,Sacado::Fad::DFad<double>> phi_e_grad=this->battery_fields->quad_fields[phi_e_index].value_grad;
	//dealii::Table<1,Sacado::Fad::DFad<double>> c_li_plus=this->battery_fields->quad_fields[c_li_plus_index].value;
	//dealii::Table<2,Sacado::Fad::DFad<double>> c_li_plus_grad=this->battery_fields->quad_fields[c_li_plus_index].value_grad;

	dealii::Table<1,double> D_li_plus=this->electricChemoFormula->D_li_plus_interface(c_li_plus_old);
	for(unsigned int q=0; q<n_q_points;q++){
    D_li_plus[q] = D_li_plus[q] / 1000;
  }
	dealii::Table<1,double> sigma_e = this->electricChemoFormula->sigma_e_interface(c_li_plus_old);

	dealii::Table<2,Sacado::Fad::DFad<double> > i_phi_e(n_q_points, dim);
	for(unsigned int q=0; q<n_q_points;q++){
    //std::cout << "----q---- "  << q << " sigma_e[q] " << sigma_e[q] << " D_li_plus[q] " << D_li_plus[q] << std::endl;
		for(unsigned int i=0; i<dim;i++){
			//i_phi_e[q][i]=-sigma_e[q]*phi_e_grad[q][i]-2*Rr*Temp/F*sigma_e[q]*(1-t_0)/c_li_plus[q]*c_li_plus_grad[q][i];
			i_phi_e[q][i]=-sigma_e[q]*phi_e_grad[q][i]-2*Rr*Temp/F*sigma_e[q]*(1-t_0)/0.001*c_li_plus_grad[q][i];
			diffu[q][i]=-D_li_plus[q]*c_li_plus_grad[q][i]+t_0/F*i_phi_e[q][i];
      //std::cout << " i_phi_e[q][i] " <<  q << " i " << i << " " << i_phi_e[q][i] << std::endl;
		}
	}
	// double M=(*params_json)["ElectroChemo"]["D_li_plus"];
	// double jn_react=(*params_json)["ElectroChemo"]["jn_react"];
	// double eps_0=1.0e-5;
	// int interface_index=this->battery_fields->active_fields_index["Diffuse_interface"];
	// diffu=table_scaling<2,Sacado::Fad::DFad<double>,double >(this->battery_fields->quad_fields[this->primiary_dof].value_grad,-M);
	
  //std::cout << " ---------- in lithium cation set diffusion reaction " << std::endl;
}

template class Lithium_cation<1>;
template class Lithium_cation<2>;
template class Lithium_cation<3>;



