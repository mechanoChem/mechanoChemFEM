/*
zhenlin wang 2020
*Electrode_potential
*/
#include "../../include/battery_components.h"


template <int dim>
void Electrolyte_potential<dim>::declare_parameters(nlohmann::json& _params){
	params_json=&_params;
}

template <int dim>
void Electrolyte_potential<dim>::set_field_and_source_term(dealii::Table<2,Sacado::Fad::DFad<double> >& field, dealii::Table<1,Sacado::Fad::DFad<double> >& source)
{
	double Rr=(*params_json)["ElectroChemo"]["Rr"];
	double t_0=(*params_json)["ElectroChemo"]["t_0"];
	double F=(*params_json)["ElectroChemo"]["F"];
	double Temp=(*params_json)["ElectroChemo"]["T_0"];
	unsigned int n_q_points= source.size(0);
	
	dealii::Table<1,Sacado::Fad::DFad<double>>sigma_e = this->electricChemoFormula->sigma_e();
	int phi_e_index=this->battery_fields->active_fields_index["Electrolyte_potential"];
	int c_li_plus_index=this->battery_fields->active_fields_index["Lithium_cation"];
	dealii::Table<2,Sacado::Fad::DFad<double>> phi_e_grad=this->battery_fields->quad_fields[phi_e_index].value_grad;
	dealii::Table<1,Sacado::Fad::DFad<double>> c_li_plus=this->battery_fields->quad_fields[c_li_plus_index].value;
	dealii::Table<2,Sacado::Fad::DFad<double>> c_li_plus_grad=this->battery_fields->quad_fields[c_li_plus_index].value_grad;
	
	for(unsigned int q=0; q<n_q_points;q++){
		for(unsigned int i=0; i<dim;i++){
			field[q][i]=-sigma_e[q]*phi_e_grad[q][i]+2*Rr*Temp/F*sigma_e[q]*(1-t_0)/c_li_plus[q]*c_li_plus_grad[q][i]; 
		}
	}
	table_scaling<1,Sacado::Fad::DFad<double> > (source,0);

  //std::cout << " ---------- in electrolyte set field and source " << std::endl;
}

template <int dim>
void Electrolyte_potential<dim>::set_field_and_source_term_interface(dealii::Table<2,Sacado::Fad::DFad<double> >& field, dealii::Table<1,Sacado::Fad::DFad<double> >& source, dealii::Table<2, Sacado::Fad::DFad<double>> &phi_e_grad, dealii::Table<2, Sacado::Fad::DFad<double>> &c_li_plus_grad, dealii::Table<1, Sacado::Fad::DFad<double>> &c_li_plus, dealii::Table<1, double> &c_li_plus_old)
{
	double Rr=(*params_json)["ElectroChemo"]["Rr"];
	double t_0=(*params_json)["ElectroChemo"]["t_0"];
	double F=(*params_json)["ElectroChemo"]["F"];
	double Temp=(*params_json)["ElectroChemo"]["T_0"];
	unsigned int n_q_points= source.size(0);
	
	dealii::Table<1,double>sigma_e = this->electricChemoFormula->sigma_e_interface(c_li_plus_old);
	int phi_e_index=this->battery_fields->active_fields_index["Electrolyte_potential"];
	int c_li_plus_index=this->battery_fields->active_fields_index["Lithium_cation"];
	//dealii::Table<2,Sacado::Fad::DFad<double>> phi_e_grad=this->battery_fields->quad_fields[phi_e_index].value_grad;
	//dealii::Table<1,Sacado::Fad::DFad<double>> c_li_plus=this->battery_fields->quad_fields[c_li_plus_index].value;
	//dealii::Table<2,Sacado::Fad::DFad<double>> c_li_plus_grad=this->battery_fields->quad_fields[c_li_plus_index].value_grad;
	
	for(unsigned int q=0; q<n_q_points;q++){
    //std::cout << "sigma_e " << q << " "  << sigma_e[q] << std::endl;
    //std::cout << "c_li_plus " << q << " "  << c_li_plus[q] << " homo " << this->battery_fields->quad_fields[c_li_plus_index].value[q]<< std::endl;
		for(unsigned int i=0; i<dim;i++){
			field[q][i]=-sigma_e[q]*phi_e_grad[q][i]+2*Rr*Temp/F*sigma_e[q]*(1-t_0)/c_li_plus[q]*c_li_plus_grad[q][i]; 
      //std::cout << "phi_e_grad " << q << " i " << i << " "  << phi_e_grad[q][i] << " homo " << this->battery_fields->quad_fields[phi_e_index].value_grad[q][i] << std::endl;
      //std::cout << "c_li_plus_grad " << q << " i " << i << " "  << c_li_plus_grad[q][i] << " homo " << this->battery_fields->quad_fields[c_li_plus_index].value_grad[q][i] << std::endl;
		}
	}
	table_scaling<1,Sacado::Fad::DFad<double> > (source,0);

  //std::cout << " ---------- in electrolyte set field and source " << std::endl;
}

template class  Electrolyte_potential<1>;
template class  Electrolyte_potential<2>;
template class  Electrolyte_potential<3>;

