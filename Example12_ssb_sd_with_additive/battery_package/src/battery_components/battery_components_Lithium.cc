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
void Lithium<dim>::set_diffusion_reaction_term(dealii::Table<2,Sacado::Fad::DFad<double> >& diffu, dealii::Table<1,Sacado::Fad::DFad<double> >& react, std::vector<double> &pressure_cell)
{
	
	double M, V;
  double k_b = 1.38e-23; // J/K
  //double T_room = 300; // K
	double Temp=(*params_json)["ElectroChemo"]["T_0"];
  // pressure 0.3GPa * V_si = 0.3e9 N/m2 * 1.6e-28
  // P * V / k_b /T_room = 4.8e-20/300/1.38e-23 =  4.8/300/1.38e-3 = 11.59
  // exp(-11.59) = 9.0e-6
  // exp(11.59) = 1.1e5
	Point<dim> center=this->battery_fields->cell_center;
	double separator_line=(*params_json)["ElectroChemo"]["separator_line"];
	int orientation=(*params_json)["ElectroChemo"]["orientation"];
	if (center[orientation]<separator_line){
		M=(*params_json)["ElectroChemo"]["D_li_neg"];
    V=(*params_json)["ElectroChemo"]["V_li_neg"];// Si-Si: 0.543nm*3=1.6e-28m^3, e-20 to offset GPa
	}
	else{
		M=(*params_json)["ElectroChemo"]["D_li_pos"];
    V=(*params_json)["ElectroChemo"]["V_li_pos"];
	}

	unsigned int n_q_points= react.size(0);

	double jn_react=(*params_json)["ElectroChemo"]["jn_react"];
	double eps_0=1.0e-5;
	int interface_index=this->battery_fields->active_fields_index["Diffuse_interface"];
	
	int phaesField_index=this->battery_fields->active_fields_index["Lithium_phaseField"];
  // in the first time step, pressure_old is not initialized
  if (pressure_cell.size() == 0) pressure_cell.resize(n_q_points);
  //std::cout << "M " << M << std::endl;
	if(phaesField_index==-1) 
  {
    for (unsigned int q=0; q<n_q_points; q++)
    {
      for (unsigned int i=0; i<dim; i++)
      {
        diffu[q][i] = this->battery_fields->quad_fields[this->primiary_dof].value_grad[q][i] * (-M) * exp(pressure_cell[q]*V/k_b/Temp);
        //std::cout << exp(pressure_cell[q]*V/k_b/Temp) << std::endl;
      }
    }
    //diffu=table_scaling<2,Sacado::Fad::DFad<double>,double >(this->battery_fields->quad_fields[this->primiary_dof].value_grad,-M);
  }
	else{
		diffu=table_scaling<2,Sacado::Fad::DFad<double>,double >(this->battery_fields->quad_fields[phaesField_index].value_grad,-M);
	}
  //std::cout << " ---------- in lithium set diffusion reaction " << std::endl;
}

template <int dim>
void Lithium<dim>::set_diffusion_reaction_term_interface(dealii::Table<2,Sacado::Fad::DFad<double> >& diffu, dealii::Table<1,Sacado::Fad::DFad<double> >& react, dealii::Table<2, Sacado::Fad::DFad<double>> &grad, std::vector<double> &pressure_cell)
{
	double M, V;
  double k_b = 1.38e-23; // J/K
  //double T_room = 300; // K
	double Temp=(*params_json)["ElectroChemo"]["T_0"];
  // pressure 0.3GPa * V_si = 0.3e9 N/m2 * 1.6e-28
  // P * V / k_b /T_room = 4.8e-20/300/1.38e-23 =  4.8/300/1.38e-3 = 11.59
  // exp(-11.59) = 9.0e-6
  // exp(11.59) = 1.1e5
	Point<dim> center=this->battery_fields->cell_center;
	double separator_line=(*params_json)["ElectroChemo"]["separator_line"];
	int orientation=(*params_json)["ElectroChemo"]["orientation"];
	if (center[orientation]<separator_line){
		M=(*params_json)["ElectroChemo"]["D_li_neg"];
    V=(*params_json)["ElectroChemo"]["V_li_neg"];// Si-Si: 0.543nm*3=1.6e-28m^3, e-20 to offset GPa
	}
	else{
		M=(*params_json)["ElectroChemo"]["D_li_pos"];
    V=(*params_json)["ElectroChemo"]["V_li_pos"];
	}
	unsigned int n_q_points= react.size(0);

	double jn_react=(*params_json)["ElectroChemo"]["jn_react"];
	double eps_0=1.0e-5;
	int interface_index=this->battery_fields->active_fields_index["Diffuse_interface"];
	
	int phaesField_index=this->battery_fields->active_fields_index["Lithium_phaseField"];
  // in the first time step, pressure_old is not initialized
  if (pressure_cell.size() == 0) pressure_cell.resize(n_q_points);

	if(phaesField_index==-1)
  {

    for (unsigned int q=0; q<n_q_points; q++)
    {
      for (unsigned int i=0; i<dim; i++)
      {
        diffu[q][i] = grad[q][i] * (-M) * exp(pressure_cell[q]*V/k_b/Temp);
        //if (abs(exp(pressure_cell[q]*V/k_b/Temp) ) < 0.1 or abs(exp(pressure_cell[q]*V/k_b/Temp) ) > 5.0)  std::cout << " electrode: " << exp(pressure_cell[q]*V/k_b/Temp) << std::endl;
      }
    }
    //diffu=table_scaling<2,Sacado::Fad::DFad<double>,double >(grad,-M);
  }
	else{
		diffu=table_scaling<2,Sacado::Fad::DFad<double>,double >(this->battery_fields->quad_fields[phaesField_index].value_grad,-M);
    std::cout << " ---------- in lithium set diffusion reaction interface, phaseField_index is not tested ----- " << std::endl;
    exit(0);
	}
  //std::cout << " ---------- in lithium set diffusion reaction interface " << std::endl;
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
	int Lithium_index=this->battery_fields->active_fields_index["Lithium"];
	if (Lithium_index==-1){std::cout<<"Lithium_phaseField is defined, but Lithium is not defined"<<std::endl; exit(-1);}
	double kappa=(*params_json)["Lithium_phaseField"]["kappa"];
	double C_alpha=(*params_json)["Lithium_phaseField"]["c_alpha"];
	double C_beta=(*params_json)["Lithium_phaseField"]["c_beta"];
	double omega=(*params_json)["Lithium_phaseField"]["omega"];
	
	field=table_scaling<2,Sacado::Fad::DFad<double>,double >(this->battery_fields->quad_fields[Lithium_index].value_grad,kappa);
	unsigned int n_q_points= source.size(0);
	dealii::Table<1,Sacado::Fad::DFad<double> > C_li(this->battery_fields->quad_fields[Lithium_index].value);
	for (unsigned int q=0;q<n_q_points;q++){
		source[q]=2*omega*(C_li[q]-C_alpha)*(C_li[q]-C_beta)*(2*C_li[q]-C_alpha-C_beta)-this->battery_fields->quad_fields[this->primiary_dof].value[q];
	}
}

template class  Lithium_phaseField<1>;
template class  Lithium_phaseField<2>;
template class  Lithium_phaseField<3>;

template <int dim>
void Diffuse_interface<dim>::declare_parameters(nlohmann::json& _params){
	params_json=&_params;
}
template class  Diffuse_interface<1>;
template class  Diffuse_interface<2>;
template class  Diffuse_interface<3>;
