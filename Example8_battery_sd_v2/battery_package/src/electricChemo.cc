#include "../include/ElectricChemo.h"

template <int dim, class T>
ElectricChemo<dim, T>::ElectricChemo(){}

template <int dim, class T>
void ElectricChemo<dim, T>::declare_parameters(nlohmann::json& _params)
{
	params_ElectricChemo_json=&_params;
	(*params_ElectricChemo_json)["ElectroChemo"]["F"]=96485.3329;
	(*params_ElectricChemo_json)["ElectroChemo"]["Rr"]=8.3144598;
	(*params_ElectricChemo_json)["ElectroChemo"]["T_0"]=340;
	(*params_ElectricChemo_json)["ElectroChemo"]["t_0"]=0.2;
	(*params_ElectricChemo_json)["ElectroChemo"]["alpha_a"]=0.5;
	(*params_ElectricChemo_json)["ElectroChemo"]["alpha_c"]=0.5;
  (*params_ElectricChemo_json)["ElectroChemo"]["k_neg"]=0.8;
	(*params_ElectricChemo_json)["ElectroChemo"]["k_pos"]=0.8;
  (*params_ElectricChemo_json)["ElectroChemo"]["c_max_neg"]=28.7e-3;
	(*params_ElectricChemo_json)["ElectroChemo"]["c_max_pos"]=37.5e-3;
	(*params_ElectricChemo_json)["ElectroChemo"]["D_li_neg"]=0;
	(*params_ElectricChemo_json)["ElectroChemo"]["D_li_pos"]=0;
	(*params_ElectricChemo_json)["ElectroChemo"]["sigma_s_neg"]=0;
	(*params_ElectricChemo_json)["ElectroChemo"]["sigma_s_pos"]=0;
	(*params_ElectricChemo_json)["ElectroChemo"]["Lithium_phaseField"]["kappa"]=0;
	(*params_ElectricChemo_json)["ElectroChemo"]["Lithium_phaseField"]["c_alpha"]=0;
	(*params_ElectricChemo_json)["ElectroChemo"]["Lithium_phaseField"]["c_beta"]=0;
	(*params_ElectricChemo_json)["ElectroChemo"]["Lithium_phaseField"]["omega"]=0;
	// (*params_ElectricChemo_json)["ElectroChemo"]["Mechanics"]["youngs_modulus"]=0;
	// (*params_ElectricChemo_json)["ElectroChemo"]["Mechanics"]["poisson_ratio"]=0;
	// (*params_ElectricChemo_json)["ElectroChemo"]["Mechanics"]["lithium_a"]=0;
	// (*params_ElectricChemo_json)["ElectroChemo"]["Mechanics"]["lithium_b"]=1;
	// (*params_ElectricChemo_json)["ElectroChemo"]["Mechanics"]["Feiga_11"]=0;
	// (*params_ElectricChemo_json)["ElectroChemo"]["Mechanics"]["Feiga_22"]=0;
	// (*params_ElectricChemo_json)["ElectroChemo"]["Mechanics"]["Feiga_33"]=0;
	// (*params_ElectricChemo_json)["ElectroChemo"]["Mechanics"]["Feigb_11"]=0;
	// (*params_ElectricChemo_json)["ElectroChemo"]["Mechanics"]["Feigb_22"]=0;
	// (*params_ElectricChemo_json)["ElectroChemo"]["Mechanics"]["Feigb_33"]=0;
}

template <int dim, class T>
void ElectricChemo<dim, T>::init(Battery_fields<dim>& _fields)
{
	battery_fields=&_fields;
	F=(*params_ElectricChemo_json)["ElectroChemo"]["F"];
	Rr=(*params_ElectricChemo_json)["ElectroChemo"]["Rr"];
	k_neg=(*params_ElectricChemo_json)["ElectroChemo"]["k_neg"];
	k_pos=(*params_ElectricChemo_json)["ElectroChemo"]["k_pos"];
	alpha_pos=(*params_ElectricChemo_json)["ElectroChemo"]["alpha_a"];	
	alpha_neg=(*params_ElectricChemo_json)["ElectroChemo"]["alpha_c"];
	c_max_neg=(*params_ElectricChemo_json)["ElectroChemo"]["c_max_neg"];
	c_max_pos=(*params_ElectricChemo_json)["ElectroChemo"]["c_max_pos"];	
}
	
/*
*type 1 Theoretical Analysis of Stresses in a Lithium Ion Cell,
Sindhuja Renganathan, Godfrey Sikha, Shriram Santhanagopalan, and Ralph E. White, Journal of The Electrochemical Society, 2010
*type 2 Multi-Domain Modeling of Lithium-Ion Batteries Encompassing Multi-Physics in Varied Length Scales,
Gi-Heon Kim, Kandler Smith, Kyu-Jin Lee, Shriram Santhanagopalan, and Ahmad Pesaran, Journal of The Electrochemical Society, 2010
*/
template <int dim, class T>
dealii::Table<1,T > ElectricChemo<dim, T>::sigma_e(int type)
{
	double Temp=(*params_ElectricChemo_json)["ElectroChemo"]["T_0"];
	int c_li_plus_index=battery_fields->active_fields_index["Lithium_cation"];
	dealii::Table<1,T > C_li_plus_q=battery_fields->quad_fields[c_li_plus_index].value;
	unsigned int n_q_points= C_li_plus_q.size(0);
	T c_li_plus;
	
	dealii::Table<1,T > _sigma_e(n_q_points);
	
	if(type==1) {
		for(unsigned int q=0; q<n_q_points;q++){
			c_li_plus=C_li_plus_q[q];
			_sigma_e[q]=c_li_plus*std::pow(-10.5+0.074*Temp-6.96e-5*Temp*Temp+668*c_li_plus-17.8*c_li_plus*Temp+0.028*c_li_plus*Temp*Temp+4.94e5*c_li_plus*c_li_plus-886*c_li_plus*c_li_plus*Temp,2.0)*1.0e8;
		}
	}
	else {std::cout<<"wrong type for Ke"<<std::endl; exit(-1);}
	return _sigma_e;
}

template <int dim, class T>
dealii::Table<1,T > ElectricChemo<dim, T>::D_li_plus(int type)
{
	double Temp=(*params_ElectricChemo_json)["ElectroChemo"]["T_0"];
	
	int c_li_plus_index=battery_fields->active_fields_index["Lithium_cation"];
	dealii::Table<1,T > C_li_plus_q=battery_fields->quad_fields[c_li_plus_index].value;
	unsigned int n_q_points= C_li_plus_q.size(0);
	
	dealii::Table<1,T > _D_li_plus(n_q_points);
	T c_li_plus;
	if(type==1) {
		for(unsigned int q=0; q<n_q_points;q++){
			c_li_plus=C_li_plus_q[q];
			_D_li_plus[q]=std::pow(10,(-4.43-54/(Temp-229-5e3*c_li_plus)-2.2e2*c_li_plus))*1.0e8;
		}
	}
	else {std::cout<<"wrong type for D_l"<<std::endl; exit(-1);}
	return _D_li_plus;
}





template <int dim, class T>
T ElectricChemo<dim, T>::formula_dUdt(T UnitC)
{
	T dUdt;
	if(UnitC<=0.2)  dUdt=0.01442*UnitC*UnitC-0.00291*UnitC-0.000138;
	else if(UnitC<=0.4) dUdt=0.00634*UnitC*UnitC*UnitC-0.006625*UnitC*UnitC+0.002635*UnitC-0.0004554;
	else if(UnitC<=0.5) dUdt=0.001059*UnitC-0.0004793;
	else if(UnitC<=0.7) dUdt=0.00025*UnitC-7.5e-5;
	else if(UnitC<=0.8) dUdt=-0.001*UnitC+0.0008;
	else if(UnitC<=0.85) dUdt=0.0333*UnitC*UnitC-0.057*UnitC+0.02427;
	else if(UnitC<=0.95) dUdt=0.002*UnitC*UnitC-0.0039*UnitC+0.00177;
  else if (UnitC<=1) dUdt=-0.0014*UnitC+0.0012;
	return dUdt;
}

template <int dim, class T>
T ElectricChemo<dim, T>::formula_Usc(T x, int domainflag)
{
	T Usc;
	if(domainflag==-1){
		double c1 =0.265697795660872, c2 =0.555104680448495, c3 =178.97991682123549, c4 =0.012357396765331, c5 =0.55727132360432;
		double c6 =0.028219099799268, c7 =0.011683704080029, c8 =0.239318990250894, c9 =0.048646992277392, c10 =0.012910660088849;
		double c11 =0.174909417938192, c12 =0.03483002163646, c13 =0.050098062010346, c14 =0.024505122678677, c15 =0.03529369961247;
		double c16 =0.011931381413342, c17 =0.12992241633878, c18 =0.019797869695897, c19 =0.152636640731331, c20 =0.030000057933125, c21 =0.022725508023415;
		
		Usc=c1+c2*exp(-c3*x)-c4*std::tanh((x-c5)/c6)-c7*std::tanh((x-c8)/c9)-c10*std::tanh((x-c11)/c12)-c13*std::tanh((x-0.99)/c14)-c15*x-c16*std::tanh((x-c17)/c18)-c19*std::tanh((x-c20)/c21);//
		//value=-0.132+1.41*std::exp(-3.52*x);
	}
	//positive electrode
	else if(domainflag==1){
	 double	b1= -0.0922859116552415, b2= -7.8680409125385697, b3= 50.072175799512607, b4= -122.28161948058685, b5= 82.985110649682696;
	 double b6= 140.2938943391359, b7= -374.73497214300698, b8= 403.2463575744942, b9= -221.19151490076541, b10= 49.3392659530526530;
	 double	b11= -0.0217591621507594, b12= -1.9006524442210881, b13= 11.726362513914014, b14= -28.784794013633256, b15= 27.542704665893613;
	 double b16= -8.6342730487746202;
		Usc = (b1+b2*x+b3*std::pow(x,2)+b4*std::pow(x,3)+b5*std::pow(x,4)+b6*std::pow(x,5)+b7*std::pow(x,6)+b8*std::pow(x,7)+b9*std::pow(x,8)+b10*std::pow(x,9))/(b11+b12*x+b13*std::pow(x,2)+b14*std::pow(x,3)+b15*std::pow(x,4)+b16*std::pow(x,5));
	}
	return Usc;
}


template <int dim, class T>
T ElectricChemo<dim, T>::formula_j0(T c_li, T c_li_plus, int domainflag)
{
	double k, alpha, c1max;
	T UnitC_surface, j0;
	if(domainflag==1) {k=k_pos;alpha=alpha_pos;c1max=c_max_pos; }
	else if(domainflag==-1) {k=k_neg;alpha=alpha_neg;c1max=c_max_neg; }
	UnitC_surface=c_li/c1max;
  j0=k*c1max*std::pow(c_li_plus,alpha)*std::pow((1-UnitC_surface),alpha)*std::pow(UnitC_surface,(1-alpha));
	return j0;
}

template <int dim, class T>
T ElectricChemo<dim, T>::formula_jn(T Temp, T c_li, T c_li_plus, T phi_s, T phi_e, int domainflag)
{
	T jn;
	if(domainflag==0)jn=0;
	else{
		double k, alpha;
		T j0, Usc, UnitC;
		if(domainflag==1) {k=k_pos; alpha=alpha_pos; UnitC=c_li/c_max_pos; }
		else if(domainflag==-1) {k=k_neg; alpha=alpha_neg; UnitC=c_li/c_max_neg; }
	
		j0=formula_j0(c_li, c_li_plus, domainflag);
		Usc=formula_Usc(UnitC, domainflag);
		
		T eta = phi_s-phi_e-Usc;
		jn=j0*(exp(alpha*F/Rr/Temp*eta)-exp(-alpha*F/Rr/Temp*eta));
	}
	return jn;
}

template class ElectricChemo<1,Sacado::Fad::DFad<double>>;
template class ElectricChemo<2,Sacado::Fad::DFad<double>>;
template class ElectricChemo<3,Sacado::Fad::DFad<double>>;
