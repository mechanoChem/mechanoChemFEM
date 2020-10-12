#include "../include/battery_components.h"

template <class T>
ElectricChemo<T>::ElectricChemo(){}

template <class T>
void ElectricChemo<T>::declare_parameters(nlohmann::json& _params)
{
	params_ElectricChemo_json=&_params;
	(*params_ElectricChemo_json)["ElectroChemo"]["F"]=96485.3329;
	(*params_ElectricChemo_json)["ElectroChemo"]["Rr"]=8.3144598;
	(*params_ElectricChemo_json)["ElectroChemo"]["alpha_a"]=0.5;
	(*params_ElectricChemo_json)["ElectroChemo"]["alpha_c"]=0.5;
  (*params_ElectricChemo_json)["ElectroChemo"]["k_neg"]=0.8;
	(*params_ElectricChemo_json)["ElectroChemo"]["k_pos"]=0.8;
  (*params_ElectricChemo_json)["ElectroChemo"]["c_max_neg"]=28.7e-3;
	(*params_ElectricChemo_json)["ElectroChemo"]["c_max_pos"]=37.5e-3;
}

template <class T>
void ElectricChemo<T>::init()
{
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
template <class T>
T ElectricChemo<T>::formula_conductivity_e(T Temp, T c_li_plus, int type)
{
	T Ke;
	if(type==1) Ke=c_li_plus*std::pow(-10.5+0.074*Temp-6.96e-5*Temp*Temp+668*c_li_plus-17.8*c_li_plus*Temp+0.028*c_li_plus*Temp*Temp+4.94e5*c_li_plus*c_li_plus-886*c_li_plus*c_li_plus*Temp,2.0)*1.0e8;
	else {std::cout<<"wrong type for Ke"<<std::endl; exit(-1);}
	return Ke;
}

template <class T>
T ElectricChemo<T>::formula_diffusivity_e(T Temp, T c_li_plus, int type)
{
	T D_l;
	if(type==1) D_l=std::pow(10,(-4.43-54/(Temp-229-5e3*c_li_plus)-2.2e2*c_li_plus))*1.0e8;
	//else if(type==2) D_l=std::pow(10,(-4.43-54/(this_Temp-299-5e3*this_c_li_plus)-2.2e2*this_c_li_plus))*1.0e8;
	else {std::cout<<"wrong type for D_l"<<std::endl; exit(-1);}
	return D_l;
}

template <class T>
T ElectricChemo<T>::formula_diffusivity_s(T Temp, T c_li, int type)
{
	T D_s;
	if(type==1) D_s=1.4523*1.0e-9*exp(68025.7/Rr*(1/318-1/Temp))*1.0e8;
	else {std::cout<<"wrong type for D_s"<<std::endl; exit(-1);}
	return D_s;
}

template <class T>
T ElectricChemo<T>::solid_particle_expansion(T unitC, int type)
{
	T gamma;
	if(type==1) gamma=1.496*std::pow(unitC,3)-1.739*unitC*unitC+1.02*unitC-0.03304*std::exp(2.972*unitC)-0.04587*tanh((unitC-0.1)/0.1)-0.003608*tanh((unitC-0.3)/0.1)+0.0214*tanh((un\
itC-0.65)/0.1);
	else {std::cout<<"wrong type for D_s"<<std::endl; exit(-1);}
	return gamma;
}


template <class T>
T ElectricChemo<T>::formula_dUdt(T UnitC)
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

template <class T>
T ElectricChemo<T>::formula_Usc(T x, int domainflag)
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


template <class T>
T ElectricChemo<T>::formula_j0(T c_li, T c_li_plus, int domainflag)
{
	double k, alpha, c1max;
	T UnitC_surface, j0;
	if(domainflag==1) {k=k_pos;alpha=alpha_pos;c1max=c_max_pos; }
	else if(domainflag==-1) {k=k_neg;alpha=alpha_neg;c1max=c_max_neg; }
	UnitC_surface=c_li/c1max;
  j0=k*c1max*std::pow(c_li_plus,alpha)*std::pow((1-UnitC_surface),alpha)*std::pow(UnitC_surface,(1-alpha));
	return j0;
}

template <class T>
T ElectricChemo<T>::formula_jn(T Temp, T c_li, T c_li_plus, T phi_s, T phi_e, int domainflag)
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

template class ElectricChemo<Sacado::Fad::DFad<double>>;
template class ElectricChemo<double>;