#include"electricChemo.h"
/*
*type 1 Theoretical Analysis of Stresses in a Lithium Ion Cell,
Sindhuja Renganathan, Godfrey Sikha, Shriram Santhanagopalan, and Ralph E. White, Journal of The Electrochemical Society, 2010
*type 2 Multi-Domain Modeling of Lithium-Ion Batteries Encompassing Multi-Physics in Varied Length Scales,
Gi-Heon Kim, Kandler Smith, Kyu-Jin Lee, Shriram Santhanagopalan, and Ahmad Pesaran, Journal of The Electrochemical Society, 2010
*/
template <class T,int dim>
T ElectricChemo<T,dim>::formula_conductivity_e(T Temp, T c_li_plus, int type)
{
	T Ke;
	if(type==1) Ke=c_li_plus*std::pow(-10.5+0.074*Temp-6.96e-5*Temp*Temp+668*c_li_plus-17.8*c_li_plus*Temp+0.028*c_li_plus*Temp*Temp+4.94e5*c_li_plus*c_li_plus-886*c_li_plus*c_li_plus*Temp,2.0)*1.0e8;
	else if(type==2) Ke=((34.5*exp(-798/Temp)*std::pow((c_li_plus*1.0e3),3)-485*exp(-1080/Temp)*std::pow((c_li_plus*1.0e3),2)+2440*exp(-1440/Temp)*(c_li_plus*1.0e3))/10)*1.0e6;
	else {std::cout<<"wrong type for Ke"<<std::endl; exit(-1);}
	return Ke;
}

template <class T,int dim>
T ElectricChemo<T,dim>::formula_diffusivity_e(T Temp, T c_li_plus, int type)
{
	T D_l;
	if(type==1) D_l=std::pow(10,(-4.43-54/(Temp-229-5e3*c_li_plus)-2.2e2*c_li_plus))*1.0e8;
	//else if(type==2) D_l=std::pow(10,(-4.43-54/(this_Temp-299-5e3*this_c_li_plus)-2.2e2*this_c_li_plus))*1.0e8;
	else {std::cout<<"wrong type for D_l"<<std::endl; exit(-1);}
	return D_l;
}

template <class T,int dim>
T ElectricChemo<T,dim>::formula_diffusivity_s(T Temp, T c_li, int type)
{
	T D_s;
	if(type==1) D_s=1.4523*1.0e-9*exp(68025.7/Rr*(1/318-1/Temp))*1.0e8;
	else {std::cout<<"wrong type for D_s"<<std::endl; exit(-1);}
	return D_s;
}

template <class T,int dim>
T ElectricChemo<T,dim>::solid_particle_expansion(T unitC, int type)
{
	T gamma;
	if(type==1) gamma=1.496*std::pow(unitC,3)-1.739*unitC*unitC+1.02*unitC-0.03304*std::exp(2.972*unitC)-0.04587*tanh((unitC-0.1)/0.1)-0.003608*tanh((unitC-0.3)/0.1)+0.0214*tanh((un\
itC-0.65)/0.1);
	else {std::cout<<"wrong type for D_s"<<std::endl; exit(-1);}
	return gamma;
}

template <class T,int dim>
T ElectricChemo<T,dim>::electrode_expansion(T unitC, int type)
{
	T gamma;
	if(type==1) gamma=0.1315*std::pow(unitC,4)-0.1989*std::pow(unitC,3)+0.06481*unitC*unitC+0.02793*unitC-0.000655;
	else {std::cout<<"wrong type for D_s"<<std::endl; exit(-1);}
	return gamma;
}


template class ElectricChemo<Sacado::Fad::DFad<double>,1>;
template class ElectricChemo<Sacado::Fad::DFad<double>,2>;
template class ElectricChemo<Sacado::Fad::DFad<double>,3>;

template class ElectricChemo<double,1>;
template class ElectricChemo<double,2>;
template class ElectricChemo<double,3>;