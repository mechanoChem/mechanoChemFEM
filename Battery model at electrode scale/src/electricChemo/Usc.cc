#include"electricChemo.h"

template <class T,int dim>
T ElectricChemo<T,dim>::formula_Usc(T x, int domainflag)
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


template class ElectricChemo<Sacado::Fad::DFad<double>,1>;
template class ElectricChemo<Sacado::Fad::DFad<double>,2>;
template class ElectricChemo<Sacado::Fad::DFad<double>,3>;

template class ElectricChemo<double,1>;
template class ElectricChemo<double,2>;
template class ElectricChemo<double,3>;