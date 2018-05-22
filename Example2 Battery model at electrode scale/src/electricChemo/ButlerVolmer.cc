#include"electricChemo.h"

template <class T,int dim>
T ElectricChemo<T,dim>::formula_j0(T c_li, T c_li_plus, int domainflag)
{
	double k, alpha, c1max;
	T UnitC_surface, j0;
	if(domainflag==1) {k=k_pos;alpha=alpha_pos;c1max=c_max_pos; }
	else if(domainflag==-1) {k=k_neg;alpha=alpha_neg;c1max=c_max_neg; }
	UnitC_surface=c_li/c1max;
  j0=k*c1max*std::pow(c_li_plus,alpha)*std::pow((1-UnitC_surface),alpha)*std::pow(UnitC_surface,(1-alpha));
	return j0;
}

template <class T,int dim>
T ElectricChemo<T,dim>::formula_jn(T Temp, T c_li, T c_li_plus, T phi_s, T phi_e, int domainflag)
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


template class ElectricChemo<Sacado::Fad::DFad<double>,1>;
template class ElectricChemo<Sacado::Fad::DFad<double>,2>;
template class ElectricChemo<Sacado::Fad::DFad<double>,3>;

template class ElectricChemo<double,1>;
template class ElectricChemo<double,2>;
template class ElectricChemo<double,3>;