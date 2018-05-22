#ifndef ElectricChemo_h
#define ElectricChemo_h

#include <deal.II/base/parameter_handler.h>
#include "supplementary/dataStruct.h"
#include <Sacado.hpp>
#include <deal.II/base/table.h>

template <class T, int dim>
class ElectricChemo
{
public:
  ElectricChemo ();
	ElectricChemo (dealii::ParameterHandler& _params);
  ~ElectricChemo();
	
	dealii::ParameterHandler* params;
	void declare_parameters();
	
	void setParametersFromHandler();
	
	/*
	*expression, formula and equations
	*/
  T formula_jn(T Temp, T c_li, T c_li_plus, T phi_s, T phi_e, int domainflag);
	T formula_j0(T c_li, T c_li_plus, int domainflag);
	T formula_Usc(T x, int domainflag);
	
	T c_li_surface_parabolic(T c_li, T jn, T D_s, T R_s_0);
	
	T formula_dUdt(T UnitC);
	
	T formula_conductivity_e(T Temp, T c_li_plus, int type);
	T formula_diffusivity_e(T Temp, T c_li_plus, int type);
	T formula_diffusivity_s(T Temp, T c_li, int type);
	
	T solid_particle_expansion(T unitC, int type);
	T electrode_expansion(T unitC, int type);

	double F, Rr;
	double k_neg, k_pos, alpha_neg, alpha_pos, c_max_neg, c_max_pos;
	
};

#endif