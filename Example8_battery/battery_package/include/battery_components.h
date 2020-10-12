/*
zhenlin wang 2019
*battery components
*/

#ifndef battery_components_h
#define battery_components_h

#include "mechanoChemFEM.h"
#include "Transportation.h"
#include "PoissonEquation.h"
#include "Mechanics.h"

template <class T>
class ElectricChemo
{
public:
  ElectricChemo ();
	
	dealii::ParameterHandler* params;
	void declare_parameters(nlohmann::json& _params);
	void init();	
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

	nlohmann::json* params_ElectricChemo_json;
	
	double F, Rr;
	double k_neg, k_pos, alpha_neg, alpha_pos, c_max_neg, c_max_pos;
	
};

template <int dim>
class Lithium:public Transportation<dim>
{
	public:
		ParameterHandler* params;
		nlohmann::json* params_json;
		void declare_parameters(nlohmann::json& _params);
		void set_diffusion_term(dealii::Table<2,Sacado::Fad::DFad<double> >& diffu);
};
/*
**
*/
template <int dim>
class Lithium_phaseField:public PoissonEquation<dim>
{
	public:
		ParameterHandler* params;
		nlohmann::json* params_json;
		void declare_parameters(nlohmann::json& _params);
		void set_field_and_source_term(dealii::Table<2,Sacado::Fad::DFad<double> >& field, dealii::Table<1,Sacado::Fad::DFad<double> >& source);
};

//****************
//****************

template <int dim>
class Lithium_cation:public Transportation<dim>
{
	public:
		ParameterHandler* params;
		nlohmann::json* params_json;
		void declare_parameters(nlohmann::json& _params);
		void set_diffusion_term(dealii::Table<2,Sacado::Fad::DFad<double> >& diffu);
};

//****************
//****************

template <int dim>
class Electrode_potential:public PoissonEquation<dim>
{
	public:
		ParameterHandler* params;
		nlohmann::json* params_json;
		void declare_parameters(nlohmann::json& _params);
		void set_field_and_source_term(dealii::Table<2,Sacado::Fad::DFad<double> >& field, dealii::Table<1,Sacado::Fad::DFad<double> >& source);
};

template <int dim>
class Electrolyte_potential:public PoissonEquation<dim>
{
	public:
		ParameterHandler* params;
		nlohmann::json* params_json;
		void declare_parameters(nlohmann::json& _params);
		void set_field_and_source_term(dealii::Table<2,Sacado::Fad::DFad<double> >& field, dealii::Table<1,Sacado::Fad::DFad<double> >& source);
};


//****************
//****************

template <int dim>
class Temperature:public Transportation<dim>
{
	public:
		ParameterHandler* params;
		nlohmann::json* params_json;
		void declare_parameters(nlohmann::json& _params);
		void set_diffusion_term(dealii::Table<2,Sacado::Fad::DFad<double> >& diffu);
};
//****************
//****************

template <int dim>
class Displacement:public Mechanics<dim>
{
	public:
		ParameterHandler* params;
		nlohmann::json* params_json;
		void declare_parameters(nlohmann::json& _params);
		void set_stress(dealii::Table<3,Sacado::Fad::DFad<double> >& F, dealii::Table<3,Sacado::Fad::DFad<double> >& P);		
};

#endif