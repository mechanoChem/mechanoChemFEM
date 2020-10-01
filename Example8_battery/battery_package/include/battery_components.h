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

template <int dim>
class Lithium:public Transportation<dim>
{
	public:
		ParameterHandler* params;
		nlohmann::json* params_json;
		void declare_parameters(ParameterHandler& _params);
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
		void declare_parameters(ParameterHandler& _params);
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
		void declare_parameters(ParameterHandler& _params);
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
		void declare_parameters(ParameterHandler& _params);
		void declare_parameters(nlohmann::json& _params);
		void set_current(dealii::Table<2,Sacado::Fad::DFad<double> >& current);
};

//****************
//****************

template <int dim>
class Electrolyte_potential:public PoissonEquation<dim>
{
	public:
		ParameterHandler* params;
		nlohmann::json* params_json;
		void declare_parameters(ParameterHandler& _params);
		void declare_parameters(nlohmann::json& _params);
		void set_current(dealii::Table<2,Sacado::Fad::DFad<double> >& current);
};

//****************
//****************

template <int dim>
class Temperature:public Transportation<dim>
{
	public:
		ParameterHandler* params;
		nlohmann::json* params_json;
		void declare_parameters(ParameterHandler& _params);
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
		void declare_parameters(ParameterHandler& _params);
		void declare_parameters(nlohmann::json& _params);
		void set_stress(dealii::Table<3,Sacado::Fad::DFad<double> >& F, dealii::Table<3,Sacado::Fad::DFad<double> >& P);		
};

#endif