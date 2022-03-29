/*
zhenlin wang 2020
*battery components
*/

#ifndef battery_components_h
#define battery_components_h

#include "mechanoChemFEM.h"
#include "Transportation.h"
#include "PoissonEquation.h"
#include "Mechanics.h"
#include "ElectricChemo.h"
/*
**
*/
template <int dim>
class Lithium:public Transportation<dim>
{
	public:
		nlohmann::json* params_json;
		void declare_parameters(nlohmann::json& _params);
		void set_diffusion_reaction_term(dealii::Table<2,Sacado::Fad::DFad<double> >& diffu, dealii::Table<1,Sacado::Fad::DFad<double> >& react, std::vector<double> &pressure_cell);
		void set_diffusion_reaction_term_interface(dealii::Table<2,Sacado::Fad::DFad<double> >& diffu, dealii::Table<1,Sacado::Fad::DFad<double> >& react, dealii::Table<2, Sacado::Fad::DFad<double>> &grad, std::vector<double> &pressure_cell);
};
/*
**
*/
template <int dim>
class Lithium_phaseField:public PoissonEquation<dim>
{
	public:
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
		nlohmann::json* params_json;
		void declare_parameters(nlohmann::json& _params);
		void set_diffusion_reaction_term(dealii::Table<2,Sacado::Fad::DFad<double> >& diffu, dealii::Table<1,Sacado::Fad::DFad<double> >& react, std::vector<double> &pressure_cell);
		void set_diffusion_reaction_term_interface(dealii::Table<2,Sacado::Fad::DFad<double> >& diffu, dealii::Table<1,Sacado::Fad::DFad<double> >& react, dealii::Table<2, Sacado::Fad::DFad<double>> &phi_e_grad, dealii::Table<2, Sacado::Fad::DFad<double>> &c_li_plus_grad, dealii::Table<1, Sacado::Fad::DFad<double>> &c_li_plus, dealii::Table<1, double> &c_li_plus_old, std::vector<double> &pressure_cell);
};

//****************
//****************

template <int dim>
class Electrode_potential:public PoissonEquation<dim>
{
	public:
		nlohmann::json* params_json;
		void declare_parameters(nlohmann::json& _params);
		void set_field_and_source_term(dealii::Table<2,Sacado::Fad::DFad<double> >& field, dealii::Table<1,Sacado::Fad::DFad<double> >& source);
		void set_field_and_source_term_interface(dealii::Table<2,Sacado::Fad::DFad<double> >& field, dealii::Table<1,Sacado::Fad::DFad<double> >& source, dealii::Table<2, Sacado::Fad::DFad<double>> &grad);
};

template <int dim>
class Electrolyte_potential:public PoissonEquation<dim>
{
	public:
		nlohmann::json* params_json;
		void declare_parameters(nlohmann::json& _params);
		void set_field_and_source_term(dealii::Table<2,Sacado::Fad::DFad<double> >& field, dealii::Table<1,Sacado::Fad::DFad<double> >& source);
    void set_field_and_source_term_interface(dealii::Table<2,Sacado::Fad::DFad<double> >& field, dealii::Table<1,Sacado::Fad::DFad<double> >& source, dealii::Table<2, Sacado::Fad::DFad<double>> &phi_e_grad, dealii::Table<2, Sacado::Fad::DFad<double>> &c_li_plus_grad, dealii::Table<1, Sacado::Fad::DFad<double>> &c_li_plus, dealii::Table<1, double> &c_li_plus_old);

};


//****************
//****************

template <int dim>
class Temperature:public Transportation<dim>
{
	public:
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
		nlohmann::json* params_json;
		void declare_parameters(nlohmann::json& _params);
		void set_stress(dealii::Table<3,Sacado::Fad::DFad<double> >& F, dealii::Table<3,Sacado::Fad::DFad<double> >& P);		
};

template <int dim>
class Diffuse_interface:public Transportation<dim>
{
	public:
		nlohmann::json* params_json;
		void declare_parameters(nlohmann::json& _params);
};

#endif