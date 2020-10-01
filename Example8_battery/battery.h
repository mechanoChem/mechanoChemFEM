/*
zhenlin wang 2019
*coupled diffusion reaction
*/
#include "mechanoChemFEM.h"
#include "battery_package/include/battery_components.h"

template <int dim>
class battery: public mechanoChemFEM<dim>
{
	public:
		battery(std::string parameter_file_Dir);
		//this is a overloaded function 
		void apply_boundary_condition();
		void apply_Neumann_boundary_condition();
		void get_residual(const typename hp::DoFHandler<dim>::active_cell_iterator &cell, const FEValues<dim>& fe_values, Table<1, Sacado::Fad::DFad<double> >& R, Table<1, Sacado::Fad::DFad<double>>& ULocal, Table<1, double >& ULocalConv);
		ParameterHandler* params;		
		nlohmann::json* params_json;
		ConstraintMatrix* constraints;
		Battery_fields battery_fields;
		Lithium<dim> lithium;
		Lithium_phaseField<dim> lithium_mu;
		Displacement<dim> displacement;
		
		
};