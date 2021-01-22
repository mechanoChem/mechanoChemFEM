/*
zhenlin wang 2019
*coupled diffusion reaction
*/
#ifndef battery_h
#include "mechanoChemFEM.h"
#include "battery_package/include/battery_components.h"

template <int dim>
class battery: public mechanoChemFEM<dim>
{
	public:
		battery(std::string parameter_file_Dir);
		~battery();
		//this is a overloaded function 
		void apply_boundary_condition();
		void apply_Neumann_boundary_condition();
		void declare_parameters();
		void get_residual(const typename hp::DoFHandler<dim>::active_cell_iterator &cell, const FEValues<dim>& fe_values, Table<1, Sacado::Fad::DFad<double> >& R, Table<1, Sacado::Fad::DFad<double>>& ULocal, Table<1, double >& ULocalConv);
		void setMultDomain();
		void output_w_domain();
		void define_battery_fields();
		void get_residual_at_diffuse_interface(const typename hp::DoFHandler<dim>::active_cell_iterator &cell, const FEValues<dim>& fe_values, Table<1, Sacado::Fad::DFad<double> >& R, Table<1, Sacado::Fad::DFad<double>>& ULocal, Table<1, double >& ULocalConv);
		
		
		ParameterHandler* params;		
		nlohmann::json* params_json;
		ConstraintMatrix* constraints;
		Battery_fields<dim> battery_fields;
		ElectricChemo<dim, Sacado::Fad::DFad<double>> electricChemoFormula;
		Lithium<dim> lithium;
		Lithium_phaseField<dim> lithium_mu;
		Lithium_cation<dim> lithium_cation;
		Electrode_potential<dim> phi_s;
		Electrolyte_potential<dim> phi_e;
		Displacement<dim> displacement;
		
		/*
		--------------------------------------------------------------------------------------------
		* FE for order parameter to define diffuse interface 
		--------------------------------------------------------------------------------------------
		*/
		// void setup_diffuse_interface();
		// std::vector<std::shared_ptr<FESystem<dim>> > fe_system_interface;
		//     hp::FECollection<dim> fe_collection_interface;
		//     hp::QCollection<dim>  q_collection_interface;
		// hp::DoFHandler<dim>*  dof_handler_interface;
		//
		// PETScWrappers::MPI::Vector diffuse_interface;
		
		
		
};

#endif