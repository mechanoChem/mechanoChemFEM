/*
zhenlin wang 2019
*coupled diffusion reaction
*/
#ifndef battery_h
#include "mechanoChemFEM.h"
#include "battery_package/include/battery_components.h"
#include "battery_package/include/SDdata.h"
#include "nodalField.h"

template <int dim>
class battery: public mechanoChemFEM<dim>
{
	public:
		battery(std::string parameter_file_Dir);
		~battery();
		//this is a overloaded function 
		void apply_boundary_condition();
		void apply_initial_condition();
		void setMultDomain();
		void make_grid();
		void apply_Neumann_boundary_condition();
		void declare_parameters();
		void get_residual(const typename hp::DoFHandler<dim>::active_cell_iterator &cell, const FEValues<dim>& fe_values, Table<1, Sacado::Fad::DFad<double> >& R, Table<1, Sacado::Fad::DFad<double>>& ULocal, Table<1, double >& ULocalConv);
		void define_battery_fields();
		void get_residual_at_diffuse_interface(const typename hp::DoFHandler<dim>::active_cell_iterator &cell, const FEValues<dim>& fe_values, Table<1, Sacado::Fad::DFad<double> >& R, Table<1, Sacado::Fad::DFad<double>>& ULocal, Table<1, double >& ULocalConv);
		void output_w_domain();
		void output_results();
		void run();
		
		
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
		Diffuse_interface<dim> diffuse_interface;
		
		void setup_diffuse_interface();		
    void identify_diffuse_interface();
    std::vector<SDdata<dim>> cell_SDdata;
		
		nodalField<dim> computedNodalField;
		double iso_value = 0.5;
		int active_particle_id=0, electrolyte_id=1, interface_id=2;

	  Vector<double> crack_id;  // has to be double for output purpose
	  Vector<double> jump_n; 
	  Vector<double> jump_m; 
	  Vector<double> jump_w; 
	  Vector<double> T_n; 
    std::vector<bool> is_new_step;
    std::vector<std::vector<double>> pressure; 
    std::vector<std::vector<double>> pressure_old; 

		PETScWrappers::MPI::Vector solution_k;
};

#endif
