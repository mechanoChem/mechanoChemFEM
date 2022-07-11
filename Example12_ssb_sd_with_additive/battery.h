/*
zhenlin wang 2019
*coupled diffusion reaction
*/
#ifndef battery_h
#define battery_h
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
		void get_residual_at_li_metal_interface(const typename hp::DoFHandler<dim>::active_cell_iterator &cell, const FEValues<dim>& fe_values, Table<1, Sacado::Fad::DFad<double> >& R, Table<1, Sacado::Fad::DFad<double>>& ULocal, Table<1, double >& ULocalConv);
		void get_residual_at_additive_interface(const typename hp::DoFHandler<dim>::active_cell_iterator &cell, const FEValues<dim>& fe_values, Table<1, Sacado::Fad::DFad<double> >& R, Table<1, Sacado::Fad::DFad<double>>& ULocal, Table<1, double >& ULocalConv);
		void output_w_domain();
		void output_results();
		void run();
    void save_sd_data();
    void load_sd_data();
		void solve_ibvp();
		
		
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
		int  electrolyte_id=0, active_particle_id=1, interface_id=2, li_metal_id=3, additive_id=4, li_metal_interface_id=5, additive_interface_id=6;
    //int anode_opt = 0; // li metal
    int anode_opt = 1; // graphite
    int cathode_opt = 0; // lco

	  Vector<double> crack_id;  // has to be double for output purpose
	  Vector<double> jump_n; 
	  Vector<double> jump_m; 
	  Vector<double> jump_w; 
	  Vector<double> pressure_gp0; 
	  Vector<double> pressure_gp1; 
	  Vector<double> pressure_gp2; 
	  Vector<double> pressure_gp3; 
	  Vector<double> jn; 
	  Vector<double> T_n; 
    std::vector<double> sim_time_info;
    std::vector<bool> is_new_step;
    std::vector<std::vector<double>> pressure; 
    std::vector<std::vector<double>> pressure_old; 
    int anode_point_dof_index = -1;
    int cathode_point_dof_index = -1;

		PETScWrappers::MPI::Vector solution_k;
    int N_GPs = 4; // # of GPs has to be 4 due to the SDs.
};

#endif
