#include"initBoundValProbs.h"

template <int dim>
void initBoundValProbs<dim>::declare_parameters()
{

	//declare problem setting
	params->enter_subsection("Problem");
	params->declare_entry("print_parameter","true",Patterns::Bool());
	params->declare_entry("dt","0",Patterns::Double() );
	params->declare_entry("totalTime","0",Patterns::Double() );
	params->declare_entry("step_load","false",Patterns::Bool());
	
	params->declare_entry("first_domain_id","0",Patterns::Integer());
	params->declare_entry("active_material_id","0",Patterns::Integer());
	params->declare_entry("electrolyte_id","0",Patterns::Integer());
	params->declare_entry("current_collector_id","0",Patterns::Integer());
	params->declare_entry("binder_id","0",Patterns::Integer());
	params->declare_entry("solid_id","0",Patterns::Integer());

	params->declare_entry("active_material_fe","0",Patterns::Integer());
	params->declare_entry("electrolyte_fe","0",Patterns::Integer());
	params->declare_entry("current_collector_fe","0",Patterns::Integer());
	params->declare_entry("binder_fe","0",Patterns::Integer());
	params->declare_entry("solid_fe","0",Patterns::Integer());
	
	//directory
	params->declare_entry("mesh","1",Patterns::FileName() );
	params->declare_entry("output_directory","1",Patterns::DirectoryName() );
	params->declare_entry("snapshot_directory","1",Patterns::DirectoryName() );
	//FEM
	params->declare_entry("volume_quadrature","0",Patterns::Double());
	params->declare_entry("face_quadrature","0",Patterns::Double());
	//applied current
	params->declare_entry("IpA","0",Patterns::Double());  
	params->leave_subsection();	
	
	//declare some useful geometry information beforehand
	params->enter_subsection("Geometry");
	params->declare_entry("X_0","0",Patterns::Double() );
	params->declare_entry("Y_0","0",Patterns::Double() );
	params->declare_entry("Z_0","0",Patterns::Double() );
	params->declare_entry("X_end","0",Patterns::Double() );
	params->declare_entry("Y_end","0",Patterns::Double() );
	params->declare_entry("Z_end","0",Patterns::Double() );
	
	params->declare_entry("electrode_Y1","0",Patterns::Double() );
	params->declare_entry("electrode_Y2","0",Patterns::Double() );
	params->declare_entry("currentCollector_Y1","0",Patterns::Double() );
	params->declare_entry("currentCollector_Y2","0",Patterns::Double() );
	
	params->declare_entry("particle_R","0",Patterns::Double() );
	params->declare_entry("particle_number","0",Patterns::Integer() );	
	params->leave_subsection();	
	
	//declare initial condition
	params->enter_subsection("Initial condition");	
	params->declare_entry("IpA","0",Patterns::Double() );
  params->declare_entry("c_li_max_neg","0",Patterns::Double() );
	params->declare_entry("c_li_max_pos","0",Patterns::Double());
	params->declare_entry("c_li_100_neg","0",Patterns::Double() );
  params->declare_entry("c_li_100_pos","0",Patterns::Double() );
	params->declare_entry("c_li_0_neg","0",Patterns::Double() );
  params->declare_entry("c_li_0_pos","0",Patterns::Double() );
	params->declare_entry("c_li_plus_ini","0",Patterns::Double() );
	params->declare_entry("T_0","0",Patterns::Double());
	params->leave_subsection();
	//declare paramter for elastiticy equations
	params->enter_subsection("Elasticity" );
	params->declare_entry("youngsModulus_Al","0",Patterns::Double() );
	params->declare_entry("youngsModulus_Cu","0",Patterns::Double() );	
	params->declare_entry("youngsModulus_binder","0",Patterns::Double());	
	params->declare_entry("youngsModulus_s_neg","0",Patterns::Double() );
	params->declare_entry("youngsModulus_s_pos","0",Patterns::Double() );
	params->declare_entry("youngsModulus_sep","0",Patterns::Double() );
	
	params->declare_entry("nu_Al","0",Patterns::Double() );
	params->declare_entry("nu_Cu","0",Patterns::Double() );
	params->declare_entry("nu_binder","0",Patterns::Double());
	params->declare_entry("nu_sep","0",Patterns::Double() );
	params->declare_entry("nu_s_neg","0",Patterns::Double() );
	params->declare_entry("nu_s_pos","0",Patterns::Double() );
	
	params->declare_entry("omega_c","0",Patterns::Double() );
	params->declare_entry("omega_t_s_neg","0",Patterns::Double() );
	params->declare_entry("omega_t_s_pos","0",Patterns::Double() );
	params->declare_entry("omega_t_sep","0",Patterns::Double() );
	params->declare_entry("omega_t_Al","0",Patterns::Double() );
	params->declare_entry("omega_t_Cu","0",Patterns::Double() );
	params->declare_entry("omega_t_binder","0",Patterns::Double());
	params->leave_subsection();	
	//declare paramter for fluid equations
	params->enter_subsection("Fluid");
	params->declare_entry("viscosity","0",Patterns::Double() );
	params->declare_entry("youngsModulus_mesh","0",Patterns::Double() );
	params->declare_entry("nu_mesh","0",Patterns::Double() );
	params->leave_subsection();	
	
	//declare paramter for electro-chemo equations
	params->enter_subsection("ElectroChemo" );
	params->declare_entry("sigma_s_neg","0",Patterns::Double() );
	params->declare_entry("sigma_s_pos","0",Patterns::Double() );
	params->declare_entry("sigma_Al","0",Patterns::Double() );
	params->declare_entry("sigma_Cu","0",Patterns::Double() );
	params->declare_entry("sigma_binder","0",Patterns::Double());
	
	params->declare_entry("t_0","0",Patterns::Double() );
  params->declare_entry("D_li_neg","0",Patterns::Double() );
	params->declare_entry("D_li_pos","0",Patterns::Double() );
	params->leave_subsection();	
	
	//declare paramter for thermal equations
	params->enter_subsection("Thermal" );
	params->declare_entry("lambda_s_neg","0",Patterns::Double() );
	params->declare_entry("lambda_s_pos","0",Patterns::Double() );
	params->declare_entry("lambda_e","0",Patterns::Double() );
	params->declare_entry("lambda_sep","0",Patterns::Double() );
	params->declare_entry("lambda_Al","0",Patterns::Double() );
	params->declare_entry("lambda_Cu","0",Patterns::Double() );
	params->declare_entry("lambda_binder","0",Patterns::Double());
	
	params->declare_entry("density_s_neg","0",Patterns::Double() );
	params->declare_entry("density_s_pos","0",Patterns::Double() );
	params->declare_entry("density_sep","0",Patterns::Double() );
	params->declare_entry("density_e","0",Patterns::Double() );
	params->declare_entry("density_Al","0",Patterns::Double() );
	params->declare_entry("density_Cu","0",Patterns::Double() );
	params->declare_entry("density_binder","0",Patterns::Double());
	
	params->declare_entry("Cp_s_neg","0",Patterns::Double() );
	params->declare_entry("Cp_s_pos","0",Patterns::Double() );
	params->declare_entry("Cp_e","0",Patterns::Double() );
	params->declare_entry("Cp_sep","0",Patterns::Double() );
	params->declare_entry("Cp_Al","0",Patterns::Double() );
	params->declare_entry("Cp_Cu","0",Patterns::Double() );
	params->declare_entry("Cp_binder","0",Patterns::Double());
	params->declare_entry("h","0",Patterns::Double() );
	params->leave_subsection();	
	
}

template class initBoundValProbs<1>;
template class initBoundValProbs<2>;
template class initBoundValProbs<3>;