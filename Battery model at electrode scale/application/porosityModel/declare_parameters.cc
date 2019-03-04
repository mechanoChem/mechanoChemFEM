#include"initBoundValProbs.h"
template <int dim>
void initBoundValProbs<dim>::declare_parameters()
{

	//declare problem setting
	params->enter_subsection("Problem");
	params->declare_entry("print_parameter","true",Patterns::Bool());
	
	params->declare_entry("snapshot_file","1",Patterns::FileName() );
	params->declare_entry("output_directory","1",Patterns::DirectoryName() );
	params->declare_entry("snapshot_directory","1",Patterns::DirectoryName() );
	params->declare_entry("resuming_from_snapshot","false",Patterns::Bool());
	
	params->declare_entry("dt","0",Patterns::Double() );
	params->declare_entry("totalTime","0",Patterns::Double() );
	params->declare_entry("step_load","false",Patterns::Bool());
	
	params->declare_entry("first_domain_id","0",Patterns::Integer());
	params->declare_entry("separator_id","0",Patterns::Integer());
	params->declare_entry("electrode_id","0",Patterns::Integer());

	params->declare_entry("separator_fe","0",Patterns::Integer());
	params->declare_entry("electrode_fe","0",Patterns::Integer());
	
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
	
	params->declare_entry("element_div_x","0",Patterns::Integer() );
	params->declare_entry("element_div_y","0",Patterns::Integer() );
	params->declare_entry("element_div_z","0",Patterns::Integer() );
	
	params->declare_entry("electrode_Y1","0",Patterns::Double() );
	params->declare_entry("electrode_Y2","0",Patterns::Double() );
	
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
	params->declare_entry("dis_top0","0",Patterns::Double());

	params->leave_subsection();
	
	//declare paramter for elastiticy equations
	params->enter_subsection("Elasticity" );
	params->declare_entry("youngsModulus_neg","0",Patterns::Double() );
	params->declare_entry("youngsModulus_sep","0",Patterns::Double() );
	params->declare_entry("youngsModulus_pos","0",Patterns::Double() );
	params->declare_entry("nu_neg","0",Patterns::Double() );
	params->declare_entry("nu_pos","0",Patterns::Double() );
	params->declare_entry("nu_sep","0",Patterns::Double() );
	
	params->declare_entry("kappa_sep","0",Patterns::Double() );
	params->declare_entry("kappa_neg","0",Patterns::Double() );
	params->declare_entry("kappa_pos","0",Patterns::Double() );
	params->declare_entry("kappa_s","0",Patterns::Double() );
	
	params->declare_entry("pl","0",Patterns::Double() );
	params->declare_entry("pb","0",Patterns::Double() );
	
	params->declare_entry("omega_neg","0",Patterns::Double() );
	params->declare_entry("omega_pos","0",Patterns::Double() );
	params->declare_entry("omega_sep","0",Patterns::Double() );
	params->declare_entry("omega_s","0",Patterns::Double() );

	params->leave_subsection();	

	//declare paramter for electro-chemo equations
	params->enter_subsection("ElectroChemo" );
	params->declare_entry("sigma_neg","0",Patterns::Double() );
	params->declare_entry("sigma_pos","0",Patterns::Double() );
	
	params->declare_entry("t_0","0",Patterns::Double() );
  params->declare_entry("D_li_neg","0",Patterns::Double() );
	params->declare_entry("D_li_pos","0",Patterns::Double() );
	
	params->declare_entry("eps_s_0_neg","0",Patterns::Double() );
	params->declare_entry("eps_s_0_pos","0",Patterns::Double() );
	params->declare_entry("eps_s_0_sep","0",Patterns::Double() );
	
	params->declare_entry("eps_l_0_neg","0",Patterns::Double() );
	params->declare_entry("eps_l_0_pos","0",Patterns::Double() );
	params->declare_entry("eps_l_0_sep","0",Patterns::Double() );
	
	params->declare_entry("eps_b_0_neg","0",Patterns::Double() );
	params->declare_entry("eps_b_0_pos","0",Patterns::Double() );
	params->declare_entry("eps_b_0_sep","0",Patterns::Double() );
	
	params->declare_entry("R_s_0_neg","0",Patterns::Double() );
	params->declare_entry("R_s_0_pos","0",Patterns::Double() );
	params->declare_entry("R_s_0_sep","0",Patterns::Double() );
	params->leave_subsection();	
	
	//declare paramter for thermal equations
	params->enter_subsection("Thermal" );
	params->declare_entry("lambda_neg","0",Patterns::Double() );
	params->declare_entry("lambda_pos","0",Patterns::Double() );
	params->declare_entry("lambda_sep","0",Patterns::Double() );
	
	params->declare_entry("density_neg","0",Patterns::Double() );
	params->declare_entry("density_pos","0",Patterns::Double() );
	params->declare_entry("density_sep","0",Patterns::Double() );
	
	params->declare_entry("Cp_s_neg","0",Patterns::Double() );
	params->declare_entry("Cp_s_pos","0",Patterns::Double() );
	params->declare_entry("Cp_sep","0",Patterns::Double() );
	
	params->declare_entry("h","0",Patterns::Double() );
	params->leave_subsection();	
}

template class initBoundValProbs<1>;
template class initBoundValProbs<2>;
template class initBoundValProbs<3>;