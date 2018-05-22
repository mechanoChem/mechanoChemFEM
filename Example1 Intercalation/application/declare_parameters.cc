#include"initBoundValProbs.h"
template <int dim>
void initBoundValProbs<dim>::declare_parameters()
{
	params->enter_subsection("Problem");
	params->declare_entry("print_parameter","true",Patterns::Bool());

	params->declare_entry("with_additional_data","true",Patterns::Bool());
	params->declare_entry("resuming_additional_data_from_snapshot","true",Patterns::Bool());
	params->declare_entry("additional_data_snapshot","1",Patterns::FileName() );
	params->declare_entry("additional_data_nodal_file","1",Patterns::FileName() );
	
	params->declare_entry("dt","0",Patterns::Double() );
	params->declare_entry("totalTime","0",Patterns::Double() );
  params->declare_entry("current_increment","0",Patterns::Integer());
  params->declare_entry("current_time","0",Patterns::Double());
  params->declare_entry("resuming_from_snapshot","false",Patterns::Bool());
	
	
	params->declare_entry("meshType","external",Patterns::Selection("external|half_hyperShell") );
	params->declare_entry("mesh","1",Patterns::FileName() );
	
	params->declare_entry("Shell_center_X","0",Patterns::Double() );
	params->declare_entry("Shell_center_Y","0",Patterns::Double() );
	params->declare_entry("Shell_center_Z","0",Patterns::Double() );
	params->declare_entry("inner_radius","0",Patterns::Double() );
	params->declare_entry("outer_radius","0",Patterns::Double() );
	params->declare_entry("initial_global_refine","0",Patterns::Integer() );
	
	params->declare_entry("a_scale","0",Patterns::Double() );
	params->declare_entry("b_scale","0",Patterns::Double() );
	params->declare_entry("c_scale","0",Patterns::Double() );
	params->declare_entry("cort_thickness","0",Patterns::Double() );
	
	params->declare_entry("snapshot_file","1",Patterns::FileName() );
	params->declare_entry("output_directory","1",Patterns::DirectoryName() );
	params->declare_entry("snapshot_directory","1",Patterns::DirectoryName() );
	
	//FEM
	params->declare_entry("volume_quadrature","0",Patterns::Integer());
	params->declare_entry("face_quadrature","0",Patterns::Integer() );
	
	params->declare_entry("first_domain_id","0",Patterns::Integer());
		
	params->declare_entry("Cortex_id","0",Patterns::Integer());
	params->declare_entry("Subcortex_id","0",Patterns::Integer());
	params->declare_entry("Ventricle_id","0",Patterns::Integer());
	
	params->declare_entry("Cortex_fe","0",Patterns::Integer());
	params->declare_entry("Subcortex_fe","0",Patterns::Integer());
	params->declare_entry("Ventricle_fe","0",Patterns::Integer());
	
	params->leave_subsection();	
	
	//declare paramters for mechanics
	params->enter_subsection("Mechanics");
	params->declare_entry("youngsModulus","0",Patterns::Double() );
	params->declare_entry("poissonRatio","0",Patterns::Double() );
	params->declare_entry("saturation_matID_Cortex","0",Patterns::Double() );
	params->declare_entry("saturation_matID_Subcortex","0",Patterns::Double() );
	
	params->declare_entry("GROWTH","Tangential",Patterns::Selection("Uniform|Isotropic|Tangential") );
	params->leave_subsection();	
	
	//declare paramters for concentrations
	params->enter_subsection("Concentration");
	params->declare_entry("inward_flux","0",Patterns::Double() );
	params->declare_entry("c1_ini","0",Patterns::Double() );
	params->declare_entry("c2_ini","0",Patterns::Double() );
	params->declare_entry("c1_ini_interface","0",Patterns::Double() );
	params->declare_entry("c2_ini_interface","0",Patterns::Double() );
	
	params->declare_entry("advection_type","fromFile",Patterns::Selection("fromFile|radial|noAdvection") );
	params->declare_entry("mobility_c1","0",Patterns::Double() );
	params->declare_entry("mobility_c2","0",Patterns::Double() );
	params->declare_entry("alpha_1","0",Patterns::Double() );
	params->declare_entry("alpha_2","0",Patterns::Double() );
	params->declare_entry("tau_1","0",Patterns::Double() );
	params->declare_entry("tau_2","0",Patterns::Double() );
	params->declare_entry("reac_10","0",Patterns::Double() );
	params->declare_entry("reac_20","0",Patterns::Double() );
	params->leave_subsection();	
	
}

template class initBoundValProbs<1>;
template class initBoundValProbs<2>;
template class initBoundValProbs<3>;
