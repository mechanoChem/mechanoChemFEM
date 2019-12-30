/*
zhenlin wang 2019
*/

#include"../../include/mechanoChemFEM.h"
template <int dim>
void mechanoChemFEM<dim>::declare_parameters_mechanoChemFEM()
{
	params->enter_subsection("Problem");
	params->declare_entry("print_parameter","true",Patterns::Bool());
	params->declare_entry("dt","0",Patterns::Double() );
	params->declare_entry("totalTime","0",Patterns::Double() );
	params->declare_entry("current_increment","0",Patterns::Integer());
	params->declare_entry("current_time","0",Patterns::Double());
	params->declare_entry("resuming_from_snapshot","false",Patterns::Bool());
	params->declare_entry("save_output","true",Patterns::Bool());
	params->declare_entry("save_snapshot","true",Patterns::Bool());
	
	params->declare_entry("off_output_index","0",Patterns::Integer() );
	params->declare_entry("skip_output","1",Patterns::Integer() );
	
	params->declare_entry("mesh","1",Patterns::FileName() );
	params->declare_entry("snapshot_file","1",Patterns::DirectoryName() );
	params->declare_entry("output_directory","1",Patterns::DirectoryName() );
	params->declare_entry("snapshot_directory","1",Patterns::DirectoryName() );
	
	//FEM
	params->declare_entry("volume_quadrature","3",Patterns::Integer());
	params->declare_entry("face_quadrature","2",Patterns::Integer() );
	params->leave_subsection();	
	
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
	params->leave_subsection();		
}

template <int dim>
void mechanoChemFEM<dim>::load_parameters()
{
	params->read_input ("parameters.prm");
	
	params->enter_subsection("Problem");
	bool printParameter=params->get_bool("print_parameter");
	params->leave_subsection();	
	if(printParameter) {
		if(this_mpi_process == 0) params->print_parameters (std::cout, ParameterHandler::Text);
	}
}

template class mechanoChemFEM<1>;
template class mechanoChemFEM<2>;
template class mechanoChemFEM<3>;
