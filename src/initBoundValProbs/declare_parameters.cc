/*
zhenlin wang 2019
*/

#include"../../include/mechanoChemFEM.h"
template <int dim>
void mechanoChemFEM<dim>::declare_parameters_mechanoChemFEM()
{
	params_mechanoChemFEM->enter_subsection("Problem");
	params_mechanoChemFEM->declare_entry("print_parameter","true",Patterns::Bool());
	params_mechanoChemFEM->declare_entry("dt","0",Patterns::Double() );
	params_mechanoChemFEM->declare_entry("totalTime","0",Patterns::Double() );
	params_mechanoChemFEM->declare_entry("current_increment","0",Patterns::Integer());
	params_mechanoChemFEM->declare_entry("current_time","0",Patterns::Double());
	params_mechanoChemFEM->declare_entry("resuming_from_snapshot","false",Patterns::Bool());
	params_mechanoChemFEM->declare_entry("save_output","true",Patterns::Bool());
	params_mechanoChemFEM->declare_entry("save_snapshot","true",Patterns::Bool());
	
	params_mechanoChemFEM->declare_entry("off_output_index","0",Patterns::Integer() );
	params_mechanoChemFEM->declare_entry("skip_output","1",Patterns::Integer() );
	
	params_mechanoChemFEM->declare_entry("mesh","1",Patterns::FileName() );
	params_mechanoChemFEM->declare_entry("snapshot_file","1",Patterns::DirectoryName() );
	params_mechanoChemFEM->declare_entry("output_directory","1",Patterns::DirectoryName() );
	params_mechanoChemFEM->declare_entry("snapshot_directory","1",Patterns::DirectoryName() );
	
	//FEM
	params_mechanoChemFEM->declare_entry("volume_quadrature","3",Patterns::Integer());
	params_mechanoChemFEM->declare_entry("face_quadrature","2",Patterns::Integer() );
	params_mechanoChemFEM->leave_subsection();	
	
	params_mechanoChemFEM->enter_subsection("Geometry");
	params_mechanoChemFEM->declare_entry("X_0","0",Patterns::Double() );
	params_mechanoChemFEM->declare_entry("Y_0","0",Patterns::Double() );
	params_mechanoChemFEM->declare_entry("Z_0","0",Patterns::Double() );
	params_mechanoChemFEM->declare_entry("X_end","0",Patterns::Double() );
	params_mechanoChemFEM->declare_entry("Y_end","0",Patterns::Double() );
	params_mechanoChemFEM->declare_entry("Z_end","0",Patterns::Double() );
	
	params_mechanoChemFEM->declare_entry("element_div_x","0",Patterns::Integer() );
	params_mechanoChemFEM->declare_entry("element_div_y","0",Patterns::Integer() );
	params_mechanoChemFEM->declare_entry("element_div_z","0",Patterns::Integer() );
	params_mechanoChemFEM->leave_subsection();		
}

template <int dim>
void mechanoChemFEM<dim>::load_parameters(std::string parametersfile)
{	
	params_mechanoChemFEM->read_input (parametersfile);
	params_mechanoChemFEM->enter_subsection("Problem");
	bool printParameter=params_mechanoChemFEM->get_bool("print_parameter");
	output_directory=params_mechanoChemFEM->get("output_directory");
	snapshot_directory=params_mechanoChemFEM->get("snapshot_directory");
	skip_output=params_mechanoChemFEM->get_integer("skip_output");
	snapfile=params_mechanoChemFEM->get("snapshot_file");
	current_dt=params_mechanoChemFEM->get_double("dt");
	total_time=params_mechanoChemFEM->get_double("totalTime");
	current_increment=params_mechanoChemFEM->get_integer("current_increment");
	current_time=params_mechanoChemFEM->get_double("current_time");
	resuming_from_snapshot=params_mechanoChemFEM->get_bool("resuming_from_snapshot");
	save_snapshot=params_mechanoChemFEM->get_bool("save_snapshot");
	save_output=params_mechanoChemFEM->get_bool("save_output");
	off_output_index=params_mechanoChemFEM->get_integer("off_output_index");

	volume_quadrature= new const QGauss<dim>(params_mechanoChemFEM->get_integer("volume_quadrature"));
	common_face_quadrature= new const QGauss<dim-1>(params_mechanoChemFEM->get_integer("face_quadrature"));
	
	params_mechanoChemFEM->leave_subsection();	

  const int dir_err1 = system(("mkdir -p " + output_directory).c_str());
  const int dir_err2 = system(("mkdir -p " + snapshot_directory).c_str());
  if (dir_err1 == -1 or dir_err2 == -1)
  {
    printf("Error creating directory!\n");
    exit(1);
  }
	
	if(printParameter) {
		if(this_mpi_process == 0) params_mechanoChemFEM->print_parameters (std::cout, ParameterHandler::Text);
	}
}

template class mechanoChemFEM<1>;
template class mechanoChemFEM<2>;
template class mechanoChemFEM<3>;
