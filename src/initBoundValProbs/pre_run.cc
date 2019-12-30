/*
zhenlin wang 2019
*/

#include"../../include/mechanoChemFEM.h"
#include <cstdlib>

template <int dim>
void mechanoChemFEM<dim>::pre_run()
{
	FEMdata_out.set_output_name(primary_variables);	
	declare_parameters_mechanoChemFEM();
	params->read_input ("parameters.prm");
	params->enter_subsection("Problem");
	bool printParameter=params->get_bool("print_parameter");
	output_directory=params->get("output_directory");
	snapshot_directory=params->get("snapshot_directory");
	skip_output=params->get_integer("skip_output");
	snapfile=params->get("snapshot_file");
	current_dt=params->get_double("dt");
	total_time=params->get_double("totalTime");
	current_increment=params->get_integer("current_increment");
	current_time=params->get_double("current_time");
	resuming_from_snapshot=params->get_bool("resuming_from_snapshot");
	save_snapshot=params->get_bool("save_snapshot");
	save_output=params->get_bool("save_output");
	off_output_index=params->get_integer("off_output_index");

	volume_quadrature= new const QGauss<dim>(params->get_integer("volume_quadrature"));
	common_face_quadrature= new const QGauss<dim-1>(params->get_integer("face_quadrature"));
	
	params->leave_subsection();	

  const int dir_err1 = system(("mkdir -p " + output_directory).c_str());
  const int dir_err2 = system(("mkdir -p " + snapshot_directory).c_str());
  if (dir_err1 == -1 or dir_err2 == -1)
  {
    printf("Error creating directory!\n");
    exit(1);
  }
	
	if(printParameter) {
		if(this_mpi_process == 0) params->print_parameters (std::cout, ParameterHandler::Text);
	}
	this->setup_FeSystem(fe_system, fe_collection, q_collection, primary_variables_dof,primary_variables,FE_support,*volume_quadrature);
}

template class mechanoChemFEM<1>;
template class mechanoChemFEM<2>;
template class mechanoChemFEM<3>;
