/*
zhenlin wang 2019
*/

#include"../../include/initBoundValProbs.h"

template <int dim>
void initBoundValProbs<dim>::pre_run()
{
	FEMdata_out.set_output_name(primary_variables);	
	declare_parameters_initBoundValProbs();
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
	off_output_index=params->get_integer("off_output_index");

	const QGauss<dim> _volume_quadrature(params->get_integer("volume_quadrature"));
	volume_quadrature=&_volume_quadrature;
	const QGauss<dim-1> _face_quadrature(params->get_integer("face_quadrature"));
	common_face_quadrature=&_face_quadrature;
	
	params->leave_subsection();	
	
	if(printParameter) {
		if(this_mpi_process == 0) params->print_parameters (std::cout, ParameterHandler::Text);
	}
	this->setup_FeSystem(fe_system, fe_collection, q_collection, primary_variables_dof,primary_variables,FE_support,*volume_quadrature);
}


template class initBoundValProbs<1>;
template class initBoundValProbs<2>;
template class initBoundValProbs<3>;
