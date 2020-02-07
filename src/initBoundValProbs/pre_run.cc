/*
zhenlin wang 2019
*/

#include"../../include/mechanoChemFEM.h"
#include <cstdlib>

template <int dim>
void mechanoChemFEM<dim>::define_primary_fields()
{
	//define_primary_fields from parameters file.
}

template <int dim>
void mechanoChemFEM<dim>::init_ibvp()
{
	FEMdata_out.set_output_name(primary_variables);	
	this->setup_FeSystem(fe_system, fe_collection, q_collection, primary_variables_dof,primary_variables,FE_support,*volume_quadrature);
	make_grid();
	refine_grid();
	setMultDomain();
  mark_boundary();
	setup_linear_system();
}

template class mechanoChemFEM<1>;
template class mechanoChemFEM<2>;
template class mechanoChemFEM<3>;
