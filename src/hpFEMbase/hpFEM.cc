#include"../../include/hpFEM.h"
template <int dim>
hpFEM<dim>::hpFEM ()
	: triangulation (Triangulation<dim>::maximum_smoothing), dof_handler (triangulation)
{

}
template <int dim>
hpFEM<dim>::~hpFEM (){dof_handler.clear ();}


template class hpFEM<1>;
template class hpFEM<2>;
template class hpFEM<3>;


