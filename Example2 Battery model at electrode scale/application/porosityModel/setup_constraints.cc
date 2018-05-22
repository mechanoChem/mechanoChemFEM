#include "initBoundValProbs.h"

template <int dim>
void initBoundValProbs<dim>::setup_constraints()
{	
	hpFEM<dim>::constraints.clear ();
	
	
	int totalDOF=this->totalDOF(primary_variables);	
  std::vector<bool> x_component (totalDOF, false); x_component[u_dof]=true; 
  std::vector<bool> y_component (totalDOF, false); y_component[u_dof+1]=true;  
	std::vector<bool> z_component (totalDOF, false); if(dim==3) z_component[u_dof+2]=true;  
	std::vector<bool> phi_s_component (totalDOF, false); phi_s_component[phi_s_dof]=true;
	
  DoFTools::make_hanging_node_constraints (hpFEM<dim>::dof_handler, hpFEM<dim>::constraints);

  VectorTools:: interpolate_boundary_values (hpFEM<dim>::dof_handler, 2, ZeroFunction<dim> (totalDOF),hpFEM<dim>::constraints, x_component);
  VectorTools:: interpolate_boundary_values (hpFEM<dim>::dof_handler, 2, ZeroFunction<dim> (totalDOF),hpFEM<dim>::constraints, y_component);
  if(dim==3) VectorTools:: interpolate_boundary_values (hpFEM<dim>::dof_handler, 2, ZeroFunction<dim> (totalDOF),hpFEM<dim>::constraints, z_component);
  VectorTools:: interpolate_boundary_values (hpFEM<dim>::dof_handler, dim+2, ZeroFunction<dim> (totalDOF),hpFEM<dim>::constraints, x_component);
  VectorTools:: interpolate_boundary_values (hpFEM<dim>::dof_handler, dim+2, ZeroFunction<dim> (totalDOF),hpFEM<dim>::constraints, y_component);
  if(dim==3) VectorTools:: interpolate_boundary_values (hpFEM<dim>::dof_handler, dim+2, ZeroFunction<dim> (totalDOF),hpFEM<dim>::constraints, z_component);
	
	VectorTools:: interpolate_boundary_values (hpFEM<dim>::dof_handler, 2, ZeroFunction<dim> (totalDOF),hpFEM<dim>::constraints,phi_s_component);
	
  hpFEM<dim>::constraints.close ();
  pcout << "   Number of active cells:       " << hpFEM<dim>::triangulation.n_active_cells() << std::endl;
  pcout << "   Number of degrees of freedom: " << hpFEM<dim>::dof_handler.n_dofs() << std::endl; 
  
}


template class initBoundValProbs<1>;
template class initBoundValProbs<2>;
template class initBoundValProbs<3>;
