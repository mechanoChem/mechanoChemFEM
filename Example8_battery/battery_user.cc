/*
zhenlin wang 2019
*coupled diffusion reaction
*/
#include "battery.h"

//set Dirichlet BC
template <int dim>
void battery<dim>::apply_boundary_condition()
{
	constraints->clear ();
	
	DoFTools::make_hanging_node_constraints (this->dof_handler, *constraints);
	
	int totalDOF=this->totalDOF(this->primary_variables);
  std::vector<bool> All_component (totalDOF, true);	
	VectorTools:: interpolate_boundary_values (this->dof_handler, dim, ZeroFunction<dim> (totalDOF),*constraints, All_component);
	
	constraints->close ();
}
template <int dim>
void battery<dim>::apply_Neumann_boundary_condition()
{
	
}

template class battery<1>;
template class battery<2>;
template class battery<3>;

template <int dim>
void InitialConditions<dim>::vector_value (const Point<dim>   &p, Vector<double>   &values) const{
	Assert (values.size() == totalDOF, ExcDimensionMismatch (values.size(), totalDOF));
	values=0;
	for(unsigned int i=0;i<primary_variables.size();i++){
		if(std::strcmp(primary_variables[i][0].c_str(),"Lithium")==0){
			values(primary_variables_dof[i])= 0.5;
			if(p[2]==0) values(primary_variables_dof[i])= 1;
		}
		//if(std::strcmp(primary_variables[i][0].c_str(),"Displacement")==0) values(primary_variables_dof[i])= 1.0;
	}
}

template class InitialConditions<1>;
template class InitialConditions<2>;
template class InitialConditions<3>;