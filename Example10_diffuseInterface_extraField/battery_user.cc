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
	VectorTools:: interpolate_boundary_values (this->dof_handler, 1, ZeroFunction<dim> (totalDOF),*constraints, All_component);
	int interface_index=battery_fields.active_fields_index["Diffuse_interface"];
	if (interface_index>-1){
		typename hp::DoFHandler<dim>::active_cell_iterator cell = this->dof_handler.begin_active(), endc=this->dof_handler.end();
		for (;cell!=endc; ++cell){
			const unsigned int dofs_per_cell = cell->get_fe().dofs_per_cell;
			std::vector<unsigned int> local_dof_indices (dofs_per_cell);
 		 	cell->get_dof_indices (local_dof_indices);
			for (unsigned int i=0; i<dofs_per_cell; ++i) {
				const unsigned int ck = cell->get_fe().system_to_component_index(i).first;
				if(ck==interface_index) constraints->add_line (local_dof_indices[i]);
			}
		}
	}
	constraints->close ();
	setup_diffuse_interface();
}
template <int dim>
void battery<dim>::apply_Neumann_boundary_condition()
{}

	
template <int dim>
void battery<dim>::setup_diffuse_interface()
{

}

template class battery<1>;
template class battery<2>;
template class battery<3>;

template <int dim>
void InitialConditions<dim>::vector_value (const Point<dim>   &p, Vector<double>   &values) const{
	Assert (values.size() == totalDOF, ExcDimensionMismatch (values.size(), totalDOF));
	values=0;
	Point<dim> origin(2,2);
	for(unsigned int i=0;i<primary_variables.size();i++){
		if(std::strcmp(primary_variables[i][0].c_str(),"Diffuse_interface")==0){
			if(p.distance(origin)<1) values(primary_variables_dof[i])=1;
			else if(p.distance(origin)<1.1) values(primary_variables_dof[i])=1-(p.distance(origin)-1)*10;
			else values(primary_variables_dof[i])=0.0;
		}
	}
	for(unsigned int i=0;i<primary_variables.size();i++){
		if(std::strcmp(primary_variables[i][0].c_str(),"Lithium")==0){
			values(primary_variables_dof[i])= 0.5+ 0.04*(static_cast <double> (rand())/(static_cast <double>(RAND_MAX))-0.5);
			values(primary_variables_dof[i])=values(primary_variables_dof[i])*values(primary_variables_dof.back());
		}
		if(std::strcmp(primary_variables[i][0].c_str(),"Lithium_cation")==0){
			values(primary_variables_dof[i])=1*(1-values(primary_variables_dof.back()));
		}
	}
}

template class InitialConditions<1>;
template class InitialConditions<2>;
template class InitialConditions<3>;