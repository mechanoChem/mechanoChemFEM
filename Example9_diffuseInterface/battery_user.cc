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
	
	constraints->close ();
}
template <int dim>
void battery<dim>::apply_Neumann_boundary_condition()
{
	
}

template <int dim>
void battery<dim>::setup_diffuse_interface()
{
	Point<dim> origin(2,2);
  typename hp::DoFHandler<dim>::active_cell_iterator cell = dof_handler_interface->begin_active(), endc=dof_handler_interface->end();
  for (;cell!=endc; ++cell){
		if (cell->subdomain_id() == this->this_mpi_process){
			Point<dim> cell_center=cell->center();
			if(cell_center.distance(origin)<1){
				const unsigned int dofs_per_cell = cell->get_fe().dofs_per_cell;
				std::vector<unsigned int> local_dof_indices (dofs_per_cell);
				cell->get_dof_indices (local_dof_indices);
				for (unsigned int i=0; i<dofs_per_cell; ++i) {
					diffuse_interface(local_dof_indices[i])=1;
				}
			}
		}		
	}
	diffuse_interface.compress(VectorOperation::insert);
	
	std::string path = this->output_directory+"diffuse_interface.vtk";
	FEMdata<dim,PETScWrappers::MPI::Vector> FEMdata_out_interface(*dof_handler_interface);
	FEMdata_out_interface.set_output_name(variables_interface);
	FEMdata_out_interface.write_vtk(diffuse_interface,path);
}
	
template <int dim>
void battery<dim>::setMultDomain()
{
	Point<dim> origin(2,2);
	for (typename Triangulation<dim>::active_cell_iterator cell = this->dof_handler.begin_active(); cell != this->dof_handler.end(); ++cell){
		Point<dim> cell_center = cell->center();
		if(cell_center.distance(origin)<1){
			cell->set_material_id(1);
		}
	}
}
	
template <int dim>
void battery<dim>::apply_initial_condition()
{
	setup_diffuse_interface_FEM();
	this->pcout << "applying initial condition\n";
	int index_lithium=battery_fields.active_fields_index["Lithium"];
	int totalDOF=this->totalDOF(this->primary_variables);
	typename hp::DoFHandler<dim>::active_cell_iterator cell_interface = dof_handler_interface->begin_active();
	typename hp::DoFHandler<dim>::active_cell_iterator cell = this->dof_handler.begin_active(), endc=this->dof_handler.end();
	for (;cell!=endc; ++cell, ++cell_interface){
		if (cell->subdomain_id() == this->this_mpi_process){
			hp::FEValues<dim> hp_fe_values (this->fe_collection, this->q_collection, update_values | update_quadrature_points);
    	hp_fe_values.reinit (cell);
			const unsigned int dofs_per_cell = cell->get_fe().dofs_per_cell;
			std::vector<unsigned int> local_dof_indices (dofs_per_cell);
			cell->get_dof_indices (local_dof_indices);
			const unsigned int dofs_per_cell_interface = cell_interface->get_fe().dofs_per_cell;
			std::vector<unsigned int> local_dof_indices_interface (dofs_per_cell_interface);
			cell_interface->get_dof_indices(local_dof_indices_interface);
			
			for (unsigned int i=0; i<dofs_per_cell_interface; ++i) {
				if(diffuse_interface(local_dof_indices_interface[i])>0.5){
					this->solution_prev(local_dof_indices[i*totalDOF+index_lithium])=0.5+ 0.04*(static_cast <double> (rand())/(static_cast <double>(RAND_MAX))-0.5);
				}
			}
		}
	}
	this->solution_prev.compress(VectorOperation::insert);
	this->solution=this->solution_prev;		
}

template class battery<1>;
template class battery<2>;
template class battery<3>;

// template <int dim>
// void InitialConditions<dim>::vector_value (const Point<dim>   &p, Vector<double>   &values) const{
// 	Assert (values.size() == totalDOF, ExcDimensionMismatch (values.size(), totalDOF));
// 	values=0;
// 	Point<dim> origin(2,2);
// 	if(p.distance(origin)<1){
// 		for(unsigned int i=0;i<primary_variables.size();i++){
// 			if(std::strcmp(primary_variables[i][0].c_str(),"Lithium")==0){
// 				values(primary_variables_dof[i])= 0.5+ 0.04*(static_cast <double> (rand())/(static_cast <double>(RAND_MAX))-0.5);;
// 			}
// 		}
// 	}
// }
//
// template class InitialConditions<1>;
// template class InitialConditions<2>;
// template class InitialConditions<3>;