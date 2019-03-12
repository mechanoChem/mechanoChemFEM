#include "initBoundValProbs.h"

template <int dim>
void initBoundValProbs<dim>::apply_advection_direction()
{ 
	params->enter_subsection("Problem");
	std::string advection_file=params->get("additional_data_nodal_file");
	params->leave_subsection();	
	std::ifstream file(advection_file.c_str(),std::ios::in);
	if (!file.is_open()) {
    std::cout << "There was a problem opening the input file!\n";
    exit(1);
	}
	std::vector< std::vector<double > > advection_array; 
	while(file)
	{
		std::vector<double > tem(2*dim);
		for (unsigned int j=0;j<dim*2;j++) {
			file >> tem[j];
		}
		advection_array.push_back(tem);
	}
	int total_active_node=advection_array.size();
	additional_data=0;
	hp::FEValues<dim> hp_fe_values (fe_collection_add, q_collection_add, update_values | update_quadrature_points  | update_JxW_values | update_gradients);
  typename hp::DoFHandler<dim>::active_cell_iterator cell = dof_handler_add->begin_active(), endc=dof_handler_add->end();
  for (;cell!=endc; ++cell){
		if(cell->material_id()==Ventricle_id) continue;
		if (cell->subdomain_id() == this_mpi_process){
    	hp_fe_values.reinit (cell);
    	const FEValues<dim> &fe_values = hp_fe_values.get_present_fe_values();
    	const unsigned int dofs_per_cell = cell->get_fe().dofs_per_cell;

			std::vector<unsigned int> local_dof_indices (dofs_per_cell);
 		 	cell->get_dof_indices (local_dof_indices);
			for (unsigned int vetex_index=0; vetex_index<GeometryInfo<dim>::vertices_per_cell; ++vetex_index) {
				Point<dim,double> vetex_coord=cell->vertex(vetex_index);
				for(unsigned int node=0;node<total_active_node;node++){
					Point<dim,double> node_coord(advection_array[node][0],advection_array[node][1],advection_array[node][2]);
					if(vetex_coord.distance(node_coord)<1.0e-3){
						for(unsigned int i=0;i<dim;i++){
							additional_data(local_dof_indices[vetex_index*dim+i])= -advection_array[node][dim+i];
						}
					}
				}	
			}
		}
	}
	additional_data.compress(VectorOperation::insert);
}
template class initBoundValProbs<1>;
template class initBoundValProbs<2>;
template class initBoundValProbs<3>;
