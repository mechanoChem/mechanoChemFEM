#include "initBoundValProbs.h"

template <int dim>
void initBoundValProbs<dim>::apply_initial_condition()
{
	params->enter_subsection("Geometry");	
	double Y_end=params->get_double("Y_end");
	params->leave_subsection();	
  typename hp::DoFHandler<dim>::active_cell_iterator cell = this->dof_handler.begin_active(), endc=this->dof_handler.end();
  FEFaceValues<dim> electrode_fe_face_values (*fe_system[electrode_fe], *common_face_quadrature, update_values | update_quadrature_points | update_JxW_values | update_normal_vectors | update_gradients);
  for (;cell!=endc; ++cell){
		if (cell->subdomain_id() == this_mpi_process){
	    const unsigned int dofs_per_cell = cell->get_fe().dofs_per_cell;
	    std::vector<unsigned int> local_dof_indices (dofs_per_cell);
	    cell->get_dof_indices (local_dof_indices);
			for (unsigned int faceID=0; faceID<2*dim; faceID++){
	    	const Point<dim> face_center_upper = cell->face(faceID)->center();
				if(face_center_upper[1]==Y_end){
	    		for (unsigned int i=0; i<dofs_per_cell; ++i) {
						electrode_fe_face_values.reinit (cell, faceID);
	        	const unsigned int ckf = electrode_fe_face_values.get_fe().system_to_component_index(i).first;
	        	if(ckf==u_dof+1) {
		  				if(current_increment>=1) solution(local_dof_indices[i])=-0.5;
		  				if(current_increment>=2) solution(local_dof_indices[i])=-0.8;
		  				if(current_increment>=3) solution(local_dof_indices[i])=-1;
		  				if(current_increment>=4) solution(local_dof_indices[i])=-1.2;
		  				if(current_increment>=5) solution(local_dof_indices[i])=-1.4;
	          	if(current_increment>=6) solution(local_dof_indices[i])=-1.6;
	          	if(current_increment>=7) solution(local_dof_indices[i])=-1.8;
	          	if(current_increment>=8) solution(local_dof_indices[i])=-2;
		  				if(current_increment>=9) solution(local_dof_indices[i])=-2.2;
		  				// if(current_increment>=0.8) U(local_dof_indices[i])=-2.4;
		  				//if(current_increment>=0.3) IpA=44;
						}
					}
				}
			}
		}
	}
	solution.compress(VectorOperation::insert);
}

template class initBoundValProbs<1>;
template class initBoundValProbs<2>;
template class initBoundValProbs<3>;