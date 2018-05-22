#include "initBoundValProbs.h"
#include "random_initial_condition.h"

template <int dim>
void initBoundValProbs<dim>::apply_initial_condition()
{ 
	pcout << "applying initial condition\n";
	params->enter_subsection("Concentration");
	double c1_ini=params->get_double("c1_ini");
	double c2_ini=params->get_double("c2_ini");
	
	double c1_ini_interface=params->get_double("c1_ini_interface");
	double c2_ini_interface=params->get_double("c2_ini_interface");
  params->leave_subsection();		
  solution_0=0;
	
	hp::FEValues<dim> hp_fe_values (fe_collection, q_collection, update_values | update_quadrature_points  | update_JxW_values | update_gradients);
	FEFaceValues<dim> fe_face_values_0 (*fe_system[Cortex_fe], *common_face_quadrature, update_values | update_quadrature_points | update_JxW_values | update_normal_vectors | update_gradients);
	FEFaceValues<dim> fe_face_values_1 (*fe_system[Subcortex_fe], *common_face_quadrature, update_values | update_quadrature_points | update_JxW_values | update_normal_vectors | update_gradients);

  typename hp::DoFHandler<dim>::active_cell_iterator cell = this->dof_handler.begin_active(), endc=this->dof_handler.end();
  for (;cell!=endc; ++cell){
		if (cell->subdomain_id() == this_mpi_process){
    	hp_fe_values.reinit (cell);
    	const FEValues<dim> &fe_values = hp_fe_values.get_present_fe_values();
    	const unsigned int dofs_per_cell = cell->get_fe().dofs_per_cell;

			std::vector<unsigned int> local_dof_indices (dofs_per_cell);
 		 	cell->get_dof_indices (local_dof_indices);
			
    	for (unsigned int i=0; i<dofs_per_cell; ++i) {
      	const unsigned int ck = fe_values.get_fe().system_to_component_index(i).first;
				if(ck==c1_dof){	
					solution_0(local_dof_indices[i])=c1_ini;//C1
		      	    	}
				
				else if(ck==c2_dof){	
					solution_0(local_dof_indices[i])=c2_ini;//C2
	    	}
			}
		}
	}
  cell = this->dof_handler.begin_active();
  for (;cell!=endc; ++cell){
		if (cell->subdomain_id() == this_mpi_process){
    	hp_fe_values.reinit (cell);
			
		  std::vector<types::global_dof_index> local_face_dof_indices0 (fe_system[Cortex_fe]->dofs_per_face);
		  std::vector<types::global_dof_index> local_face_dof_indices1 (fe_system[Subcortex_fe]->dofs_per_face);	
	    if(cell->material_id()==Cortex_id ){
	      for (unsigned int f=0; f<GeometryInfo<dim>::faces_per_cell; ++f){
					bool at_interface=false;
					//if(cell->face(f)->boundary_id()==2) at_interface=true;
					if (cell->at_boundary(f) == false){
		  			if(cell->neighbor(f)->material_id()==Ventricle_id and cell->neighbor(f)->has_children() == false){
		   			 	at_interface=true;
		 			 	}	
		  			else if (cell->neighbor(f)->has_children() == true){
		    			for (unsigned int sf=0; sf<cell->face(f)->n_children(); ++sf){
		      			if (cell->neighbor_child_on_subface(f, sf)->material_id()==Ventricle_id){
									at_interface=true;                                                              
									break;
		      			}
		    			}
		  			}	
					}
	
					if(at_interface){
		  			cell->face(f)->get_dof_indices (local_face_dof_indices0,Cortex_fe);
		  			for (unsigned int i=0; i<local_face_dof_indices0.size(); ++i){
							const unsigned int ck = fe_system[Cortex_fe]->face_system_to_component_index(i).first;
		    		 	if (ck==c1_dof){
								solution_0(local_face_dof_indices0[i])=c1_ini_interface;
		    		 	}
		    		 	if (ck==c2_dof){
								solution_0(local_face_dof_indices0[i])=c2_ini_interface;
		    		 	}
		  		 	}
				 	}		
	      }			
	    }
		
	    if(cell->material_id()==Subcortex_id){
	      for (unsigned int f=0; f<GeometryInfo<dim>::faces_per_cell; ++f){
					bool at_interface=false;
					//if(cell->face(f)->boundary_id()==2) at_interface=true;
					if (cell->at_boundary(f) == false){
		  			if(cell->neighbor(f)->material_id()==Ventricle_id and cell->neighbor(f)->has_children() == false){
		   			 	at_interface=true;
		 			 	}	
		  			else if (cell->neighbor(f)->has_children() == true){
		    			for (unsigned int sf=0; sf<cell->face(f)->n_children(); ++sf){
		      			if (cell->neighbor_child_on_subface(f, sf)->material_id()==Ventricle_id){
									at_interface=true;                                                              
									break;
		      			}
		    			}
		  			}	
					}
	
					if(at_interface){
		  			cell->face(f)->get_dof_indices (local_face_dof_indices1, Subcortex_fe);//first Fe_system is 0
		  			for (unsigned int i=0; i<local_face_dof_indices1.size(); ++i){
							const unsigned int ck = fe_system[Subcortex_fe]->face_system_to_component_index(i).first;
		    		 	if (ck==c1_dof ){
								solution_0(local_face_dof_indices1[i])=c1_ini_interface;
		    		 	}
		    		 	if (ck==c2_dof ){
								solution_0(local_face_dof_indices1[i])=c2_ini_interface;
		    		 	}
		  		 	}
				 	}		
	      }		
			}			
		}
	}
	
	solution_0.compress(VectorOperation::insert);
}


template class initBoundValProbs<1>;
template class initBoundValProbs<2>;
template class initBoundValProbs<3>;
