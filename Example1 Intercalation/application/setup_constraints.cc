#include "initBoundValProbs.h"

template <int dim>
void initBoundValProbs<dim>::setup_constraints()
{	
	hpFEM<dim>::constraints.clear ();

	
	int totalDOF=this->totalDOF(primary_variables);
  std::vector<bool> c_1_component (totalDOF, false); c_1_component[c1_dof]=true; 
	std::vector<bool> c_2_component (totalDOF, false); c_2_component[c2_dof]=true; 
  std::vector<bool> x_component (totalDOF, false); x_component[u_dof]=true; 
  std::vector<bool> y_component (totalDOF, false); y_component[u_dof+1]=true;  
  std::vector<bool> z_component (totalDOF, false); if(dim==3) z_component[u_dof+2]=true;  

  DoFTools::make_hanging_node_constraints (hpFEM<dim>::dof_handler, hpFEM<dim>::constraints);
	
	//VectorTools:: interpolate_boundary_values (hpFEM<dim>::dof_handler, 2, ZeroFunction<dim> (totalDOF),hpFEM<dim>::constraints, c_1_component);
	//VectorTools:: interpolate_boundary_values (hpFEM<dim>::dof_handler, 2, ZeroFunction<dim> (totalDOF),hpFEM<dim>::constraints, c_2_component);
	//VectorTools:: interpolate_boundary_values (hpFEM<dim>::dof_handler, 2, ZeroFunction<dim> (totalDOF),hpFEM<dim>::constraints, x_component);
	//VectorTools:: interpolate_boundary_values (hpFEM<dim>::dof_handler, 2, ZeroFunction<dim> (totalDOF),hpFEM<dim>::constraints, y_component);
	//VectorTools:: interpolate_boundary_values (hpFEM<dim>::dof_handler, 2, ZeroFunction<dim> (totalDOF),hpFEM<dim>::constraints, z_component);
  //VectorTools:: interpolate_boundary_values (hpFEM<dim>::dof_handler, 1, ZeroFunction<dim> (totalDOF),hpFEM<dim>::constraints, x_component);
  // VectorTools:: interpolate_boundary_values (hpFEM<dim>::dof_handler, 1, ZeroFunction<dim> (totalDOF),hpFEM<dim>::constraints, y_component);
  //if(dim==3) VectorTools:: interpolate_boundary_values (hpFEM<dim>::dof_handler, 1, ZeroFunction<dim> (totalDOF),hpFEM<dim>::constraints, z_component);
  
  std::vector<types::global_dof_index> local_face_dof_indices0 (fe_system[Cortex_fe]->dofs_per_face);
  std::vector<types::global_dof_index> local_face_dof_indices1 (fe_system[Subcortex_fe]->dofs_per_face);
	
  typename hp::DoFHandler<dim>::active_cell_iterator cell = this->dof_handler.begin_active(), endc=this->dof_handler.end();
  for (;cell!=endc; ++cell){
    //add constraints for dofs at interface between Cortex and Ventricle
    if(cell->material_id()==Cortex_id ){
      for (unsigned int f=0; f<GeometryInfo<dim>::faces_per_cell; ++f){
				bool at_interface=false;
				if (cell->at_boundary(f) == false){
	  			if(cell->neighbor(f)->material_id()==Ventricle_id and cell->neighbor(f)->has_children() == false){
	   			 	at_interface=true;
	 			 	}	
	  		 // check for the cell that has children cells
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
	  			cell->face(f)->get_dof_indices (local_face_dof_indices0, Cortex_fe);
	  			for (unsigned int i=0; i<local_face_dof_indices0.size(); ++i){
				  const unsigned int ck = fe_system[Cortex_fe]->face_system_to_component_index(i).first;
				  if(ck>=0 and ck<10*dim)hpFEM<dim>::constraints.add_line (local_face_dof_indices0[i]);//add constrain line for all dofs
	      		 	//constraints.set_inhomogeneity(local_face_dof_indices1[i],0);//set inhomogeneity if needed
	  		 	}
			 	}		
      }			
    }
		//add constraints for dofs at interface between Subcortex and Ventricle
    if(cell->material_id()==Subcortex_id){
      for (unsigned int f=0; f<GeometryInfo<dim>::faces_per_cell; ++f){
				bool at_interface=false;
				//Do not check boundary cell. If set_boundary_id  is not defined use somthing else 
				if (cell->at_boundary(f) == false){
	  			if(cell->neighbor(f)->material_id()==Ventricle_id and cell->neighbor(f)->has_children() == false){
	   			 	at_interface=true;
	 			 	}	
	  		 // check for the cell that has children cells
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
	  			cell->face(f)->get_dof_indices (local_face_dof_indices1, Subcortex_fe);
	  			for (unsigned int i=0; i<local_face_dof_indices1.size(); ++i){
						const unsigned int ck = fe_system[Subcortex_fe]->face_system_to_component_index(i).first;
						if(ck>=0 and ck<10*dim)hpFEM<dim>::constraints.add_line (local_face_dof_indices1[i]);//add constrain line for all dofs
	      		 	//constraints.set_inhomogeneity(local_face_dof_indices1[i],0);//set inhomogeneity if needed
	  		 	}
			 	}		
      }			
    }
  }
	
  hpFEM<dim>::constraints.close ();
  pcout << "   Number of active cells:       " << hpFEM<dim>::triangulation.n_active_cells() << std::endl;
  pcout << "   Number of degrees of freedom: " << hpFEM<dim>::dof_handler.n_dofs() << std::endl; 
  
}


template class initBoundValProbs<1>;
template class initBoundValProbs<2>;
template class initBoundValProbs<3>;
