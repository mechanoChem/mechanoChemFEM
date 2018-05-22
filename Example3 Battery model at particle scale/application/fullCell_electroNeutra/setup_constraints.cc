#include "initBoundValProbs.h"

template <int dim>
void initBoundValProbs<dim>::setup_constraints()
{	
	hpFEM<dim>::constraints.clear ();
	
	
	int totalDOF=this->totalDOF(primary_variables);
  std::vector<bool> c_1_component (totalDOF, false); c_1_component[primary_variables_dof[1] ]=true; 

	constrain_index.resize(2);
	
  std::vector<bool> x_component (totalDOF, false); x_component[u_dof]=true; 
  std::vector<bool> y_component (totalDOF, false); y_component[u_dof+1]=true;  
	std::vector<bool> z_component (totalDOF, false); if(dim==3) z_component[u_dof+2]=true;  
	
  std::vector<bool> Vx_component (totalDOF, false); Vx_component[v_dof]=true; 
  std::vector<bool> Vy_component (totalDOF, false); Vy_component[v_dof+1]=true;  
	std::vector<bool> Vz_component (totalDOF, false); if(dim==3) Vz_component[v_dof+2]=true;  
	std::vector<bool> P_component (totalDOF, false);  P_component[p_dof]=true;  
	
  std::vector<bool> c_li_component (totalDOF, false); c_li_component[c_li_dof]=true;
	std::vector<bool> phi_s_component (totalDOF, false); phi_s_component[phi_s_dof]=true;
	
  std::vector<bool> c_li_plus_component (totalDOF, false); c_li_plus_component[c_li_plus_dof]=true;
  std::vector<bool> phi_e_component (totalDOF, false); phi_e_component[phi_e_dof]=true;
	
  std::vector<bool> Um_x_component (totalDOF, false); Um_x_component[u_m_dof]=true; 
  std::vector<bool> Um_y_component (totalDOF, false); Um_y_component[u_m_dof+1]=true;  
	std::vector<bool> Um_z_component (totalDOF, false); if(dim==3) Um_z_component[u_m_dof+2]=true;  
	
  DoFTools::make_hanging_node_constraints (hpFEM<dim>::dof_handler, hpFEM<dim>::constraints);

  VectorTools:: interpolate_boundary_values (hpFEM<dim>::dof_handler, 2, ZeroFunction<dim> (totalDOF),hpFEM<dim>::constraints, x_component);
  VectorTools:: interpolate_boundary_values (hpFEM<dim>::dof_handler, 2, ZeroFunction<dim> (totalDOF),hpFEM<dim>::constraints, y_component);
  if(dim==3) VectorTools:: interpolate_boundary_values (hpFEM<dim>::dof_handler, 2, ZeroFunction<dim> (totalDOF),hpFEM<dim>::constraints, z_component);
  VectorTools:: interpolate_boundary_values (hpFEM<dim>::dof_handler, dim+2, ZeroFunction<dim> (totalDOF),hpFEM<dim>::constraints, x_component);
  VectorTools:: interpolate_boundary_values (hpFEM<dim>::dof_handler, dim+2, ZeroFunction<dim> (totalDOF),hpFEM<dim>::constraints, y_component);
  if(dim==3) VectorTools:: interpolate_boundary_values (hpFEM<dim>::dof_handler, dim+2, ZeroFunction<dim> (totalDOF),hpFEM<dim>::constraints, z_component);
	
	VectorTools:: interpolate_boundary_values (hpFEM<dim>::dof_handler, 2, ZeroFunction<dim> (totalDOF),hpFEM<dim>::constraints,phi_s_component);
	
	//side surfaces
  VectorTools:: interpolate_boundary_values (hpFEM<dim>::dof_handler, 1, ZeroFunction<dim> (totalDOF),hpFEM<dim>::constraints, x_component);
  VectorTools:: interpolate_boundary_values (hpFEM<dim>::dof_handler, dim+1, ZeroFunction<dim> (totalDOF),hpFEM<dim>::constraints, x_component);	

  //VectorTools:: interpolate_boundary_values (hpFEM<dim>::dof_handler, 1, ZeroFunction<dim> (totalDOF),hpFEM<dim>::constraints, y_component);
    //VectorTools:: interpolate_boundary_values (hpFEM<dim>::dof_handler, dim+1, ZeroFunction<dim> (totalDOF),hpFEM<dim>::constraints, y_component);
	
  VectorTools:: interpolate_boundary_values (hpFEM<dim>::dof_handler, 1, ZeroFunction<dim> (totalDOF),hpFEM<dim>::constraints, Vx_component);
  VectorTools:: interpolate_boundary_values (hpFEM<dim>::dof_handler, dim+1, ZeroFunction<dim> (totalDOF),hpFEM<dim>::constraints, Vx_component);
	
   VectorTools:: interpolate_boundary_values (hpFEM<dim>::dof_handler, 1, ZeroFunction<dim> (totalDOF),hpFEM<dim>::constraints, Um_x_component);
  VectorTools:: interpolate_boundary_values (hpFEM<dim>::dof_handler, dim+1, ZeroFunction<dim> (totalDOF),hpFEM<dim>::constraints, Um_x_component);


  
  std::vector<types::global_dof_index> local_face_dof_indices0 (fe_system[active_material_fe]->dofs_per_face);
  std::vector<types::global_dof_index> local_face_dof_indices1 (fe_system[electrolyte_fe]->dofs_per_face);
	std::vector<types::global_dof_index> local_face_dof_indices2 (fe_system[current_collector_fe]->dofs_per_face);
	std::vector<types::global_dof_index> local_face_dof_indices3 (fe_system[solid_fe]->dofs_per_face);
	
  typename hp::DoFHandler<dim>::active_cell_iterator cell = hpFEM<dim>::dof_handler.begin_active(), endc=hpFEM<dim>::dof_handler.end();
  for (;cell!=endc; ++cell){
		if(cell->material_id()==electrolyte_id){
			for (unsigned int f=0; f<GeometryInfo<dim>::faces_per_cell; ++f){
				bool interface_activeMaterial=false;
				bool interface_currentCollector=false;
				bool interface_solid=false;
				int a=cell->material_id();
				 const Point<dim> cell_center = cell->center();//face n=y-;
				if (cell->at_boundary(f) == false){
					int a= cell->neighbor(f)->material_id();
					if(cell->neighbor(f)->material_id()==active_material_id and cell->neighbor(f)->has_children() == false){
						interface_activeMaterial=true;
					}
					else if(cell->neighbor(f)->material_id()==current_collector_id and cell->neighbor(f)->has_children() == false){
							interface_currentCollector=true;
					}
					else if(cell->neighbor(f)->material_id()==solid_id and cell->neighbor(f)->has_children() == false){
						interface_solid=true;
					}
					
				  else if (cell->neighbor(f)->has_children() == true){
						for (unsigned int sf=0; sf<cell->face(f)->n_children(); ++sf){
							if (cell->neighbor_child_on_subface(f, sf)->material_id() == active_material_id){
							  interface_activeMaterial=true;                                                               
							  break;
						  }
						  if (cell->neighbor_child_on_subface(f, sf)->material_id() == current_collector_id ){
							  interface_currentCollector=true;                                                               
							  break;
						  }
						  if (cell->neighbor_child_on_subface(f, sf)->material_id() == solid_id ){
							  interface_solid=true;                                                               
							  break;
						  }
					  }
				  }
			  }
				//at interface_activeMaterial
				if (interface_activeMaterial){
					cell->face(f)->get_dof_indices (local_face_dof_indices1, electrolyte_fe);
					cell->neighbor(f)->face(cell->neighbor_of_neighbor(f))->get_dof_indices (local_face_dof_indices0, active_material_fe);
					int node_u_m_applied=0;				
					for (unsigned int i=0; i<local_face_dof_indices1.size(); ++i){
				    const unsigned int ck = fe_system[electrolyte_fe]->face_system_to_component_index(i).first-v_dof;
				    if (ck>=0 && ck<dim){
							hpFEM<dim>::constraints.add_line (local_face_dof_indices1[i]);//velocity							
							for (unsigned int j=0; j<local_face_dof_indices0.size(); ++j){
								const unsigned int ck2 = fe_system[active_material_fe]->face_system_to_component_index(j).first;
								if(ck2==ck) hpFEM<dim>::constraints.add_entry (local_face_dof_indices1[i], local_face_dof_indices0[j], 1.0/current_dt/(2*dim-2) );
							}					
						}
						const unsigned int ck1 = fe_system[electrolyte_fe]->face_system_to_component_index(i).first-u_m_dof;	
				    if (ck1>=0 && ck1<dim ){
							int dof_per_node=local_face_dof_indices0.size()/4;
							hpFEM<dim>::constraints.add_line (local_face_dof_indices1[i]);//u_mesh	
							for (unsigned int j=(node_u_m_applied/2)*dof_per_node; j<local_face_dof_indices0.size(); ++j){
								const unsigned int ck3 = fe_system[active_material_fe]->face_system_to_component_index(j).first;
								if(ck3==ck1){
									constrain_index[0].push_back(local_face_dof_indices1[i]);
									constrain_index[1].push_back(local_face_dof_indices0[j]);
									node_u_m_applied++;
									break;
								}
							}						
						}
				 	}
				}
				//at interface_currentCollector
				else if (interface_currentCollector){		
					cell->face(f)->get_dof_indices (local_face_dof_indices1, 1);
					cell->neighbor(f)->face(cell->neighbor_of_neighbor(f))->get_dof_indices (local_face_dof_indices2, 2);
					
					int node_u_m_applied=0;				
					for (unsigned int i=0; i<local_face_dof_indices1.size(); ++i){
				    const unsigned int ck = fe_system[electrolyte_fe]->face_system_to_component_index(i).first-v_dof;
				    if (ck>=0 && ck<dim){
							hpFEM<dim>::constraints.add_line (local_face_dof_indices1[i]);//velocity							
							for (unsigned int j=0; j<local_face_dof_indices2.size(); ++j){
								const unsigned int ck2 = fe_system[current_collector_fe]->face_system_to_component_index(j).first;
								if(ck2==ck) hpFEM<dim>::constraints.add_entry (local_face_dof_indices1[i], local_face_dof_indices2[j], 1.0/current_dt/(2*dim-2) );
							}					
						}
						const unsigned int ck1 = fe_system[electrolyte_fe]->face_system_to_component_index(i).first-u_m_dof;	
				    if (ck1>=0 && ck1<dim ){
							int dof_per_node=local_face_dof_indices2.size()/4;
							hpFEM<dim>::constraints.add_line (local_face_dof_indices1[i]);//u_mesh														
							for (unsigned int j=(node_u_m_applied/2)*dof_per_node; j<local_face_dof_indices2.size(); ++j){
								const unsigned int ck3 = fe_system[current_collector_fe]->face_system_to_component_index(j).first;
								if(ck3==ck1){
									constrain_index[0].push_back(local_face_dof_indices1[i]);
									constrain_index[1].push_back(local_face_dof_indices2[j]);
									node_u_m_applied++;
									break;
								}
							}						
						}	  
				 	}
				}
				//at interface_solid
				else if (interface_solid){
					cell->face(f)->get_dof_indices (local_face_dof_indices1, 1);
					cell->neighbor(f)->face(cell->neighbor_of_neighbor(f))->get_dof_indices (local_face_dof_indices3, 3);
					
					int node_u_m_applied=0;				
					for (unsigned int i=0; i<local_face_dof_indices1.size(); ++i){
				    const unsigned int ck = fe_system[electrolyte_fe]->face_system_to_component_index(i).first-v_dof;
				    if (ck>=0 && ck<dim){
							hpFEM<dim>::constraints.add_line (local_face_dof_indices1[i]);//velocity							
							for (unsigned int j=0; j<local_face_dof_indices3.size(); ++j){
								const unsigned int ck2 = fe_system[solid_fe]->face_system_to_component_index(j).first;
								if(ck2==ck) hpFEM<dim>::constraints.add_entry (local_face_dof_indices1[i], local_face_dof_indices3[j], 1.0/current_dt/(2*dim-2) );
							}					
						}
						
						const unsigned int ck1 = fe_system[electrolyte_fe]->face_system_to_component_index(i).first-u_m_dof;	
				    if (ck1>=0 && ck1<dim ){
							int dof_per_node=local_face_dof_indices3.size()/4;
							hpFEM<dim>::constraints.add_line (local_face_dof_indices1[i]);//u_mesh														
							for (unsigned int j=(node_u_m_applied/2)*dof_per_node; j<local_face_dof_indices3.size(); ++j){
								const unsigned int ck3 = fe_system[solid_fe]->face_system_to_component_index(j).first;
								if(ck3==ck1){
									constrain_index[0].push_back(local_face_dof_indices1[i]);
									constrain_index[1].push_back(local_face_dof_indices3[j]);
									node_u_m_applied++;
									break;
								}
							}						
						}	  
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
