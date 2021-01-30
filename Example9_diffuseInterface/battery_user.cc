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
void battery<dim>::identify_diffuse_interface()
{
  double iso_value = 0.5;
  int total_cell_num = this->triangulation.n_active_cells();
  cell_SDdata.resize(total_cell_num);
  typename hp::DoFHandler<dim>::active_cell_iterator cell = dof_handler_interface->begin_active(), endc=dof_handler_interface->end();
  for (;cell!=endc; ++cell){
		if (cell->subdomain_id() == this->this_mpi_process){
      int cell_id = cell->active_cell_index();
      cell_SDdata[cell_id].cell_id = cell_id;
      Point<dim> cell_center=cell->center();
			//if(cell_center.distance(origin)<1){
      const unsigned int dofs_per_cell = cell->get_fe().dofs_per_cell;
      // assume
      bool greater_flag = true;
      bool smaller_flag = true;
      std::vector<unsigned int> local_dof_indices (dofs_per_cell);
      cell->get_dof_indices (local_dof_indices);
      // 1 dof per vertex, dofs_per_cell == vertex per cell
      for (unsigned int i=0; i<dofs_per_cell; ++i) {
        if (diffuse_interface(local_dof_indices[i]) >= iso_value){
          smaller_flag = false;
        };
        if (diffuse_interface(local_dof_indices[i]) < iso_value){
          greater_flag = false;
        };
      }
      if ((not greater_flag) and (not smaller_flag))
      {
        cell_SDdata[cell_id].is_interface_element = true;

        // get the side of the local and global node number
        for (unsigned int i=0; i<dofs_per_cell; ++i) {
          if (diffuse_interface(local_dof_indices[i]) >= iso_value){
            cell_SDdata[cell_id].lnode_plus.push_back(i);
          };
          if (diffuse_interface(local_dof_indices[i]) < iso_value){
            cell_SDdata[cell_id].lnode_minus.push_back(i);
          };
        }

        std::vector<types::global_dof_index> local_face_dof_indices(fe_system_interface[0]->dofs_per_face);
        int count = 0;
        for (unsigned int f = 0; f < GeometryInfo<dim>::faces_per_cell; ++f) {
          cell->face(f)->get_dof_indices(local_face_dof_indices, 0);
          double c_1 = diffuse_interface(local_face_dof_indices[0]);
          double c_2 = diffuse_interface(local_face_dof_indices[1]);
          if (c_1 == iso_value and c_2 == iso_value){
            // for the case where the edge of element is aligned with the contour. 
            cell_SDdata[cell_id].edge1_node = cell->face(f)->vertex(0) ;
            cell_SDdata[cell_id].edge2_node = cell->face(f)->vertex(1) ;
          }
          else if ((c_1 >= iso_value and c_2 < iso_value) || (c_1 <= iso_value and c_2 > iso_value)){
            // Without equal sign between c_2 vs iso_value can prevent assigning the same node with c=iso_value to both edge1_node and edge2_node
            if (count == 0){
              cell_SDdata[cell_id].edge1_node1 = cell->face(f)->vertex(0) ;
              cell_SDdata[cell_id].edge1_node2 = cell->face(f)->vertex(1) ;
              cell_SDdata[cell_id].edge1_local_s = (c_1 - iso_value)/(c_1 - c_2);
              cell_SDdata[cell_id].edge1_node =  cell_SDdata[cell_id].edge1_node1 - cell_SDdata[cell_id].edge1_local_s * (cell_SDdata[cell_id].edge1_node1 - cell_SDdata[cell_id].edge1_node2);
            }
            else if (count == 1){
              cell_SDdata[cell_id].edge2_node1 = cell->face(f)->vertex(0) ;
              cell_SDdata[cell_id].edge2_node2 = cell->face(f)->vertex(1) ;
              cell_SDdata[cell_id].edge2_local_s = (c_1 - iso_value)/(c_1 - c_2);
              cell_SDdata[cell_id].edge2_node =  cell_SDdata[cell_id].edge2_node1 - cell_SDdata[cell_id].edge2_local_s * (cell_SDdata[cell_id].edge2_node1 - cell_SDdata[cell_id].edge2_node2);

              //std::cout << "--vertex index " << cell->face(f)->vertex_index(0) << std::endl;
            }
            count++;
          }
        }

        //// perform the calculation for interface.
        //std::cout 
          //<< " cell_id = " << cell_id 
          //<< " center = " << cell_center 
          //<< " " << cell_SDdata[cell_id].edge1_node1
          //<< " " << cell_SDdata[cell_id].edge1_node2
          //<< " " << cell_SDdata[cell_id].edge1_node
          //<< " " << cell_SDdata[cell_id].edge1_local_s
          //<< " " << cell_SDdata[cell_id].edge2_node1
          //<< " " << cell_SDdata[cell_id].edge2_node2
          //<< " " << cell_SDdata[cell_id].edge2_node
          //<< " " << cell_SDdata[cell_id].edge2_local_s
          //<< std::endl;

        double elem_length = cell_SDdata[cell_id].edge1_node.distance(cell_SDdata[cell_id].edge2_node);

        Triangulation<1> triangulation_1d;
        GridGenerator::hyper_cube	(	triangulation_1d, 0.,  elem_length);
        DoFHandler<1>      dof_handler(triangulation_1d);
        int problem_dof = 1;
        int poly_order = 1;
        int quad_order = 2;
        FESystem<1> fe (FE_Q<1>(poly_order), problem_dof);
        dof_handler.distribute_dofs (fe);

        QGauss<1>  quadrature_formula(quad_order);
        FEValues<1> fe_values_1d (fe, quadrature_formula,
                                 update_values   | update_gradients |
                                 update_quadrature_points | update_JxW_values);

        typename DoFHandler<1>::active_cell_iterator cell_1d = dof_handler.begin_active(),
                                                       endc_1d = dof_handler.end();

        double vol = 0.0;
        cell_SDdata[cell_id].shape_value_1d.reinit(2,quadrature_formula.size());
        cell_SDdata[cell_id].jxw_1d.reinit(quadrature_formula.size());

        for (; cell_1d!=endc_1d; ++cell_1d)
        {
            fe_values_1d.reinit (cell_1d);
            for (unsigned int i = 0; i < 2; ++i) {
                for (unsigned q=0; q< quadrature_formula.size(); ++q)
                {
                    //std::cout << " q " << q << std::endl;
                    vol += fe_values_1d.JxW(q);

                    cell_SDdata[cell_id].shape_value_1d(i, q) = fe_values_1d.shape_value(i, q);
                    cell_SDdata[cell_id].jxw_1d(q) = fe_values_1d.JxW(q);
                    //r_local[i] += fe_values_1d.shape_value(i, q) * dRc * fe_values_1d.JxW(q);
                } // q_point
            }
        }
        //std::cout <<  " total elem # = " << triangulation_1d.n_active_cells() << " length: " << elem_length << " vol " << vol<< std::endl;
      }
			//}
		}		
	}

  {
    typename hp::DoFHandler<dim>::active_cell_iterator cell = dof_handler_interface->begin_active(), endc=dof_handler_interface->end();
    for (;cell!=endc; ++cell){
	  	if (cell->subdomain_id() == this->this_mpi_process){
        int cell_id = cell->active_cell_index();
        if  (cell_SDdata[cell_id].is_interface_element)
        {
          Point<dim> cell_center=cell->center();
          //std::cout << " cell_id = " << cell_id << " center = " << cell_center << " has interface! " << std::endl;
        }
      }
    }
  }

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
