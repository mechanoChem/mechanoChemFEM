/*
zhenlin wang 2019
*coupled diffusion reaction
*/
#include "battery.h"
#include "nodalField.h"

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

template <int dim>
void battery<dim>::identify_diffuse_interface()
{
  int primary_dof = -1;
  int opposite_flux_dof = -1;
	if(battery_fields.active_fields_index["Diffuse_interface"]>-1) primary_dof=battery_fields.active_fields_index["Diffuse_interface"];
	if(battery_fields.active_fields_index["Lithium_cation"]>-1) opposite_flux_dof=battery_fields.active_fields_index["Lithium_cation"];
	
  std::cout << "---------- primary dof for diffusive interface ------ " << primary_dof  << " opposite dof " << opposite_flux_dof << std::endl;

  hp::FEValues<dim> hp_fe_values (this->fe_collection, this->q_collection, update_values | update_quadrature_points  | update_JxW_values | update_gradients);	

  Vector<double> localized_U(this->solution_prev);
  double iso_value = 0.5;
  int total_cell_num = this->triangulation.n_active_cells();
  cell_SDdata.resize(total_cell_num);
  typename hp::DoFHandler<dim>::active_cell_iterator cell = this->dof_handler.begin_active(), endc=this->dof_handler.end();
  for (;cell!=endc; ++cell){
		if (cell->subdomain_id() == this->this_mpi_process){

			hp_fe_values.reinit (cell);
    	const FEValues<dim> &fe_values = hp_fe_values.get_present_fe_values();

      int cell_id = cell->active_cell_index();
      cell_SDdata[cell_id].cell_id = cell_id;
      cell_SDdata[cell_id].opposite_flux_dof = opposite_flux_dof;

      Point<dim> cell_center=cell->center();
			//if(cell_center.distance(origin)<1){
      const unsigned int dofs_per_cell = cell->get_fe().dofs_per_cell;
      // assume
      bool greater_flag = true;
      bool smaller_flag = true;
      std::vector<unsigned int> local_dof_indices (dofs_per_cell);
      cell->get_dof_indices (local_dof_indices);
      std::vector<double> local_diffuse_interface;

      for (unsigned int i=0; i<dofs_per_cell; ++i) {
        unsigned int ck = fe_values.get_fe().system_to_component_index(i).first - primary_dof;
        if (ck == 0) {
          //std::cout 
            //<< "---------- primary dof for diffusive interface ------ i = " << i 
            //<< " val = " << localized_U(local_dof_indices[i])
            //<< std::endl;
          local_diffuse_interface.push_back(localized_U(local_dof_indices[i]));
        }
      }

      // 1 dof per vertex, dofs_per_cell == vertex per cell
      for (unsigned int i=0; i<local_diffuse_interface.size(); ++i) {
        if (local_diffuse_interface[i] >= iso_value){
          smaller_flag = false;
        };
        if (local_diffuse_interface[i] < iso_value){
          greater_flag = false;
        };
      }

      if (greater_flag)
      {
        // active particle = 1
        cell->set_material_id(1);
      }
      if (smaller_flag)
      {
        // electrolyte = 2
        cell->set_material_id(2);
      }

      if ((not greater_flag) and (not smaller_flag))
      {
        cell_SDdata[cell_id].is_interface_element = true;

        // get the side of the local and global node number
        for (unsigned int i=0; i<local_diffuse_interface.size(); ++i) {
          if (local_diffuse_interface[i] >= iso_value){
            cell_SDdata[cell_id].lnode_plus.push_back(i);
          };
          if (local_diffuse_interface[i] < iso_value){
            cell_SDdata[cell_id].lnode_minus.push_back(i);
          };
        }

        std::vector<types::global_dof_index> local_face_dof_indices(this->fe_system[0]->dofs_per_face);
        int count = 0;
        for (unsigned int f = 0; f < GeometryInfo<dim>::faces_per_cell; ++f) {
          cell->face(f)->get_dof_indices(local_face_dof_indices, 0);
          double c_1 = 0.0;
          double c_2 = 0.0;
          std::vector<double> local_local_diffuse_interface_face;
          for (unsigned int i = 0; i < local_face_dof_indices.size(); ++i) {
            const unsigned int ck = this->fe_system[0]->face_system_to_component_index(i).first - primary_dof;
            if (ck == 0)
              local_local_diffuse_interface_face.push_back(localized_U(local_face_dof_indices[i]));
          }

          c_1 = local_local_diffuse_interface_face[0];
          c_2 = local_local_diffuse_interface_face[1];

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
    typename hp::DoFHandler<dim>::active_cell_iterator cell = this->dof_handler.begin_active(), endc=this->dof_handler.end();
    for (;cell!=endc; ++cell){
	  	if (cell->subdomain_id() == this->this_mpi_process){
        int cell_id = cell->active_cell_index();
        if  (cell_SDdata[cell_id].is_interface_element)
        {
          Point<dim> cell_center=cell->center();
          //std::cout 
            //<< " cell_id = " << cell_id << " center = " 
            //<< cell_center << " has interface!  R = " 
            //<< std::sqrt((cell_center(0)-2)*(cell_center(0)-2) + (cell_center(1)-2)*(cell_center(1)-2))
            //<< std::endl;
        }
      }
    }
  }

}



template class battery<1>;
template class battery<2>;
template class battery<3>;

template <int dim>
void InitialConditions<dim>::vector_value (const Point<dim>   &p, Vector<double>   &values) const{
	Assert (values.size() == totalDOF, ExcDimensionMismatch (values.size(), totalDOF));
	values=0;
	double r=0.8;
	double bandwitdh=0.1;
	std::vector<std::vector<double>> origin_list={{1.5,1},{3,1},{4.5,1} };
	for(unsigned int i=0;i<primary_variables.size();i++){
		if(std::strcmp(primary_variables[i][0].c_str(),"Diffuse_interface")==0){
			bool inside_flag=false;
			bool interface_flag=false;
			int particle_index=-1;
			double current_dis=1.0e4;
			for (unsigned int ori_index=0;ori_index<origin_list.size();ori_index++){
				Point<dim> origin(origin_list[ori_index][0],origin_list[ori_index][1]);
				if(p.distance(origin)<r) {interface_flag=false; values(primary_variables_dof[i])=1; break; }
				else if(p.distance(origin)<r+bandwitdh) {
					interface_flag=true;
					if (p.distance(origin)<current_dis){
						current_dis=p.distance(origin);
						particle_index=ori_index;
					}
				}
			}
			if(interface_flag) {
				Point<dim> origin(origin_list[particle_index][0],origin_list[particle_index][1]);
				values(primary_variables_dof[i])=1-(p.distance(origin)-r)*(1/bandwitdh);
			}
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


template <int dim>
void nodalField<dim>::evaluate_vector_field(const DataPostprocessorInputs::Vector< dim > &input_data, std::vector< Vector< double >> &computed_quantities)const
{	
	const unsigned int n_q_points = computed_quantities.size();	
	double youngsModulus=(*params_json)["Mechanics"]["youngs_modulus_particle"];
	double poissonRatio=(*params_json)["Mechanics"]["poisson_ratio"];
	
	Residual<double,dim> ResidualEq;
	int lithium_index=this->battery_fields->active_fields_index["Lithium"];
	int interface_index=this->battery_fields->active_fields_index["Diffuse_interface"];
	int u_index=this->battery_fields->active_fields_index["Displacement"];
	double eps_0=1.0e-5;
	
	if(input_data.solution_values[0][interface_index]<0.5) youngsModulus=(*params_json)["Mechanics"]["youngs_modulus_electrolyte"];
	
	ResidualEq.setLameParametersByYoungsModulusPoissonRatio(youngsModulus, poissonRatio);	
	double C_a=(*params_json)["Mechanics"]["lithium_a"];
	double C_b=(*params_json)["Mechanics"]["lithium_b"];
	dealii::Table<2,double> Feiga(dim,dim);
	dealii::Table<2,double> Feigba(dim,dim);
	Feiga[0][0]=(*params_json)["Mechanics"]["Feiga_11"];
	Feigba[0][0]=(*params_json)["Mechanics"]["Feigb_11"];
	Feigba[0][0]-=(*params_json)["Mechanics"]["Feiga_11"].get<double>();
	if(dim>=2){
		Feiga[1][1]=(*params_json)["Mechanics"]["Feiga_22"];
		Feigba[1][1]=(*params_json)["Mechanics"]["Feigb_22"];
		Feigba[1][1]-=(*params_json)["Mechanics"]["Feiga_22"].get<double>();
	}
	if(dim==3){
		Feiga[2][2]=(*params_json)["Mechanics"]["Feiga_33"];
		Feigba[2][2]=(*params_json)["Mechanics"]["Feigb_33"];
		Feigba[2][2]-=(*params_json)["Mechanics"]["Feiga_33"].get<double>();
	}
	//std::cout<<"n_q_points"<<n_q_points<<std::endl;
	//std::cout<<"input_data.solution_values[q]"<<input_data.solution_values[0].size()<<std::endl;
	//std::cout<<"u_index="<<u_index<<std::endl;
	for (unsigned int q=0; q<n_q_points; ++q){
		dealii::Table<3,double > Fe(1,dim,dim);
		dealii::Table<3, double > P_stress(1,dim,dim);

		for (unsigned int i = 0; i < dim; ++i){
			for (unsigned int j = 0; j < dim; ++j){
			  Fe[0][i][j] = (i==j) + input_data.solution_gradients[q][i+u_index][j];
			}
		}
		if(input_data.solution_values[q][interface_index]>=0.5){
			double C_q=input_data.solution_values[q][lithium_index];
			dealii::Table<2,double > Feig(dim,dim);
			dealii::Table<2,double> invFeig(dim,dim);
			Feig=table_scaling<2,double,double > (Feigba, (C_q-C_a)/(C_b-C_a) );   
			Feig=table_add<2,double,double > (Feig, Feiga);
			getInverse<double,dim> (Feig,invFeig);
			for (unsigned int i=0; i<dim; ++i){
				for (unsigned int j=0; j<dim; ++j){
					for (unsigned int k=0; k<dim; ++k){
	 					Fe[0][i][j]+=Fe[0][i][k]*invFeig[k][j];
					}
				}
			}
		}
		ResidualEq.evaluateNeoHookeanStress(P_stress, Fe);
		computed_quantities[q][0]=std::sqrt(std::pow(P_stress[0][0][0],2)+std::pow(P_stress[0][1][1],2)-P_stress[0][0][0]*P_stress[0][1][1]+3*P_stress[0][0][1]*P_stress[0][0][1] );
	}	
}

template class nodalField<1>;
template class nodalField<2>;
template class nodalField<3>;
