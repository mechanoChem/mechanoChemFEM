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
	
	//int totalDOF=this->totalDOF(this->primary_variables);
    //std::vector<bool> All_component (totalDOF, true);	
	//VectorTools:: interpolate_boundary_values (this->dof_handler, 1, ZeroFunction<dim> (totalDOF),*constraints, All_component);
	//int interface_index=battery_fields.active_fields_index["Diffuse_interface"];


  // for benchmark test
	int totalDOF=this->totalDOF(this->primary_variables);
  std::vector<bool> All_component (totalDOF, false);
	if(battery_fields.active_fields_index["Electrode_potential"]>-1) All_component[battery_fields.active_fields_index["Electrode_potential"]]=true;
	if(battery_fields.active_fields_index["Electrolyte_potential"]>-1) All_component[battery_fields.active_fields_index["Electrolyte_potential"]]=true;
	
	VectorTools:: interpolate_boundary_values (this->dof_handler, 1, ZeroFunction<dim> (totalDOF),*constraints, All_component);
	
	if(battery_fields.active_fields_index["Lithium_cation"]>-1) All_component[battery_fields.active_fields_index["Lithium_cation"]]=true;
	VectorTools:: interpolate_boundary_values (this->dof_handler, 3, ZeroFunction<dim> (totalDOF),*constraints, All_component);
	int interface_index=battery_fields.active_fields_index["Diffuse_interface"];
	constraints->close ();
}
template <int dim>
void battery<dim>::apply_Neumann_boundary_condition()
{}

	
template <int dim>
void battery<dim>::setup_diffuse_interface(){}

template <int dim>
void battery<dim>::setMultDomain()
{
	std::vector<std::vector<double>> origin_list={{1.5,1},{4.5,1}};
	double r=(*params_json)["ElectroChemo"]["particle_R"];
	
  typename hp::DoFHandler<dim>::active_cell_iterator cell = this->dof_handler.begin_active(), endc=this->dof_handler.end();
  for (;cell!=endc; ++cell){
		if (cell->subdomain_id() == this->this_mpi_process){
			cell->set_material_id(electrolyte_id);
      // for rectangle
			Point<dim> center=cell->center();
			int inside_vertex=0;
			if(center[0]<=3.1 and center[0]>2.9 ) cell->set_material_id(interface_id);
			if(center[0]<=2.9) cell->set_material_id(active_particle_id);
      std::cout << "cell_id "<< int(cell->active_cell_index()) << " center " << center[0] << "\t" << center[1] << " mat_id "<< int(cell->material_id()) << std::endl;
      //------------------------------------ for circle -------------
			//int inside_vertex=0;
      //for (unsigned int vertex = 0; vertex < GeometryInfo<dim>::vertices_per_cell; ++vertex){
				//Point<dim> vertrx_point=cell->vertex(vertex);
				//for (unsigned int ori_index=0;ori_index<origin_list.size();ori_index++){
					//Point<dim> origin(origin_list[ori_index][0],origin_list[ori_index][1]);
					//if(vertrx_point.distance(origin)<r-1.0e-15) {
						//inside_vertex++;
						//break;
					//}
				//}
			//}
			//if(inside_vertex==GeometryInfo<dim>::vertices_per_cell) cell->set_material_id(active_particle_id);
			//else if (inside_vertex>0) cell->set_material_id( interface_id); 
      //------------------------------- circle---------------------

		}
	}
	this->set_active_fe_indices (this->FE_support, this->dof_handler);
}

template <int dim>
void battery<dim>::apply_initial_condition()
{
	std::vector<std::vector<double>> origin_list={{1.5,1},{4.5,1}};
	double r=(*params_json)["ElectroChemo"]["particle_R"];
	double bandwitdh=(*params_json)["ElectroChemo"]["interface_bandwitdh"];
	
	double C_li_0_neg=(*params_json)["ElectroChemo"]["C_li_0_neg"];
	double C_li_0_pos=(*params_json)["ElectroChemo"]["C_li_0_pos"];
	double C_li_plus_0=(*params_json)["ElectroChemo"]["C_li_plus_0"];
	double phi_s_0_neg=(*params_json)["ElectroChemo"]["phi_s_0_neg"];
	double phi_s_0_pos=(*params_json)["ElectroChemo"]["phi_s_0_pos"];
	double phi_e_0=(*params_json)["ElectroChemo"]["phi_e_0"];
	double separator_line=(*params_json)["ElectroChemo"]["separator_line"];
	double iso_value=(*params_json)["ElectroChemo"]["iso_value"];
  typename hp::DoFHandler<dim>::active_cell_iterator cell = this->dof_handler.begin_active(), endc=this->dof_handler.end();
  for (;cell!=endc; ++cell){
		if (cell->subdomain_id() == this->this_mpi_process){
			double C_li_0=C_li_0_neg;
			double phi_s_0=phi_s_0_neg;
			Point<dim> center=cell->center();
			if (center[1]>separator_line){
				C_li_0=C_li_0_pos;
				phi_s_0=phi_s_0_pos;
			}
    	hp::FEValues<dim> hp_fe_values (this->fe_collection, this->q_collection, update_values | update_quadrature_points);
    	hp_fe_values.reinit (cell);
    	const FEValues<dim> &fe_values = hp_fe_values.get_present_fe_values();
    	const unsigned int dofs_per_cell = cell->get_fe().dofs_per_cell;
			std::vector<unsigned int> local_dof_indices (dofs_per_cell);
			cell->get_dof_indices (local_dof_indices);
    	for (unsigned int i=0; i<dofs_per_cell; ++i) {
      	int ck = fe_values.get_fe().system_to_component_index(i).first;
				if (ck==battery_fields.active_fields_index["Lithium"]){ 
          //this->solution_prev(local_dof_indices[i])=C_li_0+0.04*(static_cast <double> (rand())/(static_cast <double>(RAND_MAX))-0.5);
          this->solution_prev(local_dof_indices[i])=C_li_0;
				}
				else if(ck==battery_fields.active_fields_index["Lithium_cation"]) this->solution_prev(local_dof_indices[i])=C_li_plus_0;
				else if(ck==battery_fields.active_fields_index["Electrode_potential"]) this->solution_prev(local_dof_indices[i])=phi_s_0;
				else if(ck==battery_fields.active_fields_index["Electrolyte_potential"]) this->solution_prev(local_dof_indices[i])=phi_e_0;
				
				if (cell->material_id()==interface_id){
          //std::cout << "---------- in interface ---------------" << std::endl;
					int vertex_id=i / (dofs_per_cell/GeometryInfo<dim>::vertices_per_cell);
					Point<dim> vertrx_point=cell->vertex(vertex_id);
					bool inside_flag=false;
					double distance=1.0e16;
					for (unsigned int ori_index=0;ori_index<origin_list.size();ori_index++){
						Point<dim> origin(origin_list[ori_index][0],origin_list[ori_index][1]);
						if (vertrx_point.distance(origin)<distance) distance=vertrx_point.distance(origin);
						if(vertrx_point.distance(origin)<=r) {inside_flag=true; break;}
					}
					if (ck==battery_fields.active_fields_index["Diffuse_interface"]) this->solution_prev(local_dof_indices[i])=(r-distance)+iso_value;

					if (ck==battery_fields.active_fields_index["Diffuse_interface"]) this->solution_prev(local_dof_indices[i])=3.0-vertrx_point[0]+iso_value;
          if (vertrx_point[0] <= 2.9 ) inside_flag = true;
					if (inside_flag){
						if (ck==battery_fields.active_fields_index["Lithium_cation"] or ck==battery_fields.active_fields_index["Electrolyte_potential"]){
							this->solution_prev(local_dof_indices[i])=0;
						}
					}
					else{
						if (ck==battery_fields.active_fields_index["Lithium"] or ck==battery_fields.active_fields_index["Electrode_potential"]){
							this->solution_prev(local_dof_indices[i])=0;
						}
					}
				}
			}
		}
	}
	this->solution_prev.compress(VectorOperation::insert);
	this->solution=this->solution_prev;		

  identify_diffuse_interface();

	constraints->clear ();
	
	//DoFTools::make_hanging_node_constraints (this->dof_handler, *constraints);
	
	//int totalDOF=this->totalDOF(this->primary_variables);
    //std::vector<bool> All_component (totalDOF, true);	
	//ectorTools:: interpolate_boundary_values (this->dof_handler, 1, ZeroFunction<dim> (totalDOF),*constraints, All_component);
	//int interface_index=battery_fields.active_fields_index["Diffuse_interface"];

  // for benchmark test
	int totalDOF=this->totalDOF(this->primary_variables);
  std::vector<bool> All_component (totalDOF, false);
	if(battery_fields.active_fields_index["Electrode_potential"]>-1) All_component[battery_fields.active_fields_index["Electrode_potential"]]=true;
	if(battery_fields.active_fields_index["Electrolyte_potential"]>-1) All_component[battery_fields.active_fields_index["Electrolyte_potential"]]=true;
	
	VectorTools:: interpolate_boundary_values (this->dof_handler, 1, ZeroFunction<dim> (totalDOF),*constraints, All_component);
	
	if(battery_fields.active_fields_index["Lithium_cation"]>-1) All_component[battery_fields.active_fields_index["Lithium_cation"]]=true;
	VectorTools:: interpolate_boundary_values (this->dof_handler, 3, ZeroFunction<dim> (totalDOF),*constraints, All_component);
	int interface_index=battery_fields.active_fields_index["Diffuse_interface"];


  {
    // remove the minus or plus side of node, not to solve
    typename hp::DoFHandler<dim>::active_cell_iterator cell = this->dof_handler.begin_active(), endc = this->dof_handler.end();
    for (; cell != endc; ++cell) {
      if (cell->subdomain_id() == this->this_mpi_process) {
        int cell_id = cell->active_cell_index();
        if (cell_SDdata[cell_id].is_interface_element) {

          const unsigned int dofs_per_cell = cell->get_fe().dofs_per_cell;
          std::vector<unsigned int> local_dof_indices(dofs_per_cell);
          cell->get_dof_indices(local_dof_indices);
          const unsigned int dofs_per_node = dofs_per_cell / 4;


          for (auto i0 : cell_SDdata[cell_id].lnode_minus)
          {
	          if(battery_fields.active_fields_index["Lithium"]>-1) 
            {
              auto globalDOF = local_dof_indices[i0*dofs_per_node + battery_fields.active_fields_index["Lithium"]];
              //std::cout << "Lithium " << i0*dofs_per_node + battery_fields.active_fields_index["Lithium"] << " i0 " << i0 << " dofs_per_node " << dofs_per_node << " globalDOF " << globalDOF << std::endl;
              constraints->add_line(globalDOF);
              constraints->set_inhomogeneity(globalDOF, 0.0);
              // initialize the xi_old to be the same value as the neighbor node
              cell_SDdata[cell_id].xi_old(0) = this->solution_prev(local_dof_indices[cell_SDdata[cell_id].lnode_plus[0]*dofs_per_node + battery_fields.active_fields_index["Lithium"]]);
              cell_SDdata[cell_id].xi_conv(0) = this->solution_prev(local_dof_indices[cell_SDdata[cell_id].lnode_plus[0]*dofs_per_node + battery_fields.active_fields_index["Lithium"]]);
              //std::cout << "xi_old Lithium " << cell_SDdata[cell_id].xi_old(0) << " plus[0] " << cell_SDdata[cell_id].lnode_plus[0] << " i0 " << cell_SDdata[cell_id].lnode_plus[0]*dofs_per_node + battery_fields.active_fields_index["Lithium"] << std::endl;
              //std::cout << this->solution_prev(local_dof_indices[cell_SDdata[cell_id].lnode_plus[0]*dofs_per_node + battery_fields.active_fields_index["Lithium"]]) << std::endl;
              //std::cout << this->solution_prev(local_dof_indices[cell_SDdata[cell_id].lnode_plus[1]*dofs_per_node + battery_fields.active_fields_index["Lithium"]]) << std::endl;
              //std::cout << this->solution_prev(local_dof_indices[cell_SDdata[cell_id].lnode_minus[0]*dofs_per_node + battery_fields.active_fields_index["Lithium"]]) << std::endl;
              //std::cout << this->solution_prev(local_dof_indices[cell_SDdata[cell_id].lnode_minus[1]*dofs_per_node + battery_fields.active_fields_index["Lithium"]]) << std::endl;
            }
            if(battery_fields.active_fields_index["Electrode_potential"]>-1)
            {
              auto globalDOF = local_dof_indices[i0*dofs_per_node + battery_fields.active_fields_index["Electrode_potential"]];
              //std::cout << "Electrode_potential " << i0*dofs_per_node + battery_fields.active_fields_index["Electrode_potential"] << " i0 " << i0 << " dofs_per_node " << dofs_per_node << " globalDOF " << globalDOF << std::endl;
              constraints->add_line(globalDOF);
              constraints->set_inhomogeneity(globalDOF, 0.0);
              cell_SDdata[cell_id].xi_old_phi_s(0) = this->solution_prev(local_dof_indices[cell_SDdata[cell_id].lnode_plus[0]*dofs_per_node + battery_fields.active_fields_index["Electrode_potential"]]);
              cell_SDdata[cell_id].xi_conv_phi_s(0) = this->solution_prev(local_dof_indices[cell_SDdata[cell_id].lnode_plus[0]*dofs_per_node + battery_fields.active_fields_index["Electrode_potential"]]);
            }

          }

          for (auto i0 : cell_SDdata[cell_id].lnode_plus)
          {
            if(battery_fields.active_fields_index["Lithium_cation"]>-1) 
            {
              auto globalDOF = local_dof_indices[i0*dofs_per_node + battery_fields.active_fields_index["Lithium_cation"]];
              //std::cout << "Lithium_cation " << i0*dofs_per_node + battery_fields.active_fields_index["Lithium_cation"] << " i0 " << i0 << " dofs_per_node " << dofs_per_node << " globalDOF " << globalDOF << std::endl;
              constraints->add_line(globalDOF);
              constraints->set_inhomogeneity(globalDOF, 0.0);
              cell_SDdata[cell_id].xi_old_c_e(0) = this->solution_prev(local_dof_indices[cell_SDdata[cell_id].lnode_minus[0]*dofs_per_node + battery_fields.active_fields_index["Lithium_cation"]]);
              cell_SDdata[cell_id].xi_conv_c_e(0) = this->solution_prev(local_dof_indices[cell_SDdata[cell_id].lnode_minus[0]*dofs_per_node + battery_fields.active_fields_index["Lithium_cation"]]);
              for (int q=0; q<4; q++) cell_SDdata[cell_id].C_Li_plus_old[q] = cell_SDdata[cell_id].xi_old_c_e(0);
            }
            if(battery_fields.active_fields_index["Electrolyte_potential"]>-1) 
            {
              auto globalDOF = local_dof_indices[i0*dofs_per_node + battery_fields.active_fields_index["Electrolyte_potential"]];
              //std::cout << "Electrolyte_potential " << i0*dofs_per_node + battery_fields.active_fields_index["Electrolyte_potential"] << " i0 " << i0 << " dofs_per_node " << dofs_per_node << " globalDOF " << globalDOF << std::endl;
              constraints->add_line(globalDOF);
              constraints->set_inhomogeneity(globalDOF, 0.0);
              cell_SDdata[cell_id].xi_old_phi_e(0) = this->solution_prev(local_dof_indices[cell_SDdata[cell_id].lnode_minus[0]*dofs_per_node + battery_fields.active_fields_index["Electrolyte_potential"]]);
              cell_SDdata[cell_id].xi_conv_phi_e(0) = this->solution_prev(local_dof_indices[cell_SDdata[cell_id].lnode_minus[0]*dofs_per_node + battery_fields.active_fields_index["Electrolyte_potential"]]);
            }
          }

        }
      }
    }
  }



	constraints->close ();


}

template class battery<1>;
template class battery<2>;
template class battery<3>;




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
