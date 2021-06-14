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
	std::vector<bool> All_component (totalDOF, false);
	double separator_line=(*params_json)["ElectroChemo"]["separator_line"];
	int orientation=(*params_json)["ElectroChemo"]["orientation"];
	if(battery_fields.active_fields_index["Displacement"]>-1) {
		for(unsigned int i=0;i<dim;i++) All_component[battery_fields.active_fields_index["Displacement"]+i]=true;
	}
	VectorTools:: interpolate_boundary_values (this->dof_handler, 1+orientation, ZeroFunction<dim> (totalDOF),*constraints, All_component);

	std::vector<types::global_dof_index> local_face_dof_indices (this->fe_system[electrolyte_id]->dofs_per_face);
  typename hp::DoFHandler<dim>::active_cell_iterator cell = this->dof_handler.begin_active(), endc=this->dof_handler.end();
  for (;cell!=endc; ++cell){
		if (cell->center()[orientation]>separator_line) continue;
		if(cell->material_id()==electrolyte_id){
			for (unsigned int f=0; f<GeometryInfo<dim>::faces_per_cell; ++f){
				if (cell->at_boundary(f) == false){
					if(cell->neighbor(f)->material_id()==active_particle_id  and cell->neighbor(f)->has_children() == false){
						cell->face(f)->get_dof_indices (local_face_dof_indices,electrolyte_id);
						for (unsigned int i=0; i<local_face_dof_indices.size(); ++i){
					    const unsigned int ck = this->fe_system[electrolyte_id]->face_system_to_component_index(i).first;
							if(ck==battery_fields.active_fields_index["Electrolyte_potential"]) constraints->add_line (local_face_dof_indices[i]);
						}
					}
				}
			}
		}
	}
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
	std::vector<std::vector<double>> origin_list_benchmark={{0,0},{0,18},{0,38},{0,56}};
	double r=(*params_json)["ElectroChemo"]["particle_R"];
	double neg_electrode_line=(*params_json)["ElectroChemo"]["neg_electrode_line"];
	double pos_electrode_line=(*params_json)["ElectroChemo"]["pos_electrode_line"];
	int orientation=(*params_json)["ElectroChemo"]["orientation"];
	this->pcout<<"setMultDomain"<<std::endl;
  typename hp::DoFHandler<dim>::active_cell_iterator cell = this->dof_handler.begin_active(), endc=this->dof_handler.end();
  for (;cell!=endc; ++cell){
		if (true){
			cell->set_material_id(electrolyte_id);
			Point<dim> center=cell->center();
			if (center[orientation]>neg_electrode_line and center[orientation]<pos_electrode_line) continue;
			for(unsigned int ori_index=0;ori_index<origin_list_benchmark.size();ori_index++){
				Point<dim> origin(origin_list_benchmark[ori_index][0],origin_list_benchmark[ori_index][1]);
				if(center.distance(origin)<r-1.0e-15) {
					cell->set_material_id(active_particle_id);
					break;
				}
			}
			
			
			// int inside_vertex=0;
			// if(center[0]<3) cell->set_material_id(active_particle_id);
			// if(center[0]>6) cell->set_material_id(active_particle_id);
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
	
	double C_li_max_neg=(*params_json)["ElectroChemo"]["c_li_max_neg"];
	double C_li_max_pos=(*params_json)["ElectroChemo"]["c_li_max_pos"];
	double C_li_100_neg=(*params_json)["ElectroChemo"]["c_li_100_neg"];
	double C_li_100_pos=(*params_json)["ElectroChemo"]["c_li_100_pos"];
	
	double C_li_plus_0=(*params_json)["ElectroChemo"]["C_li_plus_0"];
	double separator_line=(*params_json)["ElectroChemo"]["separator_line"];
	int orientation=(*params_json)["ElectroChemo"]["orientation"];
	
	double iso_value=(*params_json)["ElectroChemo"]["iso_value"];
  typename hp::DoFHandler<dim>::active_cell_iterator cell = this->dof_handler.begin_active(), endc=this->dof_handler.end();
  for (;cell!=endc; ++cell){
		if (cell->subdomain_id() == this->this_mpi_process){
			double C_li_0=C_li_100_neg*C_li_max_neg;
			Point<dim> center=cell->center();
			if (center[orientation]>separator_line){
				C_li_0=C_li_100_pos*C_li_max_pos;
			}
    	hp::FEValues<dim> hp_fe_values (this->fe_collection, this->q_collection, update_values | update_quadrature_points);
    	hp_fe_values.reinit (cell);
    	const FEValues<dim> &fe_values = hp_fe_values.get_present_fe_values();
    	const unsigned int dofs_per_cell = cell->get_fe().dofs_per_cell;
			std::vector<unsigned int> local_dof_indices (dofs_per_cell);
			cell->get_dof_indices (local_dof_indices);
    	for (unsigned int i=0; i<dofs_per_cell; ++i) {
      	int ck = fe_values.get_fe().system_to_component_index(i).first;
				if (ck==battery_fields.active_fields_index["Lithium"]) this->solution_prev(local_dof_indices[i])=C_li_0;//+0.01*(static_cast <double> (rand())/(static_cast <double>(RAND_MAX))-0.5);
				else if(ck==battery_fields.active_fields_index["Lithium_cation"]) this->solution_prev(local_dof_indices[i])=C_li_plus_0;
				else if(ck==battery_fields.active_fields_index["Electrode_potential"]){
					if(center[orientation]>separator_line) this->solution_prev(local_dof_indices[i])=electricChemoFormula.formula_Usc(C_li_100_pos,1).val();
					if(center[orientation]<separator_line) this->solution_prev(local_dof_indices[i])=electricChemoFormula.formula_Usc(C_li_100_neg,-1).val();
				}
					
				else if(ck==battery_fields.active_fields_index["Electrolyte_potential"]) this->solution_prev(local_dof_indices[i])=0;
			}
		}
	}
	this->solution_prev.compress(VectorOperation::insert);
	this->solution=this->solution_prev;		
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
	double neg_electrode_line=(*params_json)["ElectroChemo"]["neg_electrode_line"];
	double pos_electrode_line=(*params_json)["ElectroChemo"]["pos_electrode_line"];
	int orientation=(*params_json)["ElectroChemo"]["orientation"];
	
	
	Residual<double,dim> ResidualEq;
	int lithium_index=this->battery_fields->active_fields_index["Lithium"];
	int Electrode_potential_index=this->battery_fields->active_fields_index["Electrode_potential"];
	int interface_index=this->battery_fields->active_fields_index["Diffuse_interface"];
	int u_index=this->battery_fields->active_fields_index["Displacement"];
	double eps_0=1.0e-5;
	
	Point<dim> points=input_data.evaluation_points[0];
	if(input_data.solution_values[0][Electrode_potential_index]<1.0e-5){
			youngsModulus=(*params_json)["Mechanics"]["youngs_modulus_electrolyte"];
	} 
	
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
		if(input_data.solution_values[q][Electrode_potential_index]>0){
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
		//std::cout<<"P_stress[0][0][0]"<<P_stress[0][0][0]<<" P_stress[0][1][1]"<<P_stress[0][1][1]<<" P_stress[0][0][1]"<<P_stress[0][0][1]<<std::endl;
		computed_quantities[q][0]=std::sqrt(std::pow(P_stress[0][0][0],2)+std::pow(P_stress[0][1][1],2)-P_stress[0][0][0]*P_stress[0][1][1]+3*P_stress[0][0][1]*P_stress[0][0][1] );
	}	
}

template class nodalField<1>;
template class nodalField<2>;
template class nodalField<3>;