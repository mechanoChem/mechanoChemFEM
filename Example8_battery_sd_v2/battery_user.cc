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
	//int totalDOF=this->totalDOF(this->primary_variables);
  //std::vector<bool> All_component (totalDOF, false);
	//if(battery_fields.active_fields_index["Electrode_potential"]>-1) All_component[battery_fields.active_fields_index["Electrode_potential"]]=true;
	//if(battery_fields.active_fields_index["Electrolyte_potential"]>-1) All_component[battery_fields.active_fields_index["Electrolyte_potential"]]=true;
	
	//VectorTools:: interpolate_boundary_values (this->dof_handler, 1, ZeroFunction<dim> (totalDOF),*constraints, All_component);
	
	//if(battery_fields.active_fields_index["Lithium_cation"]>-1) All_component[battery_fields.active_fields_index["Lithium_cation"]]=true;
	//VectorTools:: interpolate_boundary_values (this->dof_handler, 3, ZeroFunction<dim> (totalDOF),*constraints, All_component);
	//int interface_index=battery_fields.active_fields_index["Diffuse_interface"];
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
  {
    typename hp::DoFHandler<dim>::active_cell_iterator cell = this->dof_handler.begin_active(), endc=this->dof_handler.end();
    for (;cell!=endc; ++cell){
	  	if (cell->subdomain_id() == this->this_mpi_process){
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
      }
    }
  }

  {
    // create a layer of interface element
    typename hp::DoFHandler<dim>::active_cell_iterator cell = this->dof_handler.begin_active(), endc=this->dof_handler.end();
    for (;cell!=endc; ++cell){
	  	if (cell->subdomain_id() == this->this_mpi_process){
		    if(cell->material_id()==active_particle_id){
		    	for (unsigned int f=0; f<GeometryInfo<dim>::faces_per_cell; ++f){
		    		if (cell->at_boundary(f) == false){
		    			if(cell->neighbor(f)->material_id()==electrolyte_id  and cell->neighbor(f)->has_children() == false){
	  	  			  cell->set_material_id(interface_id);
		    			}
		    		}
		    	}
		    }
      }
    }
  }

  {
    // fix a few missing interface element
    typename hp::DoFHandler<dim>::active_cell_iterator cell = this->dof_handler.begin_active(), endc=this->dof_handler.end();
    for (;cell!=endc; ++cell){
	  	if (cell->subdomain_id() == this->this_mpi_process){
		    if(cell->material_id()==active_particle_id){
	  		  Point<dim> center=cell->center();
          int _neighbor_count = 0;
		    	for (unsigned int f=0; f<GeometryInfo<dim>::faces_per_cell; ++f){
		    		if (cell->at_boundary(f) == false){
		    			if(cell->neighbor(f)->material_id()==interface_id  and cell->neighbor(f)->has_children() == false){
                _neighbor_count += 1;
		    			}
		    		}
		    	}
          if (_neighbor_count == 2)
          {
            // bottom connection corner
            if (center[orientation] >4 and  center[orientation] < 12) cell->set_material_id(interface_id);
            // top connection corner
            if (center[orientation] >42 and  center[orientation] < 50) cell->set_material_id(interface_id);
          }
		    }
      }
    }
  }

			//int inside_vertex=0;
			//if(center[0]<3) cell->set_material_id(active_particle_id);
			//if(center[0]>6) cell->set_material_id(active_particle_id);
      //if(center[0]<=3.1 and center[0]>2.9 ) cell->set_material_id(interface_id);
      //if(center[0]<=6.1 and center[0]>5.9 ) cell->set_material_id(interface_id);

			////if(center[0]<=2.9) cell->set_material_id(active_particle_id);
			////if(center[0]<=3.1 and center[0]>2.9 ) cell->set_material_id(interface_id);
			////if(center[0]<=2.9) cell->set_material_id(active_particle_id);

      //std::cout << "cell_id "<< int(cell->active_cell_index()) << " center " << center[0] << "\t" << center[1] << " mat_id "<< int(cell->material_id()) << std::endl;
      ////------------------------------------ for circle -------------
			////int inside_vertex=0;
      ////for (unsigned int vertex = 0; vertex < GeometryInfo<dim>::vertices_per_cell; ++vertex){
				////Point<dim> vertex_point=cell->vertex(vertex);
				////for (unsigned int ori_index=0;ori_index<origin_list.size();ori_index++){
					////Point<dim> origin(origin_list[ori_index][0],origin_list[ori_index][1]);
					////if(vertex_point.distance(origin)<r-1.0e-15) {
						////inside_vertex++;
						////break;
					////}
				////}
			////}
			////if(inside_vertex==GeometryInfo<dim>::vertices_per_cell) cell->set_material_id(active_particle_id);
			////else if (inside_vertex>0) cell->set_material_id( interface_id); 
      ////------------------------------- circle---------------------

	this->set_active_fe_indices (this->FE_support, this->dof_handler);
}

template <int dim>
void battery<dim>::apply_initial_condition()
{
  std::vector<std::vector<double>> origin_list={{1.5,1},{4.5,1}};
	std::vector<std::vector<double>> origin_list_benchmark={{0,0},{0,18},{0,38},{0,56}};
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
        if (cell->material_id()!=interface_id){
				  if (ck==battery_fields.active_fields_index["Lithium"]) this->solution_prev(local_dof_indices[i])=C_li_0;//+0.01*(static_cast <double> (rand())/(static_cast <double>(RAND_MAX))-0.5);
				  else if(ck==battery_fields.active_fields_index["Lithium_cation"]) this->solution_prev(local_dof_indices[i])=C_li_plus_0;
				  else if(ck==battery_fields.active_fields_index["Electrode_potential"]){
				  	if(center[orientation]>separator_line) this->solution_prev(local_dof_indices[i])=electricChemoFormula.formula_Usc(C_li_100_pos,1).val();
				  	if(center[orientation]<separator_line) this->solution_prev(local_dof_indices[i])=electricChemoFormula.formula_Usc(C_li_100_neg,-1).val();
				  }
				  	
				  else if(ck==battery_fields.active_fields_index["Electrolyte_potential"]) this->solution_prev(local_dof_indices[i])=0;
        }
        else
        {
          bool is_horizontal_element = false;
          {
            for (int i=0; i <GeometryInfo<dim>::vertices_per_cell; i++)
            {
              Point<dim> vertex_point=cell->vertex(i);
              if (vertex_point[orientation] == 18 or vertex_point[orientation] == 38 )
              {
                is_horizontal_element = true;
                break;
              }
            }
          }

          //std::cout << "---------- in interface ---------------" << std::endl;
          int vertex_id=i / (dofs_per_cell/GeometryInfo<dim>::vertices_per_cell);
          Point<dim> vertex_point=cell->vertex(vertex_id);
          bool inside_flag=false;
          double distance=1.0e16;
          for (unsigned int ori_index=0;ori_index<origin_list_benchmark.size();ori_index++){
            Point<dim> origin(origin_list_benchmark[ori_index][0],origin_list_benchmark[ori_index][1]);
            if(vertex_point.distance(origin)<distance) distance=vertex_point.distance(origin);
            //std::cout 
              //<< " i " << i 
              //<< " center " << center 
              //<< " ori_ind "<< ori_index 
              //<< " origin " << origin 
              //<< " distance " << distance 
              //<< std::endl;
            if(abs(vertex_point.distance(origin) -r) < 1.0e-5) {inside_flag=false; break;}
            if(vertex_point.distance(origin) < r - 1.0e-5) {inside_flag=true; break;}
          }

          if(abs(distance -r) < 1.0e-5) {inside_flag=false;}
          //std::cout << " " << abs(distance -r) << std::endl;


          //if (inside_flag and vertex_point[orientation] > 38 and distance > 9.0 and i%5==0)
          //{
            //std::cout 
              //<< " i " << i 
              //<< " center " << center 
              //<< " distance " << distance 
              //<< " vertex " << vertex_point
              //<< std::endl;
          //}
          // the isosurface should be very close to the edge. thus distance is chose to be the size of part of the element length dh
          double dh = 0.3;
          //std::cout << "----***---- distance " << distance << std::endl;
          if (ck==battery_fields.active_fields_index["Diffuse_interface"]) this->solution_prev(local_dof_indices[i])=(r-dh-distance)+iso_value;
          if (is_horizontal_element)
          {
            // not on the circle
            if (abs(distance - r) > 1.0e-5)
            {
              if (vertex_point[orientation] <= 18)
              {
                if (ck==battery_fields.active_fields_index["Diffuse_interface"]) this->solution_prev(local_dof_indices[i])=(18-dh-vertex_point[orientation])+iso_value;
              }
              if (vertex_point[orientation] >= 38)
              {
                if (ck==battery_fields.active_fields_index["Diffuse_interface"]) this->solution_prev(local_dof_indices[i])=(vertex_point[orientation]-38 -dh)+iso_value;
              }

              if (vertex_point[orientation] == 18 or vertex_point[orientation] == 38 )
              {
                inside_flag = false;
              }
              else
              {
                inside_flag = true;
              }
                
            }
          }
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
        } //interface
      }
    }
  }
  this->solution_prev.compress(VectorOperation::insert);
  this->solution=this->solution_prev;		

  identify_diffuse_interface();

  constraints->clear ();
  
  DoFTools::make_hanging_node_constraints (this->dof_handler, *constraints);


  {
    bool is_one_node_Electrolyte_potential_fixed = false;
    // remove the minus or plus side of node, not to solve
    typename hp::DoFHandler<dim>::active_cell_iterator cell = this->dof_handler.begin_active(), endc = this->dof_handler.end();
    for (; cell != endc; ++cell) {
      if (cell->subdomain_id() == this->this_mpi_process) {
        int cell_id = cell->active_cell_index();
        Point<dim> center=cell->center();

        const unsigned int dofs_per_cell = cell->get_fe().dofs_per_cell;
        std::vector<unsigned int> local_dof_indices(dofs_per_cell);
        cell->get_dof_indices(local_dof_indices);
        const unsigned int dofs_per_node = dofs_per_cell / 4;

        if (center[orientation] >= 27 and center[orientation] < 29  ) // for the new 3 layer geometry, special case
        {
          if (not is_one_node_Electrolyte_potential_fixed)
          {
            // add Dirichlet constraint to make some of the potential to be zero.
            int i0 = 0; // first node of one cell in the region 4 < x < 5
            auto globalDOF = local_dof_indices[i0*dofs_per_node + battery_fields.active_fields_index["Electrolyte_potential"]];
            //std::cout << "Electrolyte_potential " << i0*dofs_per_node + battery_fields.active_fields_index["Electrolyte_potential"] << " i0 " << i0 << " dofs_per_node " << dofs_per_node << " globalDOF " << globalDOF << std::endl;
            constraints->add_line(globalDOF);
            constraints->set_inhomogeneity(globalDOF, 0.0);
            is_one_node_Electrolyte_potential_fixed = true;
          }
        }

        if (cell_SDdata[cell_id].is_interface_element) {

          cell_SDdata[cell_id].reaction_rate_potential = 0.0;
          cell_SDdata[cell_id].reaction_rate_li = 0.0;

          for (auto i0 : cell_SDdata[cell_id].lnode_minus)
          {
            if(battery_fields.active_fields_index["Lithium"]>-1) 
            {
              auto globalDOF = local_dof_indices[i0*dofs_per_node + battery_fields.active_fields_index["Lithium"]];
              //std::cout << "Lithium " << i0*dofs_per_node + battery_fields.active_fields_index["Lithium"] << " i0 " << i0 << " dofs_per_node " << dofs_per_node << " globalDOF " << globalDOF << std::endl;
              constraints->add_line(globalDOF);
              constraints->set_inhomogeneity(globalDOF, 0.0);
              // initialize the xi_old to be the same value as the neighbor node
              //cell_SDdata[cell_id].xi_old(0) = this->solution_prev(local_dof_indices[cell_SDdata[cell_id].lnode_plus[0]*dofs_per_node + battery_fields.active_fields_index["Lithium"]]);
              //cell_SDdata[cell_id].xi_conv(0) = this->solution_prev(local_dof_indices[cell_SDdata[cell_id].lnode_plus[0]*dofs_per_node + battery_fields.active_fields_index["Lithium"]]);
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
              //cell_SDdata[cell_id].xi_old_phi_s(0) = this->solution_prev(local_dof_indices[cell_SDdata[cell_id].lnode_plus[0]*dofs_per_node + battery_fields.active_fields_index["Electrode_potential"]]);
              //cell_SDdata[cell_id].xi_conv_phi_s(0) = this->solution_prev(local_dof_indices[cell_SDdata[cell_id].lnode_plus[0]*dofs_per_node + battery_fields.active_fields_index["Electrode_potential"]]);
            }

          }

          int _count_plus = 0;
          cell_SDdata[cell_id].xi_old(0) = 0.0;
          cell_SDdata[cell_id].xi_conv(0) = 0.0;
          cell_SDdata[cell_id].xi_old_phi_s(0) = 0.0;
          cell_SDdata[cell_id].xi_conv_phi_s(0) = 0.0;
          for (auto i0 : cell_SDdata[cell_id].lnode_plus)
          {
              // initialize the xi_old to be the same value as the neighbor node
              cell_SDdata[cell_id].xi_old(0) += this->solution_prev(local_dof_indices[i0*dofs_per_node + battery_fields.active_fields_index["Lithium"]]);
              cell_SDdata[cell_id].xi_old_phi_s(0) += this->solution_prev(local_dof_indices[i0*dofs_per_node + battery_fields.active_fields_index["Electrode_potential"]]);
            _count_plus += 1;
          }
          cell_SDdata[cell_id].xi_old(0) = cell_SDdata[cell_id].xi_old(0)/ _count_plus;
          cell_SDdata[cell_id].xi_old_phi_s(0) = cell_SDdata[cell_id].xi_old_phi_s(0)/ _count_plus;
          cell_SDdata[cell_id].xi_conv(0) = cell_SDdata[cell_id].xi_old(0);
          cell_SDdata[cell_id].xi_conv_phi_s(0) = cell_SDdata[cell_id].xi_old_phi_s(0);

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

          int _count_minus = 0;
          cell_SDdata[cell_id].xi_old_c_e(0) = 0.0;
          cell_SDdata[cell_id].xi_conv_c_e(0) = 0.0;
          cell_SDdata[cell_id].xi_old_phi_e(0) = 0.0;
          cell_SDdata[cell_id].xi_conv_phi_e(0) = 0.0;
          for (auto i0 : cell_SDdata[cell_id].lnode_minus)
          {
              // initialize the xi_old to be the same value as the neighbor node
              cell_SDdata[cell_id].xi_old_c_e(0) += this->solution_prev(local_dof_indices[i0*dofs_per_node + battery_fields.active_fields_index["Lithium_cation"]]);
              cell_SDdata[cell_id].xi_old_phi_e(0) += this->solution_prev(local_dof_indices[i0*dofs_per_node + battery_fields.active_fields_index["Electrolyte_potential"]]);
            _count_minus += 1;
          }
          cell_SDdata[cell_id].xi_old_c_e(0) = cell_SDdata[cell_id].xi_old_c_e(0)/ _count_minus;
          cell_SDdata[cell_id].xi_old_phi_e(0) = cell_SDdata[cell_id].xi_old_phi_e(0)/ _count_minus;
          cell_SDdata[cell_id].xi_conv_c_e(0) = cell_SDdata[cell_id].xi_old_c_e(0);
          cell_SDdata[cell_id].xi_conv_phi_e(0) = cell_SDdata[cell_id].xi_old_phi_e(0);

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
