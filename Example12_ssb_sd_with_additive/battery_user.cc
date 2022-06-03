/*
zhenlin wang 2019
*coupled diffusion reaction
*/
#include "battery.h"
#include "nodalField.h"
#include <math.h>
#include <iostream>
#include <fstream>


template <int dim>
void battery<dim>::identify_diffuse_interface()
{
  int primary_dof = -1;
  int opposite_flux_dof_li = -1;
  int opposite_flux_dof_potential = -1;
    if(battery_fields.active_fields_index["Diffuse_interface"]>-1) primary_dof=battery_fields.active_fields_index["Diffuse_interface"];
    if(battery_fields.active_fields_index["Lithium_cation"]>-1) opposite_flux_dof_li=battery_fields.active_fields_index["Lithium_cation"];
    if(battery_fields.active_fields_index["Electrolyte_potential"]>-1) opposite_flux_dof_potential=battery_fields.active_fields_index["Electrolyte_potential"];

  //std::cout << "---------- primary dof for diffusive interface ------ " << primary_dof  << " opposite dof " << opposite_flux_dof_li << " "<< opposite_flux_dof_potential << std::endl;
  double iso_value=(*params_json)["ElectroChemo"]["iso_value"];

  hp::FEValues<dim> hp_fe_values (this->fe_collection, this->q_collection, update_values | update_quadrature_points  | update_JxW_values | update_gradients);

  Vector<double> localized_U(this->solution_prev);
  int total_cell_num = this->triangulation.n_active_cells();
  cell_SDdata.resize(total_cell_num);
  //std::cout << " total_cell_num " << total_cell_num << std::endl;
  typename hp::DoFHandler<dim>::active_cell_iterator cell = this->dof_handler.begin_active(), endc=this->dof_handler.end();
  for (;cell!=endc; ++cell){
        if (cell->subdomain_id() == this->this_mpi_process){
      if (cell->material_id()==interface_id or cell->material_id()==li_metal_interface_id or cell->material_id()==additive_interface_id)
      {
        unsigned int this_interface_id = cell->material_id();
                hp_fe_values.reinit (cell);
            const FEValues<dim> &fe_values = hp_fe_values.get_present_fe_values();

          int cell_id = cell->active_cell_index();
        //std::cout <<  " cell_id " << cell_id << " cell_SDdata.size() " << cell_SDdata.size() << " interface_id " << interface_id << " mpi " <<this->this_mpi_process << std::endl;
          cell_SDdata[cell_id].cell_id = cell_id;
          cell_SDdata[cell_id].opposite_flux_dof_li = opposite_flux_dof_li;
          cell_SDdata[cell_id].opposite_flux_dof_potential = opposite_flux_dof_potential;

          const unsigned int dofs_per_cell = cell->get_fe().dofs_per_cell;

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
                
        cell_SDdata[cell_id].is_interface_element = true;

        cell_SDdata[cell_id].rlocal.reinit(1);
        cell_SDdata[cell_id].rlocal(0) = 0.0;
        cell_SDdata[cell_id].xi_old.reinit(1);
        cell_SDdata[cell_id].xi_old(0) = 0.0;
        cell_SDdata[cell_id].xi_conv.reinit(1);
        cell_SDdata[cell_id].xi_conv(0) = 0.0;
        cell_SDdata[cell_id].Kcc.reinit(4,4);
        cell_SDdata[cell_id].Kcxi.reinit(4,1);
        cell_SDdata[cell_id].Kxic.reinit(1,4);
        cell_SDdata[cell_id].Kxixi_inv.reinit(1,1);

        cell_SDdata[cell_id].rlocal_c_e.reinit(1);
        cell_SDdata[cell_id].rlocal_c_e(0) = 0.0;
        cell_SDdata[cell_id].xi_old_c_e.reinit(1);
        cell_SDdata[cell_id].xi_old_c_e(0) = 0.0;
        cell_SDdata[cell_id].xi_conv_c_e.reinit(1);
        cell_SDdata[cell_id].xi_conv_c_e(0) = 0.0;
        cell_SDdata[cell_id].Kcc_c_e.reinit(4,4);
        cell_SDdata[cell_id].Kcxi_c_e.reinit(4,1);
        cell_SDdata[cell_id].Kxic_c_e.reinit(1,4);
        cell_SDdata[cell_id].Kxixi_inv_c_e.reinit(1,1);

        cell_SDdata[cell_id].C_Li_plus_old.reinit(4); // size of gps
        cell_SDdata[cell_id].C_Li_plus_new.reinit(4);

        cell_SDdata[cell_id].rlocal_phi_s.reinit(1);
        cell_SDdata[cell_id].rlocal_phi_s(0) = 0.0;
        cell_SDdata[cell_id].xi_old_phi_s.reinit(1);
        cell_SDdata[cell_id].xi_old_phi_s(0) = 0.0;
        cell_SDdata[cell_id].xi_conv_phi_s.reinit(1);
        cell_SDdata[cell_id].xi_conv_phi_s(0) = 0.0;
        cell_SDdata[cell_id].Kcc_phi_s.reinit(4,4);
        cell_SDdata[cell_id].Kcxi_phi_s.reinit(4,1);
        cell_SDdata[cell_id].Kxic_phi_s.reinit(1,4);
        cell_SDdata[cell_id].Kxixi_inv_phi_s.reinit(1,1);

        cell_SDdata[cell_id].rlocal_phi_e.reinit(1);
        cell_SDdata[cell_id].rlocal_phi_e(0) = 0.0;
        cell_SDdata[cell_id].xi_old_phi_e.reinit(1);
        cell_SDdata[cell_id].xi_old_phi_e(0) = 0.0;
        cell_SDdata[cell_id].xi_conv_phi_e.reinit(1);
        cell_SDdata[cell_id].xi_conv_phi_e(0) = 0.0;
        cell_SDdata[cell_id].Kcc_phi_e.reinit(4,4);
        cell_SDdata[cell_id].Kcxi_phi_e.reinit(4,1);
        cell_SDdata[cell_id].Kxic_phi_e.reinit(1,4);
        cell_SDdata[cell_id].Kxixi_inv_phi_e.reinit(1,1);

        cell_SDdata[cell_id].rlocal_u_sd.reinit(3);
        cell_SDdata[cell_id].rlocal_u_sd = 0.0;
        cell_SDdata[cell_id].xi_old_u_sd.reinit(3);
        cell_SDdata[cell_id].xi_old_u_sd = 0.0;
        cell_SDdata[cell_id].xi_conv_u_sd.reinit(3);
        cell_SDdata[cell_id].xi_conv_u_sd = 0.0;
        cell_SDdata[cell_id].Kuu_sd.reinit(4*dim,4*dim);
        cell_SDdata[cell_id].Kuxi_sd.reinit(4*dim,3);
        cell_SDdata[cell_id].Kxiu_sd.reinit(3,4*dim);
        cell_SDdata[cell_id].Kxixi_inv_u_sd.reinit(3,3);

        cell_SDdata[cell_id].ULocal_k.reinit(40);

        unsigned int n_q_points = fe_values.n_quadrature_points;
        for (unsigned int q = 0; q < n_q_points; ++q) {
          cell_SDdata[cell_id].area_elem += fe_values.JxW(q);
        }

        unsigned int count_larger_c = 0, count_smaller_c = 0, count_equal_c = 0;

        // get the side of the local and global node number
        for (unsigned int i=0; i<local_diffuse_interface.size(); ++i) {
          if (local_diffuse_interface[i] >= iso_value){
            cell_SDdata[cell_id].lnode_plus.push_back(i);
            cell_SDdata[cell_id].one_plus_node = cell->vertex(i);
            if (std::abs(local_diffuse_interface[i] - iso_value) < 1e-12)
            {
              count_equal_c += 1;
            }
            else
            {
              count_larger_c += 1;
            }
          };
          if (local_diffuse_interface[i] < iso_value){
            cell_SDdata[cell_id].lnode_minus.push_back(i);
            //std::cout << " interface " << cell_id << " lnode_minus " << i << std::endl;
            count_smaller_c += 1;
          };
          //std::cout << " --- ** -- " << i << std::endl;
        }

        std::vector<types::global_dof_index> local_face_dof_indices(this->fe_system[this_interface_id]->dofs_per_face);
        int count = 0;
        for (unsigned int f = 0; f < GeometryInfo<dim>::faces_per_cell; ++f) {
          cell->face(f)->get_dof_indices(local_face_dof_indices, this_interface_id);
          double c_1 = 0.0;
          double c_2 = 0.0;
          std::vector<double> local_local_diffuse_interface_face;
          for (unsigned int i = 0; i < local_face_dof_indices.size(); ++i) {
            const unsigned int ck = this->fe_system[this_interface_id]->face_system_to_component_index(i).first - primary_dof;
            if (ck == 0) local_local_diffuse_interface_face.push_back(localized_U(local_face_dof_indices[i]));
          }

          c_1 = local_local_diffuse_interface_face[0];
          c_2 = local_local_diffuse_interface_face[1];
          //std::cout << " count " << count << " c_1 " << c_1 << " c_2 " << c_2
            //<< " larger " << count_larger_c
            //<< " smaller " << count_smaller_c
            //<< " equal " << count_equal_c
            //<< std::endl;
          // slightly perturb the iso_value to avoid node cut
          if (count_larger_c == 0 and count_equal_c > 0)
          {
            if (std::abs(c_1 - iso_value) < 1e-12) c_1 = iso_value + iso_value * 0.001;
            if (std::abs(c_2 - iso_value) < 1e-12) c_2 = iso_value + iso_value * 0.001;
          }
          else if (count_smaller_c == 0 and count_equal_c > 0)
          {
            if (std::abs(c_1 - iso_value) < 1e-12) c_1 = iso_value - iso_value * 0.001;
            if (std::abs(c_2 - iso_value) < 1e-12) c_2 = iso_value - iso_value * 0.001;
          }
          else
          {
            if (count_equal_c > 0)
            {
              // all set to larger value
              if (std::abs(c_1 - iso_value) < 1e-12) c_1 = iso_value + iso_value * 0.001;
              if (std::abs(c_2 - iso_value) < 1e-12) c_2 = iso_value + iso_value * 0.001;
            }
          }

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
            }
            count++;
          }
        }
        double elem_length = cell_SDdata[cell_id].edge1_node.distance(cell_SDdata[cell_id].edge2_node);
        //std::cout << "elem_length " << elem_length << " " << cell_SDdata[cell_id].edge1_node << " "<< cell_SDdata[cell_id].edge2_node << std::endl;
        if (elem_length < 0.05) elem_length = 0.05;
        cell_SDdata[cell_id].interface_length = elem_length;

        double dx = cell_SDdata[cell_id].edge1_node[0] - cell_SDdata[cell_id].edge2_node[0];
        double dy = cell_SDdata[cell_id].edge1_node[1] - cell_SDdata[cell_id].edge2_node[1];
        double mid_x = 0.5*(cell_SDdata[cell_id].edge1_node[0] + cell_SDdata[cell_id].edge2_node[0]);
        double mid_y = 0.5*(cell_SDdata[cell_id].edge1_node[1] + cell_SDdata[cell_id].edge2_node[1]);
        //std::cout << "----- p1 ---- " << cell_SDdata[cell_id].edge1_node <<  " p2 " << cell_SDdata[cell_id].edge2_node << " dx " << dx <<  " dy " << dy  << std::endl;
        // two possible normal directions
        //std::cout << "----- normal ---- " << -dy <<  " " << dx << " or " << dy <<  " " << -dx  <<  " plus_node " << cell_SDdata[cell_id].one_plus_node<< std::endl;

        // correct outward normal for the plus region
        if ( (-dy * (cell_SDdata[cell_id].one_plus_node[0]-mid_x) + dx * (cell_SDdata[cell_id].one_plus_node[1] -mid_y)) < 0)
        {
            cell_SDdata[cell_id].crk_n[0] = -dy / sqrt(dy*dy+dx*dx);
            cell_SDdata[cell_id].crk_n[1] = dx / sqrt(dy*dy+dx*dx);
        }
        else
        {
            cell_SDdata[cell_id].crk_n[0] = dy / sqrt(dy*dy+dx*dx);
            cell_SDdata[cell_id].crk_n[1] = -dx / sqrt(dy*dy+dx*dx);
        }
        //std::cout << "----- final normal ---- " << cell_SDdata[cell_id].crk_n[0] <<  " " << cell_SDdata[cell_id].crk_n[1]  << " length :" << elem_length << std::endl;

        /// update the area_elem for the actual sizes
        /// should not do the following. As the crack length is smaller if the cutting region is changed. This is reflected in the local residual function.
        if ( cell_SDdata[cell_id].lnode_plus.size() == 1)
        {
          // https://www.mathopenref.com/coordtrianglearea.html
          Point<dim, double> A = cell->vertex(cell_SDdata[cell_id].lnode_plus[0]);
          Point<dim, double> B = cell_SDdata[cell_id].edge1_node;
          Point<dim, double> C = cell_SDdata[cell_id].edge2_node;
          double area = 0.5 * std::abs(A[0]*(B[1]-C[1]) + B[0]*(C[1]-A[1]) +  C[0]*(A[1]-B[1]));
          //std::cout << "----- plus node size :" << cell_SDdata[cell_id].lnode_plus.size() << " new area " << area << " old area: " << cell_SDdata[cell_id].area_elem << std::endl;
          //std::cout << " A " << A << std::endl;
          //std::cout << " B " << B << std::endl;
          //std::cout << " C " << C << std::endl;
          if (abs(area) < 1.e-12) area = 1.0e-10;
          cell_SDdata[cell_id].computed_area = area;
        }

        if ( cell_SDdata[cell_id].lnode_plus.size() == 3)
        {
          // https://www.mathopenref.com/coordtrianglearea.html
          Point<dim, double> A = cell->vertex(cell_SDdata[cell_id].lnode_minus[0]);
          Point<dim, double> B = cell_SDdata[cell_id].edge1_node;
          Point<dim, double> C = cell_SDdata[cell_id].edge2_node;
          double area = 0.5 * std::abs(A[0]*(B[1]-C[1]) + B[0]*(C[1]-A[1]) +  C[0]*(A[1]-B[1]));
          //std::cout << "----- plus node size :" << cell_SDdata[cell_id].lnode_plus.size() << " new area " << cell_SDdata[cell_id].area_elem - area << " old area: " << cell_SDdata[cell_id].area_elem << std::endl;
          if (abs(area) < 1.e-12) area = 1.0e-10;
          cell_SDdata[cell_id].computed_area = cell_SDdata[cell_id].area_elem - area;
        }

        if ( cell_SDdata[cell_id].lnode_plus.size() == 2)
        {
          // https://www.mathopenref.com/coordtrianglearea.html
          double area_1 =0, area_2 =0, area_3 =0, area_4 = 0;
          {
            Point<dim, double> A = cell->vertex(cell_SDdata[cell_id].lnode_plus[0]);
            Point<dim, double> B = cell_SDdata[cell_id].edge1_node;
            Point<dim, double> C = cell_SDdata[cell_id].edge2_node;
            area_1 = 0.5 * std::abs(A[0]*(B[1]-C[1]) + B[0]*(C[1]-A[1]) +  C[0]*(A[1]-B[1]));
            if (abs(area_1) < 1.e-12) area_1 = 1.0e-10;
          }

          {
            Point<dim, double> A = cell->vertex(cell_SDdata[cell_id].lnode_plus[1]);
            Point<dim, double> B = cell_SDdata[cell_id].edge1_node;
            Point<dim, double> C = cell_SDdata[cell_id].edge2_node;
            area_2 = 0.5 * std::abs(A[0]*(B[1]-C[1]) + B[0]*(C[1]-A[1]) +  C[0]*(A[1]-B[1]));
            if (abs(area_2) < 1.e-12) area_2 = 1.0e-10;
          }

          {
            Point<dim, double> A = cell->vertex(cell_SDdata[cell_id].lnode_plus[1]);
            Point<dim, double> B = cell->vertex(cell_SDdata[cell_id].lnode_plus[0]);
            Point<dim, double> C = cell_SDdata[cell_id].edge2_node;
            area_3 = 0.5 * std::abs(A[0]*(B[1]-C[1]) + B[0]*(C[1]-A[1]) +  C[0]*(A[1]-B[1]));
            if (abs(area_3) < 1.e-12) area_3 = 1.0e-10;
          }
          {
            Point<dim, double> A = cell->vertex(cell_SDdata[cell_id].lnode_plus[1]);
            Point<dim, double> B = cell->vertex(cell_SDdata[cell_id].lnode_plus[0]);
            Point<dim, double> C = cell_SDdata[cell_id].edge1_node;
            area_4 = 0.5 * std::abs(A[0]*(B[1]-C[1]) + B[0]*(C[1]-A[1]) +  C[0]*(A[1]-B[1]));
            if (abs(area_4) < 1.e-12) area_4 = 1.0e-10;
          }
          double area = 0.5 * (area_1 + area_2 + area_3 + area_4 );
          //std::cout << "----- plus node size :" << cell_SDdata[cell_id].lnode_plus.size() << " new area " << area << " old area: " << cell_SDdata[cell_id].area_elem << std::endl;
          cell_SDdata[cell_id].computed_area = area;
        }

        //
        //

        Triangulation<1> triangulation_1d;
        GridGenerator::hyper_cube    (    triangulation_1d, 0.,  elem_length);
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
            for (unsigned q=0; q< quadrature_formula.size(); ++q)
            {
              //std::cout << " q " << q << std::endl;
              vol += fe_values_1d.JxW(q);
              for (unsigned int i = 0; i < 2; ++i) {
                cell_SDdata[cell_id].shape_value_1d(i, q) = fe_values_1d.shape_value(i, q);
                cell_SDdata[cell_id].jxw_1d(q) = fe_values_1d.JxW(q);
                //r_local[i] += fe_values_1d.shape_value(i, q) * dRc * fe_values_1d.JxW(q);
              } // q_point
            }
        } // cell_1d
        //std::cout <<  " total elem # = " << triangulation_1d.n_active_cells() << " length: " << elem_length << " vol " << vol<< std::endl;
        //for (auto p0: cell_SDdata[cell_id].lnode_plus) std::cout << " plus node: " << p0 << std::endl;
        //for (auto p0: cell_SDdata[cell_id].lnode_minus) std::cout << " minus node: " << p0 << std::endl;
        //std::cout << " crk_n: " << cell_SDdata[cell_id].crk_n[0] << "\t" << cell_SDdata[cell_id].crk_n[1] << std::endl;
      } // interface id
        }        // this process
    } // cell

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


//set Dirichlet BC
template <int dim>
void battery<dim>::apply_boundary_condition()
{
  //std::cout << "apply BCs" << std::endl;
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

    std::cout << " setMultDomain, domain id might be changed!!!" << std::endl;
  std::vector<std::vector<double>> origin_list_json;
  (*params_json)["ElectroChemo"]["Origin list" ]["value"].get_to(origin_list_json);
  //for (auto i0 : origin_list_json)
    //for (auto j0 : i0)
      //std::cout << " j0 " << j0 << std::endl;

  double neg_electrode_line=(*params_json)["ElectroChemo"]["neg_electrode_line"];
  double pos_electrode_line=(*params_json)["ElectroChemo"]["pos_electrode_line"];
  int orientation=(*params_json)["ElectroChemo"]["orientation"];

  double iso_value= (*params_json)["ElectroChemo"]["iso_value"];
  iso_value = -1.0 * iso_value;

  // local refinement mesh

  // assign values to the diffuse interface
  typename hp::DoFHandler<dim>::active_cell_iterator cell = this->dof_handler.begin_active(), endc=this->dof_handler.end();
  for (;cell!=endc; ++cell){
    //
    //if (cell->subdomain_id() == this->this_mpi_process)
    {
      // try to detect if on the neg or pos line
      bool all_greater = true;
      bool all_smaller = true;
      bool on_neg_line = false;
      bool on_pos_line = false;
      Point<dim> center=cell->center();

      for (unsigned int i = 0; i < GeometryInfo<dim>::vertices_per_cell; ++i) {
        unsigned int vertex_id = i;
        if (cell->vertex(vertex_id)[orientation] >= neg_electrode_line) all_smaller = false;
        if (cell->vertex(vertex_id)[orientation] < neg_electrode_line) all_greater = false;
      }
      if ((not all_greater) and (not all_smaller)) on_neg_line = true;
      all_greater = true;
      all_smaller = true;

      for (unsigned int i = 0; i < GeometryInfo<dim>::vertices_per_cell; ++i) {
        unsigned int vertex_id = i;
        if (cell->vertex(vertex_id)[orientation] >= pos_electrode_line) all_smaller = false;
        if (cell->vertex(vertex_id)[orientation] < pos_electrode_line) all_greater = false;
      }
      if ((not all_greater) and (not all_smaller)) on_pos_line = true;
      all_greater = true;
      all_smaller = true;


      for (unsigned int i = 0; i < GeometryInfo<dim>::vertices_per_cell; ++i) {
        unsigned int vertex_id = i;

        // assign value based on multiple particle locations
        double Val = 1.0e12;
        for (unsigned int o1=0;o1<origin_list_json.size();o1++){
          double x0 = origin_list_json[o1][0];
          double y0 = origin_list_json[o1][1];
          double a = origin_list_json[o1][2];
          double b = origin_list_json[o1][3];
          double alpha = origin_list_json[o1][4]/180.0*3.14159265359;

          double cosa=cos(alpha);
          double sina=sin(alpha);
          double _x = cell->vertex(vertex_id)[0];
          double _y = cell->vertex(vertex_id)[1];
          double _val=0.5*(((_x-x0)*cosa + (_y-y0)*sina)*((_x-x0)*cosa + (_y-y0)*sina)/a/a + ((_x-x0)*sina - (_y-y0)*cosa)*((_x-x0)*sina - (_y-y0)*cosa)/b/b);
          if (_val < Val) Val = _val;
        }
        //if (Val < 0.6) std::cout << " Val" << Val << std::endl;
        //double val = iso_value + origin_list_json[origin_index][2] - distance;

        //if ((not on_neg_line) and (not on_pos_line) and (cell->vertex(vertex_id)[orientation] > neg_electrode_line and cell->vertex(vertex_id)[orientation] < pos_electrode_line))
        //{
          //val = 0.0;
          ////std::cout << "val=0.0" << std::endl;
        //}

        ////std::cout <<  " val = " <<  val << std::endl;
        //if (on_neg_line and val > iso_value)
        //{
          //val = iso_value + neg_electrode_line - cell->vertex(vertex_id)[orientation] ;
        ////std::cout <<  " neg val = " <<  val << std::endl;
        //}
        //if (on_pos_line and val > iso_value)
        //{
          //val = iso_value + cell->vertex(vertex_id)[orientation] - pos_electrode_line ;
        ////std::cout <<  " pos val = " <<  val << std::endl;
        //}

        double val = Val;
        if (val >= iso_value)
        {
          all_smaller = false;
        }
        if (val < iso_value)
        {
          all_greater = false;
        }
      }  // per_cell
      //std::cout
        //<< " all_smaller " << all_smaller
        //<< " all_greater " << all_greater
        //<< " not and not " << ((not all_greater) and (not all_smaller))
        //<< std::endl;
      if (all_smaller) cell->set_material_id(active_particle_id);
      if (all_greater) cell->set_material_id(electrolyte_id);
      if ((not all_greater) and (not all_smaller)) cell->set_material_id(interface_id);
      //if (cell->material_id() == 2) std::cout << " mat_id " << cell->material_id() << " center " << center << std::endl;

      //if (center[orientation] < neg_electrode_line) cell->set_material_id(li_metal_id); // disable li metal id 2022-05-19
      //if (on_neg_line) cell->set_material_id(li_metal_interface_id);// disable li metal id 2022-05-19

    } // this_mpi_process
  }

  {
      std::ofstream myfile;
      myfile.open ("dealii_mat.csv");
        myfile << "center_x,center_y,cell_id,mat_id\n";
    typename hp::DoFHandler<dim>::active_cell_iterator cell = this->dof_handler.begin_active(), endc=this->dof_handler.end();
    for (;cell!=endc; ++cell){
      //if (cell->subdomain_id() == this->this_mpi_process)
      { // wrong
        Point<dim> center=cell->center();
        int mat_id = cell->material_id();
        int cell_id = cell->active_cell_index();

        myfile << center[0] << "," << center[1] << "," << cell_id << ", " << mat_id << "\n";
      }
    }
      myfile.close();
  }
  {
        std::map<int, int> cell_mat_map;
        std::ifstream ifile("dealii_updated_mat_id.csv");
        std::string one_cell_id;
        std::string one_mat_id;

        if (ifile) {
          while ( !ifile.eof () ) {
                ifile >> one_cell_id;
                ifile >> one_mat_id;
                //std::cout << " one_cell_id " << one_cell_id << " mat_id " << one_mat_id << std::endl;
                cell_mat_map[std::stoi(one_cell_id)] = std::stoi(one_mat_id);
          }
          typename hp::DoFHandler<dim>::active_cell_iterator cell = this->dof_handler.begin_active(), endc=this->dof_handler.end();
          for (;cell!=endc; ++cell){
            //if (cell->subdomain_id() == this->this_mpi_process)
            { // wrong
              int cell_id = cell->active_cell_index();
              int mat_id = cell->material_id();
              int new_mat_id = cell_mat_map[cell_id];
              if (mat_id != new_mat_id and mat_id == electrolyte_id)
              {
                  cell->set_material_id(new_mat_id);
                  //std::cout << " assigned binder mat_id: old mat id " << mat_id << " new mat id " << new_mat_id << std::endl;
              }
              if (mat_id != new_mat_id and mat_id == interface_id)
              {
                  cell->set_material_id(new_mat_id);
                  //std::cout << " assigned binder mat_id: old mat id " << mat_id << " new mat id " << new_mat_id << std::endl;
              }

              //if (cell->material_id() == 4 or cell->material_id() == 6) {cell->set_material_id(1);} // for debugging purpose
              //if (cell->material_id() == 6) {cell->set_material_id(1);} // for debugging purpose
            }
          }
        } // ifile
  }


  {
      double separator_line=(*params_json)["ElectroChemo"]["separator_line"];
      typename hp::DoFHandler<dim>::active_cell_iterator cell = this->dof_handler.begin_active(), endc=this->dof_handler.end();
      for (;cell!=endc; ++cell){
          {
              int mat_id = cell->material_id();
       	      Point<dim> center=cell->center();
              if (center[orientation] < separator_line)
              {
                  if (mat_id == active_particle_id) cell->set_material_id(li_metal_id);
                  if (mat_id == interface_id)  cell->set_material_id(li_metal_interface_id);
          //std::cout << " assigned binder mat_id: old mat id " << mat_id << " new mat id " << new_mat_id << std::endl;
              }
             }
      }
  }


    this->pcout<<"setMultDomain"<<std::endl;
    //double r=(*params_json)["ElectroChemo"]["particle_R"];
    //double neg_electrode_line=(*params_json)["ElectroChemo"]["neg_electrode_line"];
    //double pos_electrode_line=(*params_json)["ElectroChemo"]["pos_electrode_line"];
    //int orientation=(*params_json)["ElectroChemo"]["orientation"];
  //{
    //typename hp::DoFHandler<dim>::active_cell_iterator cell = this->dof_handler.begin_active(), endc=this->dof_handler.end();
    //for (;cell!=endc; ++cell){
            //if (cell->subdomain_id() == this->this_mpi_process){
                //cell->set_material_id(electrolyte_id);
                //Point<dim> center=cell->center();
                //if (center[orientation]>neg_electrode_line and center[orientation]<pos_electrode_line) continue;
                //for(unsigned int ori_index=0;ori_index<origin_list_benchmark.size();ori_index++){
                    //Point<dim> origin(origin_list_benchmark[ori_index][0],origin_list_benchmark[ori_index][1]);
                    //if(center.distance(origin)<r-1.0e-15) {
                        //cell->set_material_id(active_particle_id);
                        //break;
                    //}
                //}
      //}
    //}
  //}

  //{
    //// create a layer of interface element
    //typename hp::DoFHandler<dim>::active_cell_iterator cell = this->dof_handler.begin_active(), endc=this->dof_handler.end();
    //for (;cell!=endc; ++cell){
            //if (cell->subdomain_id() == this->this_mpi_process){
                //if(cell->material_id()==active_particle_id){
                    //for (unsigned int f=0; f<GeometryInfo<dim>::faces_per_cell; ++f){
                        //if (cell->at_boundary(f) == false){
                            //if(cell->neighbor(f)->material_id()==electrolyte_id  and cell->neighbor(f)->has_children() == false){
                                //cell->set_material_id(interface_id);
                            //}
                        //}
                    //}
                //}
      //}
    //}
  //}

  //{
    //// fix a few missing interface element
    //typename hp::DoFHandler<dim>::active_cell_iterator cell = this->dof_handler.begin_active(), endc=this->dof_handler.end();
    //for (;cell!=endc; ++cell){
            //if (cell->subdomain_id() == this->this_mpi_process){
                //if(cell->material_id()==active_particle_id){
                    //Point<dim> center=cell->center();
          //int _neighbor_count = 0;
                    //for (unsigned int f=0; f<GeometryInfo<dim>::faces_per_cell; ++f){
                        //if (cell->at_boundary(f) == false){
                            //if(cell->neighbor(f)->material_id()==interface_id  and cell->neighbor(f)->has_children() == false){
                //_neighbor_count += 1;
                            //}
                        //}
                    //}
          //if (_neighbor_count == 2)
          //{
            //// bottom connection corner
            //if (center[orientation] >4 and  center[orientation] < 12) cell->set_material_id(interface_id);
            //// top connection corner
            //if (center[orientation] >42 and  center[orientation] < 50) cell->set_material_id(interface_id);
          //}
                //}
      //}
    //}
  //}

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

    int input_electrolyte_id=(*params_json)["Mechanics"]["electrolyte_id"];
    int input_active_particle_id=(*params_json)["Mechanics"]["active_particle_id"];
    int input_li_metal_id=(*params_json)["Mechanics"]["li_metal_id"];
  if (input_electrolyte_id != electrolyte_id)
  {
    std::cout << " input_electrolyte_id = " << input_electrolyte_id << " != electrolyte_id in battery.h = " << electrolyte_id  << " stop! "<< std::endl;
    exit(0);
  }
  if (input_active_particle_id != active_particle_id)
  {
    std::cout << " input_active_particle_id = " << input_active_particle_id << " != active_particle_id in battery.h = " << active_particle_id  << " stop! "<< std::endl;
    exit(0);
  }
  if (input_li_metal_id != li_metal_id)
  {
    std::cout << " input_li_metal_id = " << input_li_metal_id << " != li_metal_id in battery.h = " << li_metal_id  << " stop! "<< std::endl;
    exit(0);
  }

}

template <int dim>
void battery<dim>::apply_initial_condition()
{
  //std::cout<< this->this_mpi_process << "   Number of degrees of freedom: " << hpFEM<dim>::dof_handler.n_dofs() << std::endl;
  //std::cout<< this->this_mpi_process <<" apply init:  sol_prev size = " << this->solution_prev.size() <<std::endl;
  std::vector<std::vector<double>> origin_list_json;
  (*params_json)["ElectroChemo"]["Origin list" ]["value"].get_to(origin_list_json);

  double r=(*params_json)["ElectroChemo"]["particle_R"];
  double bandwitdh=(*params_json)["ElectroChemo"]["interface_bandwitdh"];
  
  double C_li_max_neg=(*params_json)["ElectroChemo"]["c_li_max_neg"];
  double C_li_max_pos=(*params_json)["ElectroChemo"]["c_li_max_pos"];
  double C_li_100_neg=(*params_json)["ElectroChemo"]["c_li_100_neg"];
  double C_li_100_pos=(*params_json)["ElectroChemo"]["c_li_100_pos"];
  
  double C_li_plus_0=(*params_json)["ElectroChemo"]["C_li_plus_0"];
  double separator_line=(*params_json)["ElectroChemo"]["separator_line"];
  int orientation=(*params_json)["ElectroChemo"]["orientation"];
  double neg_electrode_line=(*params_json)["ElectroChemo"]["neg_electrode_line"];
  double pos_electrode_line=(*params_json)["ElectroChemo"]["pos_electrode_line"];
  
  double iso_value=(*params_json)["ElectroChemo"]["iso_value"];

  typename hp::DoFHandler<dim>::active_cell_iterator cell = this->dof_handler.begin_active(), endc=this->dof_handler.end();
  for (;cell!=endc; ++cell){
    if (cell->subdomain_id() == this->this_mpi_process){
      //std::cout << " this->this_mpi_process: " << this->this_mpi_process << std::endl;
      double C_li_0=C_li_100_neg*C_li_max_neg;
      Point<dim> center=cell->center();
      if (center[orientation]>separator_line){ C_li_0=C_li_100_pos*C_li_max_pos;}
      hp::FEValues<dim> hp_fe_values (this->fe_collection, this->q_collection, update_values | update_quadrature_points);
      hp_fe_values.reinit (cell);
      const FEValues<dim> &fe_values = hp_fe_values.get_present_fe_values();
      const unsigned int dofs_per_cell = cell->get_fe().dofs_per_cell;
      std::vector<unsigned int> local_dof_indices (dofs_per_cell);
      cell->get_dof_indices (local_dof_indices);


      bool on_neg_line = false;
      bool on_pos_line = false;
      if (cell->material_id()==interface_id){
        bool all_greater = true;
        bool all_smaller = true;
        for (unsigned int i = 0; i < GeometryInfo<dim>::vertices_per_cell; ++i) {
          unsigned int vertex_id = i;
          if (cell->vertex(vertex_id)[orientation] >= neg_electrode_line) all_smaller = false;
          if (cell->vertex(vertex_id)[orientation] < neg_electrode_line) all_greater = false;
        }
        if ((not all_greater) and (not all_smaller)) on_neg_line = true;
        all_greater = true;
        all_smaller = true;

        for (unsigned int i = 0; i < GeometryInfo<dim>::vertices_per_cell; ++i) {
          unsigned int vertex_id = i;
          if (cell->vertex(vertex_id)[orientation] >= pos_electrode_line) all_smaller = false;
          if (cell->vertex(vertex_id)[orientation] < pos_electrode_line) all_greater = false;
        }
        if ((not all_greater) and (not all_smaller)) on_pos_line = true;
        all_greater = true;
        all_smaller = true;
      }

      for (unsigned int i=0; i<dofs_per_cell; ++i) {
        int ck = fe_values.get_fe().system_to_component_index(i).first;
        if (cell->material_id()==active_particle_id){
          if (ck==battery_fields.active_fields_index["Lithium"]) this->solution_prev(local_dof_indices[i])=C_li_0;//+0.01*(static_cast <double> (rand())/(static_cast <double>(RAND_MAX))-0.5);
          else if(ck==battery_fields.active_fields_index["Electrode_potential"]){
            //if(center[orientation]<separator_line) this->solution_prev(local_dof_indices[i])=electricChemoFormula.formula_Usc(C_li_100_neg,-1).val();
            //if(center[orientation]>separator_line)
              this->solution_prev(local_dof_indices[i])=electricChemoFormula.formula_Usc(C_li_100_pos,1).val();
          }
        }
        else if (cell->material_id()==electrolyte_id){
          if(ck==battery_fields.active_fields_index["Lithium_cation"]) this->solution_prev(local_dof_indices[i])=C_li_plus_0;
          else if(ck==battery_fields.active_fields_index["Electrolyte_potential"]) this->solution_prev(local_dof_indices[i])=1e-8;
        }
        else if (cell->material_id()==li_metal_id){
          if (ck==battery_fields.active_fields_index["Lithium"]) this->solution_prev(local_dof_indices[i])=C_li_0;//+0.01*(static_cast <double> (rand())/(static_cast <double>(RAND_MAX))-0.5);
          if(ck==battery_fields.active_fields_index["Electrode_potential"]){
            this->solution_prev(local_dof_indices[i])=electricChemoFormula.formula_Usc(C_li_100_neg,-1).val();
          }
        } // no potential, no lithium
        else if (cell->material_id()==additive_id){
          if(ck==battery_fields.active_fields_index["Electrode_potential"]){
            if(center[orientation]>separator_line) this->solution_prev(local_dof_indices[i])=electricChemoFormula.formula_Usc(C_li_100_pos,1).val();
          }
        }
        else if (cell->material_id()==li_metal_interface_id)
        {
          if (anode_opt == 0) //li metal
          {
            int vertex_id=i / (dofs_per_cell/GeometryInfo<dim>::vertices_per_cell);

            double val = - iso_value +  cell->vertex(vertex_id)[orientation] - neg_electrode_line;
            val = -1.0 * val;
            //std::cout << " val " << val << " iso_value " << iso_value << " cell_vertex " << cell->vertex(vertex_id)[orientation]  << " neg_line " << neg_electrode_line << std::endl;
            if (ck==battery_fields.active_fields_index["Diffuse_interface"]) this->solution_prev(local_dof_indices[i])= val;
            bool inside_li_metal_flag=false;
            if (val >= iso_value){inside_li_metal_flag=true;}

            if (inside_li_metal_flag){ // Li metal
              if (ck==battery_fields.active_fields_index["Lithium_cation"] or ck==battery_fields.active_fields_index["Electrolyte_potential"]){
                this->solution_prev(local_dof_indices[i])=0;
              }

              if (ck==battery_fields.active_fields_index["Lithium"]) this->solution_prev(local_dof_indices[i])=C_li_0;//+0.01*(static_cast <double> (rand())/(static_cast <double>(RAND_MAX))-0.5);
              if (ck==battery_fields.active_fields_index["Electrode_potential"]){
                //if(center[orientation]>separator_line) this->solution_prev(local_dof_indices[i])=electricChemoFormula.formula_Usc(C_li_100_pos,1).val();
                //if(center[orientation]<separator_line)
                this->solution_prev(local_dof_indices[i])=electricChemoFormula.formula_Usc(C_li_100_neg,-1).val();
              }
            }
            else{ // solid electrolyte
              if (ck==battery_fields.active_fields_index["Lithium"] or ck==battery_fields.active_fields_index["Electrode_potential"]){
                this->solution_prev(local_dof_indices[i])=0;
              }
              if (ck==battery_fields.active_fields_index["Lithium_cation"]) this->solution_prev(local_dof_indices[i])=C_li_plus_0;
              if (ck==battery_fields.active_fields_index["Electrolyte_potential"]) this->solution_prev(local_dof_indices[i])=1e-8;
            }
          }
          else if (anode_opt == 1) // graphite
          {
            //std::cout << "---------- in interface ---------------" << std::endl;
            int vertex_id=i / (dofs_per_cell/GeometryInfo<dim>::vertices_per_cell);
            Point<dim> vertex_point=cell->vertex(vertex_id);

            // assign value based on multiple particle locations
            double Val = 1.0e12;
            for (unsigned int o1=0;o1<origin_list_json.size();o1++){
              double x0 = origin_list_json[o1][0];
              double y0 = origin_list_json[o1][1];
              double a = origin_list_json[o1][2];
              double b = origin_list_json[o1][3];
              double alpha = origin_list_json[o1][4]/180.0*3.14159265359;

              double cosa=cos(alpha);
              double sina=sin(alpha);
              double _x = cell->vertex(vertex_id)[0];
              double _y = cell->vertex(vertex_id)[1];
              double _val=0.5*(((_x-x0)*cosa + (_y-y0)*sina)*((_x-x0)*cosa + (_y-y0)*sina)/a/a + ((_x-x0)*sina - (_y-y0)*cosa)*((_x-x0)*sina - (_y-y0)*cosa)/b/b);
              if (_val < Val) Val = _val;
            }
            double val = Val;

            // the -1.0 here is used to reuse the original implementation.
            val = -1.0 * val;
            if (ck==battery_fields.active_fields_index["Diffuse_interface"]) this->solution_prev(local_dof_indices[i]) = val;

            bool inside_particle_flag=false;

            if (val >= iso_value){inside_particle_flag=true;}

            if (inside_particle_flag){
              if (ck==battery_fields.active_fields_index["Lithium_cation"] or ck==battery_fields.active_fields_index["Electrolyte_potential"]){
                this->solution_prev(local_dof_indices[i])=0;
                //std::cout << " center " << center << " i " << i << " Lithium_cation &  Electrolyte_potential = 0 " << std::endl;
              }

              if (ck==battery_fields.active_fields_index["Lithium"]) this->solution_prev(local_dof_indices[i])=C_li_0;//+0.01*(static_cast <double> (rand())/(static_cast <double>(RAND_MAX))-0.5);
              if (ck==battery_fields.active_fields_index["Electrode_potential"]){
                //if(center[orientation]>separator_line)
                this->solution_prev(local_dof_indices[i])=electricChemoFormula.formula_Usc(C_li_100_neg,-1).val();
                //if(center[orientation]<separator_line) this->solution_prev(local_dof_indices[i])=electricChemoFormula.formula_Usc(C_li_100_neg,-1).val();
              }
            }
            else{
              if (ck==battery_fields.active_fields_index["Lithium"] or ck==battery_fields.active_fields_index["Electrode_potential"]){
                this->solution_prev(local_dof_indices[i])=0;
                //std::cout << " center " << center << " i " << i << " Lithium &  Electrode_potential = 0 " << std::endl;
              }
              if (ck==battery_fields.active_fields_index["Lithium_cation"]) this->solution_prev(local_dof_indices[i])=C_li_plus_0;
              if (ck==battery_fields.active_fields_index["Electrolyte_potential"]) this->solution_prev(local_dof_indices[i])=1e-8;
            }

          }
          else
          {
            std::cout << " anode option = " << anode_opt << " is not known! " << std::endl;
            exit(0);
          }


          //if (ck==battery_fields.active_fields_index["Lithium"] and this->this_mpi_process == 1) std::cout << this->this_mpi_process << " x_metal " << cell->vertex(vertex_id)[orientation] << " flag " << inside_li_metal_flag << " y= "<< cell->vertex(vertex_id)[1] << " val " << this->solution_prev(local_dof_indices[i]) << std::endl;

        }
        else if (cell->material_id()==additive_interface_id){

          int vertex_id=i / (dofs_per_cell/GeometryInfo<dim>::vertices_per_cell);
          Point<dim> vertex_point=cell->vertex(vertex_id);

          bool inside_flag=false;

          // assign value based on multiple particle locations
          double Val = 1.0e12;
          for (unsigned int o1=0;o1<origin_list_json.size();o1++){
            double x0 = origin_list_json[o1][0];
            double y0 = origin_list_json[o1][1];
            double a = origin_list_json[o1][2];
            double b = origin_list_json[o1][3];
            double alpha = origin_list_json[o1][4]/180.0*3.14159265359;

            double cosa=cos(alpha);
            double sina=sin(alpha);
            double _x = cell->vertex(vertex_id)[0];
            double _y = cell->vertex(vertex_id)[1];
            double _val=0.5*(((_x-x0)*cosa + (_y-y0)*sina)*((_x-x0)*cosa + (_y-y0)*sina)/a/a + ((_x-x0)*sina - (_y-y0)*cosa)*((_x-x0)*sina - (_y-y0)*cosa)/b/b);
            if (_val < Val) Val = _val;
          }
          double val = Val;
          val = -1.0 * val;

          if (ck==battery_fields.active_fields_index["Diffuse_interface"]) this->solution_prev(local_dof_indices[i])= val;// reverse the sign to make sure we can use the old code

          bool inside_additive_flag=false;
          if (val < iso_value){inside_additive_flag=true;}
          // inside and outside are all the same for the electrode potential
          if (ck==battery_fields.active_fields_index["Electrode_potential"]){

            if(center[orientation]>separator_line) this->solution_prev(local_dof_indices[i])=electricChemoFormula.formula_Usc(C_li_100_pos,1).val();
            if(center[orientation]<=separator_line) this->solution_prev(local_dof_indices[i])=electricChemoFormula.formula_Usc(C_li_100_neg,-1).val();
          }
          if (inside_additive_flag){ // additive
            if (ck==battery_fields.active_fields_index["Lithium"]) this->solution_prev(local_dof_indices[i])=0.0;
          }
          else{ // active particle
            if (ck==battery_fields.active_fields_index["Lithium"]) this->solution_prev(local_dof_indices[i])=C_li_0;//+0.01*(static_cast <double> (rand())/(static_cast <double>(RAND_MAX))-0.5);
          }
        }
        else if (cell->material_id()==interface_id) // deal with the interface for cathode
        {
          //std::cout << "---------- in interface ---------------" << std::endl;
          int vertex_id=i / (dofs_per_cell/GeometryInfo<dim>::vertices_per_cell);
          Point<dim> vertex_point=cell->vertex(vertex_id);

          // assign value based on multiple particle locations
          double Val = 1.0e12;
          for (unsigned int o1=0;o1<origin_list_json.size();o1++){
            double x0 = origin_list_json[o1][0];
            double y0 = origin_list_json[o1][1];
            double a = origin_list_json[o1][2];
            double b = origin_list_json[o1][3];
            double alpha = origin_list_json[o1][4]/180.0*3.14159265359;

            double cosa=cos(alpha);
            double sina=sin(alpha);
            double _x = cell->vertex(vertex_id)[0];
            double _y = cell->vertex(vertex_id)[1];
            double _val=0.5*(((_x-x0)*cosa + (_y-y0)*sina)*((_x-x0)*cosa + (_y-y0)*sina)/a/a + ((_x-x0)*sina - (_y-y0)*cosa)*((_x-x0)*sina - (_y-y0)*cosa)/b/b);
            if (_val < Val) Val = _val;
          }
          double val = Val;

          // the -1.0 here is used to reuse the original implementation.
          val = -1.0 * val;
          if (ck==battery_fields.active_fields_index["Diffuse_interface"]) this->solution_prev(local_dof_indices[i]) = val;

          bool inside_particle_flag=false;

          if (val >= iso_value){inside_particle_flag=true;}

          if (inside_particle_flag){
            if (ck==battery_fields.active_fields_index["Lithium_cation"] or ck==battery_fields.active_fields_index["Electrolyte_potential"]){
              this->solution_prev(local_dof_indices[i])=0;
              //std::cout << " center " << center << " i " << i << " Lithium_cation &  Electrolyte_potential = 0 " << std::endl;
            }

            if (ck==battery_fields.active_fields_index["Lithium"]) this->solution_prev(local_dof_indices[i])=C_li_0;//+0.01*(static_cast <double> (rand())/(static_cast <double>(RAND_MAX))-0.5);
            if (ck==battery_fields.active_fields_index["Electrode_potential"]){
              //if(center[orientation]>separator_line)
              this->solution_prev(local_dof_indices[i])=electricChemoFormula.formula_Usc(C_li_100_pos,1).val();
              //if(center[orientation]<separator_line) this->solution_prev(local_dof_indices[i])=electricChemoFormula.formula_Usc(C_li_100_neg,-1).val();
            }
          }
          else{
            if (ck==battery_fields.active_fields_index["Lithium"] or ck==battery_fields.active_fields_index["Electrode_potential"]){
              this->solution_prev(local_dof_indices[i])=0;
              //std::cout << " center " << center << " i " << i << " Lithium &  Electrode_potential = 0 " << std::endl;
            }
            if (ck==battery_fields.active_fields_index["Lithium_cation"]) this->solution_prev(local_dof_indices[i])=C_li_plus_0;
            if (ck==battery_fields.active_fields_index["Electrolyte_potential"]) this->solution_prev(local_dof_indices[i])=1e-8;
          }
        } //interface
      }
    }
  }



  this->solution_prev.compress(VectorOperation::insert);
  this->solution=this->solution_prev;

  //std::cout << "size of solution: " << this->solution.size() << std::endl;
  std::cout << " identify diffuse interface (start) " << std::endl;
  identify_diffuse_interface();
  std::cout << " identify diffuse interface (end) " << this->this_mpi_process << std::endl;

  constraints->clear ();
  
  DoFTools::make_hanging_node_constraints (this->dof_handler, *constraints);

  {

    Vector<double> localized_U(this->solution_prev);
    double separator_line=(*params_json)["ElectroChemo"]["separator_line"];
        double X_0=(*params_json)["Geometry"]["x_min"];
        double Y_0=(*params_json)["Geometry"]["y_min"];
        double Z_0=(*params_json)["Geometry"]["z_min"];
    
        double X_end=(*params_json)["Geometry"]["x_max"];
        double Y_end=(*params_json)["Geometry"]["y_max"];
        double Z_end=(*params_json)["Geometry"]["z_max"];
    
    int num_elem_x =(*params_json)["Geometry"]["num_elem_x"];
    int num_elem_y =(*params_json)["Geometry"]["num_elem_y"];
    int num_elem_z =(*params_json)["Geometry"]["num_elem_z"];

    double dh_x = (X_end - X_0)/num_elem_x;
    double dh_y = (Y_end - Y_0)/num_elem_y;
    double dh_z = (Z_end - Z_0)/num_elem_z;

    hp::FEValues<dim> hp_fe_values (this->fe_collection, this->q_collection, update_values | update_quadrature_points);
    bool is_one_node_Electrolyte_potential_fixed = false;
    // remove the minus or plus side of node, not to solve
    typename hp::DoFHandler<dim>::active_cell_iterator cell = this->dof_handler.begin_active(), endc = this->dof_handler.end();
    for (; cell != endc; ++cell) {
      if (cell->subdomain_id() == this->this_mpi_process) {
        int cell_id = cell->active_cell_index();
        int mat_id = cell->material_id();
        Point<dim> center=cell->center();
          hp_fe_values.reinit (cell);
          const FEValues<dim> &fe_values = hp_fe_values.get_present_fe_values();


        const unsigned int dofs_per_cell = cell->get_fe().dofs_per_cell;
        std::vector<unsigned int> local_dof_indices;
        local_dof_indices.resize(0);
        local_dof_indices.resize(dofs_per_cell);
        cell->get_dof_indices(local_dof_indices);
        const unsigned int dofs_per_node = dofs_per_cell / 4;

        {
          int _v_id = -1;
            int DOF_Electrolyte_potential = battery_fields.active_fields_index["Electrolyte_potential"];
          if (not is_one_node_Electrolyte_potential_fixed)
          {
            for (unsigned int i=0; i<dofs_per_cell; ++i) {
              const unsigned int ck = fe_values.get_fe().system_to_component_index(i).first;
              //std::cout << " ck " << ck << std::endl;
              if (ck==DOF_Electrolyte_potential) _v_id += 1;
              //if (ck==DOF_Displacement or ck== DOF_Displacement+1) std::cout << " cell_id " << cell_id << " ck " << ck << " _v_id " << _v_id << " vertex " << cell->vertex(_v_id) << std::endl;
              if (ck==DOF_Electrolyte_potential)
              {
                Point<dim> vertex_point=cell->vertex(_v_id);
                //std::cout << " vertex_point " << vertex_point << " X_0 " << X_0 << " Y_0 " << Y_0 << std::endl;
                //double dh = 0.0;
                //if (orientation == 0) dh = dh_x;
                //if (orientation == 1) dh = dh_y;
                //if (orientation == 2) dh = dh_z;

                if (
                    ( (orientation == 0) and (vertex_point[0] >= separator_line-0.5*dh_x and vertex_point[0] < separator_line+0.5*dh_x and abs(vertex_point[1] - Y_0) < 1e-5) )
                    or
                    ( (orientation == 0) and (vertex_point[0] >= separator_line-0.5*dh_x and vertex_point[0] < separator_line+0.5*dh_x and abs(vertex_point[1] - Y_end) < 1e-5) )
                    or
                    ( (orientation == 1) and (vertex_point[1] >= separator_line-0.5*dh_y and vertex_point[1] < separator_line+0.5*dh_y and abs(vertex_point[0] - X_0) < 1e-5) )
                    or
                    ( (orientation == 1) and (vertex_point[1] >= separator_line-0.5*dh_y and vertex_point[1] < separator_line+0.5*dh_y and abs(vertex_point[0] - X_end) < 1e-5) )
                    )
                 //if ( (vertex_point[orientation] >= separator_line-dh and vertex_point[orientation] < separator_line+dh) )
                //if (abs(vertex_point[1] - Y_0) < 1e-5 and abs(vertex_point[0] - X_0) < 1e-5)
                {
                  std::cout << "electrolyte potential fixed:  vertex_point " << vertex_point << " ck " << ck << " X_0 " << X_0 << " X_end " << X_end << " Y_0 " << Y_0 << " Y_end " << Y_end << " _v_id " << _v_id << std::endl;
                  auto globalDOF = local_dof_indices[i];
                  constraints->add_line(globalDOF);
                  constraints->set_inhomogeneity(globalDOF, 0.0);
                  //is_one_node_Electrolyte_potential_fixed = true;

                  //break;
                }
              }
            }
          }
        }




      //---------------debug-------------------
      //-------------- fix all displacement
      for (unsigned int i=0; i<dofs_per_cell; ++i) {
        const unsigned int ck = fe_values.get_fe().system_to_component_index(i).first;
        //std::cout << cell_id << " ck " << ck << std::endl;
        if (ck==battery_fields.active_fields_index["Diffuse_interface"])
        {
          auto globalDOF = local_dof_indices[i];
          constraints->add_line(globalDOF);
          constraints->set_inhomogeneity(globalDOF, 0.0);
          //std::cout << " fix " << ck << std::endl;
        }

        //if (cell_SDdata[cell_id].is_interface_element and mat_id == additive_interface_id) {
        //}

        auto globalDOF = local_dof_indices[i];
        if (globalDOF == 31930)
        {
          std::cout << i << " globalDOF " << local_dof_indices[i] << " is interface " << cell_SDdata[cell_id].is_interface_element << " cell_id " << cell_id << " center " << center << " mat_id " << mat_id << std::endl;
        }

        //if (ck==battery_fields.active_fields_index["Electrolyte_potential"])
        //{
          //auto globalDOF = local_dof_indices[i];
          //constraints->add_line(globalDOF);
          //constraints->set_inhomogeneity(globalDOF, 0.0);
          //std::cout << " fix " << ck << std::endl;
        //}
        //if (ck==battery_fields.active_fields_index["Electrode_potential"])
        //{
          //auto globalDOF = local_dof_indices[i];
          //constraints->add_line(globalDOF);
          //constraints->set_inhomogeneity(globalDOF, 0.0);
          //std::cout << " fix " << ck << std::endl;
        //}
        //if (ck==battery_fields.active_fields_index["Lithium"])
        //{
          //auto globalDOF = local_dof_indices[i];
          //constraints->add_line(globalDOF);
          //constraints->set_inhomogeneity(globalDOF, 0.0);
          //std::cout << " fix " << ck << std::endl;
        //}
        //if (ck==battery_fields.active_fields_index["Lithium_cation"])
        //{
          //auto globalDOF = local_dof_indices[i];
          //constraints->add_line(globalDOF);
          //constraints->set_inhomogeneity(globalDOF, 0.0);
          //std::cout << " fix " << ck << std::endl;
        //}
                //int DOF_Displacement = battery_fields.active_fields_index["Displacement"];
        //if (ck==DOF_Displacement or ck== DOF_Displacement+1)
        //{
          //auto globalDOF = local_dof_indices[i];
          //constraints->add_line(globalDOF);
          //constraints->set_inhomogeneity(globalDOF, 0.0);
          //std::cout << " fix " << ck << std::endl;
        //}
      }
      //---------------debug-------------------


        int _v_id = -1;
          int DOF_Displacement = battery_fields.active_fields_index["Displacement"];
        for (unsigned int i=0; i<dofs_per_cell; ++i) {
          const unsigned int ck = fe_values.get_fe().system_to_component_index(i).first;
          //std::cout << " ck " << ck << std::endl;
          if (ck==DOF_Displacement) _v_id += 1;
          //if (ck==DOF_Displacement or ck== DOF_Displacement+1) std::cout << " cell_id " << cell_id << " ck " << ck << " _v_id " << _v_id << " vertex " << cell->vertex(_v_id) << std::endl;
          if (ck==DOF_Displacement or ck== DOF_Displacement+1)
          {

            Point<dim> vertex_point=cell->vertex(_v_id);
            //std::cout
              //<< " v[0] " << vertex_point[0]
              //<< " v[1] " << vertex_point[1]
              //<< " ?1 " << (vertex_point[1] == Y_0)
              //<< " ?2 " << (vertex_point[1] == Y_end)
              //<< " ?3 " << (vertex_point[0] == X_0)
              //<< " ?4 " << (vertex_point[0] == X_end)
              //<< std::endl;
            // 1
            //if (vertex_point[1] == 0.0 or vertex_point[1] == 56.0 or vertex_point[0] == -15.0 or vertex_point[0] == 15.0) // fix x & y at all edges
            // 2
            //if (ck== DOF_Displacement+1 and (vertex_point[1] == 0.0 ) // fix y displacement
                //or
                //ck== DOF_Displacement and (vertex_point[0] == -15.0 )) // fix x displacement
            // 3
            //if ( (ck== DOF_Displacement+1 and (vertex_point[1] == 0.0 or vertex_point[1] == 56.0)) // fix y displacement
                //or
                // (ck== DOF_Displacement and (vertex_point[0] == -15.0 or vertex_point[0] == 15.0))) // fix x displacement
            // 4
            //if (vertex_point[1] == Y_0 or vertex_point[1] == Y_end or vertex_point[0] == X_0 or vertex_point[0] == X_end) // fix x & y at all edges

            if ((ck== DOF_Displacement+1 and (abs(vertex_point[1] - Y_0) < 1e-5 or abs(vertex_point[1] - Y_end) < 1e-5)) // fix y displacement
                or
               (ck== DOF_Displacement and ( abs(vertex_point[0] - X_0) < 1e-5 or abs(vertex_point[0] - X_end) < 1e-5))) // fix x displacement
            {
              //std::cout << " vertex_point " << vertex_point << " ck " << ck << " X_0 " << X_0 << " X_end " << X_end << " Y_0 " << Y_0 << " Y_end " << Y_end << " _v_id " << _v_id << std::endl;
              auto globalDOF = local_dof_indices[i];
              constraints->add_line(globalDOF);
              constraints->set_inhomogeneity(globalDOF, 0.0);
            }
          }
        }


        if (cell_SDdata[cell_id].is_interface_element and (mat_id == li_metal_interface_id or mat_id == interface_id)) {
          // for interface element, the following globalDOF is working fine, as each dof is defined.

          cell_SDdata[cell_id].reaction_rate_potential = 0.0;
          cell_SDdata[cell_id].reaction_rate_li = 0.0;

          for (auto i0 : cell_SDdata[cell_id].lnode_minus) // minus node list: lithium and phi_s is fixed = 0
          {
            if(battery_fields.active_fields_index["Lithium"]>-1)
            {
              auto globalDOF = local_dof_indices[i0*dofs_per_node + battery_fields.active_fields_index["Lithium"]];
              constraints->add_line(globalDOF);
              constraints->set_inhomogeneity(globalDOF, 0.0);
            }
            if(battery_fields.active_fields_index["Electrode_potential"]>-1)
            {
              auto globalDOF = local_dof_indices[i0*dofs_per_node + battery_fields.active_fields_index["Electrode_potential"]];
              constraints->add_line(globalDOF);
              constraints->set_inhomogeneity(globalDOF, 0.0);
            }
          }
          // plus node list: assign correct value of Lithium and electrode potential
          int _count_plus = 0;
          cell_SDdata[cell_id].xi_old(0) = 0.0;
          cell_SDdata[cell_id].xi_conv(0) = 0.0;
          cell_SDdata[cell_id].xi_old_phi_s(0) = 0.0;
          cell_SDdata[cell_id].xi_conv_phi_s(0) = 0.0;
          for (auto i0 : cell_SDdata[cell_id].lnode_plus)
          {
              // initialize the xi_old to be the same value as the neighbor node
              cell_SDdata[cell_id].xi_old(0) += localized_U(local_dof_indices[i0*dofs_per_node + battery_fields.active_fields_index["Lithium"]]);
              cell_SDdata[cell_id].xi_old_phi_s(0) += localized_U(local_dof_indices[i0*dofs_per_node + battery_fields.active_fields_index["Electrode_potential"]]);
            _count_plus += 1;
          }
          cell_SDdata[cell_id].xi_old(0) = cell_SDdata[cell_id].xi_old(0)/ _count_plus;
          cell_SDdata[cell_id].xi_old_phi_s(0) = cell_SDdata[cell_id].xi_old_phi_s(0)/ _count_plus;
          cell_SDdata[cell_id].xi_conv(0) = cell_SDdata[cell_id].xi_old(0);
          cell_SDdata[cell_id].xi_conv_phi_s(0) = cell_SDdata[cell_id].xi_old_phi_s(0);

          for (auto i0 : cell_SDdata[cell_id].lnode_plus) // plus node list: Lithium cation and phi_e is fixed
          {
            if(battery_fields.active_fields_index["Lithium_cation"]>-1)
            {
              auto globalDOF = local_dof_indices[i0*dofs_per_node + battery_fields.active_fields_index["Lithium_cation"]];
              constraints->add_line(globalDOF);
              constraints->set_inhomogeneity(globalDOF, 0.0);
            }

            if(battery_fields.active_fields_index["Electrolyte_potential"]>-1)
            {
              auto globalDOF = local_dof_indices[i0*dofs_per_node + battery_fields.active_fields_index["Electrolyte_potential"]];
              constraints->add_line(globalDOF);
              constraints->set_inhomogeneity(globalDOF, 0.0);
            }
          }

          int _count_minus = 0;
          cell_SDdata[cell_id].xi_old_c_e(0) = 0.0;
          cell_SDdata[cell_id].xi_conv_c_e(0) = 0.0;
          cell_SDdata[cell_id].xi_old_phi_e(0) = 0.0;
          cell_SDdata[cell_id].xi_conv_phi_e(0) = 0.0;
          for (auto i0 : cell_SDdata[cell_id].lnode_minus)
          {
              cell_SDdata[cell_id].xi_old_c_e(0) += localized_U(local_dof_indices[i0*dofs_per_node + battery_fields.active_fields_index["Lithium_cation"]]);
              cell_SDdata[cell_id].xi_old_phi_e(0) += localized_U(local_dof_indices[i0*dofs_per_node + battery_fields.active_fields_index["Electrolyte_potential"]]);
            _count_minus += 1;
          }
          cell_SDdata[cell_id].xi_old_c_e(0) = cell_SDdata[cell_id].xi_old_c_e(0)/ _count_minus;
          cell_SDdata[cell_id].xi_old_phi_e(0) = cell_SDdata[cell_id].xi_old_phi_e(0)/ _count_minus;
          cell_SDdata[cell_id].xi_conv_c_e(0) = cell_SDdata[cell_id].xi_old_c_e(0);
          cell_SDdata[cell_id].xi_conv_phi_e(0) = cell_SDdata[cell_id].xi_old_phi_e(0);
          for (int q=0; q<4; q++) cell_SDdata[cell_id].C_Li_plus_old[q] = cell_SDdata[cell_id].xi_old_c_e(0);

        } // is interface element

        if (cell_SDdata[cell_id].is_interface_element and mat_id == additive_interface_id) {
          // for interface element, the following globalDOF is working fine, as each dof is defined.
          cell_SDdata[cell_id].reaction_rate_potential = 0.0;
          cell_SDdata[cell_id].reaction_rate_li = 0.0;

          int _count_minus = 0;
          cell_SDdata[cell_id].xi_old(0) = 0.0;
          cell_SDdata[cell_id].xi_conv(0) = 0.0;
          cell_SDdata[cell_id].xi_old_phi_s(0) = 0.0;
          cell_SDdata[cell_id].xi_conv_phi_s(0) = 0.0;
          for (auto i0 : cell_SDdata[cell_id].lnode_minus)
          {
              // initialize the xi_old to be the same value as the neighbor node
              cell_SDdata[cell_id].xi_old(0) += localized_U(local_dof_indices[i0*dofs_per_node + battery_fields.active_fields_index["Lithium"]]);
              cell_SDdata[cell_id].xi_old_phi_s(0) += localized_U(local_dof_indices[i0*dofs_per_node + battery_fields.active_fields_index["Electrode_potential"]]);
            _count_minus += 1;
          }
          cell_SDdata[cell_id].xi_old(0) = cell_SDdata[cell_id].xi_old(0)/ _count_minus;
          cell_SDdata[cell_id].xi_old_phi_s(0) = cell_SDdata[cell_id].xi_old_phi_s(0)/ _count_minus;
          cell_SDdata[cell_id].xi_conv(0) = cell_SDdata[cell_id].xi_old(0);
          cell_SDdata[cell_id].xi_conv_phi_s(0) = cell_SDdata[cell_id].xi_old_phi_s(0);


          for (auto i0 : cell_SDdata[cell_id].lnode_minus)
          {
            if(battery_fields.active_fields_index["Lithium"]>-1)
            {
              auto globalDOF = local_dof_indices[i0*dofs_per_node + battery_fields.active_fields_index["Lithium"]];
              constraints->add_line(globalDOF);
              constraints->set_inhomogeneity(globalDOF, 0.0);
            }
          }

          cell_SDdata[cell_id].xi_old_c_e(0) = 0.0;
          cell_SDdata[cell_id].xi_conv_c_e(0) = 0.0;
          cell_SDdata[cell_id].xi_old_phi_e(0) = 0.0;
          cell_SDdata[cell_id].xi_conv_phi_e(0) = 0.0;
          for (int q=0; q<4; q++) cell_SDdata[cell_id].C_Li_plus_old[q] = cell_SDdata[cell_id].xi_old_c_e(0);

        } // is interface element
      }
    }
  }

  //std::cout<< this->this_mpi_process <<" before sleep "<<std::endl;
//https://dealii.org/7.3.0/doxygen/deal.II/classPETScWrappers_1_1MPI_1_1Vector.html
  //std::this_thread::sleep_for(std::chrono::seconds(5));
  //std::cout<< this->this_mpi_process <<" after sleep "<<std::endl;
    //int totalDOF=this->totalDOF(this->primary_variables);
    //std::vector<bool> All_component (totalDOF, false);
    //if(battery_fields.active_fields_index["Displacement"]>-1) {
        //for(unsigned int i=0;i<dim;i++) All_component[battery_fields.active_fields_index["Displacement"]+i]=true;
    //}
    //VectorTools:: interpolate_boundary_values (this->dof_handler, 1+orientation, ZeroFunction<dim> (totalDOF),*constraints, All_component);

  constraints->close ();
  std::cout<< this->this_mpi_process <<" end of apply init"<<std::endl;

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
