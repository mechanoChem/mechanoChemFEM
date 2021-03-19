#include "battery.h"


template <int dim>
void battery<dim>::get_residual_at_diffuse_interface(const typename hp::DoFHandler<dim>::active_cell_iterator &cell, const FEValues<dim>& fe_values, Table<1, Sacado::Fad::DFad<double> >& R, Table<1, Sacado::Fad::DFad<double>>& ULocal, Table<1, double >& ULocalConv)
{	
  int cell_id = cell->active_cell_index();
  unsigned int dofs_per_cell= fe_values.dofs_per_cell;
  unsigned int n_q_points= fe_values.n_quadrature_points;

  //std::cout << "compute residual at the interface " << cell_id << std::endl;
  //
	double reaction_rate=(*params_json)["ElectroChemo"]["jn_react"];
	double F=(*params_json)["ElectroChemo"]["F"];

  cell_SDdata[cell_id].reaction_rate_potential = reaction_rate*F;
  cell_SDdata[cell_id].reaction_rate_li = reaction_rate;
  //cell_SDdata[cell_id].reaction_rate_potential = 0.0;
  //cell_SDdata[cell_id].reaction_rate_li = 0.0;
  std::cout << " reaction rate: " << cell_SDdata[cell_id].reaction_rate_li << " potential " << cell_SDdata[cell_id].reaction_rate_potential << std::endl;

  // update reaction rate at the interface 
	double tem=(*params_json)["ElectroChemo"]["jn_react"];
	double fliptime=(*params_json)["ElectroChemo"]["flip_time"];
	int Li_index=battery_fields.active_fields_index["Lithium"];
	int Li_plus_index=battery_fields.active_fields_index["Lithium_cation"];



	int DOF_Lithium =  battery_fields.active_fields_index["Lithium"];
	int DOF_Lithium_phaseField = battery_fields.active_fields_index["Lithium_phaseField"];
	int DOF_Lithium_cation = battery_fields.active_fields_index["Lithium_cation"];
	int DOF_Electrode_potential = battery_fields.active_fields_index["Electrode_potential"];
	int DOF_Electrolyte_potential = battery_fields.active_fields_index["Electrolyte_potential"];
	int DOF_Displacement = battery_fields.active_fields_index["Displacement"];
	int DOF_Diffuse_interface = battery_fields.active_fields_index["Diffuse_interface"];

  std::vector<int> this_dof_local_index_Lithium;
  std::vector<int> this_dof_local_index_Lithium_phaseField;
  std::vector<int> this_dof_local_index_Lithium_cation;
  std::vector<int> this_dof_local_index_Electrode_potential;
  std::vector<int> this_dof_local_index_Electrolyte_potential;
  std::vector<int> this_dof_local_index_Displacement;
  std::vector<int> this_dof_local_index_Diffuse_interface;
  {
    //evaluate Residual: need to be modified later to include the jump of concentration
    for (unsigned int i=0; i<dofs_per_cell; ++i) {
      const unsigned int ck = fe_values.get_fe().system_to_component_index(i).first;
      if (ck==DOF_Lithium) this_dof_local_index_Lithium.push_back(i);
      if (ck==DOF_Lithium_phaseField) this_dof_local_index_Lithium_phaseField.push_back(i);
      if (ck==DOF_Lithium_cation) this_dof_local_index_Lithium_cation.push_back(i);
      if (ck==DOF_Electrode_potential) this_dof_local_index_Electrode_potential.push_back(i);
      if (ck==DOF_Electrolyte_potential) this_dof_local_index_Electrolyte_potential.push_back(i);
      if (ck==DOF_Displacement) this_dof_local_index_Displacement.push_back(i);
      if (ck==DOF_Diffuse_interface) this_dof_local_index_Diffuse_interface.push_back(i);
    }	
  }
  //std::cout << "--b0-1--" << std::endl;

  //std::cout 
    //<< " " << DOF_Lithium
    //<< " " << DOF_Lithium_phaseField
    //<< " " << DOF_Lithium_cation
    //<< " " << DOF_Electrode_potential
    //<< " " << DOF_Electrolyte_potential
    //<< " " << DOF_Displacement
    //<< " " << DOF_Diffuse_interface
    //<< std::endl;


  std::vector<unsigned int> local_dof_indices(dofs_per_cell);
  cell->get_dof_indices(local_dof_indices);


  // will remove later
  dealii::Table<1, double> coeff(n_q_points);
  for (unsigned int q = 0; q < n_q_points; q++) {
    coeff[q] = 1.0;
  }
		
  //std::cout << "--b0-0 (primary_dof)--: " << primiary_dof << " " << this->primiary_dof  << std::endl;
  //call residual functions
  //ResidualEq->residualForDiff_ReacEq(fe_values,primiary_dof, R,battery_fields->quad_fields[primiary_dof].value, battery_fields->quad_fields[primiary_dof].value_conv, diffu,react);

  // add AssertThrow message here later to replace 999.
  int ind_Lithium = 0;
  int ind_Lithium_phaseField = 999;
  int ind_Lithium_cation = 1;
  int ind_Electrode_potential = 2;
  int ind_Electrolyte_potential = 3;
  int ind_Displacement = 999;
  int ind_Diffuse_interface = 999;

  //Vector<double> dC_k1;
  Vector<double> dC_k1_Lithium;
  Vector<double> dC_k1_Lithium_phaseField;
  Vector<double> dC_k1_Lithium_cation;
  Vector<double> dC_k1_Electrode_potential;
  Vector<double> dC_k1_Electrolyte_potential;
  Vector<double> dC_k1_Displacement;
  Vector<double> dC_k1_Diffuse_interface;

  dC_k1_Lithium.reinit(4);
  dC_k1_Lithium_phaseField.reinit(4);
  dC_k1_Lithium_cation.reinit(4);
  dC_k1_Electrode_potential.reinit(4);
  dC_k1_Electrolyte_potential.reinit(4);
  dC_k1_Displacement.reinit(4);
  dC_k1_Diffuse_interface.reinit(4);

  ////std::cout << "--b0-2--" << std::endl;

  {
    int _i = -1;
    for (unsigned int i = 0; i < dofs_per_cell; ++i) {
      const unsigned int ck = fe_values.get_fe().system_to_component_index(i).first;
      if (ck == 0) _i += 1;
      if (ck==DOF_Lithium) dC_k1_Lithium[_i] = (ULocal[i].val()  - cell_SDdata[cell_id].ULocal_k[i]);
      if (ck==DOF_Lithium_phaseField) dC_k1_Lithium_phaseField[_i] = (ULocal[i].val()  - cell_SDdata[cell_id].ULocal_k[i]);
      if (ck==DOF_Lithium_cation) dC_k1_Lithium_cation[_i] = (ULocal[i].val()  - cell_SDdata[cell_id].ULocal_k[i]);
      if (ck==DOF_Electrode_potential) dC_k1_Electrode_potential[_i] = (ULocal[i].val()  - cell_SDdata[cell_id].ULocal_k[i]);
      if (ck==DOF_Electrolyte_potential) dC_k1_Electrolyte_potential[_i] = (ULocal[i].val()  - cell_SDdata[cell_id].ULocal_k[i]);
      if (ck==DOF_Displacement) dC_k1_Displacement[_i] = (ULocal[i].val()  - cell_SDdata[cell_id].ULocal_k[i]);
      if (ck==DOF_Diffuse_interface) dC_k1_Diffuse_interface[_i] = (ULocal[i].val()  - cell_SDdata[cell_id].ULocal_k[i]);
    }
    //std::cout << "--size-1--" <<  dC_k1_Lithium.size() << std::endl;
  }

  //Vector<double> dxi_k1;
  Vector<double> dxi_k1_Lithium;
  Vector<double> dxi_k1_Lithium_phaseField;
  Vector<double> dxi_k1_Lithium_cation;
  Vector<double> dxi_k1_Electrode_potential;
  Vector<double> dxi_k1_Electrolyte_potential;
  Vector<double> dxi_k1_Displacement;
  Vector<double> dxi_k1_Diffuse_interface;

  //dxi_k1.reinit(1);
  dxi_k1_Lithium.reinit(1);
  dxi_k1_Lithium_phaseField.reinit(1);
  dxi_k1_Lithium_cation.reinit(1);
  dxi_k1_Electrode_potential.reinit(1);
  dxi_k1_Electrolyte_potential.reinit(1);
  dxi_k1_Displacement.reinit(1);
  dxi_k1_Diffuse_interface.reinit(1);

  //Table<1, Sacado::Fad::DFad<double>> xi_0(1);  // define sacado xi_0 for stiffness calculation.
  Table<1, Sacado::Fad::DFad<double>> xi_0_Lithium(1);
  Table<1, Sacado::Fad::DFad<double>> xi_0_Lithium_phaseField(1);
  Table<1, Sacado::Fad::DFad<double>> xi_0_Lithium_cation(1);
  Table<1, Sacado::Fad::DFad<double>> xi_0_Electrode_potential(1);
  Table<1, Sacado::Fad::DFad<double>> xi_0_Electrolyte_potential(1);
  Table<1, Sacado::Fad::DFad<double>> xi_0_Displacement(1);
  Table<1, Sacado::Fad::DFad<double>> xi_0_Diffuse_interface(1);

  //if (primiary_dof != cell_SDdata[cell_id].opposite_flux_dof_li)
  //{ // for electrode
    cell_SDdata[cell_id].Kxic.vmult(dxi_k1_Lithium, dC_k1_Lithium);
    cell_SDdata[cell_id].rlocal -= dxi_k1_Lithium;
    cell_SDdata[cell_id].rlocal[0] += cell_SDdata[cell_id].reaction_rate_li * cell_SDdata[cell_id].interface_length;
    cell_SDdata[cell_id].Kxixi_inv.vmult(dxi_k1_Lithium, cell_SDdata[cell_id].rlocal);
    xi_0_Lithium[0] = cell_SDdata[cell_id].xi_old(0) + dxi_k1_Lithium(0);  
    cell_SDdata[cell_id].xi_old(0) = xi_0_Lithium[0].val();
    //std::cout << "--a0-0--"  << std::endl;
  //}
  //else
  //{// for electrolyte
    ////std::cout << "--a0-1--" << std::endl;
    cell_SDdata[cell_id].Kxic_c_e.vmult(dxi_k1_Lithium_cation, dC_k1_Lithium_cation);
    cell_SDdata[cell_id].rlocal_c_e -= dxi_k1_Lithium_cation;
    cell_SDdata[cell_id].rlocal_c_e[0] += (-1 *cell_SDdata[cell_id].reaction_rate_li) * cell_SDdata[cell_id].interface_length; // reaction rate li direction should not change
    cell_SDdata[cell_id].Kxixi_inv_c_e.vmult(dxi_k1_Lithium_cation, cell_SDdata[cell_id].rlocal_c_e);
    xi_0_Lithium_cation[0] = cell_SDdata[cell_id].xi_old_c_e(0) + dxi_k1_Lithium_cation(0);  
    cell_SDdata[cell_id].xi_old_c_e(0) = xi_0_Lithium_cation[0].val();
  //}


    cell_SDdata[cell_id].Kxic_phi_s.vmult(dxi_k1_Electrode_potential, dC_k1_Electrode_potential);
    cell_SDdata[cell_id].rlocal_phi_s -= dxi_k1_Electrode_potential;
    cell_SDdata[cell_id].rlocal_phi_s[0] += cell_SDdata[cell_id].reaction_rate_potential * cell_SDdata[cell_id].interface_length;
    cell_SDdata[cell_id].Kxixi_inv_phi_s.vmult(dxi_k1_Electrode_potential, cell_SDdata[cell_id].rlocal_phi_s);
    xi_0_Electrode_potential[0] = cell_SDdata[cell_id].xi_old_phi_s(0) + dxi_k1_Electrode_potential(0);  
    cell_SDdata[cell_id].xi_old_phi_s(0) = xi_0_Electrode_potential[0].val();

    //std::cout << "--a0-1--" << std::endl;
    cell_SDdata[cell_id].Kxic_phi_e.vmult(dxi_k1_Electrolyte_potential, dC_k1_Electrolyte_potential);
    cell_SDdata[cell_id].rlocal_phi_e -= dxi_k1_Electrolyte_potential;
    cell_SDdata[cell_id].rlocal_phi_e[0] += (-1 * cell_SDdata[cell_id].reaction_rate_potential) * cell_SDdata[cell_id].interface_length; // reaction rate li direction should not change
    cell_SDdata[cell_id].Kxixi_inv_phi_e.vmult(dxi_k1_Electrolyte_potential, cell_SDdata[cell_id].rlocal_phi_e);
    xi_0_Electrolyte_potential[0] = cell_SDdata[cell_id].xi_old_phi_e(0) + dxi_k1_Electrolyte_potential(0);  
    cell_SDdata[cell_id].xi_old_phi_e(0) = xi_0_Electrolyte_potential[0].val();

  //xi_0[0].diff(0, 1);

  const unsigned int total_local_xi_dof = 4;
  Table<1, Sacado::Fad::DFad<double> > ULocal_xi(dofs_per_cell + total_local_xi_dof);
  for (unsigned int i = 0; i < dofs_per_cell; ++i) {
    ULocal_xi[i]=ULocal[i].val();
  }
  ULocal_xi[dofs_per_cell + ind_Lithium] = xi_0_Lithium[0].val();
  ULocal_xi[dofs_per_cell + ind_Lithium_cation] = xi_0_Lithium_cation[0].val();
  ULocal_xi[dofs_per_cell + ind_Electrode_potential] = xi_0_Electrode_potential[0].val();
  ULocal_xi[dofs_per_cell + ind_Electrolyte_potential] = xi_0_Electrolyte_potential[0].val();
  for (unsigned int i = 0; i < dofs_per_cell + total_local_xi_dof; ++i) {
    ULocal_xi[i].diff (i, dofs_per_cell + total_local_xi_dof);
  }


	battery_fields.update_fields(cell, fe_values, ULocal_xi, ULocalConv);
  //double D_1 = 1.0;
  ////evaluate primary fields
  /// \todo add new interface related methods for the following in the battery components to address the tilde part.
  dealii::Table<2,Sacado::Fad::DFad<double> > diffu_Lithium(n_q_points, dim);
  dealii::Table<1,Sacado::Fad::DFad<double> > react_Lithium(n_q_points);
  dealii::Table<2,Sacado::Fad::DFad<double> > diffu_Lithium_cation(n_q_points, dim);
  dealii::Table<1,Sacado::Fad::DFad<double> > react_Lithium_cation(n_q_points);
  lithium.set_diffusion_reaction_term(diffu_Lithium, react_Lithium);
  lithium_cation.set_diffusion_reaction_term(diffu_Lithium_cation, react_Lithium_cation);

	dealii::Table<2,Sacado::Fad::DFad<double> > field_Electrode_potential(n_q_points,dim);
	dealii::Table<1,Sacado::Fad::DFad<double> > source_Electrode_potential(n_q_points);
	dealii::Table<2,Sacado::Fad::DFad<double> > field_Electrolyte_potential(n_q_points,dim);
	dealii::Table<1,Sacado::Fad::DFad<double> > source_Electrolyte_potential(n_q_points);
  phi_s.set_field_and_source_term(field_Electrode_potential, source_Electrode_potential);
  phi_e.set_field_and_source_term(field_Electrolyte_potential, source_Electrolyte_potential);


  ////std::cout << "--a1--" << xi_0[0]<< std::endl;

  Table<1, Sacado::Fad::DFad<double>> Rcc_Lithium(4);
  Table<1, Sacado::Fad::DFad<double>> Rcxi_Lithium(4);
  Table<1, Sacado::Fad::DFad<double>> rxic_Lithium(1);
  Table<1, Sacado::Fad::DFad<double>> rxixi_Lithium(1);
  Table<1, Sacado::Fad::DFad<double>> rr_Lithium(1);
  rxic_Lithium[0] = 0.0;
  rxixi_Lithium[0] = 0.0;
  rr_Lithium[0] = 0.0;
  for (unsigned int i = 0; i < 4; ++i) {
    Rcc_Lithium[i] = 0.0;
    Rcxi_Lithium[i] = 0.0;
  }

  Table<1, Sacado::Fad::DFad<double>> Rcc_Lithium_cation(4);
  Table<1, Sacado::Fad::DFad<double>> Rcxi_Lithium_cation(4);
  Table<1, Sacado::Fad::DFad<double>> rxic_Lithium_cation(1);
  Table<1, Sacado::Fad::DFad<double>> rxixi_Lithium_cation(1);
  Table<1, Sacado::Fad::DFad<double>> rr_Lithium_cation(1);
  rxic_Lithium_cation[0] = 0.0;
  rxixi_Lithium_cation[0] = 0.0;
  rr_Lithium_cation[0] = 0.0;
  for (unsigned int i = 0; i < 4; ++i) {
    Rcc_Lithium_cation[i] = 0.0;
    Rcxi_Lithium_cation[i] = 0.0;
  }

  Table<1, Sacado::Fad::DFad<double>> Rcc_Electrode_potential(4);
  Table<1, Sacado::Fad::DFad<double>> Rcxi_Electrode_potential(4);
  Table<1, Sacado::Fad::DFad<double>> rxic_Electrode_potential(1);
  Table<1, Sacado::Fad::DFad<double>> rxixi_Electrode_potential(1);
  Table<1, Sacado::Fad::DFad<double>> rr_Electrode_potential(1);
  rxic_Electrode_potential[0] = 0.0;
  rxixi_Electrode_potential[0] = 0.0;
  rr_Electrode_potential[0] = 0.0;
  for (unsigned int i = 0; i < 4; ++i) {
    Rcc_Electrode_potential[i] = 0.0;
    Rcxi_Electrode_potential[i] = 0.0;
  }

  Table<1, Sacado::Fad::DFad<double>> Rcc_Electrolyte_potential(4);
  Table<1, Sacado::Fad::DFad<double>> Rcxi_Electrolyte_potential(4);
  Table<1, Sacado::Fad::DFad<double>> rxic_Electrolyte_potential(1);
  Table<1, Sacado::Fad::DFad<double>> rxixi_Electrolyte_potential(1);
  Table<1, Sacado::Fad::DFad<double>> rr_Electrolyte_potential(1);
  rxic_Electrolyte_potential[0] = 0.0;
  rxixi_Electrolyte_potential[0] = 0.0;
  rr_Electrolyte_potential[0] = 0.0;
  for (unsigned int i = 0; i < 4; ++i) {
    Rcc_Electrolyte_potential[i] = 0.0;
    Rcxi_Electrolyte_potential[i] = 0.0;
  }

  dealii::Table<2, Sacado::Fad::DFad<double>> c_1_tilde_grad_Lithium(n_q_points, dim);
  dealii::Table<1, Sacado::Fad::DFad<double>> c_1_tilde_Lithium(n_q_points);
  dealii::Table<1, double> c_1_tilde_conv_Lithium(n_q_points);

  for (unsigned int q = 0; q < n_q_points; ++q) {
      c_1_tilde_Lithium[q] = 0.0;
      c_1_tilde_conv_Lithium[q] = 0.0;
    for (unsigned int j = 0; j < dim; j++) {
      c_1_tilde_grad_Lithium[q][j] = 0.0;
    }
  }

  dealii::Table<2, Sacado::Fad::DFad<double>> c_1_tilde_grad_Lithium_cation(n_q_points, dim);
  dealii::Table<1, Sacado::Fad::DFad<double>> c_1_tilde_Lithium_cation(n_q_points);
  dealii::Table<1, double> c_1_tilde_conv_Lithium_cation(n_q_points);

  for (unsigned int q = 0; q < n_q_points; ++q) {
      c_1_tilde_Lithium_cation[q] = 0.0;
      c_1_tilde_conv_Lithium_cation[q] = 0.0;
    for (unsigned int j = 0; j < dim; j++) {
      c_1_tilde_grad_Lithium_cation[q][j] = 0.0;
    }
  }

  dealii::Table<2, Sacado::Fad::DFad<double>> c_1_tilde_grad_Electrode_potential(n_q_points, dim);
  dealii::Table<1, Sacado::Fad::DFad<double>> c_1_tilde_Electrode_potential(n_q_points);
  dealii::Table<1, double> c_1_tilde_conv_Electrode_potential(n_q_points);

  for (unsigned int q = 0; q < n_q_points; ++q) {
      c_1_tilde_Electrode_potential[q] = 0.0;
      c_1_tilde_conv_Electrode_potential[q] = 0.0;
    for (unsigned int j = 0; j < dim; j++) {
      c_1_tilde_grad_Electrode_potential[q][j] = 0.0;
    }
  }

  dealii::Table<2, Sacado::Fad::DFad<double>> c_1_tilde_grad_Electrolyte_potential(n_q_points, dim);
  dealii::Table<1, Sacado::Fad::DFad<double>> c_1_tilde_Electrolyte_potential(n_q_points);
  dealii::Table<1, double> c_1_tilde_conv_Electrolyte_potential(n_q_points);

  for (unsigned int q = 0; q < n_q_points; ++q) {
      c_1_tilde_Electrolyte_potential[q] = 0.0;
      c_1_tilde_conv_Electrolyte_potential[q] = 0.0;
    for (unsigned int j = 0; j < dim; j++) {
      c_1_tilde_grad_Electrolyte_potential[q][j] = 0.0;
    }
  }

    
  //if (primiary_dof != cell_SDdata[cell_id].opposite_flux_dof_li)
  //{ // for electrode
    std::vector<double> Ms_list;
    std::vector<double> Ms_list_opposite;
    if (n_q_points != 4)
    {
      std::cout << "Please make sure total bulk quadrature point == 4, exiting ...." << std::endl;
      exit(0);
    }
    for (unsigned int q = 0; q < n_q_points; ++q) {
      Ms_list.push_back(0.0);
      Ms_list_opposite.push_back(0.0);
    }
    for (unsigned i = 0; i < cell_SDdata[cell_id].lnode_plus.size(); ++i) {
      Ms_list[cell_SDdata[cell_id].lnode_plus[i]] = 1.0;
    }

    for (unsigned i = 0; i < cell_SDdata[cell_id].lnode_minus.size(); ++i) {
      Ms_list_opposite[cell_SDdata[cell_id].lnode_minus[i]] = 1.0;
    }

    double dummy_area = 0.0;
    double dummy_area_opposite = 0.0;
    for (unsigned int q = 0; q < n_q_points; ++q) {
      //double Ms = 1.0;
      //if (fe_values.quadrature_point(q)[0] < 0.5){
        //Ms = 0.0;
      //}
      double Ms = Ms_list[q];
      dummy_area += Ms * fe_values.JxW(q); // the actually int area

      for (unsigned i = 0; i < cell_SDdata[cell_id].lnode_plus.size(); ++i) {
        int plus_node = this_dof_local_index_Lithium[cell_SDdata[cell_id].lnode_plus[i]]; // same as the Electrode_potential
        Ms -= fe_values.shape_value(plus_node, q);
        for (unsigned int j = 0; j < dim; j++) {
          c_1_tilde_grad_Lithium[q][j] -= fe_values.shape_grad(plus_node, q)[j] * ULocal_xi[dofs_per_cell+ ind_Lithium];  
          c_1_tilde_grad_Electrode_potential[q][j] -= fe_values.shape_grad(plus_node, q)[j] * ULocal_xi[dofs_per_cell+ ind_Electrode_potential];  
        }
      }
      c_1_tilde_Lithium[q] = Ms * ULocal_xi[dofs_per_cell+ ind_Lithium];
      c_1_tilde_conv_Lithium[q] = Ms * cell_SDdata[cell_id].xi_conv[0];

      c_1_tilde_Electrode_potential[q] = Ms * ULocal_xi[dofs_per_cell+ ind_Electrode_potential];
      c_1_tilde_conv_Electrode_potential[q] = Ms * cell_SDdata[cell_id].xi_conv_phi_s[0];

      double Ms_opposite = Ms_list_opposite[q];
      dummy_area_opposite += Ms_opposite * fe_values.JxW(q); // the actually int area
      for (unsigned i = 0; i < cell_SDdata[cell_id].lnode_minus.size(); ++i) {
        int minus_node = this_dof_local_index_Lithium_cation[cell_SDdata[cell_id].lnode_minus[i]]; // same as the electrolyte potential
        Ms_opposite -= fe_values.shape_value(minus_node, q);
        for (unsigned int j = 0; j < dim; j++) {
          c_1_tilde_grad_Lithium_cation[q][j] -= fe_values.shape_grad(minus_node, q)[j] * ULocal_xi[dofs_per_cell+ ind_Lithium_cation];  
          c_1_tilde_grad_Electrolyte_potential[q][j] -= fe_values.shape_grad(minus_node, q)[j] * ULocal_xi[dofs_per_cell+ ind_Electrolyte_potential];  
        }
      }
      c_1_tilde_Lithium_cation[q] = Ms_opposite * ULocal_xi[dofs_per_cell+ ind_Lithium_cation];
      c_1_tilde_conv_Lithium_cation[q] = Ms_opposite * cell_SDdata[cell_id].xi_conv_c_e[0];

      c_1_tilde_Electrolyte_potential[q] = Ms_opposite * ULocal_xi[dofs_per_cell+ ind_Electrolyte_potential];
      c_1_tilde_conv_Electrolyte_potential[q] = Ms_opposite * cell_SDdata[cell_id].xi_conv_phi_e[0];
    }



    rr_Lithium[0] = - cell_SDdata[cell_id].reaction_rate_li * cell_SDdata[cell_id].interface_length;
    rr_Lithium_cation[0] = - ( -1.0 * cell_SDdata[cell_id].reaction_rate_li) * cell_SDdata[cell_id].interface_length;
    rr_Electrode_potential[0] = - cell_SDdata[cell_id].reaction_rate_potential * cell_SDdata[cell_id].interface_length;
    rr_Electrolyte_potential[0] = - (-1.0 * cell_SDdata[cell_id].reaction_rate_potential) * cell_SDdata[cell_id].interface_length;

    ////std::cout << "--a2--" << rr[0] << std::endl;

    for (unsigned int q = 0; q < n_q_points; ++q) {
          rxixi_Lithium[0] += cell_SDdata[cell_id].interface_length / dummy_area * (cell_SDdata[cell_id].computed_area /dummy_area)  * ( c_1_tilde_grad_Lithium[q][0] * cell_SDdata[cell_id].crk_n[0] +  c_1_tilde_grad_Lithium[q][1] * cell_SDdata[cell_id].crk_n[1]) * fe_values.JxW(q); 
          rxixi_Lithium_cation[0] += cell_SDdata[cell_id].interface_length / dummy_area_opposite * (cell_SDdata[cell_id].computed_area /dummy_area_opposite)  * ( - c_1_tilde_grad_Lithium_cation[q][0] * cell_SDdata[cell_id].crk_n[0] -  c_1_tilde_grad_Lithium_cation[q][1] * cell_SDdata[cell_id].crk_n[1]) * fe_values.JxW(q); // crk_n direction is reversed
          rxixi_Electrode_potential[0] += cell_SDdata[cell_id].interface_length / dummy_area * (cell_SDdata[cell_id].computed_area /dummy_area)  * (  c_1_tilde_grad_Electrode_potential[q][0] * cell_SDdata[cell_id].crk_n[0] +  c_1_tilde_grad_Electrode_potential[q][1] * cell_SDdata[cell_id].crk_n[1]) * fe_values.JxW(q); // crk_n direction is reversed
          rxixi_Electrolyte_potential[0] += cell_SDdata[cell_id].interface_length / dummy_area_opposite * (cell_SDdata[cell_id].computed_area /dummy_area_opposite)  * (- c_1_tilde_grad_Electrolyte_potential[q][0] * cell_SDdata[cell_id].crk_n[0] - c_1_tilde_grad_Electrolyte_potential[q][1] * cell_SDdata[cell_id].crk_n[1]) * fe_values.JxW(q); 
     }

    //std::cout << "--a2-1--" << rxixi[0]  << std::endl;

    for (unsigned int q = 0; q < n_q_points; ++q) {
        rxic_Lithium[0] += cell_SDdata[cell_id].interface_length /  dummy_area * (cell_SDdata[cell_id].computed_area /dummy_area)  * (diffu_Lithium[q][0] * cell_SDdata[cell_id].crk_n[0] + diffu_Lithium[q][1] * cell_SDdata[cell_id].crk_n[1]) * fe_values.JxW(q); 
        rxic_Lithium_cation[0] += cell_SDdata[cell_id].interface_length /  dummy_area_opposite * (cell_SDdata[cell_id].computed_area /dummy_area_opposite)  * (- diffu_Lithium_cation[q][0] * cell_SDdata[cell_id].crk_n[0] - diffu_Lithium_cation[q][1] * cell_SDdata[cell_id].crk_n[1]) * fe_values.JxW(q); 

        rxic_Electrode_potential[0] += cell_SDdata[cell_id].interface_length /  dummy_area * (cell_SDdata[cell_id].computed_area /dummy_area)  * (field_Electrode_potential[q][0] * cell_SDdata[cell_id].crk_n[0] + field_Electrode_potential[q][1] * cell_SDdata[cell_id].crk_n[1]) * fe_values.JxW(q); 
        rxic_Electrolyte_potential[0] += cell_SDdata[cell_id].interface_length /  dummy_area_opposite * (cell_SDdata[cell_id].computed_area /dummy_area_opposite)  * (- field_Electrolyte_potential[q][0] * cell_SDdata[cell_id].crk_n[0] - field_Electrolyte_potential[q][1] * cell_SDdata[cell_id].crk_n[1]) * fe_values.JxW(q); 
    }

    rr_Lithium[0] += rxixi_Lithium[0] + rxic_Lithium[0];
    rr_Lithium_cation[0] += rxixi_Lithium_cation[0] + rxic_Lithium_cation[0];
    rr_Electrode_potential[0] += rxixi_Electrode_potential[0] + rxic_Electrode_potential[0];
    rr_Electrolyte_potential[0] += rxixi_Electrolyte_potential[0] + rxic_Electrolyte_potential[0];

    //std::cout << "--rxic a3-- [0]" << rxic_Lithium_cation[0] << std::endl;
    //std::cout << "--rr a3-- [0]" << rr_Lithium_cation[0] << std::endl;
    //
    //std::cout << "--rxic a3-- " << rxic_Electrode_potential[0] << std::endl;
    //std::cout << "--rr a3-- " << rr_Electrode_potential[0] << std::endl;



    dealii::Table<2,Sacado::Fad::DFad<double> > diffu_Lithium_tilde(n_q_points, dim);
    dealii::Table<1,Sacado::Fad::DFad<double> > react_Lithium_tilde(n_q_points);
    lithium.set_diffusion_reaction_term_interface(diffu_Lithium_tilde, react_Lithium_tilde, c_1_tilde_grad_Lithium);
    //lithium.set_diffusion_reaction_term_interface(diffu_Lithium_tilde, react_Lithium_tilde, battery_fields.quad_fields[DOF_Electrode_potential].value_grad);

    for (unsigned int q = 0; q < n_q_points; ++q) {
      for (unsigned i = 0; i < cell_SDdata[cell_id].lnode_plus.size(); ++i) {
        int plus_node = this_dof_local_index_Lithium[cell_SDdata[cell_id].lnode_plus[i]];
          //std::cout << "--a2-3--" << i << " " << q << " " << plus_node << " " << cell_SDdata[cell_id].lnode_plus[i] << std::endl;

        Rcxi_Lithium[cell_SDdata[cell_id].lnode_plus[i]] +=  fe_values.shape_value(plus_node, q)*((c_1_tilde_Lithium[q]- c_1_tilde_conv_Lithium [q])/this->current_dt)*fe_values.JxW(q);
        for (unsigned int j = 0; j < dim; j++) {
          /// \todo introduce the flux tilde based on the existing, not use D_1
          /// \todo check if - or + sign
          Rcxi_Lithium[cell_SDdata[cell_id].lnode_plus[i]] += diffu_Lithium_tilde[q][j] * fe_values.shape_grad(plus_node, q)[j] * fe_values.JxW(q);  
          //std::cout << "--a2-4--" << plus_node << " " << cell_SDdata[cell_id].lnode_plus[i]<< std::endl;
        }
      }
    }

    dealii::Table<2,Sacado::Fad::DFad<double> > diffu_Lithium_cation_tilde(n_q_points, dim);
    dealii::Table<1,Sacado::Fad::DFad<double> > react_Lithium_cation_tilde(n_q_points);
    //lithium_cation.set_diffusion_reaction_term_interface(diffu_Lithium_cation_tilde, react_Lithium_cation_tilde, c_1_tilde_grad_Electrolyte_potential, c_1_tilde_grad_Lithium_cation, c_1_tilde_Lithium_cation, battery_fields.quad_fields[DOF_Lithium_cation].value_conv);
    lithium_cation.set_diffusion_reaction_term_interface(diffu_Lithium_cation_tilde, react_Lithium_cation_tilde, battery_fields.quad_fields[DOF_Electrolyte_potential].value_grad, battery_fields.quad_fields[DOF_Lithium_cation].value_grad, battery_fields.quad_fields[DOF_Lithium_cation].value, battery_fields.quad_fields[DOF_Lithium_cation].value_conv);
    // carefully rewrite the code this part to make sure everything is correct.
    for (unsigned int q = 0; q < n_q_points; ++q) {
      for (unsigned i = 0; i < cell_SDdata[cell_id].lnode_minus.size(); ++i) {
        int minus_node = this_dof_local_index_Lithium_cation[cell_SDdata[cell_id].lnode_minus[i]];
          //std::cout << "--a2-3--" << i << " " << q << " " << minus_node << " " << cell_SDdata[cell_id].lnode_minus[i] << std::endl;
        Rcxi_Lithium_cation[cell_SDdata[cell_id].lnode_minus[i]] +=  fe_values.shape_value(minus_node, q)*((c_1_tilde_Lithium_cation[q]- c_1_tilde_conv_Lithium_cation[q])/this->current_dt)*fe_values.JxW(q);
        for (unsigned int j = 0; j < dim; j++) {
          Rcxi_Lithium_cation[cell_SDdata[cell_id].lnode_minus[i]] += diffu_Lithium_cation_tilde[q][j] * fe_values.shape_grad(minus_node, q)[j] * fe_values.JxW(q);  
          //std::cout << "--a2-4--" << minus_node << " " << cell_SDdata[cell_id].lnode_minus[i]<< std::endl;
        }
      }
    }


	  dealii::Table<2,Sacado::Fad::DFad<double> > field_Electrode_potential_tilde(n_q_points,dim);
	  dealii::Table<1,Sacado::Fad::DFad<double> > source_Electrode_potential_tilde(n_q_points);
    phi_s.set_field_and_source_term_interface(field_Electrode_potential_tilde, source_Electrode_potential_tilde, c_1_tilde_grad_Electrode_potential);
    for (unsigned int q = 0; q < n_q_points; ++q) {
      for (unsigned i = 0; i < cell_SDdata[cell_id].lnode_plus.size(); ++i) {
        int plus_node = this_dof_local_index_Electrode_potential[cell_SDdata[cell_id].lnode_plus[i]];
          //std::cout << "--a2-3--" << i << " " << q << " " << plus_node << " " << cell_SDdata[cell_id].lnode_plus[i] << std::endl;
        for (unsigned int j = 0; j < dim; j++) {
          Rcxi_Electrode_potential[cell_SDdata[cell_id].lnode_plus[i]] += field_Electrode_potential_tilde[q][j] * fe_values.shape_grad(plus_node, q)[j] * fe_values.JxW(q);  
          //std::cout << "--a2-4--" << plus_node << " " << cell_SDdata[cell_id].lnode_plus[i]<< std::endl;
          //std::cout << "--a2-4-field: " << q << " " << field_Electrode_potential_tilde[q][j] << " tidle " << c_1_tilde_grad_Electrode_potential[q][j]  << std::endl;
      
        }
      }
    }

	  dealii::Table<2,Sacado::Fad::DFad<double> > field_Electrolyte_potential_tilde(n_q_points,dim);
	  dealii::Table<1,Sacado::Fad::DFad<double> > source_Electrolyte_potential_tilde(n_q_points);
    //phi_e.set_field_and_source_term_interface(field_Electrolyte_potential_tilde, source_Electrolyte_potential_tilde, c_1_tilde_grad_Electrolyte_potential, c_1_tilde_grad_Lithium_cation, c_1_tilde_Lithium_cation, battery_fields.quad_fields[DOF_Lithium_cation].value_conv);
    phi_e.set_field_and_source_term_interface(field_Electrolyte_potential_tilde, source_Electrolyte_potential_tilde, battery_fields.quad_fields[DOF_Electrolyte_potential].value_grad, battery_fields.quad_fields[DOF_Lithium_cation].value_grad, battery_fields.quad_fields[DOF_Lithium_cation].value, battery_fields.quad_fields[DOF_Lithium_cation].value_conv);
    for (unsigned int q = 0; q < n_q_points; ++q) {
      for (unsigned i = 0; i < cell_SDdata[cell_id].lnode_minus.size(); ++i) {
        int minus_node = this_dof_local_index_Electrolyte_potential[cell_SDdata[cell_id].lnode_minus[i]];
          //std::cout << "--a2-3--" << i << " " << q << " " << minus_node << " " << cell_SDdata[cell_id].lnode_minus[i] << std::endl;
        for (unsigned int j = 0; j < dim; j++) {
          Rcxi_Electrolyte_potential[cell_SDdata[cell_id].lnode_minus[i]] += field_Electrolyte_potential_tilde[q][j] * fe_values.shape_grad(minus_node, q)[j] * fe_values.JxW(q);  
          //std::cout << "--a2-4--" << minus_node << " " << cell_SDdata[cell_id].lnode_minus[i]<< std::endl;
          //std::cout << "--a2-4--" << minus_node << " q " << q << " j " << j << " " << field_Electrolyte_potential_tilde[q][j]<< std::endl;
        }
      }
    }
    ////std::cout << "--a3--" << primiary_dof << " rr[0]" << rr[0] << std::endl;
    //std::cout << "--Rcxi a3-- [0]" << Rcxi_Electrolyte_potential[0] << std::endl;
    //std::cout << "--Rcxi a3-- [1]" << Rcxi_Electrolyte_potential[1] << std::endl;
    //std::cout << "--Rcxi a3-- [2]" << Rcxi_Electrolyte_potential[2] << std::endl;
    //std::cout << "--Rcxi a3-- [3]" << Rcxi_Electrolyte_potential[3] << std::endl;
    //std::cout << "--Rcxi a3-- [0]" << Rcxi_Lithium_cation[0] << std::endl;
    //std::cout << "--Rcxi a3-- [1]" << Rcxi_Lithium_cation[1] << std::endl;
    //std::cout << "--Rcxi a3-- [2]" << Rcxi_Lithium_cation[2] << std::endl;
    //std::cout << "--Rcxi a3-- [3]" << Rcxi_Lithium_cation[3] << std::endl;
    //std::cout << "--Rcxi a3-- [0]" << Rcxi_Electrode_potential[0] << std::endl;
    //std::cout << "--Rcxi a3-- [1]" << Rcxi_Electrode_potential[1] << std::endl;
    //std::cout << "--Rcxi a3-- [2]" << Rcxi_Electrode_potential[2] << std::endl;
    //std::cout << "--Rcxi a3-- [3]" << Rcxi_Electrode_potential[3] << std::endl;
    //
    //

    FullMatrix<double> Kxixi_Lithium;
    Kxixi_Lithium.reinit(1, 1);
    FullMatrix<double> Kxic_Lithium;
    Kxic_Lithium.reinit(1, 4);
    FullMatrix<double> Kcxi_Lithium;
    Kcxi_Lithium.reinit(4, 1);

    FullMatrix<double> Kxixi_Lithium_cation;
    Kxixi_Lithium_cation.reinit(1, 1);
    FullMatrix<double> Kxic_Lithium_cation;
    Kxic_Lithium_cation.reinit(1, 4);
    FullMatrix<double> Kcxi_Lithium_cation;
    Kcxi_Lithium_cation.reinit(4, 1);

    FullMatrix<double> Kxixi_Electrode_potential;
    Kxixi_Electrode_potential.reinit(1, 1);
    FullMatrix<double> Kxic_Electrode_potential;
    Kxic_Electrode_potential.reinit(1, 4);
    FullMatrix<double> Kcxi_Electrode_potential;
    Kcxi_Electrode_potential.reinit(4, 1);

    FullMatrix<double> Kxixi_Electrolyte_potential;
    Kxixi_Electrolyte_potential.reinit(1, 1);
    FullMatrix<double> Kxic_Electrolyte_potential;
    Kxic_Electrolyte_potential.reinit(1, 4);
    FullMatrix<double> Kcxi_Electrolyte_potential;
    Kcxi_Electrolyte_potential.reinit(4, 1);

    Kxixi_Lithium(0, 0) = rxixi_Lithium[0].dx(dofs_per_cell+ ind_Lithium);
    Kxixi_Lithium_cation(0, 0) = rxixi_Lithium_cation[0].dx(dofs_per_cell+ ind_Lithium_cation);
    Kxixi_Electrode_potential(0, 0) = rxixi_Electrode_potential[0].dx(dofs_per_cell+ ind_Electrode_potential);
    Kxixi_Electrolyte_potential(0, 0) = rxixi_Electrolyte_potential[0].dx(dofs_per_cell+ ind_Electrolyte_potential);

    {
      unsigned int _i = 0;
      for (unsigned int i=0; i< this_dof_local_index_Lithium.size(); i++)
      {
        Kxic_Lithium(0, _i) = rxic_Lithium[0].dx(this_dof_local_index_Lithium[i]);
        Kcxi_Lithium(_i, 0) = Rcxi_Lithium[_i].dx(dofs_per_cell + ind_Lithium);
        Rcc_Lithium[_i] = ULocal_xi[0] * 0.0;
        //std::cout << "rxic" << rxic_Lithium[0] << " ind " << this_dof_local_index_Lithium[i] << " Kxic " << Kxic_Lithium(0, _i) << std::endl;
        //std::cout << "Rcxi" << Rcxi_Lithium[_i] << " ind " << this_dof_local_index_Lithium[i] << " Kcxi " << Kcxi_Lithium(_i, 0) << std::endl;
        //std::cout << "Rcc" << Rcc_Lithium[_i] << std::endl;
        _i += 1;
      }
    }

    {
      unsigned int _i = 0;
      for (unsigned int i=0; i< this_dof_local_index_Lithium_cation.size(); i++)
      {
        Kxic_Lithium_cation(0, _i) = rxic_Lithium_cation[0].dx(this_dof_local_index_Lithium_cation[i]);
        Kcxi_Lithium_cation(_i, 0) = Rcxi_Lithium_cation[_i].dx(dofs_per_cell + ind_Lithium_cation);
        Rcc_Lithium_cation[_i] = ULocal_xi[0] * 0.0;
        _i += 1;
      }
    }

    {
      unsigned int _i = 0;
      for (unsigned int i=0; i< this_dof_local_index_Electrode_potential.size(); i++)
      {
        Kxic_Electrode_potential(0, _i) = rxic_Electrode_potential[0].dx(this_dof_local_index_Electrode_potential[i]);
        Kcxi_Electrode_potential(_i, 0) = Rcxi_Electrode_potential[_i].dx(dofs_per_cell + ind_Electrode_potential);
        Rcc_Electrode_potential[_i] = ULocal_xi[0] * 0.0;
        _i += 1;
      }
    }

    {
      unsigned int _i = 0;
      for (unsigned int i=0; i< this_dof_local_index_Electrolyte_potential.size(); i++)
      {
        Kxic_Electrolyte_potential(0, _i) = rxic_Electrolyte_potential[0].dx(this_dof_local_index_Electrolyte_potential[i]);
        Kcxi_Electrolyte_potential(_i, 0) = Rcxi_Electrolyte_potential[_i].dx(dofs_per_cell + ind_Electrolyte_potential);
        Rcc_Electrolyte_potential[_i] = ULocal_xi[0] * 0.0;
        _i += 1;
      }
    }

  //---------- compute Rcc-------------------

    for (unsigned int q = 0; q < n_q_points; ++q) {
      for (unsigned i = 0; i < cell_SDdata[cell_id].lnode_plus.size(); ++i) {
        int plus_node = this_dof_local_index_Lithium[cell_SDdata[cell_id].lnode_plus[i]];
          Rcc_Lithium[cell_SDdata[cell_id].lnode_plus[i]] +=  fe_values.shape_value(plus_node, q)*((battery_fields.quad_fields[DOF_Lithium].value[q]-battery_fields.quad_fields[DOF_Lithium].value_conv[q])/this->current_dt)*fe_values.JxW(q);
        for (unsigned int j = 0; j < dim; j++) {
            Rcc_Lithium[cell_SDdata[cell_id].lnode_plus[i]] += -fe_values.shape_grad(plus_node, q)[j]*diffu_Lithium[q][j]*fe_values.JxW(q); // oscillation
        }
      }
    }

    for (unsigned int q = 0; q < n_q_points; ++q) {
      for (unsigned i = 0; i < cell_SDdata[cell_id].lnode_plus.size(); ++i) {
        int plus_node = this_dof_local_index_Electrode_potential[cell_SDdata[cell_id].lnode_plus[i]];
        for (unsigned int j = 0; j < dim; j++) {
            Rcc_Electrode_potential[cell_SDdata[cell_id].lnode_plus[i]] += -fe_values.shape_grad(plus_node, q)[j]*field_Electrode_potential[q][j]*fe_values.JxW(q); // oscillation
        }
      }
    }

    for (unsigned int q = 0; q < n_q_points; ++q) {
      for (unsigned i = 0; i < cell_SDdata[cell_id].lnode_minus.size(); ++i) {
        int minus_node = this_dof_local_index_Lithium_cation[cell_SDdata[cell_id].lnode_minus[i]];
          Rcc_Lithium_cation[cell_SDdata[cell_id].lnode_minus[i]] +=  fe_values.shape_value(minus_node, q)*((battery_fields.quad_fields[DOF_Lithium_cation].value[q]-battery_fields.quad_fields[DOF_Lithium_cation].value_conv[q])/this->current_dt)*fe_values.JxW(q);
        for (unsigned int j = 0; j < dim; j++) {
            Rcc_Lithium_cation[cell_SDdata[cell_id].lnode_minus[i]] += -fe_values.shape_grad(minus_node, q)[j]*diffu_Lithium_cation[q][j]*fe_values.JxW(q); // oscillation
        }
      }
    }

    for (unsigned int q = 0; q < n_q_points; ++q) {
      for (unsigned i = 0; i < cell_SDdata[cell_id].lnode_minus.size(); ++i) {
        int minus_node = this_dof_local_index_Electrolyte_potential[cell_SDdata[cell_id].lnode_minus[i]];
        for (unsigned int j = 0; j < dim; j++) {
            Rcc_Electrolyte_potential[cell_SDdata[cell_id].lnode_minus[i]] += -fe_values.shape_grad(minus_node, q)[j]*field_Electrolyte_potential[q][j]*fe_values.JxW(q); // oscillation
        }
      }
    }

    //std::cout << "--Rcc a4-- [0]" << Rcc_Lithium_cation[0] << std::endl;
    //std::cout << "--Rcc a4-- [1]" << Rcc_Lithium_cation[1] << std::endl;
    //std::cout << "--Rcc a4-- [2]" << Rcc_Lithium_cation[2] << std::endl;
    //std::cout << "--Rcc a4-- [3]" << Rcc_Lithium_cation[3] << std::endl;

    //std::cout << "--Rcc a4-- [0]" << Rcc_Electrode_potential[0] << std::endl;
    //std::cout << "--Rcc a4-- [1]" << Rcc_Electrode_potential[1] << std::endl;
    //std::cout << "--Rcc a4-- [2]" << Rcc_Electrode_potential[2] << std::endl;
    //std::cout << "--Rcc a4-- [3]" << Rcc_Electrode_potential[3] << std::endl;

    //std::cout << "--Rcc a4-- [0]" << Rcc_Electrolyte_potential[0] << std::endl;
    //std::cout << "--Rcc a4-- [1]" << Rcc_Electrolyte_potential[1] << std::endl;
    //std::cout << "--Rcc a4-- [2]" << Rcc_Electrolyte_potential[2] << std::endl;
    //std::cout << "--Rcc a4-- [3]" << Rcc_Electrolyte_potential[3] << std::endl;
    //std::cout << "--a4-- [0]" << Rcc_Lithium[0] << std::endl;
    //std::cout << "--a4-- [1]" << Rcc_Lithium[1] << std::endl;
    //std::cout << "--a4-- [2]" << Rcc_Lithium[2] << std::endl;
    //std::cout << "--a4-- [3]" << Rcc_Lithium[3] << std::endl;
    //std::cout << "--a4-- [0]" << Rcc_Lithium_cation[0] << std::endl;
    //std::cout << "--a4-- [1]" << Rcc_Lithium_cation[1] << std::endl;
    //std::cout << "--a4-- [2]" << Rcc_Lithium_cation[2] << std::endl;
    //std::cout << "--a4-- [3]" << Rcc_Lithium_cation[3] << std::endl;

    int _i = -1;
    for (unsigned int i = 0; i < dofs_per_cell; ++i) {
      const unsigned int ck = fe_values.get_fe().system_to_component_index(i).first;
      if (ck == 0) _i++;
      if (ck == DOF_Lithium ) {
        //std::cout << "R[i] " << R[i] << std::endl;
        R[i] = (Rcc_Lithium[_i] + Rcxi_Lithium[_i] - Kcxi_Lithium(_i,0) / Kxixi_Lithium(0,0) * rr_Lithium[0]) * (cell_SDdata[cell_id].computed_area /dummy_area) ; 
        //std::cout << "Rcc " << Rcc[_i] << " " << Rcxi[_i] << " " << Kcxi(_i,0) << " " << Kxixi(0,0) * rr[0] << std::endl;
        //std::cout << "--a4-1-- " << i <<" " << R[i] << cell_SDdata[cell_id].computed_area /dummy_area<< std::endl;
      }
      if (ck == DOF_Electrode_potential ) {
        R[i] = (Rcc_Electrode_potential[_i] + Rcxi_Electrode_potential[_i] - Kcxi_Electrode_potential(_i,0) / Kxixi_Electrode_potential(0,0) * rr_Electrode_potential[0]) * (cell_SDdata[cell_id].computed_area /dummy_area) ; 
      }

      if (ck == DOF_Lithium_cation ) {
        R[i] = (Rcc_Lithium_cation[_i] + Rcxi_Lithium_cation[_i] - Kcxi_Lithium_cation(_i,0) / Kxixi_Lithium_cation(0,0) * rr_Lithium_cation[0]) * (cell_SDdata[cell_id].computed_area /dummy_area_opposite); 
      }

      if (ck == DOF_Electrolyte_potential ) {
        //std::cout << "R[i] " << R[i] << std::endl;
        R[i] = (Rcc_Electrolyte_potential[_i] + Rcxi_Electrolyte_potential[_i] - Kcxi_Electrolyte_potential(_i,0) / Kxixi_Electrolyte_potential(0,0) * rr_Electrolyte_potential[0]) * (cell_SDdata[cell_id].computed_area /dummy_area_opposite); 
        //std::cout << "Rcc " << Rcc_Electrolyte_potential[_i] << " " << Rcxi_Electrolyte_potential[_i] << " " << Kcxi_Electrolyte_potential(_i,0) << " " << Kxixi_Electrolyte_potential(0,0) * rr_Electrolyte_potential[0] << std::endl;
        //std::cout << "--a4-1-- " << i <<" " << R[i] << cell_SDdata[cell_id].computed_area /dummy_area_opposite<< std::endl;
      }
    }

    cell_SDdata[cell_id].Kxic = Kxic_Lithium;
    cell_SDdata[cell_id].rlocal[0] = rr_Lithium[0].val();
    cell_SDdata[cell_id].Kxixi_inv(0,0) = 1.0/Kxixi_Lithium(0,0);

    cell_SDdata[cell_id].Kxic_c_e = Kxic_Lithium_cation;
    cell_SDdata[cell_id].rlocal_c_e[0] = rr_Lithium_cation[0].val();
    cell_SDdata[cell_id].Kxixi_inv_c_e(0,0) = 1.0/Kxixi_Lithium_cation(0,0);

    cell_SDdata[cell_id].Kxic_phi_s = Kxic_Electrode_potential;
    cell_SDdata[cell_id].rlocal_phi_s[0] = rr_Electrode_potential[0].val();
    cell_SDdata[cell_id].Kxixi_inv_phi_s(0,0) = 1.0/Kxixi_Electrode_potential(0,0);

    cell_SDdata[cell_id].Kxic_phi_e = Kxic_Electrolyte_potential;
    cell_SDdata[cell_id].rlocal_phi_e[0] = rr_Electrolyte_potential[0].val();
    cell_SDdata[cell_id].Kxixi_inv_phi_e(0,0) = 1.0/Kxixi_Electrolyte_potential(0,0);



  // not going to work in this two fields
	if(battery_fields.active_fields_index["Diffuse_interface"]>-1) diffuse_interface.r_get_residual(fe_values, R, ULocal, ULocalConv);
	if(battery_fields.active_fields_index["Lithium_phaseField"]>-1) lithium_mu.r_get_residual(fe_values, R, ULocal, ULocalConv);
  //if(battery_fields.active_fields_index["Diffuse_interface"]>-1)   std::cout << "!!!!!!!!!! In diffusive interface !!!!!" << std::endl;

  //// the following should still work, as Lithium and Electrode_potential are not coupled
  //if(battery_fields.active_fields_index["Lithium"]>-1) lithium.r_get_residual_with_interface(cell, fe_values, R, ULocal, ULocalConv, cell_SDdata);
  //if(battery_fields.active_fields_index["Electrode_potential"]>-1) phi_s.r_get_residual_with_interface(cell, fe_values, R, ULocal, ULocalConv, cell_SDdata);

  //// the following needs to be replaced.
  //if(battery_fields.active_fields_index["Lithium_cation"]>-1) lithium_cation.r_get_residual_with_interface(cell, fe_values, R, ULocal, ULocalConv, cell_SDdata);
  //if(battery_fields.active_fields_index["Electrolyte_potential"]>-1) phi_e.r_get_residual_with_interface(cell, fe_values, R, ULocal, ULocalConv, cell_SDdata);


    for (unsigned int i = 0; i < dofs_per_cell; ++i) {
      const unsigned int ck = fe_values.get_fe().system_to_component_index(i).first;
      if (ck == DOF_Lithium ) {
        //std::cout << "R_Lithium[i] " << R[i] << std::endl;
      }
      if (ck == DOF_Electrode_potential ) {
        //std::cout << "R_Electrode_potential[i] " << R[i] << std::endl;
      }

      if (ck == DOF_Lithium_cation ) {
        //std::cout << "R_Lithium_cation[i] " << R[i] << std::endl;
      }

      if (ck == DOF_Electrolyte_potential ) {
        //std::cout << "R_Electrolyte_potential[i] " << R[i] << std::endl;
      }
      if (R[i].val() != R[i].val()) exit(-1);
    }

    //std::cout <<  " Xi[Lithium] " << ULocal_xi[dofs_per_cell+ ind_Lithium] << std::endl;
    //std::cout <<  " Xi[Lithium_cation] " << ULocal_xi[dofs_per_cell+ ind_Lithium_cation] << std::endl;
    //std::cout <<  " Xi[Electrode_potential] " << ULocal_xi[dofs_per_cell+ ind_Electrode_potential] << std::endl;
    //std::cout <<  " Xi[Electrolyte_potential] " << ULocal_xi[dofs_per_cell+ ind_Electrolyte_potential] << std::endl;

  // save previous iteration solution for SD
	for (unsigned int i=0; i<dofs_per_cell; ++i){
    cell_SDdata[cell_id].ULocal_k[i] = ULocal[i].val(); 
	}
}

template class battery<1>;
template class battery<2>;
template class battery<3>;
