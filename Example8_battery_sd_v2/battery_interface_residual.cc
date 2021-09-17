#include "battery.h"


template <int dim>
void battery<dim>::get_residual_at_diffuse_interface(const typename hp::DoFHandler<dim>::active_cell_iterator &cell, const FEValues<dim>& fe_values, Table<1, Sacado::Fad::DFad<double> >& R, Table<1, Sacado::Fad::DFad<double>>& ULocal, Table<1, double >& ULocalConv)
{	
	double separator_line=(*params_json)["ElectroChemo"]["separator_line"];
	int orientation=(*params_json)["ElectroChemo"]["orientation"];
	Point<dim> center=cell->center();
	int domainflag=-1;
	if (center[orientation]>separator_line){
		domainflag=1;
	}

  bool disable_wd = false;
  if ( int((*params_json)["Mechanics"]["EnableWeak"]) == 0 )
  {
    disable_wd = true;
  }

  // update reaction rate at the interface 
	double Temp=(*params_json)["ElectroChemo"]["T_0"];

  int cell_id = cell->active_cell_index();
  unsigned int dofs_per_cell= fe_values.dofs_per_cell;
  unsigned int n_q_points= fe_values.n_quadrature_points;


  //if (center[0] < 5) std::cout << "\n--------------" << std::endl;
  //std::cout << "compute residual at the interface " << cell_id  << " center: " << center << std::endl;
  //
	double reaction_rate=(*params_json)["ElectroChemo"]["jn_react"];
	double F=(*params_json)["ElectroChemo"]["F"];

  //cell_SDdata[cell_id].reaction_rate_potential = 0.0;
  //cell_SDdata[cell_id].reaction_rate_li = 0.0;
  //std::cout << "*reaction rate (before) : " << cell_SDdata[cell_id].reaction_rate_li.val() << " potential (rate) " << cell_SDdata[cell_id].reaction_rate_potential.val() << std::endl;

  // update reaction rate at the interface 
	double tem=(*params_json)["ElectroChemo"]["jn_react"];
	double fliptime=(*params_json)["ElectroChemo"]["flip_time"];
	double cod_max=(*params_json)["ElectroChemo"]["cod_max"];
	double cod_coef_C0=(*params_json)["ElectroChemo"]["cod_coef_C0"];

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
      if (ck>=DOF_Displacement and ck <DOF_Displacement+dim) this_dof_local_index_Displacement.push_back(i);
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
  int ind_Displacement_sd_1 = 4;
  int ind_Displacement_sd_2 = 5;
  int ind_Displacement_wd = 6;
  int ind_Diffuse_interface = 999;

  std::vector<double> N_reference;
  std::vector<double> M_reference;
  N_reference.resize(0);
  M_reference.resize(0);
  N_reference.push_back(cell_SDdata[cell_id].crk_n[0]);
  N_reference.push_back(cell_SDdata[cell_id].crk_n[1]);
  M_reference.push_back(-cell_SDdata[cell_id].crk_n[1]);
  M_reference.push_back(cell_SDdata[cell_id].crk_n[0]);

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
  dC_k1_Displacement.reinit(4*dim);
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
      if (ck==DOF_Diffuse_interface) dC_k1_Diffuse_interface[_i] = (ULocal[i].val()  - cell_SDdata[cell_id].ULocal_k[i]);
    }
    //std::cout << "--size-1--" <<  dC_k1_Lithium.size() << std::endl;
  }

  {
    int _i = -1;
    for (unsigned int i = 0; i < dofs_per_cell; ++i) {
      const unsigned int ck = fe_values.get_fe().system_to_component_index(i).first;
      if (ck>=DOF_Displacement and ck <DOF_Displacement+dim ) {
        _i += 1;
        dC_k1_Displacement[_i] = (ULocal[i].val()  - cell_SDdata[cell_id].ULocal_k[i]);
      }
    }
  }


    //std::cout <<  " Xi_old[Lithium] " << cell_SDdata[cell_id].xi_old(0) << std::endl;
    //std::cout <<  " Xi_old[Lithium_cation] " << cell_SDdata[cell_id].xi_old_c_e(0) << std::endl;
    //std::cout <<  " Xi_old[Electrode_potential] " << cell_SDdata[cell_id].xi_old_phi_s(0) << std::endl;
    //std::cout <<  " Xi_old[Electrolyte_potential] " << cell_SDdata[cell_id].xi_old_phi_e(0) << std::endl;

  //Vector<double> dxi_k1;
  Vector<double> dxi_k1_Lithium;
  Vector<double> dxi_k1_Lithium_phaseField;
  Vector<double> dxi_k1_Lithium_cation;
  Vector<double> dxi_k1_Electrode_potential;
  Vector<double> dxi_k1_Electrolyte_potential;
  Vector<double> dxi_k1_Displacement_sd;
  Vector<double> dxi_k1_Diffuse_interface;

  //dxi_k1.reinit(1);
  dxi_k1_Lithium.reinit(1);
  dxi_k1_Lithium_phaseField.reinit(1);
  dxi_k1_Lithium_cation.reinit(1);
  dxi_k1_Electrode_potential.reinit(1);
  dxi_k1_Electrolyte_potential.reinit(1);
  dxi_k1_Displacement_sd.reinit(3);
  dxi_k1_Diffuse_interface.reinit(1);

  //Table<1, Sacado::Fad::DFad<double>> xi_0(1);  // define sacado xi_0 for stiffness calculation.
  Table<1, Sacado::Fad::DFad<double>> xi_0_Lithium(1);
  Table<1, Sacado::Fad::DFad<double>> xi_0_Lithium_phaseField(1);
  Table<1, Sacado::Fad::DFad<double>> xi_0_Lithium_cation(1);
  Table<1, Sacado::Fad::DFad<double>> xi_0_Electrode_potential(1);
  Table<1, Sacado::Fad::DFad<double>> xi_0_Electrolyte_potential(1);
  Table<1, Sacado::Fad::DFad<double>> xi_0_Displacement_sd(3);
  Table<1, Sacado::Fad::DFad<double>> xi_0_Diffuse_interface(1);

  //if (primiary_dof != cell_SDdata[cell_id].opposite_flux_dof_li)
  //{ // for electrode
    cell_SDdata[cell_id].rlocal[0] = 0.0;
    cell_SDdata[cell_id].Kxic.vmult(dxi_k1_Lithium, dC_k1_Lithium);
    cell_SDdata[cell_id].rlocal -= dxi_k1_Lithium;
    //cell_SDdata[cell_id].rlocal[0] += cell_SDdata[cell_id].reaction_rate_li.val() * cell_SDdata[cell_id].interface_length;
    cell_SDdata[cell_id].Kxixi_inv.vmult(dxi_k1_Lithium, cell_SDdata[cell_id].rlocal);
    xi_0_Lithium[0] = cell_SDdata[cell_id].xi_old(0) + dxi_k1_Lithium(0);  
    cell_SDdata[cell_id].xi_old(0) = xi_0_Lithium[0].val();
    //std::cout << "--a0-0--"  << std::endl;
  //}
  //else
  //{// for electrolyte
    ////std::cout << "--a0-1--" << std::endl;
    cell_SDdata[cell_id].rlocal_c_e[0] = 0.0;
    cell_SDdata[cell_id].Kxic_c_e.vmult(dxi_k1_Lithium_cation, dC_k1_Lithium_cation);
    cell_SDdata[cell_id].rlocal_c_e -= dxi_k1_Lithium_cation;
    //cell_SDdata[cell_id].rlocal_c_e[0] += (-1 *cell_SDdata[cell_id].reaction_rate_li.val()) * cell_SDdata[cell_id].interface_length; // reaction rate li direction should not change
    cell_SDdata[cell_id].Kxixi_inv_c_e.vmult(dxi_k1_Lithium_cation, cell_SDdata[cell_id].rlocal_c_e);
    xi_0_Lithium_cation[0] = cell_SDdata[cell_id].xi_old_c_e(0) + dxi_k1_Lithium_cation(0);  
    cell_SDdata[cell_id].xi_old_c_e(0) = xi_0_Lithium_cation[0].val();
  //}


    cell_SDdata[cell_id].rlocal_phi_s[0] = 0.0;
    cell_SDdata[cell_id].Kxic_phi_s.vmult(dxi_k1_Electrode_potential, dC_k1_Electrode_potential);
    cell_SDdata[cell_id].rlocal_phi_s -= dxi_k1_Electrode_potential;
    //cell_SDdata[cell_id].rlocal_phi_s[0] += cell_SDdata[cell_id].reaction_rate_potential.val() * cell_SDdata[cell_id].interface_length;
    cell_SDdata[cell_id].Kxixi_inv_phi_s.vmult(dxi_k1_Electrode_potential, cell_SDdata[cell_id].rlocal_phi_s);
    xi_0_Electrode_potential[0] = cell_SDdata[cell_id].xi_old_phi_s(0) + dxi_k1_Electrode_potential(0);  
    cell_SDdata[cell_id].xi_old_phi_s(0) = xi_0_Electrode_potential[0].val();

    //std::cout << "xi_0_Electrode_potential[0] " << xi_0_Electrode_potential[0] << std::endl;

    //std::cout << "--a0-1--" << std::endl;
    cell_SDdata[cell_id].rlocal_phi_e[0] = 0.0;
    cell_SDdata[cell_id].Kxic_phi_e.vmult(dxi_k1_Electrolyte_potential, dC_k1_Electrolyte_potential);
    cell_SDdata[cell_id].rlocal_phi_e -= dxi_k1_Electrolyte_potential;
    //cell_SDdata[cell_id].rlocal_phi_e[0] += (-1 * cell_SDdata[cell_id].reaction_rate_potential.val()) * cell_SDdata[cell_id].interface_length; // reaction rate li direction should not change
    cell_SDdata[cell_id].Kxixi_inv_phi_e.vmult(dxi_k1_Electrolyte_potential, cell_SDdata[cell_id].rlocal_phi_e);
    xi_0_Electrolyte_potential[0] = cell_SDdata[cell_id].xi_old_phi_e(0) + dxi_k1_Electrolyte_potential(0);  
    cell_SDdata[cell_id].xi_old_phi_e(0) = xi_0_Electrolyte_potential[0].val();

    //std::cout << "xi_0_Electrolyte_potential[0] " << xi_0_Electrolyte_potential[0] << std::endl;
  //std::cout << " delta xi: Li " << dxi_k1_Lithium[0] << "\t Li_plus " << dxi_k1_Lithium_cation[0] << "\t phi_s " << dxi_k1_Electrode_potential[0] << "\t phi_e " << dxi_k1_Electrolyte_potential[0] << std::endl;
  //xi_0[0].diff(0, 1);

  dxi_k1_Displacement_sd = 0.0;
  for (int i=0; i<4*dim; ++i)
  {
    // Kxic * delta c^{k+1}
    dxi_k1_Displacement_sd[0] += cell_SDdata[cell_id].Kxiu_sd(0,i) * dC_k1_Displacement[i];
    dxi_k1_Displacement_sd[1] += cell_SDdata[cell_id].Kxiu_sd(1,i) * dC_k1_Displacement[i];
    dxi_k1_Displacement_sd[2] += cell_SDdata[cell_id].Kxiu_sd(2,i) * dC_k1_Displacement[i];
  }
  cell_SDdata[cell_id].rlocal_u_sd[0] -= dxi_k1_Displacement_sd[0] ;
  cell_SDdata[cell_id].rlocal_u_sd[1] -= dxi_k1_Displacement_sd[1] ;
  cell_SDdata[cell_id].rlocal_u_sd[2] -= dxi_k1_Displacement_sd[2] ;

  if (cell_SDdata[cell_id].is_fractured)
  {
    dxi_k1_Displacement_sd[0] = cell_SDdata[cell_id].Kxixi_inv_u_sd(0,0) * cell_SDdata[cell_id].rlocal_u_sd(0) + cell_SDdata[cell_id].Kxixi_inv_u_sd(0,1) * cell_SDdata[cell_id].rlocal_u_sd(1) + cell_SDdata[cell_id].Kxixi_inv_u_sd(0,2) * cell_SDdata[cell_id].rlocal_u_sd(2);
    dxi_k1_Displacement_sd[1] = cell_SDdata[cell_id].Kxixi_inv_u_sd(1,0) * cell_SDdata[cell_id].rlocal_u_sd(0) + cell_SDdata[cell_id].Kxixi_inv_u_sd(1,1) * cell_SDdata[cell_id].rlocal_u_sd(1) + cell_SDdata[cell_id].Kxixi_inv_u_sd(1,2) * cell_SDdata[cell_id].rlocal_u_sd(2);
    dxi_k1_Displacement_sd[2] = cell_SDdata[cell_id].Kxixi_inv_u_sd(2,0) * cell_SDdata[cell_id].rlocal_u_sd(0) + cell_SDdata[cell_id].Kxixi_inv_u_sd(2,1) * cell_SDdata[cell_id].rlocal_u_sd(1) + cell_SDdata[cell_id].Kxixi_inv_u_sd(2,2) * cell_SDdata[cell_id].rlocal_u_sd(2);
  }
  else
  {
    dxi_k1_Displacement_sd[0] = cell_SDdata[cell_id].Kxixi_inv_u_sd(0,0) * cell_SDdata[cell_id].rlocal_u_sd(0) + cell_SDdata[cell_id].Kxixi_inv_u_sd(0,1) * cell_SDdata[cell_id].rlocal_u_sd(1);
    dxi_k1_Displacement_sd[1] = cell_SDdata[cell_id].Kxixi_inv_u_sd(1,0) * cell_SDdata[cell_id].rlocal_u_sd(0) + cell_SDdata[cell_id].Kxixi_inv_u_sd(1,1) * cell_SDdata[cell_id].rlocal_u_sd(1);
    dxi_k1_Displacement_sd[2] = cell_SDdata[cell_id].Kxixi_inv_u_sd(2,2) * cell_SDdata[cell_id].rlocal_u_sd(2);
  }


  xi_0_Displacement_sd[0] = cell_SDdata[cell_id].xi_old_u_sd(0) + dxi_k1_Displacement_sd(0);  
  xi_0_Displacement_sd[1] = cell_SDdata[cell_id].xi_old_u_sd(1) + dxi_k1_Displacement_sd(1);  
  xi_0_Displacement_sd[2] = cell_SDdata[cell_id].xi_old_u_sd(2) + dxi_k1_Displacement_sd(2);  

  cell_SDdata[cell_id].xi_old_u_sd(0) = xi_0_Displacement_sd[0].val();
  cell_SDdata[cell_id].xi_old_u_sd(1) = xi_0_Displacement_sd[1].val();
  cell_SDdata[cell_id].xi_old_u_sd(2) = xi_0_Displacement_sd[2].val();

  if(disable_wd)
  {
    xi_0_Displacement_sd[2] = 0.0;
  }

  const unsigned int total_local_xi_dof = 7;
  Table<1, Sacado::Fad::DFad<double> > ULocal_xi(dofs_per_cell + total_local_xi_dof);
  for (unsigned int i = 0; i < dofs_per_cell; ++i) {
    ULocal_xi[i]=ULocal[i].val();
  }
  ULocal_xi[dofs_per_cell + ind_Lithium] = xi_0_Lithium[0].val();
  ULocal_xi[dofs_per_cell + ind_Lithium_cation] = xi_0_Lithium_cation[0].val();
  ULocal_xi[dofs_per_cell + ind_Electrode_potential] = xi_0_Electrode_potential[0].val();
  ULocal_xi[dofs_per_cell + ind_Electrolyte_potential] = xi_0_Electrolyte_potential[0].val();
  ULocal_xi[dofs_per_cell + ind_Displacement_sd_1] = xi_0_Displacement_sd[0].val();
  ULocal_xi[dofs_per_cell + ind_Displacement_sd_2] = xi_0_Displacement_sd[1].val();
  ULocal_xi[dofs_per_cell + ind_Displacement_wd] = xi_0_Displacement_sd[2].val();
  for (unsigned int i = 0; i < dofs_per_cell + total_local_xi_dof; ++i) {
    ULocal_xi[i].diff (i, dofs_per_cell + total_local_xi_dof);
  }


  // XZ: I think the defmap should belong to the following function call as well.
	battery_fields.update_fields(cell, fe_values, ULocal_xi, ULocalConv);
  //double D_1 = 1.0;
  ////evaluate primary fields
  /// \todo add new interface related methods for the following in the battery components to address the tilde part.
  dealii::Table<2,Sacado::Fad::DFad<double> > diffu_Lithium(n_q_points, dim);
  dealii::Table<1,Sacado::Fad::DFad<double> > react_Lithium(n_q_points);
  dealii::Table<2,Sacado::Fad::DFad<double> > diffu_Lithium_cation(n_q_points, dim);
  dealii::Table<1,Sacado::Fad::DFad<double> > react_Lithium_cation(n_q_points);

	dealii::Table<2,Sacado::Fad::DFad<double> > field_Electrode_potential(n_q_points,dim);
	dealii::Table<1,Sacado::Fad::DFad<double> > source_Electrode_potential(n_q_points);
	dealii::Table<2,Sacado::Fad::DFad<double> > field_Electrolyte_potential(n_q_points,dim);
	dealii::Table<1,Sacado::Fad::DFad<double> > source_Electrolyte_potential(n_q_points);


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

  Table<1, Sacado::Fad::DFad<double>> Ruu_Displacement_sd(4*dim);
  Table<1, Sacado::Fad::DFad<double>> Ruxi_Displacement_sd(4*dim);
  Table<1, Sacado::Fad::DFad<double>> rxiu_Displacement_sd(3);
  Table<1, Sacado::Fad::DFad<double>> rr_u(3);
  Table<1, Sacado::Fad::DFad<double>> rxixi_Displacement_sd(3);
  Table<1, Sacado::Fad::DFad<double>> rr_Displacement_sd(3);
  rxiu_Displacement_sd[0] = 0.0;
  rxiu_Displacement_sd[1] = 0.0;
  rxiu_Displacement_sd[2] = 0.0;
  rxixi_Displacement_sd[0] = 0.0;
  rxixi_Displacement_sd[1] = 0.0;
  rxixi_Displacement_sd[2] = 0.0;
  rr_Displacement_sd[0] = 0.0;
  rr_Displacement_sd[1] = 0.0;
  rr_Displacement_sd[2] = 0.0;
  rr_u[0] = 0.0;
  rr_u[1] = 0.0;
  rr_u[2] = 0.0;
  for (unsigned int i = 0; i < 4*dim; ++i) {
    Ruu_Displacement_sd[i] = 0.0;
    Ruxi_Displacement_sd[i] = 0.0;
  }

  double element_edge_length = 0.0;
  element_edge_length = cell->vertex(0).distance(cell->vertex(1));
  double m_plus = cell_SDdata[cell_id].computed_area;
  double m_minus = cell_SDdata[cell_id].area_elem - cell_SDdata[cell_id].computed_area;
  double l_minus = element_edge_length;
  double l_plus = - element_edge_length * m_plus / m_minus;


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

  dealii::Table<3, Sacado::Fad::DFad<double>> F_tilde_wd(n_q_points, dim, dim);
  dealii::Table<3, Sacado::Fad::DFad<double>> F_tilde_sd(n_q_points, dim, dim);
  dealii::Table<2, Sacado::Fad::DFad<double>> F_hat_sd(dim, dim);
  dealii::Table<2, Sacado::Fad::DFad<double>> F_hat_sd_wd(dim, dim);
  dealii::Table<2, Sacado::Fad::DFad<double>> F_const_wd(dim, dim);
  dealii::Table<2, Sacado::Fad::DFad<double>> T_at_gp(n_q_points, dim);
  dealii::Table<1, Sacado::Fad::DFad<double>> TN_at_gp(n_q_points);
  dealii::Table<1, Sacado::Fad::DFad<double>> TM_at_gp(n_q_points);

	deformationMap<Sacado::Fad::DFad<double>, dim> defMap(n_q_points); 
	getDeformationMap<Sacado::Fad::DFad<double>, dim>(fe_values, DOF_Displacement, ULocal_xi, defMap);

  for (unsigned int q = 0; q < n_q_points; ++q) {
    for (unsigned int j = 0; j < dim; j++) {
      for (unsigned int k = 0; k < dim; k++) {
        F_tilde_wd[q][j][k] = 0.0;
        F_tilde_sd[q][j][k] = 0.0;
        F_hat_sd[j][k] = 0.0;
        F_hat_sd_wd[j][k] = 0.0;
        F_const_wd[j][k] = 0.0;
        T_at_gp[q][j] = 0.0;
        TN_at_gp[q] = 0.0;
        TM_at_gp[q] = 0.0;
      }
    }
  }

  for (unsigned int q = 0; q < n_q_points; ++q) {
    for (unsigned int c = 0; c < dim; c++) {
      for (unsigned int d = 0; d < dim; d++) {
        F_hat_sd[c][d] += defMap.F[q][c][d] * fe_values.JxW(q);
        //std::cout << " F " << defMap.F[q][c][d] << std::endl;
      }
    }
  }
  for (unsigned int c = 0; c < dim; c++) {
    for (unsigned int d = 0; d < dim; d++) {
      F_hat_sd[c][d] = F_hat_sd[c][d] / cell_SDdata[cell_id].area_elem;
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
      //std::cout << "plus node: " <<  cell_SDdata[cell_id].lnode_plus[i] << std::endl;
    }

    for (unsigned i = 0; i < cell_SDdata[cell_id].lnode_minus.size(); ++i) {
      Ms_list_opposite[cell_SDdata[cell_id].lnode_minus[i]] = 1.0;
      //std::cout << "minus node: " <<  cell_SDdata[cell_id].lnode_minus[i] << std::endl;
    }


    double dummy_area = 0.0;
    double dummy_area_opposite = 0.0;
    std::vector<Sacado::Fad::DFad<double>> J_jump; 
    J_jump.resize(2);
    J_jump[0] = ULocal_xi[dofs_per_cell+ ind_Displacement_sd_1] * N_reference[0] +  ULocal_xi[dofs_per_cell+ ind_Displacement_sd_2] * M_reference[0];
    J_jump[1] = ULocal_xi[dofs_per_cell+ ind_Displacement_sd_1] * N_reference[1] +  ULocal_xi[dofs_per_cell+ ind_Displacement_sd_2] * M_reference[1];
    for (unsigned int q = 0; q < n_q_points; ++q) {
      //double Ms = 1.0;
      //if (fe_values.quadrature_point(q)[0] < 0.5){
        //Ms = 0.0;
      //}
      double Ms = Ms_list[q];
      dummy_area += Ms * fe_values.JxW(q); // the actually int area

      std::vector<double> N_tilde;
      N_tilde.resize(2);

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

      //std::cout <<  "C_Li new tilde " << battery_fields.quad_fields[DOF_Lithium].value[q].val() + c_1_tilde_Lithium[q].val() << std::endl;
      //std::cout <<  "Phi_s new tilde " << battery_fields.quad_fields[DOF_Electrode_potential].value[q].val() + Ms * ULocal_xi[dofs_per_cell+ ind_Electrode_potential].val() << std::endl;

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

      cell_SDdata[cell_id].C_Li_plus_new[q] = battery_fields.quad_fields[DOF_Lithium_cation].value[q].val() + c_1_tilde_Lithium_cation[q].val();

      //std::cout <<  "C_Li_plus_new tilde " << cell_SDdata[cell_id].C_Li_plus_new[q] << std::endl;
      //std::cout <<  "Phi_e new tilde " << battery_fields.quad_fields[DOF_Electrolyte_potential].value[q].val() + Ms * ULocal_xi[dofs_per_cell+ ind_Electrolyte_potential].val() << std::endl;
      //
      for (unsigned i = 0; i < cell_SDdata[cell_id].lnode_minus.size(); ++i) {
        int minus_node = this_dof_local_index_Lithium[cell_SDdata[cell_id].lnode_minus[i]]; // same as the Electrode_potential
        for (unsigned int j = 0; j < dim; j++) {
          N_tilde[j] += fe_values.shape_grad(minus_node, q)[j];
        }
      }
      {
        F_tilde_sd[q][0][0] = - J_jump[0] * N_tilde[0];  // x,1, x,2
        F_tilde_sd[q][0][1] = - J_jump[0] * N_tilde[1];  // x,1, x,2
        F_tilde_sd[q][1][0] = - J_jump[1] * N_tilde[0];  // y,1, y,2
        F_tilde_sd[q][1][1] = - J_jump[1] * N_tilde[1];  // y,1, y,2
      }
    }
	  dealii::Table<1,double > C_Li_plus_old(n_q_points);
    for (unsigned int q=0; q<n_q_points; ++q) C_Li_plus_old[q] = Ms_list_opposite[q] * cell_SDdata[cell_id].C_Li_plus_old[q];


    // for updating the reaction rate, will need to carefully revise this to use the interface quantity
    Sacado::Fad::DFad<double> c_li_ave=0.0, c_li_plus_ave=0.0, phi_s_ave=0.0, phi_e_ave=0.0;
    {
      int _count_plus = 0;
      for (unsigned i = 0; i < cell_SDdata[cell_id].lnode_plus.size(); ++i) {
        int plus_node_c_li = this_dof_local_index_Lithium[cell_SDdata[cell_id].lnode_plus[i]]; 
        int plus_node_phi_s = this_dof_local_index_Electrode_potential[cell_SDdata[cell_id].lnode_plus[i]];
        c_li_ave += ULocal_xi[plus_node_c_li];
        phi_s_ave += ULocal_xi[plus_node_phi_s];
        //std::cout << " i " << i << " c_li " << ULocal_xi[plus_node_c_li].val() << " phi_s " << ULocal_xi[plus_node_phi_s].val() << std::endl;
        _count_plus += 1;
      }
      c_li_ave = c_li_ave /_count_plus;
      phi_s_ave = phi_s_ave /_count_plus;

      int _count_minus = 0;
      for (unsigned i = 0; i < cell_SDdata[cell_id].lnode_minus.size(); ++i) {
        int minus_node_c_li_plus = this_dof_local_index_Lithium_cation[cell_SDdata[cell_id].lnode_minus[i]]; 
        int minus_node_phi_e = this_dof_local_index_Electrolyte_potential[cell_SDdata[cell_id].lnode_minus[i]];
        c_li_plus_ave += ULocal_xi[minus_node_c_li_plus];
        phi_e_ave += ULocal_xi[minus_node_phi_e];
        //std::cout << " i " << i << " c_li_plus " << ULocal_xi[minus_node_c_li_plus].val() << " phi_e " << ULocal_xi[minus_node_phi_e].val() << std::endl;
        _count_minus += 1;
      }
      c_li_plus_ave = c_li_plus_ave /_count_minus;
      phi_e_ave = phi_e_ave /_count_minus;
    }
    // for updating the reaction rate, use the interface quantity
    Sacado::Fad::DFad<double> c_li_tld=0.0, c_li_plus_tld=0.0, phi_s_tld=0.0, phi_e_tld=0.0;
    {
      int _count_plus = 0;
      for (unsigned i = 0; i < cell_SDdata[cell_id].lnode_plus.size(); ++i) {
        int q = cell_SDdata[cell_id].lnode_plus[i];
        int plus_node_c_li = this_dof_local_index_Lithium[cell_SDdata[cell_id].lnode_plus[i]]; 
        int plus_node_phi_s = this_dof_local_index_Electrode_potential[cell_SDdata[cell_id].lnode_plus[i]];

        c_li_tld += battery_fields.quad_fields[DOF_Lithium].value[q] + c_1_tilde_Lithium[q];
        phi_s_tld += battery_fields.quad_fields[DOF_Electrode_potential].value[q] + c_1_tilde_Electrode_potential[q];
        //std::cout << " i " << i << " c_li " << ULocal_xi[plus_node_c_li].val() << " phi_s " << ULocal_xi[plus_node_phi_s].val() << std::endl;
        _count_plus += 1;
      }
      c_li_tld = c_li_tld /_count_plus;
      phi_s_tld = phi_s_tld /_count_plus;

      int _count_minus = 0;
      for (unsigned i = 0; i < cell_SDdata[cell_id].lnode_minus.size(); ++i) {
        int q = cell_SDdata[cell_id].lnode_minus[i];
        int minus_node_c_li_plus = this_dof_local_index_Lithium_cation[cell_SDdata[cell_id].lnode_minus[i]]; 
        int minus_node_phi_e = this_dof_local_index_Electrolyte_potential[cell_SDdata[cell_id].lnode_minus[i]];
        c_li_plus_tld += battery_fields.quad_fields[DOF_Lithium_cation].value[q] + c_1_tilde_Lithium_cation[q];
        phi_e_tld += battery_fields.quad_fields[DOF_Electrolyte_potential].value[q] + c_1_tilde_Electrolyte_potential[q];
        //std::cout << " i " << i << " c_li_plus " << ULocal_xi[minus_node_c_li_plus].val() << " phi_e " << ULocal_xi[minus_node_phi_e].val() << std::endl;
        _count_minus += 1;
      }
      c_li_plus_tld = c_li_plus_tld /_count_minus;
      phi_e_tld = phi_e_tld /_count_minus;
    }

    //std::cout 
      //<< "\tc_li_ave      \t" << c_li_ave.val()
      //<< "\tphi_s_ave     \t" << phi_s_ave.val()
      //<< "\tc_li_plus_ave \t" << c_li_plus_ave.val()
      //<< "\tphi_e_ave     \t" << phi_e_ave.val()
      //<< std::endl;
    //std::cout 
      //<< "\tc_li_tld      \t" << c_li_tld.val()
      //<< "\tphi_s_tld     \t" << phi_s_tld.val()
      //<< "\tc_li_plus_tld \t" << c_li_plus_tld.val()
      //<< "\tphi_e_tld     \t" << phi_e_tld.val()
      //<< std::endl;

    Sacado::Fad::DFad<double> jn = 0.0; 
		//jn = electricChemoFormula.formula_jn(Temp, c_li_ave, c_li_plus_ave, phi_s_ave, phi_e_ave, domainflag);
		jn = electricChemoFormula.formula_jn(Temp, c_li_tld, c_li_plus_tld, phi_s_tld, phi_e_tld, domainflag);


   if (jn.val() != jn.val()) 
   {
    std::cout 
      << " jn " << jn.val()
      << " Temp " << Temp
      << " c_li_tld " << c_li_tld.val()
      << " c_li_plus_tld " << c_li_plus_tld.val()
      << " phi_s_tld " << phi_s_tld.val()
      << " phi_e_tld " << phi_e_tld.val()
      << " domainflag " << domainflag
      << std::endl;
   }

    //std::cout << " jn " << jn << std::endl;

	  //double cod_coef_C0=(*params_json)["ElectroChemo"]["cod_coef_C0"];
    //double u_sd_0_max = 0.1; // a value needs to be specified --> cod_max: crack opening displacement

    jn = jn * 1.0 / (1.0 + exp(cod_coef_C0 * (cell_SDdata[cell_id].xi_conv_u_sd[0] - 0.5* cod_max )));

    //jn = 0.00001;

    //if (1.0 / (1.0 + exp(150 * (cell_SDdata[cell_id].xi_conv_u_sd[0] - u_sd_0_max ))) < 0.95)
    //{
      //std::cout << " f(d) " << 1.0 / (1.0 + exp(150 * (cell_SDdata[cell_id].xi_conv_u_sd[0] - u_sd_0_max ))) << std::endl;
    //}
    //if (cell_SDdata[cell_id].xi_conv_u_sd[0] > 0) std::cout << "jump_n: " << cell_SDdata[cell_id].xi_conv_u_sd[0] << std::endl;
    cell_SDdata[cell_id].reaction_rate_potential = jn*F;
    cell_SDdata[cell_id].reaction_rate_li = jn;
    //std::cout << "*reaction rate (new) : " << cell_SDdata[cell_id].reaction_rate_li.val() << " potential (rate) " << cell_SDdata[cell_id].reaction_rate_potential.val() << std::endl;


  for (unsigned int q = 0; q < n_q_points; ++q) {
    for (unsigned int c = 0; c < dim; c++) {
      for (unsigned int d = 0; d < dim; d++) {
        F_hat_sd_wd[c][d] += (defMap.F[q][c][d] + F_tilde_sd[q][c][d])* fe_values.JxW(q);
      }
    }
  }

  for (unsigned int c = 0; c < dim; c++) {
    for (unsigned int d = 0; d < dim; d++) {
      F_hat_sd_wd[c][d] = F_hat_sd_wd[c][d] / cell_SDdata[cell_id].area_elem;
    }
  }

  for (unsigned int c = 0; c < dim; c++) {
    for (unsigned int d = 0; d < dim; d++) {
      F_const_wd[c][d] = N_reference[c] * N_reference[d] + M_reference[c] * M_reference[d] ; // won't work, not linear independent
    }
  }

  for (unsigned int q = 0; q < n_q_points; ++q) {
      for (unsigned int i = 0; i < dim; i++) {
        for (unsigned int j = 0; j < dim; j++) {
          // plus side
          if (Ms_list[q] > 0.0){
            F_tilde_wd[q][i][j] = F_hat_sd_wd[i][j] / l_plus;
 
          }
          else
          {
            F_tilde_wd[q][i][j] = F_hat_sd_wd[i][j] / l_minus;
          }
        }
      }
  }


    if (disable_wd)
    {
      for (unsigned int q = 0; q < n_q_points; ++q) {
        for (unsigned int i = 0; i < dim; i++) {
          for (unsigned int j = 0; j < dim; j++) {
            F_tilde_wd[q][i][j] = 0.0 * F_tilde_wd[q][i][j];
          }
        }
      }
    }




	  dealii::Table<2,Sacado::Fad::DFad<double> > field_Electrolyte_potential_tilde(n_q_points,dim);
	  dealii::Table<1,Sacado::Fad::DFad<double> > source_Electrolyte_potential_tilde(n_q_points);
    //phi_e.set_field_and_source_term_interface(field_Electrolyte_potential_tilde, source_Electrolyte_potential_tilde, c_1_tilde_grad_Electrolyte_potential, c_1_tilde_grad_Lithium_cation, c_1_tilde_Lithium_cation, battery_fields.quad_fields[DOF_Lithium_cation].value_conv);
    //phi_e.set_field_and_source_term_interface(field_Electrolyte_potential_tilde, source_Electrolyte_potential_tilde, battery_fields.quad_fields[DOF_Electrolyte_potential].value_grad, battery_fields.quad_fields[DOF_Lithium_cation].value_grad, battery_fields.quad_fields[DOF_Lithium_cation].value, battery_fields.quad_fields[DOF_Lithium_cation].value_conv);
    //std::cout <<  "---------------- homo field---------------" << std::endl;
    phi_e.set_field_and_source_term_interface(field_Electrolyte_potential, source_Electrolyte_potential, battery_fields.quad_fields[DOF_Electrolyte_potential].value_grad, battery_fields.quad_fields[DOF_Lithium_cation].value_grad, battery_fields.quad_fields[DOF_Lithium_cation].value, C_Li_plus_old);
    //std::cout <<  "---------------- tilde field---------------" << std::endl;
    phi_e.set_field_and_source_term_interface(field_Electrolyte_potential_tilde, source_Electrolyte_potential_tilde, c_1_tilde_grad_Electrolyte_potential, c_1_tilde_grad_Lithium_cation, c_1_tilde_Lithium_cation, C_Li_plus_old);

    dealii::Table<2,Sacado::Fad::DFad<double> > diffu_Lithium_cation_tilde(n_q_points, dim);
    dealii::Table<1,Sacado::Fad::DFad<double> > react_Lithium_cation_tilde(n_q_points);
    //lithium_cation.set_diffusion_reaction_term_interface(diffu_Lithium_cation_tilde, react_Lithium_cation_tilde, c_1_tilde_grad_Electrolyte_potential, c_1_tilde_grad_Lithium_cation, c_1_tilde_Lithium_cation, battery_fields.quad_fields[DOF_Lithium_cation].value_conv);
    lithium_cation.set_diffusion_reaction_term_interface(diffu_Lithium_cation_tilde, react_Lithium_cation_tilde, c_1_tilde_grad_Electrolyte_potential, c_1_tilde_grad_Lithium_cation, c_1_tilde_Lithium_cation, C_Li_plus_old, pressure_old[cell_id]);
    // carefully rewrite the code this part to make sure everything is correct.
    //lithium_cation.set_diffusion_reaction_term(diffu_Lithium_cation, react_Lithium_cation);
    lithium_cation.set_diffusion_reaction_term_interface(diffu_Lithium_cation, react_Lithium_cation, battery_fields.quad_fields[DOF_Electrolyte_potential].value_grad, battery_fields.quad_fields[DOF_Lithium_cation].value_grad, battery_fields.quad_fields[DOF_Lithium_cation].value, C_Li_plus_old, pressure_old[cell_id]);
    //phi_e.set_field_and_source_term(field_Electrolyte_potential, source_Electrolyte_potential);


	  dealii::Table<2,Sacado::Fad::DFad<double> > field_Electrode_potential_tilde(n_q_points,dim);
	  dealii::Table<1,Sacado::Fad::DFad<double> > source_Electrode_potential_tilde(n_q_points);
    phi_s.set_field_and_source_term_interface(field_Electrode_potential, source_Electrode_potential,battery_fields.quad_fields[DOF_Electrode_potential].value_grad);
    phi_s.set_field_and_source_term_interface(field_Electrode_potential_tilde, source_Electrode_potential_tilde, c_1_tilde_grad_Electrode_potential);


    dealii::Table<2,Sacado::Fad::DFad<double> > diffu_Lithium_tilde(n_q_points, dim);
    dealii::Table<1,Sacado::Fad::DFad<double> > react_Lithium_tilde(n_q_points);
    lithium.set_diffusion_reaction_term_interface(diffu_Lithium, react_Lithium, battery_fields.quad_fields[DOF_Lithium].value_grad, pressure_old[cell_id]);
    lithium.set_diffusion_reaction_term_interface(diffu_Lithium_tilde, react_Lithium_tilde, c_1_tilde_grad_Lithium, pressure_old[cell_id]);


    rr_Lithium[0] = - cell_SDdata[cell_id].reaction_rate_li * cell_SDdata[cell_id].interface_length;
    rr_Lithium_cation[0] = - ( -1.0 * cell_SDdata[cell_id].reaction_rate_li) * cell_SDdata[cell_id].interface_length;
    rr_Electrode_potential[0] = - cell_SDdata[cell_id].reaction_rate_potential * cell_SDdata[cell_id].interface_length;
    rr_Electrolyte_potential[0] = - (-1.0 * cell_SDdata[cell_id].reaction_rate_potential) * cell_SDdata[cell_id].interface_length;
    //std::cout << "rr_Electrode_potential[0] " << rr_Electrode_potential[0] << std::endl;
    //std::cout << "rr_Electrolyte_potential[0] " << rr_Electrolyte_potential[0] << std::endl;
    //std::cout << "rr_Lithium[0] " << rr_Lithium[0] << std::endl;

    ////std::cout << "--a2--" << rr[0] << std::endl;
    for (unsigned int q = 0; q < n_q_points; ++q) {
    //std::cout << "-------q-------" << q << std::endl;
    //std::cout << "diffu_lithium [0] " << diffu_Lithium[q][0] << " tilde " << diffu_Lithium_tilde[q][0] << std::endl;
    //std::cout << "diffu_lithium [1] " << diffu_Lithium[q][1] << " tilde " << diffu_Lithium_tilde[q][1] << std::endl;
    //std::cout << "grad_lithium [0]  " << battery_fields.quad_fields[DOF_Lithium].value_grad[q][0] << " tilde " << c_1_tilde_grad_Lithium[q][0] << std::endl;
    //std::cout << "grad_lithium [1]  " << battery_fields.quad_fields[DOF_Lithium].value_grad[q][1] << " tilde " << c_1_tilde_grad_Lithium[q][1] << std::endl;
    //std::cout << "diffu_Lithium_cation [0] " << diffu_Lithium_cation[q][0] << " tilde " << diffu_Lithium_cation_tilde[q][0] << std::endl;
    //std::cout << "diffu_Lithium_cation [1] " << diffu_Lithium_cation[q][1] << " tilde " << diffu_Lithium_cation_tilde[q][1] << std::endl;
    //std::cout << "field_Electrolyte_potential [0] " << field_Electrolyte_potential[q][0] << " tilde " << field_Electrolyte_potential_tilde[q][0] << std::endl;
    //std::cout << "field_Electrolyte_potential [1] " << field_Electrolyte_potential[q][1] << " tilde " << field_Electrolyte_potential_tilde[q][1] << std::endl;

    }

    //std::cout << " !!!! dummy_area  !!!"  << dummy_area << " computed_area "  <<  cell_SDdata[cell_id].computed_area << " dummy_opposite " << dummy_area_opposite << " elem_area " << cell_SDdata[cell_id].area_elem << std::endl; 

    for (unsigned int q = 0; q < n_q_points; ++q) {
          //rxixi_Lithium[0] += cell_SDdata[cell_id].interface_length / dummy_area * (cell_SDdata[cell_id].computed_area /dummy_area)  * ( diffu_Lithium_tilde[q][0] * cell_SDdata[cell_id].crk_n[0] +  diffu_Lithium_tilde[q][1] * cell_SDdata[cell_id].crk_n[1]) * fe_values.JxW(q); 
          rxixi_Lithium[0] += cell_SDdata[cell_id].interface_length / cell_SDdata[cell_id].area_elem  * ( diffu_Lithium_tilde[q][0] * cell_SDdata[cell_id].crk_n[0] +  diffu_Lithium_tilde[q][1] * cell_SDdata[cell_id].crk_n[1]) * fe_values.JxW(q); 
          //rxixi_Lithium_cation[0] += cell_SDdata[cell_id].interface_length / dummy_area_opposite * (cell_SDdata[cell_id].computed_area /dummy_area_opposite)  * ( - diffu_Lithium_cation_tilde[q][0] * cell_SDdata[cell_id].crk_n[0] -  diffu_Lithium_cation_tilde[q][1] * cell_SDdata[cell_id].crk_n[1]) * fe_values.JxW(q); // crk_n direction is reversed
          rxixi_Lithium_cation[0] += cell_SDdata[cell_id].interface_length / cell_SDdata[cell_id].area_elem  * ( - diffu_Lithium_cation_tilde[q][0] * cell_SDdata[cell_id].crk_n[0] -  diffu_Lithium_cation_tilde[q][1] * cell_SDdata[cell_id].crk_n[1]) * fe_values.JxW(q); // crk_n direction is reversed
          //rxixi_Electrode_potential[0] += cell_SDdata[cell_id].interface_length / dummy_area * (cell_SDdata[cell_id].computed_area /dummy_area)  * (  field_Electrode_potential_tilde[q][0] * cell_SDdata[cell_id].crk_n[0] +  field_Electrode_potential_tilde[q][1] * cell_SDdata[cell_id].crk_n[1]) * fe_values.JxW(q); // crk_n direction is reversed
          rxixi_Electrode_potential[0] += cell_SDdata[cell_id].interface_length / cell_SDdata[cell_id].area_elem * (  field_Electrode_potential_tilde[q][0] * cell_SDdata[cell_id].crk_n[0] +  field_Electrode_potential_tilde[q][1] * cell_SDdata[cell_id].crk_n[1]) * fe_values.JxW(q); // crk_n direction is reversed
          //rxixi_Electrolyte_potential[0] += cell_SDdata[cell_id].interface_length / dummy_area_opposite * (cell_SDdata[cell_id].computed_area /dummy_area_opposite)  * (- field_Electrolyte_potential_tilde[q][0] * cell_SDdata[cell_id].crk_n[0] - field_Electrolyte_potential_tilde[q][1] * cell_SDdata[cell_id].crk_n[1]) * fe_values.JxW(q); 
          rxixi_Electrolyte_potential[0] += cell_SDdata[cell_id].interface_length / cell_SDdata[cell_id].area_elem  * (- field_Electrolyte_potential_tilde[q][0] * cell_SDdata[cell_id].crk_n[0] - field_Electrolyte_potential_tilde[q][1] * cell_SDdata[cell_id].crk_n[1]) * fe_values.JxW(q); 
          //std::cout << " len " << cell_SDdata[cell_id].interface_length 
            //<< " area "  << cell_SDdata[cell_id].area_elem  
            //<< " tilde "  << (- field_Electrolyte_potential_tilde[q][0] * cell_SDdata[cell_id].crk_n[0] - field_Electrolyte_potential_tilde[q][1] * cell_SDdata[cell_id].crk_n[1])
            //<< " crk[0] " << cell_SDdata[cell_id].crk_n[0]
            //<< " crk[1] " << cell_SDdata[cell_id].crk_n[1]
            //<< " tilde[0] " << field_Electrolyte_potential_tilde[q][0]
            //<< " tilde[1] " << field_Electrolyte_potential_tilde[q][1]
            //<< " jwx " << fe_values.JxW(q) 
            //<< std::endl;
     }
    //std::cout << "rxixi_Electrode_potential[0] " << rxixi_Electrode_potential[0] << std::endl;
    //std::cout << "rxixi_Electrolyte_potential[0] " << rxixi_Electrolyte_potential[0] << std::endl;
    //std::cout << "rxixi_Lithium_cation[0] " << rxixi_Lithium_cation[0] << std::endl;
    //std::cout << "rxixi_Lithium[0] " << rxixi_Lithium[0] << std::endl;

    //std::cout << "--a2-1--" << rxixi[0]  << std::endl;

    for (unsigned int q = 0; q < n_q_points; ++q) {
        //rxic_Lithium[0] += cell_SDdata[cell_id].interface_length /  dummy_area * (cell_SDdata[cell_id].computed_area /dummy_area)  * (diffu_Lithium[q][0] * cell_SDdata[cell_id].crk_n[0] + diffu_Lithium[q][1] * cell_SDdata[cell_id].crk_n[1]) * fe_values.JxW(q); 
        rxic_Lithium[0] += cell_SDdata[cell_id].interface_length /  cell_SDdata[cell_id].area_elem * (diffu_Lithium[q][0] * cell_SDdata[cell_id].crk_n[0] + diffu_Lithium[q][1] * cell_SDdata[cell_id].crk_n[1]) * fe_values.JxW(q); 
        //rxic_Lithium_cation[0] += cell_SDdata[cell_id].interface_length /  dummy_area_opposite * (cell_SDdata[cell_id].computed_area /dummy_area_opposite)  * (- diffu_Lithium_cation[q][0] * cell_SDdata[cell_id].crk_n[0] - diffu_Lithium_cation[q][1] * cell_SDdata[cell_id].crk_n[1]) * fe_values.JxW(q); 
        rxic_Lithium_cation[0] += cell_SDdata[cell_id].interface_length /  cell_SDdata[cell_id].area_elem * (- diffu_Lithium_cation[q][0] * cell_SDdata[cell_id].crk_n[0] - diffu_Lithium_cation[q][1] * cell_SDdata[cell_id].crk_n[1]) * fe_values.JxW(q); 

        //rxic_Electrode_potential[0] += cell_SDdata[cell_id].interface_length /  dummy_area * (cell_SDdata[cell_id].computed_area /dummy_area)  * (field_Electrode_potential[q][0] * cell_SDdata[cell_id].crk_n[0] + field_Electrode_potential[q][1] * cell_SDdata[cell_id].crk_n[1]) * fe_values.JxW(q); 
        rxic_Electrode_potential[0] += cell_SDdata[cell_id].interface_length /  cell_SDdata[cell_id].area_elem * (field_Electrode_potential[q][0] * cell_SDdata[cell_id].crk_n[0] + field_Electrode_potential[q][1] * cell_SDdata[cell_id].crk_n[1]) * fe_values.JxW(q); 
        //rxic_Electrolyte_potential[0] += cell_SDdata[cell_id].interface_length /  dummy_area_opposite * (cell_SDdata[cell_id].computed_area /dummy_area_opposite)  * (- field_Electrolyte_potential[q][0] * cell_SDdata[cell_id].crk_n[0] - field_Electrolyte_potential[q][1] * cell_SDdata[cell_id].crk_n[1]) * fe_values.JxW(q); 
        rxic_Electrolyte_potential[0] += cell_SDdata[cell_id].interface_length /  cell_SDdata[cell_id].area_elem * (- field_Electrolyte_potential[q][0] * cell_SDdata[cell_id].crk_n[0] - field_Electrolyte_potential[q][1] * cell_SDdata[cell_id].crk_n[1]) * fe_values.JxW(q); 
    }

    //std::cout << "rxic_Electrode_potential[0] " << rxic_Electrode_potential[0] << std::endl;
    //std::cout << "rxic_Electrolyte_potential[0] " << rxic_Electrolyte_potential[0] << std::endl;
    //std::cout << "rxic_Lithium_cation[0] " << rxic_Lithium_cation[0] << std::endl;
    //std::cout << "rxic_Lithium[0] " << rxic_Lithium[0] << std::endl;

    rr_Lithium[0] += rxixi_Lithium[0] + rxic_Lithium[0];
    rr_Lithium_cation[0] += rxixi_Lithium_cation[0] + rxic_Lithium_cation[0];
    rr_Electrode_potential[0] += rxixi_Electrode_potential[0] + rxic_Electrode_potential[0];
    rr_Electrolyte_potential[0] += rxixi_Electrolyte_potential[0] + rxic_Electrolyte_potential[0];

    //std::cout << "rr_Electrode_potential[0] " << rr_Electrode_potential[0] << std::endl;
    //std::cout << "rr_Electrolyte_potential[0] " << rr_Electrolyte_potential[0] << std::endl;
    //std::cout << "rr_Lithium_cation[0] " << rr_Lithium_cation[0] << std::endl;
    //std::cout << "rr_Lithium[0] " << rr_Lithium[0] << std::endl;
    //std::cout << "--rxic a3-- [0]" << rxic_Lithium_cation[0] << std::endl;
    //std::cout << "--rr a3-- [0]" << rr_Lithium_cation[0] << std::endl;
    //
    //std::cout << "--rxic a3-- " << rxic_Electrode_potential[0] << std::endl;
    //std::cout << "--rr a3-- " << rr_Electrode_potential[0] << std::endl;



    for (unsigned int q = 0; q < n_q_points; ++q) {
      for (unsigned i = 0; i < cell_SDdata[cell_id].lnode_plus.size(); ++i) {
        int plus_node = this_dof_local_index_Lithium[cell_SDdata[cell_id].lnode_plus[i]];
          //std::cout << "--a2-3--" << i << " " << q << " " << plus_node << " " << cell_SDdata[cell_id].lnode_plus[i] << std::endl;

        Rcxi_Lithium[cell_SDdata[cell_id].lnode_plus[i]] +=  fe_values.shape_value(plus_node, q)*((c_1_tilde_Lithium[q]- c_1_tilde_conv_Lithium [q])/this->current_dt)*fe_values.JxW(q);
        //std::cout << "Lithium tilde "<< c_1_tilde_Lithium[q] << " Lithium converge tilde " << c_1_tilde_conv_Lithium [q] << std::endl;
        for (unsigned int j = 0; j < dim; j++) {
          /// \todo introduce the flux tilde based on the existing, not use D_1
          /// \todo check if - or + sign
          Rcxi_Lithium[cell_SDdata[cell_id].lnode_plus[i]] += - diffu_Lithium_tilde[q][j] * fe_values.shape_grad(plus_node, q)[j] * fe_values.JxW(q);  
          //std::cout << "--a2-4--" << plus_node << " " << cell_SDdata[cell_id].lnode_plus[i]<< std::endl;
        }
      }
    }

    for (unsigned int q = 0; q < n_q_points; ++q) {
      for (unsigned i = 0; i < cell_SDdata[cell_id].lnode_minus.size(); ++i) {
        int minus_node = this_dof_local_index_Lithium_cation[cell_SDdata[cell_id].lnode_minus[i]];
          //std::cout << "--a2-3--" << i << " " << q << " " << minus_node << " " << cell_SDdata[cell_id].lnode_minus[i] << std::endl;
        Rcxi_Lithium_cation[cell_SDdata[cell_id].lnode_minus[i]] +=  fe_values.shape_value(minus_node, q)*((c_1_tilde_Lithium_cation[q]- c_1_tilde_conv_Lithium_cation[q])/this->current_dt)*fe_values.JxW(q);
        //std::cout << "Lithium cation tilde "<< c_1_tilde_Lithium_cation[q] << " Lithium cation converge tilde " << c_1_tilde_conv_Lithium_cation[q] << std::endl;
        for (unsigned int j = 0; j < dim; j++) {
          Rcxi_Lithium_cation[cell_SDdata[cell_id].lnode_minus[i]] += - diffu_Lithium_cation_tilde[q][j] * fe_values.shape_grad(minus_node, q)[j] * fe_values.JxW(q);  
          //std::cout << "--a2-4--" << minus_node << " " << cell_SDdata[cell_id].lnode_minus[i]<< std::endl;
        }
      }
    }


    for (unsigned int q = 0; q < n_q_points; ++q) {
      for (unsigned i = 0; i < cell_SDdata[cell_id].lnode_plus.size(); ++i) {
        int plus_node = this_dof_local_index_Electrode_potential[cell_SDdata[cell_id].lnode_plus[i]];
          //std::cout << "--a2-3--" << i << " " << q << " " << plus_node << " " << cell_SDdata[cell_id].lnode_plus[i] << std::endl;
        for (unsigned int j = 0; j < dim; j++) {
          Rcxi_Electrode_potential[cell_SDdata[cell_id].lnode_plus[i]] += - field_Electrode_potential_tilde[q][j] * fe_values.shape_grad(plus_node, q)[j] * fe_values.JxW(q);  
          //std::cout << "--a2-4--" << plus_node << " " << cell_SDdata[cell_id].lnode_plus[i]<< std::endl;
          //std::cout << "--a2-4-field: " << q << " " << field_Electrode_potential_tilde[q][j] << " tidle " << c_1_tilde_grad_Electrode_potential[q][j]  << std::endl;
      
        }
      }
    }

    //std::cout << "Rcxi_Electrode_potential[cell_SDdata[cell_id].lnode_plus[0]] " << Rcxi_Electrode_potential[cell_SDdata[cell_id].lnode_plus[0]] << std::endl;
    //std::cout << "Rcxi_Electrode_potential[cell_SDdata[cell_id].lnode_plus[1]] " << Rcxi_Electrode_potential[cell_SDdata[cell_id].lnode_plus[1]] << std::endl;

    for (unsigned int q = 0; q < n_q_points; ++q) {
      for (unsigned i = 0; i < cell_SDdata[cell_id].lnode_minus.size(); ++i) {
        int minus_node = this_dof_local_index_Electrolyte_potential[cell_SDdata[cell_id].lnode_minus[i]];
          //std::cout << "--a2-3--" << i << " " << q << " " << minus_node << " " << cell_SDdata[cell_id].lnode_minus[i] << std::endl;
        for (unsigned int j = 0; j < dim; j++) {
          Rcxi_Electrolyte_potential[cell_SDdata[cell_id].lnode_minus[i]] += - field_Electrolyte_potential_tilde[q][j] * fe_values.shape_grad(minus_node, q)[j] * fe_values.JxW(q);  
          //std::cout << "--a2-4--" << minus_node << " " << cell_SDdata[cell_id].lnode_minus[i]<< std::endl;
          //std::cout << "--a2-4--" << minus_node << " q " << q << " j " << j << " " << field_Electrolyte_potential_tilde[q][j]<< std::endl;
        }
      }
    }
    //std::cout << "Rcxi_Electrolyte_potential[cell_SDdata[cell_id].lnode_minus[0]] " << Rcxi_Electrolyte_potential[cell_SDdata[cell_id].lnode_minus[0]] << std::endl;
    //std::cout << "Rcxi_Electrolyte_potential[cell_SDdata[cell_id].lnode_minus[1]] " << Rcxi_Electrolyte_potential[cell_SDdata[cell_id].lnode_minus[1]] << std::endl;
    //std::cout << "Rcxi_Lithium_cation[cell_SDdata[cell_id].lnode_minus[0]] " << Rcxi_Lithium_cation[cell_SDdata[cell_id].lnode_minus[0]] << std::endl;
    //std::cout << "Rcxi_Lithium_cation[cell_SDdata[cell_id].lnode_minus[1]] " << Rcxi_Lithium_cation[cell_SDdata[cell_id].lnode_minus[1]] << std::endl;
    //std::cout << "Rcxi_Lithium[cell_SDdata[cell_id].lnode_plus[0]] " << Rcxi_Lithium[cell_SDdata[cell_id].lnode_plus[0]] << std::endl;
    //std::cout << "Rcxi_Lithium[cell_SDdata[cell_id].lnode_plus[1]] " << Rcxi_Lithium[cell_SDdata[cell_id].lnode_plus[1]] << std::endl;

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

	  dealii::Table<3, Sacado::Fad::DFad<double> > P(n_q_points,dim,dim);
    {
      double mu = displacement.ResidualEq->mu; 
      double lambda = displacement.ResidualEq->lambda; 

      //double E10 = this->MechMatData["Young's modulus"] * this->MatData["RightMatScale"];
      //double nu = this->MechMatData["Poisson's ratio"];
      //double mu10 = E10/2.0/(1.0+nu);
      //double lambda10 = (E10*nu)/(1.0+nu)/(1.0-2.0*nu);

		  dealii::Table<3,Sacado::Fad::DFad<double> > Fe(n_q_points,dim,dim);

      double swell_ratio = (*params_json)["Mechanics"]["SwellRatio"];
      double C_li_100_neg=(*params_json)["ElectroChemo"]["c_li_100_neg"];
      double C_li_100_pos=(*params_json)["ElectroChemo"]["c_li_100_pos"];
      double C_li_max_neg=(*params_json)["ElectroChemo"]["c_li_max_neg"];
      double C_li_max_pos=(*params_json)["ElectroChemo"]["c_li_max_pos"];

      double C_0 = 0.0;
		  C_0=C_li_100_neg * C_li_max_pos;
		  Point<dim> center=(* battery_fields.current_cell)->center();
		  if (center[orientation]>separator_line){
		  	C_0=C_li_100_pos * C_li_max_pos;
		  }

      for (unsigned int q = 0; q < n_q_points; ++q) {
        Sacado::Fad::DFad<double> c_li_tld = 0.0;
        c_li_tld = battery_fields.quad_fields[DOF_Lithium].value[q] + c_1_tilde_Lithium[q];
        //std::cout << " c_li_tld  = " << c_li_tld.val() << " " << Ms_list[q] << " " << Ms_list_opposite[q] << std::endl;

				dealii::Table<2,Sacado::Fad::DFad<double> > Feig(dim,dim);
				dealii::Table<2,Sacado::Fad::DFad<double>> invFeig(dim,dim);
        Feig[0][0]=(c_li_tld - C_0) * swell_ratio * Ms_list[q] + 1.0; 
        if(dim>=2){
          Feig[1][1]=(c_li_tld - C_0) * swell_ratio * Ms_list[q] + 1.0; 
        }
        if(dim==3){
          Feig[2][2]=(c_li_tld - C_0) * swell_ratio * Ms_list[q] + 1.0; 
        }

				getInverse<Sacado::Fad::DFad<double>,dim> (Feig,invFeig);
        // Fe = F inv(Fg)
        for (unsigned int i = 0; i < dim; ++i) {
          for (unsigned int j = 0; j < dim; ++j) {
            for (unsigned int k = 0; k < dim; ++k) {
              Fe[q][i][j] += (defMap.F[q][i][k] + F_tilde_sd[q][i][k] + ULocal_xi[dofs_per_cell+ ind_Displacement_wd] * F_tilde_wd[q][i][k]) * invFeig[k][j];
            }
            //Fe[q][i][j] = (defMap.F[q][i][j] + F_tilde_sd[q][i][j] + ULocal_xi[dofs_per_cell+ ind_Displacement_wd] * F_tilde_wd[q][i][j]); //  no swelling
            //Fe[q][i][j] = defMap.F[q][i][j]; //  no swelling
          }
        }

		    Table<2, Sacado::Fad::DFad<double> > C (dim, dim);
		    Table<2, Sacado::Fad::DFad<double> > S(dim, dim); 
		    for (unsigned int i=0; i<dim; ++i){
		    	 for (unsigned int j=0; j<dim; ++j){
		    		 C[i][j] = 0.0;
		    		 for (unsigned int k=0; k<dim; ++k){
		    			 C[i][j] += Fe[q][k][i]*Fe[q][k][j];
		    		 }
		    	 }
		    }
    //std::cout << " residual C at the interface " << std::endl;
		     
  	    Sacado::Fad::DFad<double> detC = 0.0;
  	    Table<2, Sacado::Fad::DFad<double>> invC(dim, dim); 
		    
  	    for (unsigned int i = 0; i < dim; ++i){
        	for (unsigned int j = 0; j < dim; ++j)
          {
		    		invC[i][j] = 0.0;
		    	}
		    }
  	    getInverse<Sacado::Fad::DFad<double>, dim>(C, invC, detC);
  	    for (unsigned int i = 0; i < dim; ++i){
		    	for (unsigned int j = 0; j < dim; ++j){
		    		S[i][j] += 0.5*lambda*detC*invC[i][j] - (0.5*lambda + mu)*invC[i][j] + mu*(i==j);
            //std::cout << " S " << i << j << " " << S[i][j] << std::endl;
        	}
 		    }
    //std::cout << " residual S at the interface " << std::endl;

        // P
        for (unsigned int i = 0; i < dim; ++i) {
          for (unsigned int j = 0; j < dim; ++j) {
            P[q][i][j] = 0.0;
            for (unsigned int k = 0; k < dim; ++k) {
              // option: F or with F_tidle, which one to use
              P[q][i][j] += (defMap.F[q][i][k] + F_tilde_sd[q][i][k] + ULocal_xi[dofs_per_cell+ ind_Displacement_wd] * F_tilde_wd[q][i][k] ) * S[k][j];
              //P[q][i][j] += defMap.F[q][i][k] * S[k][j];

            }
            //std::cout << " P " << q << i << j << " " << P[q][i][j] << std::endl;
          }
        }
      } // q
    } // code block {}

	dealii::Table<2, Sacado::Fad::DFad<double> > sigma(dim,dim);
  pressure[cell_id].resize(n_q_points);
  for (unsigned int q = 0; q < n_q_points; ++q){
		for(unsigned int i=0;i<dim;i++){
			for (unsigned int j=0; j<dim; ++j){
        sigma[i][j] = 0.0;
				for (unsigned int k=0; k<dim; ++k){
          //sigma[i][j]+=P[q][i][k]*defMap.F[q][j][k]/defMap.detF[q]; // may need to use the next line
          sigma[i][j]+=P[q][i][k]*(defMap.F[q][j][k] + F_tilde_sd[q][j][k] + ULocal_xi[dofs_per_cell+ ind_Displacement_wd] * F_tilde_wd[q][j][k])/defMap.detF[q]; // 
				} // k
			} // j
		} //i
    double _pressure = 0.0;
		for(unsigned int i=0;i<dim;i++){
      _pressure += 1.0/3.0 * sigma[i][i].val();
    }
    pressure[cell_id][q] = _pressure;
    //std::cout << "pressure: " << _pressure << " q " << q << std::endl;
  }

    //std::cout << " residual P at the interface " << std::endl;

  {
    for (unsigned int q = 0; q < n_q_points; ++q) {
      int _i = -1;
      for (unsigned i = 0; i < dofs_per_cell; ++i) {
        const unsigned int ck = fe_values.get_fe().system_to_component_index(i).first;
        if ( ck >= DOF_Displacement and ck < DOF_Displacement + dim){
          unsigned int c = ck - DOF_Displacement;
          _i += 1;
          for (unsigned int d = 0; d < dim; d++) Ruu_Displacement_sd[_i] += P[q][c][d] * fe_values.shape_grad_component(i, q, ck)[d] * fe_values.JxW(q);
        //std::cout << " residual R_uu at the interface " << _i 
          //<< " i " << i 
          //<< " ck " << ck 
          //<< " c " << c 
          //<< " "  << Ruu_Displacement_sd[_i]
          //<< " R0 "  << P[q][c][0] * fe_values.shape_grad_component(i, q, c)[0] * fe_values.JxW(q)
          //<< " R1 "  << P[q][c][1] * fe_values.shape_grad_component(i, q, c)[1] * fe_values.JxW(q)
          //<< " P[0] " << P[q][c][0] 
          //<< " P[1] " << P[q][c][1] 
          //<< fe_values.shape_grad_component(i, q, ck)[0]
          //<< fe_values.shape_grad_component(i, q, ck)[1]
          //<< std::endl;
        }
      }
    }
  }
    //std::cout << " residual R_uu at the interface " << std::endl;

  //N = [cell_SDdata[cell_id].crk_n[0], cell_SDdata[cell_id].crk_n[1]]
  //M = [-cell_SDdata[cell_id].crk_n[1], cell_SDdata[cell_id].crk_n[0]]
  //n = N
  //m = M
  dealii::Table<1,Sacado::Fad::DFad<double> > n_spatial(dim);
  dealii::Table<1,Sacado::Fad::DFad<double> > m_spatial(dim);
  dealii::Table<1,Sacado::Fad::DFad<double> > T_gamma(dim);

  for (unsigned int q = 0; q < n_q_points; ++q) {
    for (unsigned int i = 0; i < dim; i++) {
      for (unsigned int j = 0; j < dim; j++) {
        T_at_gp[q][i] += P[q][i][j] * N_reference[j];
      }
    }
  }
  for (unsigned int q = 0; q < n_q_points; ++q) {
    for (unsigned int i = 0; i < dim; i++) {
      TN_at_gp[q] += T_at_gp[q][i] * N_reference[i];
      TM_at_gp[q] += T_at_gp[q][i] * M_reference[i];
    }
  }




  for (unsigned int i = 0; i < dim; i++) {
    n_spatial[i] = 0.0;
    m_spatial[i] = 0.0;
    T_gamma[i] = 0.0;
  }
  for (unsigned int i = 0; i < dim; i++) {
      for (unsigned int j = 0; j < dim; j++) {
        n_spatial[i] += F_hat_sd[i][j] * N_reference[j];
        m_spatial[i] += F_hat_sd[i][j] * M_reference[j];
      }
  }

  // should be updated to use the correct spatial form
  for (unsigned int i = 0; i < dim; i++) {
    //n_spatial[i] = n_spatial[i]/sqrt(n_spatial[0]*n_spatial[0] + n_spatial[1]*n_spatial[1]);
    //m_spatial[i] = m_spatial[i]/sqrt(m_spatial[0]*m_spatial[0] + m_spatial[1]*m_spatial[1]);
    n_spatial[i] = N_reference[i];
    m_spatial[i] = M_reference[i];
  }


  double h_e = cell_SDdata[cell_id].area_elem / cell_SDdata[cell_id].interface_length;
  for (unsigned int q = 0; q < n_q_points; ++q) {
    for (unsigned int i = 0; i < dim; i++) {
      for (unsigned int j = 0; j < dim; j++) {
        if (cell_SDdata[cell_id].is_fractured)
        {
          //- (-1/h_e) = 1/h_e
          rr_u[0] -= (1.0/h_e) * n_spatial[i] * N_reference[j] * P[q][i][j] * fe_values.JxW(q); 
          rr_u[1] -= (1.0/h_e) * m_spatial[i] * N_reference[j] * P[q][i][j] * fe_values.JxW(q); 
        }
      }
    }
  }

  for (unsigned int q = 0; q < n_q_points; ++q) {
    for (unsigned int i = 0; i < dim; i++) {
      for (unsigned int j = 0; j < dim; j++) {
        rr_u[2] += F_tilde_wd[q][i][j] * P[q][i][j] * fe_values.JxW(q); 
      }
    }
  }


  // compute T_gamma

  if (cell_SDdata[cell_id].is_fractured)
  {
    //T_gamma[0] = this->MatData["MaximumStressN"] + this->MatData["LinearSoftenModulusN"] * ULocal_xi[dofs_per_cell + ind_Displacement_sd_1]; // this is wrong, normal direction is not accounted for!
    //T_gamma[0] = this->MatData["MaximumStressN"] + this->MatData["LinearSoftenModulusN"] * J_jump[0];
    //T_gamma[0] = this->MatData["MaximumStressN"] + this->MatData["LinearSoftenModulusN"] * ULocal_xi[dofs_per_cell + ind_Displacement_sd_1]; // this is actually correct, jump is always non-negative
    //T_gamma[0] = cell_SDdata[cell_id].max_Tn + this->MatData["LinearSoftenModulusN"] * ULocal_xi[dofs_per_cell + ind_Displacement_sd_1]; // this is actually correct, jump is always non-negative

    if (ULocal_xi[dofs_per_cell + ind_Displacement_sd_1] < 0.0)
    {
      // regularized stiffness
      T_gamma[0] = (cell_SDdata[cell_id].max_Tn + 1e11 * ULocal_xi[dofs_per_cell + ind_Displacement_sd_1]) * cell_SDdata[cell_id].interface_length; 
      this->T_n[cell_id] = (cell_SDdata[cell_id].max_Tn + 1e11 * ULocal_xi[dofs_per_cell + ind_Displacement_sd_1]).val();
    }
    else
    {
      T_gamma[0] = (cell_SDdata[cell_id].max_Tn + (*displacement.params_json)["Mechanics"]["LinearSoftenModulusN"] * ULocal_xi[dofs_per_cell + ind_Displacement_sd_1] ) * cell_SDdata[cell_id].interface_length; // this is actually correct, jump is always non-negative
      this->T_n[cell_id] = (cell_SDdata[cell_id].max_Tn + (*displacement.params_json)["Mechanics"]["LinearSoftenModulusN"] * ULocal_xi[dofs_per_cell + ind_Displacement_sd_1] ).val();
      if (T_gamma[0] < 0.0)
      {
        //std::cout << " T_gamma[0] " << T_gamma[0] << std::endl;
        T_gamma[0] = 1.0e-6 * T_gamma[0];
        this->T_n[cell_id] = 0.0;
      }
    }
    //std::cout << " xi " << ULocal_xi[dofs_per_cell + ind_Displacement_sd_1].val() 
      //<< " T_gamma[0] " << T_gamma[0].val() 
      //<< " P000 " << P[0][0][0].val() 
      //<< " P011 " << P[0][1][1].val() 
      //<< " P100 " << P[1][0][0].val() 
      //<< " P111 " << P[1][1][1].val() 
      //<< " P200 " << P[2][0][0].val() 
      //<< " P211 " << P[2][1][1].val() 
      //<< " P300 " << P[3][0][0].val() 
      //<< " P311 " << P[3][1][1].val()
      //<< std::endl;



    //T_gamma[1] = this->MatData["MaximumStressM"] + this->MatData["LinearSoftenModulusM"] * J_jump[1];
    T_gamma[1] = 0.0 * ULocal_xi[dofs_per_cell + ind_Displacement_sd_2];
    if (T_gamma[1] < 0.0)
    {
      T_gamma[1] = 0.0 * T_gamma[1];
    }

    // 1 gauss point integration for constants over the crack surface, weight is the length of the element
    {
      //for (unsigned int i = 0; i < dim; i++) {
        //rr_u[0] -= N_reference[i] * T_gamma[i] * cell_SDdata[cell_id].interface_length;
        //rr_u[1] -= M_reference[i] * T_gamma[i] * cell_SDdata[cell_id].interface_length;
      //}
      // for constant jump mode, this could be [1 0; 0 1], not [N; M].
      rr_u[0] += T_gamma[0] * cell_SDdata[cell_id].interface_length;
      rr_u[1] += T_gamma[1] * cell_SDdata[cell_id].interface_length;
    }
  }

  FullMatrix<double> Kxixi_Displacement_sd;
  FullMatrix<double> Kxixi_Displacement_sd_inv;
  Kxixi_Displacement_sd.reinit(3, 3);
  Kxixi_Displacement_sd_inv.reinit(3, 3);
  FullMatrix<double> Kxiu_Displacement_sd;
  Kxiu_Displacement_sd.reinit(3, 4*dim);
  FullMatrix<double> Kuxi_Displacement_sd;
  Kuxi_Displacement_sd.reinit(4*dim, 3);

  if (cell_SDdata[cell_id].is_fractured)
  {
    Kxixi_Displacement_sd(0, 0) = -rr_u[0].dx(dofs_per_cell+ ind_Displacement_sd_1);
    Kxixi_Displacement_sd(0, 1) = -rr_u[0].dx(dofs_per_cell+ ind_Displacement_sd_2);
    Kxixi_Displacement_sd(0, 2) = -rr_u[0].dx(dofs_per_cell+ ind_Displacement_wd);
    Kxixi_Displacement_sd(1, 0) = -rr_u[1].dx(dofs_per_cell+ ind_Displacement_sd_1);
    Kxixi_Displacement_sd(1, 1) = -rr_u[1].dx(dofs_per_cell+ ind_Displacement_sd_2);
    Kxixi_Displacement_sd(1, 2) = -rr_u[1].dx(dofs_per_cell+ ind_Displacement_wd);
    Kxixi_Displacement_sd(2, 0) = -rr_u[2].dx(dofs_per_cell+ ind_Displacement_sd_1);
    Kxixi_Displacement_sd(2, 1) = -rr_u[2].dx(dofs_per_cell+ ind_Displacement_sd_2);
    Kxixi_Displacement_sd(2, 2) = -rr_u[2].dx(dofs_per_cell+ ind_Displacement_wd);

    // because rr_u[1] is not enabled
    Kxixi_Displacement_sd(1, 1) =  Kxixi_Displacement_sd(0, 0);

    if (disable_wd)
    {
      Kxixi_Displacement_sd(2, 2) = 1.0;
    }

    Kxixi_Displacement_sd_inv.invert(Kxixi_Displacement_sd);
  }
  else
  {
    Kxixi_Displacement_sd(0, 0) = 1;
    Kxixi_Displacement_sd(0, 1) = 0;
    Kxixi_Displacement_sd(0, 2) = 0;
    Kxixi_Displacement_sd(1, 0) = 0;
    Kxixi_Displacement_sd(1, 1) = 1;
    Kxixi_Displacement_sd(1, 2) = 0;
    Kxixi_Displacement_sd(2, 0) = 0;
    Kxixi_Displacement_sd(2, 1) = 0;
    Kxixi_Displacement_sd(2, 2) = -rr_u[2].dx(dofs_per_cell+ ind_Displacement_wd);
    if (disable_wd)
    {
      Kxixi_Displacement_sd(2, 2) = 1.0;
    }

    Kxixi_Displacement_sd_inv.invert(Kxixi_Displacement_sd);
  }

  if (cell_SDdata[cell_id].is_fractured)
  {
    unsigned int _i = 0;
    for (unsigned int i=0; i< this_dof_local_index_Displacement.size(); i++)
    {
      Kxiu_Displacement_sd(0, _i) = - rr_u[0].dx(this_dof_local_index_Displacement[i]);
      Kxiu_Displacement_sd(1, _i) = - rr_u[1].dx(this_dof_local_index_Displacement[i]);
      Kxiu_Displacement_sd(2, _i) = - rr_u[2].dx(this_dof_local_index_Displacement[i]);
      Kuxi_Displacement_sd(_i, 0) = - Ruu_Displacement_sd[_i].dx(dofs_per_cell + ind_Displacement_sd_1);
      Kuxi_Displacement_sd(_i, 1) = - Ruu_Displacement_sd[_i].dx(dofs_per_cell + ind_Displacement_sd_2);
      Kuxi_Displacement_sd(_i, 2) = - Ruu_Displacement_sd[_i].dx(dofs_per_cell + ind_Displacement_wd);
      _i += 1;
    }
  }

    //std::cout << " residual K_u at the interface " << std::endl;

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
    //std::cout << " rxixi_Electrolyte_potential[0] = " <<  rxixi_Electrolyte_potential[0] << std::endl;

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

    //std::cout << "Rcc_Electrode_potential[cell_SDdata[cell_id].lnode_plus[0]] " << Rcc_Electrode_potential[cell_SDdata[cell_id].lnode_plus[0]] << std::endl;
    //std::cout << "Rcc_Electrode_potential[cell_SDdata[cell_id].lnode_plus[1]] " << Rcc_Electrode_potential[cell_SDdata[cell_id].lnode_plus[1]] << std::endl;

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
    //std::cout << "Rcc_Electrolyte_potential[cell_SDdata[cell_id].lnode_minus[0]] " << Rcc_Electrolyte_potential[cell_SDdata[cell_id].lnode_minus[0]] << std::endl;
    //std::cout << "Rcc_Electrolyte_potential[cell_SDdata[cell_id].lnode_minus[1]] " << Rcc_Electrolyte_potential[cell_SDdata[cell_id].lnode_minus[1]] << std::endl;
    //std::cout << "Rcc_Lithium_cation[cell_SDdata[cell_id].lnode_minus[0]] " << Rcc_Lithium_cation[cell_SDdata[cell_id].lnode_minus[0]] << std::endl;
    //std::cout << "Rcc_Lithium_cation[cell_SDdata[cell_id].lnode_minus[1]] " << Rcc_Lithium_cation[cell_SDdata[cell_id].lnode_minus[1]] << std::endl;
    //std::cout << "Rcc_Lithium[cell_SDdata[cell_id].lnode_plus[0]] " << Rcc_Lithium[cell_SDdata[cell_id].lnode_plus[0]] << std::endl;
    //std::cout << "Rcc_Lithium[cell_SDdata[cell_id].lnode_plus[1]] " << Rcc_Lithium[cell_SDdata[cell_id].lnode_plus[1]] << std::endl;

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
    int _ii = -1;
    for (unsigned int i = 0; i < dofs_per_cell; ++i) {
      const unsigned int ck = fe_values.get_fe().system_to_component_index(i).first;
      if (ck == 0) _i++;
      if (ck == DOF_Lithium ) {
        //std::cout << "R[i] " << R[i] << std::endl;
        //R[i] = (Rcc_Lithium[_i] + Rcxi_Lithium[_i] - Kcxi_Lithium(_i,0) / Kxixi_Lithium(0,0) * rr_Lithium[0]) * (cell_SDdata[cell_id].computed_area /dummy_area) ; 
        R[i] = (Rcc_Lithium[_i] + Rcxi_Lithium[_i] - Kcxi_Lithium(_i,0) / Kxixi_Lithium(0,0) * rr_Lithium[0])  ; 
        //std::cout << "R[i] (Lithium) " << R[i] << std::endl;
        //std::cout << "Rcc " << Rcc_Lithium[_i] << " " << Rcxi_Lithium[_i] << " " << Kcxi_Lithium(_i,0) << " Kxixi " << Kxixi_Lithium(0,0) << " rr " << rr_Lithium[0] << std::endl;
        //std::cout << "--a4-1-- " << i <<" " << R[i] << cell_SDdata[cell_id].computed_area /dummy_area<< std::endl;
      }
      if (ck == DOF_Electrode_potential ) {
        //R[i] = (Rcc_Electrode_potential[_i] + Rcxi_Electrode_potential[_i] - Kcxi_Electrode_potential(_i,0) / Kxixi_Electrode_potential(0,0) * rr_Electrode_potential[0]) * (cell_SDdata[cell_id].computed_area /dummy_area) ; 
        R[i] = (Rcc_Electrode_potential[_i] + Rcxi_Electrode_potential[_i] - Kcxi_Electrode_potential(_i,0) / Kxixi_Electrode_potential(0,0) * rr_Electrode_potential[0]) ; 

       //std::cout << "R[i] (Electrode_potential) " << R[i] << std::endl;
      }

      if (ck == DOF_Lithium_cation ) {
        //R[i] = (Rcc_Lithium_cation[_i] + Rcxi_Lithium_cation[_i] - Kcxi_Lithium_cation(_i,0) / Kxixi_Lithium_cation(0,0) * rr_Lithium_cation[0]) * (cell_SDdata[cell_id].computed_area /dummy_area_opposite); 
        R[i] = (Rcc_Lithium_cation[_i] + Rcxi_Lithium_cation[_i] - Kcxi_Lithium_cation(_i,0) / Kxixi_Lithium_cation(0,0) * rr_Lithium_cation[0]) ; 
       //std::cout << "R[i] (Lithium_cation) " << R[i] << std::endl;
      }

      if (ck == DOF_Electrolyte_potential ) {
        //std::cout << "R[i] " << R[i] << std::endl;
        //R[i] = (Rcc_Electrolyte_potential[_i] + Rcxi_Electrolyte_potential[_i] - Kcxi_Electrolyte_potential(_i,0) / Kxixi_Electrolyte_potential(0,0) * rr_Electrolyte_potential[0]) * (cell_SDdata[cell_id].computed_area /dummy_area_opposite); 
        R[i] = (Rcc_Electrolyte_potential[_i] + Rcxi_Electrolyte_potential[_i] - Kcxi_Electrolyte_potential(_i,0) / Kxixi_Electrolyte_potential(0,0) * rr_Electrolyte_potential[0]) ; 
       //std::cout << "R[i] (Electrolyte_potential) " << R[i] << std::endl;
        //std::cout << "Rcc " << Rcc_Electrolyte_potential[_i] << " Rcxi " << Rcxi_Electrolyte_potential[_i] << " Kcxi " << Kcxi_Electrolyte_potential(_i,0) << " Kxixi " << Kxixi_Electrolyte_potential(0,0) * rr_Electrolyte_potential[0] << std::endl;
        //std::cout << "--a4-1-- " << i <<" " << R[i] << cell_SDdata[cell_id].computed_area /dummy_area_opposite<< std::endl;
      }

     if ( ck >= DOF_Displacement and ck < DOF_Displacement + dim){
       _ii++;
        R[i] = Ruu_Displacement_sd[_ii];
        Sacado::Fad::DFad<double> _tmp = 0;
        //std::cout << " ii " << _ii << std::endl;
        //if (cell_SDdata[cell_id].is_fractured)
        {
          if (disable_wd)
          {
           _tmp = Kuxi_Displacement_sd(_ii,0) * Kxixi_Displacement_sd_inv(0,0) * rr_u[0] + Kuxi_Displacement_sd(_ii, 1) * Kxixi_Displacement_sd_inv(1,1) * rr_u[1]; 
          }
          else
          {
           _tmp = Kuxi_Displacement_sd(_ii,0) * Kxixi_Displacement_sd_inv(0,0) * rr_u[0] + Kuxi_Displacement_sd(_ii, 1) * Kxixi_Displacement_sd_inv(1,1) * rr_u[1] + Kuxi_Displacement_sd(_ii, 2) * Kxixi_Displacement_sd_inv(2,2) * rr_u[2]; 
          }
        }
        R[i] = R[i] - _tmp;

        //std::cout << "R[i] (Displacement) " << R[i] << std::endl;
     }
     //std::cout << " R[i] " << i << " _i " << _i << " _ii " << _ii << " R[i]" << R[i]<< std::endl;
    }
    //exit(-1);

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

    cell_SDdata[cell_id].Kxiu_sd = Kxiu_Displacement_sd;
    cell_SDdata[cell_id].rlocal_u_sd[0] = rr_u[0].val();
    cell_SDdata[cell_id].rlocal_u_sd[1] = rr_u[1].val();
    cell_SDdata[cell_id].rlocal_u_sd[2] = rr_u[2].val();
    for (unsigned int j = 0; j < 3; j++) {
    for (unsigned int k = 0; k < 3; k++) {
      cell_SDdata[cell_id].Kxixi_inv_u_sd(j,k) = Kxixi_Displacement_sd_inv(j,k);
    }
    }

      this->jump_n[cell_id] = ULocal_xi[dofs_per_cell + ind_Displacement_sd_1].val() ;
      this->jump_m[cell_id] = ULocal_xi[dofs_per_cell + ind_Displacement_sd_2].val() ;
      this->jump_w[cell_id] = ULocal_xi[dofs_per_cell + ind_Displacement_wd].val() ;

  for (unsigned int q = 0; q < n_q_points; ++q) {
    if (not cell_SDdata[cell_id].is_fractured)
    {
      if (TN_at_gp[q] > this->T_n[cell_id])  this->T_n[cell_id] = TN_at_gp[q].val();
    }
      //std::cout << " TN_at_gp[q] " << q << " " << TN_at_gp[q].val() << " TM_at_gp[q] " << TM_at_gp[q].val() << " loc = " << fe_values.quadrature_point(q) << " is_new " << this->is_new_step[cell_id] << std::endl;
    //if (TN_at_gp[q] > this->MatData["MaximumStressN"] and not cell_SDdata[cell_id].is_fractured)
    if (TN_at_gp[q] > (*displacement.params_json)["Mechanics"]["MaximumStressN"] and not cell_SDdata[cell_id].is_fractured and this->is_new_step[cell_id])
    {
      std::cout << "-------------------- crack detected ---------------------------"
        << " TN_at_gp[q] " << q << " " << TN_at_gp[q].val() << " TM_at_gp[q] " << TM_at_gp[q].val() << " loc = " << fe_values.quadrature_point(q) << std::endl;
      cell_SDdata[cell_id].is_fractured = true;
      cell_SDdata[cell_id].max_Tn = TN_at_gp[q].val();
      this->crack_id[cell_id] = cell_id ;
    }
  }
  //if (cell_SDdata[cell_id].is_fractured)  OutData.Data["crack_id"].push_back(cell_id);



  // not going to work in this two fields
  std::vector<std::vector<double>> dummy_pressure;
	if(battery_fields.active_fields_index["Diffuse_interface"]>-1) diffuse_interface.r_get_residual(fe_values, R, ULocal, ULocalConv, dummy_pressure);
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
      if (ck >= DOF_Displacement and ck < DOF_Displacement + dim) {
        //std::cout << "R_Displacement[i] " << R[i] << " " << dofs_per_cell << std::endl;
      }
      if (R[i].val() != R[i].val()) 
      {
        std::cout << "compute residual at the interface " << cell_id  << " center: " << center << std::endl;
        std::cout << " i " << i << " R[i] " << R[i] << std::endl;
    for (unsigned int q = 0; q < n_q_points; ++q) {
    std::cout << "-------q-------" << q << std::endl;
    std::cout << "diffu_lithium [0] " << diffu_Lithium[q][0] << " tilde " << diffu_Lithium_tilde[q][0] << std::endl;
    std::cout << "diffu_lithium [1] " << diffu_Lithium[q][1] << " tilde " << diffu_Lithium_tilde[q][1] << std::endl;
    std::cout << "grad_lithium [0]  " << battery_fields.quad_fields[DOF_Lithium].value_grad[q][0] << " tilde " << c_1_tilde_grad_Lithium[q][0] << std::endl;
    std::cout << "grad_lithium [1]  " << battery_fields.quad_fields[DOF_Lithium].value_grad[q][1] << " tilde " << c_1_tilde_grad_Lithium[q][1] << std::endl;
    std::cout << "diffu_Lithium_cation [0] " << diffu_Lithium_cation[q][0] << " tilde " << diffu_Lithium_cation_tilde[q][0] << std::endl;
    std::cout << "diffu_Lithium_cation [1] " << diffu_Lithium_cation[q][1] << " tilde " << diffu_Lithium_cation_tilde[q][1] << std::endl;
    std::cout << "field_Electrolyte_potential [0] " << field_Electrolyte_potential[q][0] << " tilde " << field_Electrolyte_potential_tilde[q][0] << std::endl;
    std::cout << "field_Electrolyte_potential [1] " << field_Electrolyte_potential[q][1] << " tilde " << field_Electrolyte_potential_tilde[q][1] << std::endl;

    }

        exit(-1);
      }
    }
    //std::cout << " residual done at the interface " << std::endl;

    //std::cout <<  " Xi[Lithium] " << ULocal_xi[dofs_per_cell+ ind_Lithium].val() << std::endl;
    //std::cout <<  " Xi[Lithium_cation] " << ULocal_xi[dofs_per_cell+ ind_Lithium_cation].val() << std::endl;
    //std::cout <<  " Xi[Electrode_potential] " << ULocal_xi[dofs_per_cell+ ind_Electrode_potential].val() << std::endl;
    //std::cout <<  " Xi[Electrolyte_potential] " << ULocal_xi[dofs_per_cell+ ind_Electrolyte_potential].val() << std::endl;

  //if (cell_SDdata[cell_id].is_fractured)
    //std::cout << "cell_id " << cell_id 
      //<< " xi_u_0 " << xi_0_Displacement_sd[0].val() 
      //<< " xi_u_1 " << xi_0_Displacement_sd[1].val() 
      //<< " TN[0] " << TN_at_gp[0].val() 
      //<< " 1 " << TN_at_gp[1].val()
      //<< " 2 " << TN_at_gp[2].val()
      //<< " 3 " << TN_at_gp[3].val()
      //<< std::endl;

  // save previous iteration solution for SD
	for (unsigned int i=0; i<dofs_per_cell; ++i){
    cell_SDdata[cell_id].ULocal_k[i] = ULocal[i].val(); 
	}
  if (this->is_new_step[cell_id]) this->is_new_step[cell_id] = false;
}

template class battery<1>;
template class battery<2>;
template class battery<3>;
