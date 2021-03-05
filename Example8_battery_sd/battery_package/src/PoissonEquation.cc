/*
zhenlin wang 2019
*module transportation
*/
#include "../include/PoissonEquation.h"

template <int dim>
PoissonEquation<dim>::PoissonEquation(){}

template <int dim>
PoissonEquation<dim>::PoissonEquation(Battery_fields<dim>& _fields, Residual<Sacado::Fad::DFad<double>,dim>& _ResidualEq, int _primiary_dof)
{
	battery_fields=&_fields;
	ResidualEq=&_ResidualEq;
	primiary_dof=_primiary_dof;
}

template <int dim>
void PoissonEquation<dim>::set_up_fields(Battery_fields<dim>& _fields, Residual<Sacado::Fad::DFad<double>,dim>& _ResidualEq, int _primiary_dof)
{
	battery_fields=&_fields;
	ResidualEq=&_ResidualEq;
	primiary_dof=_primiary_dof;
}

template <int dim>
void PoissonEquation<dim>::set_up_fields(Battery_fields<dim>& _fields, ElectricChemo<dim,Sacado::Fad::DFad<double>>& _electricChemoFormula, Residual<Sacado::Fad::DFad<double>,dim>& _ResidualEq, int _primiary_dof)
{
	battery_fields=&_fields;
	electricChemoFormula=&_electricChemoFormula;
	ResidualEq=&_ResidualEq;
	primiary_dof=_primiary_dof;
}

template <int dim>
void PoissonEquation<dim>::r_get_residual(const FEValues<dim>& fe_values, Table<1, Sacado::Fad::DFad<double> >& R, Table<1, Sacado::Fad::DFad<double>>& ULocal, Table<1, double >& ULocalConv)
{
//evaluate primary fields
	unsigned int n_q_points= fe_values.n_quadrature_points;
	
	dealii::Table<2,Sacado::Fad::DFad<double> > field(n_q_points,dim);
	dealii::Table<1,Sacado::Fad::DFad<double> > source(n_q_points);
	set_field_and_source_term(field,source);
	
	//call residual functions
	ResidualEq->residualForPoissonEq(fe_values, primiary_dof, R, field, source);
	
}

template <int dim>
void PoissonEquation<dim>::r_get_residual_with_interface(const typename hp::DoFHandler<dim>::active_cell_iterator &cell, const FEValues<dim>& fe_values, Table<1, Sacado::Fad::DFad<double> >& R, Table<1, Sacado::Fad::DFad<double>>& ULocal, Table<1, double >& ULocalConv, std::vector<SDdata<dim>> &cell_SDdata)
{
  int cell_id = cell->active_cell_index();
  unsigned int dofs_per_cell= fe_values.dofs_per_cell;
  std::vector<unsigned int> local_dof_indices(dofs_per_cell);
  cell->get_dof_indices(local_dof_indices);

  double D_1 = 1.0;
//evaluate primary fields
	unsigned int n_q_points= fe_values.n_quadrature_points;

	dealii::Table<2,Sacado::Fad::DFad<double> > field(n_q_points,dim);
	dealii::Table<1,Sacado::Fad::DFad<double> > source(n_q_points);
	source=table_scaling<1,Sacado::Fad::DFad<double> > (source,0);
  // Question: what's the diffusion coefficient for this part?

  dealii::Table<1, double> coeff(n_q_points);
  for (unsigned int q = 0; q < n_q_points; q++) {
    coeff[q] = 1.0;
    for (unsigned int i=0; i<dim; ++i) field[q][i] = D_1 * battery_fields->quad_fields[primiary_dof].value_grad[q][i]; // is there a minus sign needed?
		//field=battery_fields->quad_fields[primiary_dof].value_grad;
  }
		
  //std::cout << "--b0-0 (primary_dof)--: " << primiary_dof << " " << this->primiary_dof  << std::endl;
	//call residual functions
  int DOF = primiary_dof;
  std::vector<int> this_dof_local_index;
  {
    unsigned int dofs_per_cell= fe_values.dofs_per_cell;
    unsigned int n_q_points= fe_values.n_quadrature_points;
    
    //evaluate Residual: need to be modified later to include the jump of concentration
    for (unsigned int i=0; i<dofs_per_cell; ++i) {
      const unsigned int ck = fe_values.get_fe().system_to_component_index(i).first-DOF;
      if (ck==0) this_dof_local_index.push_back(i);
    }	
  }
  //std::cout << "--b0-1--" << std::endl;

  Vector<double> dC_k1;
  dC_k1.reinit(4); // 4 = this_dof_local_index.size()
  //std::cout << "--b0-2--" << std::endl;

  {
    unsigned int _i = 0;
    for (unsigned int i = 0; i < dofs_per_cell; ++i) {
      const unsigned int ck = fe_values.get_fe().system_to_component_index(i).first - DOF;
      if (ck == 0) {
        dC_k1(_i) = ULocal[i].val()  - cell_SDdata[cell_id].ULocal_k[i];
        _i ++;
      }
    }
  }

  Vector<double> dxi_k1;
  dxi_k1.reinit(1);
  Table<1, Sacado::Fad::DFad<double>> xi_0(1);  // define sacado xi_0 for stiffness calculation.
  if (primiary_dof != cell_SDdata[cell_id].opposite_flux_dof_potential)
  {
    cell_SDdata[cell_id].Kxic_phi_s.vmult(dxi_k1, dC_k1);
    cell_SDdata[cell_id].rlocal_phi_s -= dxi_k1;
    cell_SDdata[cell_id].rlocal_phi_s[0] += cell_SDdata[cell_id].reaction_rate_potential * cell_SDdata[cell_id].interface_length;
    cell_SDdata[cell_id].Kxixi_inv_phi_s.vmult(dxi_k1, cell_SDdata[cell_id].rlocal_phi_s);
    xi_0[0] = cell_SDdata[cell_id].xi_old_phi_s(0) + dxi_k1(0);  
    cell_SDdata[cell_id].xi_old_phi_s(0) = xi_0[0].val();
    //std::cout << "--a0-0--"  << std::endl;
  }
  else
  {
    //std::cout << "--a0-1--" << std::endl;
    cell_SDdata[cell_id].Kxic_phi_e.vmult(dxi_k1, dC_k1);
    cell_SDdata[cell_id].rlocal_phi_e -= dxi_k1;
    cell_SDdata[cell_id].rlocal_phi_e[0] += cell_SDdata[cell_id].reaction_rate_potential * cell_SDdata[cell_id].interface_length; // reaction rate li direction should not change
    cell_SDdata[cell_id].Kxixi_inv_phi_e.vmult(dxi_k1, cell_SDdata[cell_id].rlocal_phi_e);
    xi_0[0] = cell_SDdata[cell_id].xi_old_phi_e(0) + dxi_k1(0);  
    cell_SDdata[cell_id].xi_old_phi_e(0) = xi_0[0].val();
  }
  xi_0[0].diff(0, 1);

  Table<1, Sacado::Fad::DFad<double> > ULocal_xi(dofs_per_cell + 1);
  for (unsigned int i = 0; i < dofs_per_cell; ++i) {
	  ULocal_xi[i]=ULocal[i].val();
  }
  ULocal_xi[dofs_per_cell] = xi_0[0].val();
  for (unsigned int i = 0; i < dofs_per_cell + 1; ++i) {
	  ULocal_xi[i].diff (i, dofs_per_cell + 1);
  }

  //std::cout << "--a1--" << xi_0[0]<< std::endl;

  Table<1, Sacado::Fad::DFad<double>> Rcc(4);
  Table<1, Sacado::Fad::DFad<double>> Rcxi(4);
  Table<1, Sacado::Fad::DFad<double>> rxic(1);
  Table<1, Sacado::Fad::DFad<double>> rxixi(1);
  Table<1, Sacado::Fad::DFad<double>> rr(1);
  rxic[0] = 0.0;
  rxixi[0] = 0.0;
  rr[0] = 0.0;
  for (unsigned int i = 0; i < 4; ++i) {
    Rcc[i] = 0.0;
    Rcxi[i] = 0.0;
  }

  dealii::Table<2, Sacado::Fad::DFad<double>> c_1_tilde_grad(n_q_points, dim);
  dealii::Table<1, Sacado::Fad::DFad<double>> c_1_tilde(n_q_points);
  dealii::Table<1, double> c_1_tilde_conv(n_q_points);

  for (unsigned int q = 0; q < n_q_points; ++q) {
      c_1_tilde[q] = 0.0;
      c_1_tilde_conv[q] = 0.0;
    for (unsigned int j = 0; j < dim; j++) {
      c_1_tilde_grad[q][j] = 0.0;
    }
  }
    
  if (primiary_dof != cell_SDdata[cell_id].opposite_flux_dof_potential)
  {
    std::vector<double> Ms_list;
    for (unsigned int q = 0; q < n_q_points; ++q) {
      Ms_list.push_back(0.0);
    }
    for (unsigned i = 0; i < cell_SDdata[cell_id].lnode_plus.size(); ++i) {
      Ms_list[cell_SDdata[cell_id].lnode_plus[i]] = 1.0;
    }

    double dummy_area = 0.0;
    for (unsigned int q = 0; q < n_q_points; ++q) {
      //double Ms = 1.0;
      //if (fe_values.quadrature_point(q)[0] < 0.5){
        //Ms = 0.0;
      //}
      double Ms = Ms_list[q];
      dummy_area += Ms * fe_values.JxW(q); // the actually int area
      for (unsigned i = 0; i < cell_SDdata[cell_id].lnode_plus.size(); ++i) {
        int plus_node = this_dof_local_index[cell_SDdata[cell_id].lnode_plus[i]];
        Ms -= fe_values.shape_value(plus_node, q);
        for (unsigned int j = 0; j < dim; j++) {
          c_1_tilde_grad[q][j] -= fe_values.shape_grad(plus_node, q)[j] * ULocal_xi[dofs_per_cell];  
        }
      }
      c_1_tilde[q] = Ms * ULocal_xi[dofs_per_cell];
      c_1_tilde_conv[q] = Ms * cell_SDdata[cell_id].xi_conv_phi_s[0];
    }
    rr[0] = - cell_SDdata[cell_id].reaction_rate_potential * cell_SDdata[cell_id].interface_length;

    //std::cout << "--a2--" << rr[0] << std::endl;

    for (unsigned int q = 0; q < n_q_points; ++q) {
          rxixi[0] += cell_SDdata[cell_id].interface_length / dummy_area * (cell_SDdata[cell_id].computed_area /dummy_area)  * ( c_1_tilde_grad[q][0] * cell_SDdata[cell_id].crk_n[0] +  c_1_tilde_grad[q][1] * cell_SDdata[cell_id].crk_n[1]) * fe_values.JxW(q); 
     }
    //std::cout << "--a2-1--" << rxixi[0]  << std::endl;

    for (unsigned int q = 0; q < n_q_points; ++q) {
        rxic[0] += cell_SDdata[cell_id].interface_length /  dummy_area * (cell_SDdata[cell_id].computed_area /dummy_area)  * (field[q][0] * cell_SDdata[cell_id].crk_n[0] + field[q][1] * cell_SDdata[cell_id].crk_n[1]) * fe_values.JxW(q); 
    }
    //std::cout << "--a2-2--" << rxic[0] << std::endl;

    rr[0] += rxixi[0] + rxic[0];

    for (unsigned int q = 0; q < n_q_points; ++q) {
      for (unsigned i = 0; i < cell_SDdata[cell_id].lnode_plus.size(); ++i) {
        int plus_node = this_dof_local_index[cell_SDdata[cell_id].lnode_plus[i]];
          //std::cout << "--a2-3--" << i << " " << q << " " << plus_node << " " << cell_SDdata[cell_id].lnode_plus[i] << std::endl;
        for (unsigned int j = 0; j < dim; j++) {
          Rcxi[cell_SDdata[cell_id].lnode_plus[i]] += - coeff[q] * D_1 * c_1_tilde_grad[q][j] * fe_values.shape_grad(plus_node, q)[j] * fe_values.JxW(q);  
          //std::cout << "--a2-4--" << plus_node << " " << cell_SDdata[cell_id].lnode_plus[i]<< std::endl;
        }
      }
    }

    //std::cout << "--a3--" << primiary_dof << " rr[0]" << rr[0] << std::endl;

    FullMatrix<double> Kxixi;
    Kxixi.reinit(1, 1);
    FullMatrix<double> Kxic;
    Kxic.reinit(1, 4);
    FullMatrix<double> Kcxi;
    Kcxi.reinit(4, 1);

    Kxixi(0, 0) = rxixi[0].dx(dofs_per_cell);
    Kxic(0, 0) = rxic[0].dx(0);
    Kxic(0, 1) = rxic[0].dx(1);
    Kxic(0, 2) = rxic[0].dx(2);
    Kxic(0, 3) = rxic[0].dx(3);
    Kcxi(0, 0) = Rcxi[0].dx(dofs_per_cell);
    Kcxi(1, 0) = Rcxi[1].dx(dofs_per_cell);
    Kcxi(2, 0) = Rcxi[2].dx(dofs_per_cell);
    Kcxi(3, 0) = Rcxi[3].dx(dofs_per_cell);

    for (unsigned int i = 0; i < 4; ++i) {
      Rcc[i] = ULocal_xi[0] * 0.0;
    }
    

    for (unsigned int q = 0; q < n_q_points; ++q) {
      for (unsigned i = 0; i < cell_SDdata[cell_id].lnode_plus.size(); ++i) {
        int plus_node = this_dof_local_index[cell_SDdata[cell_id].lnode_plus[i]];
        for (unsigned int j = 0; j < dim; j++) {
            Rcc[cell_SDdata[cell_id].lnode_plus[i]] += -fe_values.shape_grad(plus_node, q)[j]*field[q][j]*fe_values.JxW(q); // oscillation
        }
      }
    }

    //std::cout << "--a4-- [0]" << Rcc[0] << std::endl;
    //std::cout << "--a4-- [1]" << Rcc[1] << std::endl;
    //std::cout << "--a4-- [2]" << Rcc[2] << std::endl;
    //std::cout << "--a4-- [3]" << Rcc[3] << std::endl;
    //for (unsigned int i = 0; i < dofs_per_cell; ++i) {
      //R[i] = (Rcc[i] + Rcxi[i] - Kcxi(i,0) / Kxixi(0,0) * rr[0]) * (cell_SDdata[cell_id].computed_area /dummy_area) ; 
    //}

    unsigned int _i = 0;
    for (unsigned int i = 0; i < dofs_per_cell; ++i) {
      const unsigned int ck = fe_values.get_fe().system_to_component_index(i).first - DOF;
      if (ck == 0) {
        //std::cout << "R[i] " << R[i] << std::endl;
        R[i] = (Rcc[_i] + Rcxi[_i] - Kcxi(_i,0) / Kxixi(0,0) * rr[0]) * (cell_SDdata[cell_id].computed_area /dummy_area) ; 
        //std::cout << "Rcc " << Rcc[_i] << " " << Rcxi[_i] << " " << Kcxi(_i,0) << " " << Kxixi(0,0) * rr[0] << std::endl;
        //std::cout << "--a4-1-- " << i <<" " << R[i] << cell_SDdata[cell_id].computed_area /dummy_area<< std::endl;
        _i ++;
      }
    }

    cell_SDdata[cell_id].Kxic_phi_s = Kxic;
    cell_SDdata[cell_id].rlocal_phi_s[0] = rr[0].val();
    cell_SDdata[cell_id].Kxixi_inv_phi_s(0,0) = 1.0/Kxixi(0,0);
    //std::cout << "--a4-2-- " << std::endl;

  }
  else
  {
    std::vector<double> Ms_list;
    for (unsigned int q = 0; q < n_q_points; ++q) {
      Ms_list.push_back(0.0);
    }
    for (unsigned i = 0; i < cell_SDdata[cell_id].lnode_minus.size(); ++i) {
      Ms_list[cell_SDdata[cell_id].lnode_minus[i]] = 1.0;
    }

    double dummy_area = 0.0;
    for (unsigned int q = 0; q < n_q_points; ++q) {
      //double Ms = 1.0;
      //if (fe_values.quadrature_point(q)[0] < 0.5){
        //Ms = 0.0;
      //}
      double Ms = Ms_list[q];
      dummy_area += Ms * fe_values.JxW(q); // the actually int area
      for (unsigned i = 0; i < cell_SDdata[cell_id].lnode_minus.size(); ++i) {
        int minus_node = this_dof_local_index[cell_SDdata[cell_id].lnode_minus[i]];
        Ms -= fe_values.shape_value(minus_node, q);
        for (unsigned int j = 0; j < dim; j++) {
          c_1_tilde_grad[q][j] -= fe_values.shape_grad(minus_node, q)[j] * ULocal_xi[dofs_per_cell];  
        }
      }
      c_1_tilde[q] = Ms * ULocal_xi[dofs_per_cell];
      c_1_tilde_conv[q] = Ms * cell_SDdata[cell_id].xi_conv_phi_e[0];
    }
    rr[0] = - cell_SDdata[cell_id].reaction_rate_potential * cell_SDdata[cell_id].interface_length;

    //std::cout << "--a2--" << rr[0] << std::endl;

    for (unsigned int q = 0; q < n_q_points; ++q) {
          rxixi[0] += cell_SDdata[cell_id].interface_length / dummy_area * (cell_SDdata[cell_id].computed_area /dummy_area)  * ( - c_1_tilde_grad[q][0] * cell_SDdata[cell_id].crk_n[0] -  c_1_tilde_grad[q][1] * cell_SDdata[cell_id].crk_n[1]) * fe_values.JxW(q); // crk_n direction is reversed
     }
    //std::cout << "--a2-1--" << rxixi[0]  << std::endl;

    for (unsigned int q = 0; q < n_q_points; ++q) {
        rxic[0] += cell_SDdata[cell_id].interface_length /  dummy_area * (cell_SDdata[cell_id].computed_area /dummy_area)  * (- field[q][0] * cell_SDdata[cell_id].crk_n[0] - field[q][1] * cell_SDdata[cell_id].crk_n[1]) * fe_values.JxW(q); 
    }
    //std::cout << "--a2-2--" << rxic[0] << std::endl;

    rr[0] += rxixi[0] + rxic[0];

    for (unsigned int q = 0; q < n_q_points; ++q) {
      for (unsigned i = 0; i < cell_SDdata[cell_id].lnode_minus.size(); ++i) {
        int minus_node = this_dof_local_index[cell_SDdata[cell_id].lnode_minus[i]];
          //std::cout << "--a2-3--" << i << " " << q << " " << minus_node << " " << cell_SDdata[cell_id].lnode_minus[i] << std::endl;
        for (unsigned int j = 0; j < dim; j++) {
          Rcxi[cell_SDdata[cell_id].lnode_minus[i]] += - coeff[q] * D_1 * c_1_tilde_grad[q][j] * fe_values.shape_grad(minus_node, q)[j] * fe_values.JxW(q);  
          //std::cout << "--a2-4--" << minus_node << " " << cell_SDdata[cell_id].lnode_minus[i]<< std::endl;
        }
      }
    }

    //std::cout << "--a3--" << primiary_dof << " rr[0]" << rr[0] << std::endl;

    FullMatrix<double> Kxixi;
    Kxixi.reinit(1, 1);
    FullMatrix<double> Kxic;
    Kxic.reinit(1, 4);
    FullMatrix<double> Kcxi;
    Kcxi.reinit(4, 1);

    Kxixi(0, 0) = rxixi[0].dx(dofs_per_cell);
    Kxic(0, 0) = rxic[0].dx(0);
    Kxic(0, 1) = rxic[0].dx(1);
    Kxic(0, 2) = rxic[0].dx(2);
    Kxic(0, 3) = rxic[0].dx(3);
    Kcxi(0, 0) = Rcxi[0].dx(dofs_per_cell);
    Kcxi(1, 0) = Rcxi[1].dx(dofs_per_cell);
    Kcxi(2, 0) = Rcxi[2].dx(dofs_per_cell);
    Kcxi(3, 0) = Rcxi[3].dx(dofs_per_cell);

    for (unsigned int i = 0; i < 4; ++i) {
      Rcc[i] = ULocal_xi[0] * 0.0;
    }
    

    for (unsigned int q = 0; q < n_q_points; ++q) {
      for (unsigned i = 0; i < cell_SDdata[cell_id].lnode_minus.size(); ++i) {
        int minus_node = this_dof_local_index[cell_SDdata[cell_id].lnode_minus[i]];
        for (unsigned int j = 0; j < dim; j++) {
            Rcc[cell_SDdata[cell_id].lnode_minus[i]] += -fe_values.shape_grad(minus_node, q)[j]*field[q][j]*fe_values.JxW(q); // oscillation
        }
      }
    }

    //std::cout << "--a4-- [0]" << Rcc[0] << std::endl;
    //std::cout << "--a4-- [1]" << Rcc[1] << std::endl;
    //std::cout << "--a4-- [2]" << Rcc[2] << std::endl;
    //std::cout << "--a4-- [3]" << Rcc[3] << std::endl;
    //for (unsigned int i = 0; i < dofs_per_cell; ++i) {
      //R[i] = (Rcc[i] + Rcxi[i] - Kcxi(i,0) / Kxixi(0,0) * rr[0]) * (cell_SDdata[cell_id].computed_area /dummy_area) ; 
    //}

    unsigned int _i = 0;
    for (unsigned int i = 0; i < dofs_per_cell; ++i) {
      const unsigned int ck = fe_values.get_fe().system_to_component_index(i).first - DOF;
      if (ck == 0) {
        //std::cout << "R[i] " << R[i] << std::endl;
        R[i] = (Rcc[_i] + Rcxi[_i] - Kcxi(_i,0) / Kxixi(0,0) * rr[0]) * (cell_SDdata[cell_id].computed_area /dummy_area) ; 
        //std::cout << "Rcc " << Rcc[_i] << " " << Rcxi[_i] << " " << Kcxi(_i,0) << " " << Kxixi(0,0) * rr[0] << std::endl;
        //std::cout << "--a4-1 R[i]-- " << i <<" " << R[i] << cell_SDdata[cell_id].computed_area /dummy_area<< std::endl;
        _i ++;
      }
    }

    cell_SDdata[cell_id].Kxic_phi_e = Kxic;
    cell_SDdata[cell_id].rlocal_phi_e[0] = rr[0].val();
    cell_SDdata[cell_id].Kxixi_inv_phi_e(0,0) = 1.0/Kxixi(0,0);
    //std::cout << "--a4-2-- " << std::endl;
  }
    //std::cout << "--a4-4-- " << std::endl;
	
}

// should here a minus sign for the field?

template <int dim>
void PoissonEquation<dim>::set_field_and_source_term(dealii::Table<2,Sacado::Fad::DFad<double> >& field, dealii::Table<1,Sacado::Fad::DFad<double> >& source)
{
	field=battery_fields->quad_fields[primiary_dof].value_grad;
	source=table_scaling<1,Sacado::Fad::DFad<double> > (source,0);
}

template class PoissonEquation<1>;
template class PoissonEquation<2>;
template class PoissonEquation<3>;
