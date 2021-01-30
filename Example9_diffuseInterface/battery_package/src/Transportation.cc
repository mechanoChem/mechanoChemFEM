/*
zhenlin wang 2020
*module transportation
*/
#include "../include/Transportation.h"

template <int dim>
Transportation<dim>::Transportation(){}

template <int dim>
Transportation<dim>::Transportation(Battery_fields<dim>& _fields, Residual<Sacado::Fad::DFad<double>,dim>& _ResidualEq)
{
	battery_fields=&_fields;
	ResidualEq=&_ResidualEq;
}

template <int dim>
Transportation<dim>::Transportation(Battery_fields<dim>& _fields, Residual<Sacado::Fad::DFad<double>,dim>& _ResidualEq, int _primiary_dof)
{
	battery_fields=&_fields;
	ResidualEq=&_ResidualEq;
	primiary_dof=_primiary_dof;
}

template <int dim>
void Transportation<dim>::set_up_fields(Battery_fields<dim>& _battery_fields, Residual<Sacado::Fad::DFad<double>,dim>& _ResidualEq, int _primiary_dof)
{
	battery_fields=&_battery_fields;
	ResidualEq=&_ResidualEq;
	primiary_dof=_primiary_dof;
}

template <int dim>
void Transportation<dim>::set_up_fields(Battery_fields<dim>& _battery_fields, ElectricChemo<dim,Sacado::Fad::DFad<double>>& _electricChemoFormula, Residual<Sacado::Fad::DFad<double>,dim>& _ResidualEq, int _primiary_dof)
{
	battery_fields=&_battery_fields;
	electricChemoFormula=&_electricChemoFormula;
	ResidualEq=&_ResidualEq;
	primiary_dof=_primiary_dof;
}


template <int dim>
void Transportation<dim>::set_primiary_dof(int _primiary_dof)
{
	primiary_dof=_primiary_dof;
}

template <int dim>
void Transportation<dim>::r_get_residual(const FEValues<dim>& fe_values, Table<1, Sacado::Fad::DFad<double> >& R, Table<1, Sacado::Fad::DFad<double>>& ULocal, Table<1, double >& ULocalConv)
{
//evaluate primary fields
	unsigned int n_q_points= fe_values.n_quadrature_points;
	dealii::Table<1,Sacado::Fad::DFad<double> > react(n_q_points);
	dealii::Table<2,Sacado::Fad::DFad<double> > diffu(n_q_points, dim);
	diffu=table_scaling<2,Sacado::Fad::DFad<double>,double >(this->battery_fields->quad_fields[this->primiary_dof].value_grad,0);
	react=table_scaling<1,Sacado::Fad::DFad<double>,double >(react,0);
	set_diffusion_reaction_term(diffu,react);
		
	//call residual functions
	ResidualEq->residualForDiff_ReacEq(fe_values,primiary_dof, R,battery_fields->quad_fields[primiary_dof].value, battery_fields->quad_fields[primiary_dof].value_conv, diffu,react);
	
}

template <int dim>
void Transportation<dim>::r_get_residual_with_interface(const typename hp::DoFHandler<dim>::active_cell_iterator &cell, const FEValues<dim>& fe_values, Table<1, Sacado::Fad::DFad<double> >& R, Table<1, Sacado::Fad::DFad<double>>& ULocal, Table<1, double >& ULocalConv, std::vector<SDdata<dim>> &cell_SDdata)
{
//evaluate primary fields
	unsigned int n_q_points= fe_values.n_quadrature_points;
	dealii::Table<1,Sacado::Fad::DFad<double> > react(n_q_points);
	dealii::Table<2,Sacado::Fad::DFad<double> > diffu(n_q_points, dim);
	diffu=table_scaling<2,Sacado::Fad::DFad<double>,double >(this->battery_fields->quad_fields[this->primiary_dof].value_grad,0);
	react=table_scaling<1,Sacado::Fad::DFad<double>,double >(react,0);
	set_diffusion_reaction_term(diffu,react);
		
	//call residual functions
	//ResidualEq->residualForDiff_ReacEq(fe_values,primiary_dof, R,battery_fields->quad_fields[primiary_dof].value, battery_fields->quad_fields[primiary_dof].value_conv, diffu,react);
  int DOF = primiary_dof;
  std::vector<int> this_dof_local_index;
  {
    unsigned int dofs_per_cell= fe_values.dofs_per_cell;
    unsigned int n_q_points= fe_values.n_quadrature_points;
    
    //evaluate Residual: need to be modified later to include the jump of concentration
    for (unsigned int i=0; i<dofs_per_cell; ++i) {
      const unsigned int ck = fe_values.get_fe().system_to_component_index(i).first-DOF;
      if (ck==0) this_dof_local_index.push_back(i);
      for (unsigned int q=0; q<n_q_points; ++q){
        if (ck==0){
    			R[i] +=  fe_values.shape_value(i, q)*((battery_fields->quad_fields[primiary_dof].value[q]-battery_fields->quad_fields[primiary_dof].value_conv[q])/ResidualEq->dt)*fe_values.JxW(q);
    			for (unsigned int j = 0; j < dim; j++){
     			 	R[i] += -fe_values.shape_grad(i, q)[j]*diffu[q][j]*fe_values.JxW(q);
    			}
    			R[i] +=  -fe_values.shape_value(i, q)*react[q]*fe_values.JxW(q);
        }
      }
    }	
  }

  // surface integration for 2D element with a line interface
  int cell_id = cell->active_cell_index();
  //std::cout <<  " with interface residual, cell_id = " << cell_id << std::endl;

  Vector<double> r_local(2);
  double vol = 0.0;
  for (unsigned int i = 0; i < 2; ++i) {
      for (unsigned q=0; q< cell_SDdata[cell_id].jxw_1d.size(); ++q)
      {
          //std::cout << " q " << q << std::endl;
          vol += cell_SDdata[cell_id].jxw_1d(q);
          r_local[i] -= cell_SDdata[cell_id].shape_value_1d(i, q) * cell_SDdata[cell_id].reaction_rate * cell_SDdata[cell_id].jxw_1d(q);
      } // q_point
  }
  double total_r_local = r_local[0]  + r_local[1] ;
  double r_for_each_R = total_r_local / cell_SDdata[cell_id].lnode_plus.size();
  //std::cout 
      //<< " r_l1: " << r_local[0] 
      //<< " r_l2: " << r_local[1] 
      //<< " vol: " << vol
      //<< std::endl;
  for (unsigned i=0; i< cell_SDdata[cell_id].lnode_plus.size(); ++i)
  {
      int plus_node = this_dof_local_index[cell_SDdata[cell_id].lnode_plus[i]];
      //std::cout 
          //<< " plus_node: " << cell_SDdata[cell_id].lnode_plus[i] 
          //<< " dof_local_index " << plus_node
          //<< " R[plus_nod] "  << R[plus_node].val()
          //<< " r_for_each_R "  << r_for_each_R
          //<< std::endl;

      R[plus_node] += r_for_each_R;
  }

  //for (unsigned i=0; i< cell_SDdata[cell_id].lnode_minus.size(); ++i)
  //{
      //int minus_node = this_dof_local_index[cell_SDdata[cell_id].lnode_minus[i]];
    //std::cout 
        //<< " minus_node: " << cell_SDdata[cell_id].lnode_minus[i] 
        //<< " dof_local_index " << minus_node
        //<< std::endl;
  //}
	
}

template <int dim>
void Transportation<dim>::set_diffusion_reaction_term(dealii::Table<2,Sacado::Fad::DFad<double> >& diffu, dealii::Table<1,Sacado::Fad::DFad<double> >& react)
{
	unsigned int n_q_points= react.size(0);
	for (unsigned int q=0; q<n_q_points; ++q) {
		react[q]=0;
		for (unsigned int i=0; i<dim; ++i) diffu[q][i]=-battery_fields->quad_fields[primiary_dof].value_grad[q][i];
	}
}


template class Transportation<1>;
template class Transportation<2>;
template class Transportation<3>;
