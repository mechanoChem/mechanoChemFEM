/*
zhenlin wang 2019
*/

#include"../../include/mechanoChemFEM.h"

template <int dim>
void mechanoChemFEM<dim>::updateLinearSystem()
{
	this->reinitLinearSystem();
	ini_updateLinearSystem();
	
  hp::FEValues<dim> hp_fe_values (fe_collection, q_collection, update_values | update_quadrature_points  | update_JxW_values | update_gradients);	
  FullMatrix<double> local_matrix;
  Vector<double>            local_rhs;
  std::vector<types::global_dof_index> local_dof_indices;
	
  //loop over cells
  dealii::Vector<double> localized_U(solution);
  dealii::Vector<double> localized_Un(solution_prev);
	
	ResidualEq.dt=current_dt;
	
  typename hp::DoFHandler<dim>::active_cell_iterator cell = hpFEM<dim>::dof_handler.begin_active(), endc=hpFEM<dim>::dof_handler.end();
  for (;cell!=endc; ++cell){
		if (cell->subdomain_id() == this_mpi_process){	
			hp_fe_values.reinit (cell);
			const unsigned int dofs_per_cell = cell->get_fe().dofs_per_cell;
    	const FEValues<dim> &fe_values = hp_fe_values.get_present_fe_values();
			
    	local_matrix.reinit (dofs_per_cell, dofs_per_cell);
    	local_rhs.reinit (dofs_per_cell);
    	local_dof_indices.resize (dofs_per_cell);
			cell->get_dof_indices (local_dof_indices);
			
    	//AD variables
    	Table<1, Sacado::Fad::DFad<double> > ULocal(dofs_per_cell);
			Table<1, double > ULocalConv(dofs_per_cell);
	  	for (unsigned int i=0; i<dofs_per_cell; ++i){
				if (std::abs(localized_U(local_dof_indices[i]))<1.0e-16) ULocal[i]=0.0;
				else{ULocal[i]=localized_U(local_dof_indices[i]);}
				ULocal[i].diff (i, dofs_per_cell);
				ULocalConv[i]= localized_Un(local_dof_indices[i]);
	  	}
		
    	Table<1, Sacado::Fad::DFad<double> > R(dofs_per_cell); 
			for(unsigned int i=0;i<dofs_per_cell;i++) R[i]=0.0;
								
			get_residual(cell,fe_values,R, ULocal, ULocalConv);
		
    	//Residual(R) and Jacobian(R')		
    	for (unsigned int i=0; i<dofs_per_cell; ++i) {
      	for (unsigned int j=0; j<dofs_per_cell; ++j){
					// R' by AD
					local_matrix(i,j)= R[i].dx(j);
      	}
      	//R
      	local_rhs(i) = -R[i].val(); 
   	  }
			this->distribute_local_to_global(local_matrix, local_rhs, local_dof_indices);
		}
	}
	this->LinearSystemCompressAdd();
}

template class mechanoChemFEM<1>;
template class mechanoChemFEM<2>;
template class mechanoChemFEM<3>;
