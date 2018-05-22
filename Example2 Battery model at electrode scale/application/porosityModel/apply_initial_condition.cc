#include "initBoundValProbs.h"

template <int dim>
void initBoundValProbs<dim>::apply_initial_condition()
{ 
	pcout << "applying initial condition\n";
		
	params->enter_subsection("Initial condition");
  double c_li_max_neg=params->get_double("c_li_max_neg");
	double c_li_max_pos=params->get_double("c_li_max_pos");
	double c_li_100_neg=params->get_double("c_li_100_neg");
  double c_li_100_pos=params->get_double("c_li_100_pos");
	double c_li_plus_ini=params->get_double("c_li_plus_ini");
	double T_ini=params->get_double("T_0");
  params->leave_subsection();		
  
	solution_prev=0;
  typename hp::DoFHandler<dim>::active_cell_iterator cell = this->dof_handler.begin_active(), endc=this->dof_handler.end();
  for (;cell!=endc; ++cell){
		if (cell->subdomain_id() == this_mpi_process){
    	hp::FEValues<dim> hp_fe_values (fe_collection, q_collection, update_values | update_quadrature_points  | update_JxW_values | update_gradients);
    	hp_fe_values.reinit (cell);
			const Point<dim> cell_center = cell->center();
    	const FEValues<dim> &fe_values = hp_fe_values.get_present_fe_values();
    	const unsigned int dofs_per_cell = cell->get_fe().dofs_per_cell;

			std::vector<unsigned int> local_dof_indices (dofs_per_cell);
 		 	cell->get_dof_indices (local_dof_indices);
    	for (unsigned int i=0; i<dofs_per_cell; ++i) {
      	const unsigned int ck = fe_values.get_fe().system_to_component_index(i).first;
				if(ck==c_li_plus_dof){	
					solution_prev(local_dof_indices[i])=c_li_plus_ini;//C_li_plus
	    	}
				if(ck==phi_e_dof){	
					solution_prev(local_dof_indices[i])=-electricChemoFormula->formula_Usc(c_li_100_neg,-1).val();//
				}
				
				if(ck==T_dof){	
					solution_prev(local_dof_indices[i])=T_ini;//
				}
			
	      if(cell_center[1]<=electrode_Y1){
			  	if(ck==c_li_dof){
						solution_prev(local_dof_indices[i])=c_li_100_neg*c_li_max_neg;//C_li
			  	}
		  	}
	    	if(cell_center[1]>=electrode_Y2 ){
			  	if(ck==c_li_dof){
						solution_prev(local_dof_indices[i])=c_li_100_pos*c_li_max_pos;//C_li
			 	  }
					if(ck==phi_s_dof){
						solution_prev(local_dof_indices[i])=electricChemoFormula->formula_Usc(c_li_100_pos,1).val()-electricChemoFormula->formula_Usc(c_li_100_neg,-1).val();//Phi_s
						//solution_prev(local_dof_indices[i])=0.0;
					}
	      }
			}
		}
	}
	solution_prev.compress(VectorOperation::insert);
	solution=solution_prev;
}

template class initBoundValProbs<1>;
template class initBoundValProbs<2>;
template class initBoundValProbs<3>;
