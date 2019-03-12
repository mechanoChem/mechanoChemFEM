#include "initBoundValProbs.h"

template <int dim>
void initBoundValProbs<dim>::refine_grid (){
	//nothing need here currently
  pcout << "Further refinement in progress" << std::endl;
  bool furtherRefine=true;
	params->enter_subsection("Problem");
		double inner_radius=params->get_double("inner_radius");
		double outer_radius=params->get_double("outer_radius");
		double a_scale=params->get_double("a_scale");
		double b_scale=params->get_double("b_scale");
		double c_scale=params->get_double("c_scale");
		double cort_thickness=params->get_double("cort_thickness");
	params->leave_subsection();	
	double finestCellWidth = 0.06*(outer_radius-inner_radius);
	
	while (furtherRefine){
		furtherRefine=false;
	   typename hp::DoFHandler<dim>::active_cell_iterator cell = this->dof_handler.begin_active(), endc = this->dof_handler.end();
		 for (;cell!=endc; ++cell){
			 Point<dim> cellcenter = cell->center();
			 if (std::sqrt(cell->center().square()) >= 0.2*inner_radius + 0.8*outer_radius){
			 	  double cellWidth=std::pow(cell->measure(), 1.0/dim);
			 	  if (cellWidth>finestCellWidth){
			 	    cell->set_refine_flag();
			 	    if (cellWidth>3*finestCellWidth){
			 	      furtherRefine=true;
			 	    }
			 	  }
			 	}
			 /*
			 if ((std::pow(cellcenter(0),2))/(std::pow((a_scale*outer_radius),2))+(std::pow(cellcenter(1),2))/(std::pow((b_scale*outer_radius),2))+std::pow(cellcenter(2),2)/(std::pow((c_scale*outer_radius),2)) >= 1-cort_thickness){
				 double cellWidth=std::pow(cell->measure(), 1.0/dim);
				 if (cellWidth>finestCellWidth){
					 cell->set_refine_flag();
					 if (cellWidth>3*finestCellWidth){
						 furtherRefine=true;
					 }
				 }
			 }
			 */
		 }
	   this->triangulation.execute_coarsening_and_refinement ();
	 }
	 pcout << "Further refinement DONE" << std::endl;
}

template class initBoundValProbs<1>;
template class initBoundValProbs<2>;
template class initBoundValProbs<3>;
