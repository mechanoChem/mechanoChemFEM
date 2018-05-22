#include "initBoundValProbs.h"

/*
*interface_electrolyte_activeMaterial
*/
template <int dim>
void initBoundValProbs<dim>::assemble_interface_electrolyte_activeMaterial_term(const typename hp::DoFHandler<dim>::active_cell_iterator &cell, unsigned int f,
																					 					 dealii::Table<1, Sacado::Fad::DFad<double> > &ULocal,
																					 					 Table<1, double > &ULocalConv,PETScWrappers::Vector &localized_U,
																					 					 dealii::Table<1, Sacado::Fad::DFad<double> >& R)
{
	params->enter_subsection("ElectroChemo");
	double F=params->get_double("F");
	double Rr=params->get_double("Rr");
	params->leave_subsection();	
	
	hp::FEValues<dim> hp_fe_values (fe_collection, q_collection, update_values | update_quadrature_points  | update_JxW_values | update_gradients);
  FEFaceValues<dim> activeMaterial_fe_face_values (*fe_system[active_material_fe], *common_face_quadrature, update_values | update_quadrature_points | update_JxW_values | update_normal_vectors | update_gradients);
  FEFaceValues<dim> electrolyte_fe_face_values (*fe_system[electrolyte_fe], *common_face_quadrature, update_values | update_quadrature_points | update_JxW_values | update_normal_vectors | update_gradients); 
	
  const unsigned int activeMaterial_dofs_per_cell = fe_system[active_material_fe]->dofs_per_cell;
  const unsigned int electrolyte_dofs_per_cell = fe_system[electrolyte_fe]->dofs_per_cell;
	
	hp_fe_values.reinit (cell);
	const FEValues<dim> &electrolyte_fe_values = hp_fe_values.get_present_fe_values();
	hp_fe_values.reinit (cell->neighbor(f));
	const FEValues<dim> &activeMaterial_fe_values = hp_fe_values.get_present_fe_values();
	
	electrolyte_fe_face_values.reinit(cell, f);
	activeMaterial_fe_face_values.reinit(cell->neighbor(f), cell->neighbor_of_neighbor(f));
	std::vector<types::global_dof_index> electrolyte_neighbor_dof_indices(activeMaterial_dofs_per_cell);
	cell->neighbor(f)->get_dof_indices (electrolyte_neighbor_dof_indices);
	
	//value of activeMaterial cell
	Table<1, Sacado::Fad::DFad<double> > ULocal_activeMaterial(activeMaterial_dofs_per_cell);
  for (unsigned int i=0; i<activeMaterial_dofs_per_cell; ++i){
		if (std::abs(localized_U(electrolyte_neighbor_dof_indices[i]))<1.0e-16) ULocal_activeMaterial[i]=0.0;
		else{ULocal_activeMaterial[i]=localized_U(electrolyte_neighbor_dof_indices[i]);}
  }
	/*
	//if using previous time step Please check localized_U should be localized_Un!
	Table<1, Sacado::Fad::DFad<double> > UnLocal(electrolyte_dofs_per_cell);
	for (unsigned int i=0; i<electrolyte_dofs_per_cell; ++i){
		UnLocal[i]=ULocalConv[i];
	}
	*/
  const unsigned int n_face_quadrature_points = activeMaterial_fe_face_values.n_quadrature_points;
	deformationMap<Sacado::Fad::DFad<double>, dim> defMap_mesh_surface(n_face_quadrature_points); 
	getDeformationMap<Sacado::Fad::DFad<double>, dim>(electrolyte_fe_values, electrolyte_fe_face_values, u_m_dof, ULocal, defMap_mesh_surface);
	
	
  dealii::Table<1,Sacado::Fad::DFad<double> > T(n_face_quadrature_points), c_li(n_face_quadrature_points),c_e(n_face_quadrature_points), c_li_plus(n_face_quadrature_points), phi_s(n_face_quadrature_points), phi_e(n_face_quadrature_points);
	dealii::Table<1,Sacado::Fad::DFad<double> > j_plus_interface(n_face_quadrature_points);
  dealii::Table<1,Sacado::Fad::DFad<double> > Grad_phi(n_face_quadrature_points);
	
	evaluateScalarFunction<Sacado::Fad::DFad<double>,dim>(electrolyte_fe_values, electrolyte_fe_face_values, c_li_plus_dof, ULocal, c_li_plus);
	evaluateScalarFunction<Sacado::Fad::DFad<double>,dim>(electrolyte_fe_values, electrolyte_fe_face_values, phi_e_dof, ULocal, phi_e);

  evaluateScalarFunction<Sacado::Fad::DFad<double>,dim>(activeMaterial_fe_values, activeMaterial_fe_face_values, T_dof, ULocal_activeMaterial, T);
	evaluateScalarFunction<Sacado::Fad::DFad<double>,dim>(activeMaterial_fe_values, activeMaterial_fe_face_values, c_li_dof, ULocal_activeMaterial, c_li);
	evaluateScalarFunction<Sacado::Fad::DFad<double>,dim>(activeMaterial_fe_values, activeMaterial_fe_face_values, phi_s_dof, ULocal_activeMaterial, phi_s);

	int domainflag;
	int Ni;
	const Point<dim> X = activeMaterial_fe_face_values.quadrature_point(0);
	if(X[1]<electrode_Y1+3){
	  Ni=-1;
		 domainflag=-1;
	 }
	else if(X[1]>electrode_Y2-3){
	  Ni=1;
		 domainflag=1;
	}
	
	for(unsigned int q=0;q<n_face_quadrature_points;q++){
		Sacado::Fad::DFad<double> jn=electricChemoFormula->formula_jn(T[q], c_li[q], c_li_plus[q], phi_s[q], phi_e[q], domainflag);
		j_plus_interface[q]=-jn;
		Grad_phi[q]=-jn*F;
	}

	ResidualEq->residualForNeummanBC(electrolyte_fe_values, electrolyte_fe_face_values, c_li_plus_dof, R, defMap_mesh_surface, j_plus_interface);
	ResidualEq->residualForNeummanBC(electrolyte_fe_values, electrolyte_fe_face_values, phi_e_dof, R, defMap_mesh_surface, Grad_phi);			
}

template class initBoundValProbs<1>;
template class initBoundValProbs<2>;
template class initBoundValProbs<3>;
