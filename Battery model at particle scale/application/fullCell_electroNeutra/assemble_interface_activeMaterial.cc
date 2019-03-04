#include "initBoundValProbs.h"

template <int dim>
void initBoundValProbs<dim>::assemble_interface_activeMaterial_electrolyte_term(const typename hp::DoFHandler<dim>::active_cell_iterator &cell, unsigned int f,
																					 					 dealii::Table<1, Sacado::Fad::DFad<double> > &ULocal,
																					 					 Table<1, double > &ULocalConv,PETScWrappers::Vector &localized_U,
																					 					 dealii::Table<1, Sacado::Fad::DFad<double> >& R)
{	
	params->enter_subsection("Initial condition");
	double c_li_max_neg=params->get_double("c_li_max_neg");
	double c_li_max_pos=params->get_double("c_li_max_pos");
	params->leave_subsection();	
		
	params->enter_subsection("Fluid");
	double viscosity=params->get_double("viscosity");
	params->leave_subsection();	
	
	params->enter_subsection("ElectroChemo");
	double F=params->get_double("F");
	double Rr=params->get_double("Rr");
	params->leave_subsection();	
	
	
  const unsigned int activeMaterial_dofs_per_cell = fe_system[active_material_fe]->dofs_per_cell;
  const unsigned int electrolyte_dofs_per_cell = fe_system[electrolyte_fe]->dofs_per_cell;
	
	hp::FEValues<dim> hp_fe_values (fe_collection, q_collection, update_values | update_quadrature_points  | update_JxW_values | update_gradients);
  FEFaceValues<dim> activeMaterial_fe_face_values (*fe_system[active_material_fe], *common_face_quadrature, update_values | update_quadrature_points | update_JxW_values | update_normal_vectors | update_gradients);
  FEFaceValues<dim> electrolyte_fe_face_values (*fe_system[electrolyte_fe], *common_face_quadrature, update_values | update_quadrature_points | update_JxW_values | update_normal_vectors | update_gradients); 
	
	hp_fe_values.reinit (cell);
	const FEValues<dim> &activeMaterial_fe_values = hp_fe_values.get_present_fe_values();
	hp_fe_values.reinit (cell->neighbor(f));
	const FEValues<dim> &electrolyte_fe_values = hp_fe_values.get_present_fe_values();
	
	activeMaterial_fe_face_values.reinit(cell, f);
	electrolyte_fe_face_values.reinit(cell->neighbor(f), cell->neighbor_of_neighbor(f));
	std::vector<types::global_dof_index> activeMaterial_neighbor_dof_indices (electrolyte_dofs_per_cell);
	cell->neighbor(f)->get_dof_indices (activeMaterial_neighbor_dof_indices);
	
	//value of electrolyte cell
	Table<1, Sacado::Fad::DFad<double> > ULocal_electrolyte(electrolyte_dofs_per_cell);
  for (unsigned int i=0; i<electrolyte_dofs_per_cell; ++i){
		if (std::abs(localized_U(activeMaterial_neighbor_dof_indices[i]))<1.0e-16) ULocal_electrolyte[i]=0.0;
		else{ULocal_electrolyte[i]=localized_U(activeMaterial_neighbor_dof_indices[i]);}
  }

	/*
	*the following feilds are on surface
	*/
  const unsigned int n_face_quadrature_points = activeMaterial_fe_face_values.n_quadrature_points;
	
	deformationMap<Sacado::Fad::DFad<double>, dim> defMap_surface(n_face_quadrature_points); 
	getDeformationMap<Sacado::Fad::DFad<double>, dim>(activeMaterial_fe_values, activeMaterial_fe_face_values, u_dof, ULocal, defMap_surface);
	
	deformationMap<Sacado::Fad::DFad<double>, dim> defMap_mesh_surface(n_face_quadrature_points); 
	getDeformationMap<Sacado::Fad::DFad<double>, dim>(electrolyte_fe_values, electrolyte_fe_face_values, u_m_dof, ULocal_electrolyte, defMap_mesh_surface);
	
	
  dealii::Table<1,Sacado::Fad::DFad<double> > T(n_face_quadrature_points), Pressure(n_face_quadrature_points), c_li(n_face_quadrature_points),  c_e(n_face_quadrature_points), c_li_plus(n_face_quadrature_points),phi_s(n_face_quadrature_points), phi_e(n_face_quadrature_points);
	dealii::Table<1,Sacado::Fad::DFad<double> > I_interface(n_face_quadrature_points), j_interface(n_face_quadrature_points),j_q_interface(n_face_quadrature_points);
  dealii::Table<1,Sacado::Fad::DFad<double> > Grad_phi(n_face_quadrature_points),Q(n_face_quadrature_points);
  dealii::Table<2,Sacado::Fad::DFad<double> > Phis_grad(n_face_quadrature_points,dim),stress(n_face_quadrature_points,dim);
	dealii::Table<3,Sacado::Fad::DFad<double> > velocity_grad(n_face_quadrature_points, dim,dim),P_stoke(n_face_quadrature_points,dim,dim);
	// initialize variables
	
	//evaluateScalarFunction<double,dim>(fe_values, 3, ULocalConv, c_li_conv);
	evaluateVectorFunctionGradient<Sacado::Fad::DFad<double>,dim>(electrolyte_fe_values, electrolyte_fe_face_values, v_dof, ULocal_electrolyte, velocity_grad,defMap_mesh_surface);
	evaluateScalarFunction<Sacado::Fad::DFad<double>,dim>(electrolyte_fe_values, electrolyte_fe_face_values, p_dof, ULocal_electrolyte, Pressure);
	evaluateScalarFunction<Sacado::Fad::DFad<double>,dim>(electrolyte_fe_values, electrolyte_fe_face_values, c_li_plus_dof, ULocal_electrolyte, c_li_plus);
	evaluateScalarFunction<Sacado::Fad::DFad<double>,dim>(electrolyte_fe_values, electrolyte_fe_face_values, phi_e_dof, ULocal_electrolyte, phi_e);
	
	evaluateScalarFunction<Sacado::Fad::DFad<double>,dim>(activeMaterial_fe_values, activeMaterial_fe_face_values, T_dof, ULocal, T);
	evaluateScalarFunction<Sacado::Fad::DFad<double>,dim>(activeMaterial_fe_values, activeMaterial_fe_face_values, c_li_dof, ULocal, c_li);
	evaluateScalarFunction<Sacado::Fad::DFad<double>,dim>(activeMaterial_fe_values, activeMaterial_fe_face_values, phi_s_dof, ULocal, phi_s);
	int Ni;
	int domainflag;
	double c_li_max;
	double k;
	const Point<dim> X = activeMaterial_fe_face_values.quadrature_point(0);
	if(X[1]<electrode_Y1+3){
	  Ni=1;
		 domainflag=-1;
		 c_li_max=c_li_max_neg;
	 }
	else if(X[1]>electrode_Y2-3){
	   Ni=-1;
		 domainflag=1;
		 c_li_max=c_li_max_pos;
	}
	
	for(unsigned int q=0;q<n_face_quadrature_points;q++){
		Sacado::Fad::DFad<double> usc=electricChemoFormula->formula_Usc(c_li[q]/c_li_max,domainflag);
		Sacado::Fad::DFad<double> jn=electricChemoFormula->formula_jn(T[q], c_li[q], c_li_plus[q], phi_s[q], phi_e[q], domainflag);
		I_interface[q]=jn*F;
		j_interface[q]=jn;
		j_q_interface[q]=-jn;
		Grad_phi[q]=I_interface[q];
		if(std::abs(jn*F)>3000) Q[q]=0;
	       else Q[q]=-(F*jn*(phi_s[q]-phi_e[q]-usc)+F*jn*T[q]*electricChemoFormula->formula_dUdt(c_li[q]/c_li_max));
		//Stokes's stress constitutive equation
		for (unsigned int i=0; i<dim; ++i){
	  	for (unsigned int j=0; j<dim; ++j){
				P_stoke[q][i][j]=viscosity*1.0e-9*(velocity_grad[q][i][j]+velocity_grad[q][j][i]);
			}
			P_stoke[q][i][i]+=-Pressure[q]*1.0e-9;
		}
		const Tensor<1,dim> normal=activeMaterial_fe_face_values.normal_vector(q);
		for (unsigned int i=0; i<dim; ++i){
			stress[q][i]=0;
			for (unsigned int j=0; j<dim; ++j){
				stress[q][i]+=P_stoke[q][i][j]*normal[j];
			}
		}		
	}
				
	ResidualEq->residualForNeummanBC(activeMaterial_fe_values, activeMaterial_fe_face_values, u_dof, R, stress);	
	ResidualEq->residualForNeummanBC(activeMaterial_fe_values, activeMaterial_fe_face_values, c_li_dof, R, defMap_surface, j_interface);
	ResidualEq->residualForNeummanBC(activeMaterial_fe_values, activeMaterial_fe_face_values, phi_s_dof, R, defMap_surface, Grad_phi);
	ResidualEq->residualForNeummanBC(activeMaterial_fe_values, activeMaterial_fe_face_values, T_dof, R, defMap_surface, Q);
}

template class initBoundValProbs<1>;
template class initBoundValProbs<2>;
template class initBoundValProbs<3>;
