#include "initBoundValProbs.h"

template <int dim>
void initBoundValProbs<dim>::assemble_interface_solid_electrolyte_term(const typename hp::DoFHandler<dim>::active_cell_iterator &cell, unsigned int f,
																					 					 dealii::Table<1, Sacado::Fad::DFad<double> > &ULocal,
																					 					 Table<1, double > &ULocalConv,PETScWrappers::Vector &localized_U,
																					 					 dealii::Table<1, Sacado::Fad::DFad<double> >& R)
{	
	params->enter_subsection("Fluid");
	double viscosity=params->get_double("viscosity");
	params->leave_subsection();	
	
  const unsigned int solid_dofs_per_cell = fe_system[solid_fe]->dofs_per_cell;
  const unsigned int electrolyte_dofs_per_cell = fe_system[electrolyte_fe]->dofs_per_cell;
	
	hp::FEValues<dim> hp_fe_values (fe_collection, q_collection, update_values | update_quadrature_points  | update_JxW_values | update_gradients);
  FEFaceValues<dim> solid_fe_face_values (*fe_system[solid_fe], *common_face_quadrature, update_values | update_quadrature_points | update_JxW_values | update_normal_vectors | update_gradients);
  FEFaceValues<dim> electrolyte_fe_face_values (*fe_system[electrolyte_fe], *common_face_quadrature, update_values | update_quadrature_points | update_JxW_values | update_normal_vectors | update_gradients); 
	
	hp_fe_values.reinit (cell);
	const FEValues<dim> &solid_fe_values = hp_fe_values.get_present_fe_values();
	hp_fe_values.reinit (cell->neighbor(f));
	const FEValues<dim> &electrolyte_fe_values = hp_fe_values.get_present_fe_values();
	
	solid_fe_face_values.reinit(cell, f);
	electrolyte_fe_face_values.reinit(cell->neighbor(f), cell->neighbor_of_neighbor(f));
	std::vector<types::global_dof_index> solid_neighbor_dof_indices (electrolyte_dofs_per_cell);
	cell->neighbor(f)->get_dof_indices (solid_neighbor_dof_indices);
	
	//value of electrolyte cell
	Table<1, Sacado::Fad::DFad<double> > ULocal_electrolyte(electrolyte_dofs_per_cell);
  for (unsigned int i=0; i<electrolyte_dofs_per_cell; ++i){
		if (std::abs(localized_U(solid_neighbor_dof_indices[i]))<1.0e-16) ULocal_electrolyte[i]=0.0;
		else{ULocal_electrolyte[i]=localized_U(solid_neighbor_dof_indices[i]);}
  }

	/*
	*the following feilds are on surface
	*/
  const unsigned int n_face_quadrature_points = solid_fe_face_values.n_quadrature_points;
	deformationMap<Sacado::Fad::DFad<double>, dim> defMap_mesh_surface(n_face_quadrature_points); 
	getDeformationMap<Sacado::Fad::DFad<double>, dim>(electrolyte_fe_values, electrolyte_fe_face_values, u_m_dof, ULocal_electrolyte, defMap_mesh_surface);
		
  dealii::Table<1,Sacado::Fad::DFad<double> > Pressure(n_face_quadrature_points);
  dealii::Table<2,Sacado::Fad::DFad<double> > stress(n_face_quadrature_points,dim);
	dealii::Table<3,Sacado::Fad::DFad<double> > velocity_grad(n_face_quadrature_points, dim,dim),P_stoke(n_face_quadrature_points,dim,dim);
	// initialize variables
	
	//evaluateScalarFunction<double,dim>(fe_values, 3, ULocalConv, c_li_conv);
	evaluateVectorFunctionGradient<Sacado::Fad::DFad<double>,dim>(electrolyte_fe_values, electrolyte_fe_face_values, v_dof, ULocal_electrolyte, velocity_grad,defMap_mesh_surface);
	evaluateScalarFunction<Sacado::Fad::DFad<double>,dim>(electrolyte_fe_values, electrolyte_fe_face_values, p_dof, ULocal_electrolyte, Pressure);
	
	for(unsigned int q=0;q<n_face_quadrature_points;q++){
		//Stokes's stress constitutive equation
		for (unsigned int i=0; i<dim; ++i){
	  	for (unsigned int j=0; j<dim; ++j){
				P_stoke[q][i][j]=viscosity*1.0e-9*(velocity_grad[q][i][j]+velocity_grad[q][j][i]);
			}
			P_stoke[q][i][i]+=-Pressure[q]*1.0e-9;
		}
		const Tensor<1,dim> normal=solid_fe_face_values.normal_vector(q);
		for (unsigned int i=0; i<dim; ++i){
			stress[q][i]=0;
			for (unsigned int j=0; j<dim; ++j){
				stress[q][i]+=P_stoke[q][i][j]*normal[j];
			}
		}		
	}			
	ResidualEq->residualForNeummanBC(solid_fe_values, solid_fe_face_values, u_dof, R, stress);	
}


template <int dim>
void initBoundValProbs<dim>::assemble_interface_currentCollector_electrolyte_term(const typename hp::DoFHandler<dim>::active_cell_iterator &cell, unsigned int f,
																					 					 dealii::Table<1, Sacado::Fad::DFad<double> > &ULocal,
																					 					 Table<1, double > &ULocalConv,PETScWrappers::Vector &localized_U,
																					 					 dealii::Table<1, Sacado::Fad::DFad<double> >& R)
{
	
	params->enter_subsection("Fluid");
	double viscosity=params->get_double("viscosity");
	params->leave_subsection();	
	
  const unsigned int currentCollector_dofs_per_cell = fe_system[current_collector_fe]->dofs_per_cell;
  const unsigned int electrolyte_dofs_per_cell = fe_system[electrolyte_fe]->dofs_per_cell;
	
	hp::FEValues<dim> hp_fe_values (fe_collection, q_collection, update_values | update_quadrature_points  | update_JxW_values | update_gradients);
  FEFaceValues<dim> currentCollector_fe_face_values (*fe_system[current_collector_fe], *common_face_quadrature, update_values | update_quadrature_points | update_JxW_values | update_normal_vectors | update_gradients);
  FEFaceValues<dim> electrolyte_fe_face_values (*fe_system[electrolyte_fe], *common_face_quadrature, update_values | update_quadrature_points | update_JxW_values | update_normal_vectors | update_gradients); 
	
	hp_fe_values.reinit (cell);
	const FEValues<dim> &currentCollector_fe_values = hp_fe_values.get_present_fe_values();
	hp_fe_values.reinit (cell->neighbor(f));
	const FEValues<dim> &electrolyte_fe_values = hp_fe_values.get_present_fe_values();
	
	currentCollector_fe_face_values.reinit(cell, f);
	electrolyte_fe_face_values.reinit(cell->neighbor(f), cell->neighbor_of_neighbor(f));
	std::vector<types::global_dof_index> currentCollector_neighbor_dof_indices (electrolyte_dofs_per_cell);
	cell->neighbor(f)->get_dof_indices (currentCollector_neighbor_dof_indices);
	
	//value of electrolyte cell
	Table<1, Sacado::Fad::DFad<double> > ULocal_electrolyte(electrolyte_dofs_per_cell);
  for (unsigned int i=0; i<electrolyte_dofs_per_cell; ++i){
		if (std::abs(localized_U(currentCollector_neighbor_dof_indices[i]))<1.0e-16) ULocal_electrolyte[i]=0.0;
		else{ULocal_electrolyte[i]=localized_U(currentCollector_neighbor_dof_indices[i]);}
  }

	/*
	*the following feilds are on surface
	*/
  const unsigned int n_face_quadrature_points = currentCollector_fe_face_values.n_quadrature_points;
	
	deformationMap<Sacado::Fad::DFad<double>, dim> defMap_mesh_surface(n_face_quadrature_points); 
	getDeformationMap<Sacado::Fad::DFad<double>, dim>(electrolyte_fe_values, electrolyte_fe_face_values, u_m_dof, ULocal_electrolyte, defMap_mesh_surface);
	
  dealii::Table<1,Sacado::Fad::DFad<double> > Pressure(n_face_quadrature_points);
  dealii::Table<2,Sacado::Fad::DFad<double> > stress(n_face_quadrature_points,dim);
	dealii::Table<3,Sacado::Fad::DFad<double> > velocity_grad(n_face_quadrature_points, dim,dim),P_stoke(n_face_quadrature_points,dim,dim);
	// initialize variables
	
	//evaluateScalarFunction<double,dim>(fe_values, 3, ULocalConv, c_li_conv);
	evaluateVectorFunctionGradient<Sacado::Fad::DFad<double>,dim>(electrolyte_fe_values, electrolyte_fe_face_values, v_dof, ULocal_electrolyte, velocity_grad,defMap_mesh_surface);
	evaluateScalarFunction<Sacado::Fad::DFad<double>,dim>(electrolyte_fe_values, electrolyte_fe_face_values, p_dof, ULocal_electrolyte, Pressure);
	
	for(unsigned int q=0;q<n_face_quadrature_points;q++){
		//Stokes's stress constitutive equation
		for (unsigned int i=0; i<dim; ++i){
	  	for (unsigned int j=0; j<dim; ++j){
				P_stoke[q][i][j]=viscosity*1.0e-9*(velocity_grad[q][i][j]+velocity_grad[q][j][i]);
			}
			P_stoke[q][i][i]+=-Pressure[q]*1.0e-9;
		}
		const Tensor<1,dim> normal=currentCollector_fe_face_values.normal_vector(q);
		for (unsigned int i=0; i<dim; ++i){
			stress[q][i]=0;
			for (unsigned int j=0; j<dim; ++j){
				stress[q][i]+=P_stoke[q][i][j]*normal[j];
			}
		}		
	}			
	ResidualEq->residualForNeummanBC(currentCollector_fe_values, currentCollector_fe_face_values, u_dof, R, stress);	
}

template class initBoundValProbs<1>;
template class initBoundValProbs<2>;
template class initBoundValProbs<3>;