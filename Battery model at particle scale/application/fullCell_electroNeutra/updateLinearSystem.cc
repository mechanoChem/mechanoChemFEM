#include "initBoundValProbs.h"

template <int dim>
void initBoundValProbs<dim>::updateLinearSystem()
{	
	params->enter_subsection("Geometry");
	double currentCollector_Y1=params->get_double("currentCollector_Y1");
	double currentCollector_Y2=params->get_double("currentCollector_Y2");
	int particle_num=params->get_integer("particle_number");
	params->leave_subsection();		
	
	params->enter_subsection("Initial condition");
  double c_li_max_neg=params->get_double("c_li_max_neg");
	double c_li_max_pos=params->get_double("c_li_max_pos");
	double c_li_100_neg=params->get_double("c_li_100_neg");
  double c_li_100_pos=params->get_double("c_li_100_pos");
	double c_li_plus_ini=params->get_double("c_li_plus_ini");
	double T_ini=params->get_double("T_0");
  params->leave_subsection();	
	
	params->enter_subsection("Elasticity");
	double youngsModulus_Al=params->get_double("youngsModulus_Al");
	double youngsModulus_Cu=params->get_double("youngsModulus_Cu");	
	double youngsModulus_binder=params->get_double("youngsModulus_binder");	
	double youngsModulus_s_neg=params->get_double("youngsModulus_s_neg");
	double youngsModulus_s_pos=params->get_double("youngsModulus_s_pos");
	double youngsModulus_sep=params->get_double("youngsModulus_sep");
	
	double nu_Al=params->get_double("nu_Al");
	double nu_Cu=params->get_double("nu_Cu");
	double nu_binder=params->get_double("nu_binder");
	double nu_sep=params->get_double("nu_sep");
	double nu_s_neg=params->get_double("nu_s_neg");
	double nu_s_pos=params->get_double("nu_s_pos");
	
	double omega_t_s_neg=params->get_double("omega_t_s_neg");
	double omega_t_s_pos=params->get_double("omega_t_s_pos");
	double omega_t_sep=params->get_double("omega_t_sep");
	double omega_t_Al=params->get_double("omega_t_Al");
	double omega_t_Cu=params->get_double("omega_t_Cu");
	double omega_t_binder=params->get_double("omega_t_binder");
	params->leave_subsection();	
	
	params->enter_subsection("Fluid");
	double viscosity=params->get_double("viscosity");
	double youngsModulus_mesh=params->get_double("youngsModulus_mesh");
	double nu_mesh=params->get_double("nu_mesh"); //cite 4
	
	params->leave_subsection();	
	
	params->enter_subsection("ElectroChemo");
	double F=params->get_double("F");
	double Rr=params->get_double("Rr");
	double sigma_s_neg=params->get_double("sigma_s_neg");
	double sigma_s_pos=params->get_double("sigma_s_pos");
	double sigma_binder=params->get_double("sigma_binder");
	double sigma_Al=params->get_double("sigma_Al");
	double sigma_Cu=params->get_double("sigma_Cu");
	double t_0=params->get_double("t_0");
  double D_li_neg=params->get_double("D_li_neg");
	double D_li_pos=params->get_double("D_li_pos");
	params->leave_subsection();	
	
	params->enter_subsection("Thermal");
	double lambda_s_neg=params->get_double("lambda_s_neg");
	double lambda_s_pos=params->get_double("lambda_s_pos");
	double lambda_e=params->get_double("lambda_e");
	double lambda_sep=params->get_double("lambda_sep");
	double lambda_Al=params->get_double("lambda_Al");
	double lambda_Cu=params->get_double("lambda_Cu");
	double lambda_binder=params->get_double("lambda_binder");
	
	double density_s_neg=params->get_double("density_s_neg");
	double density_s_pos=params->get_double("density_s_pos");
	double density_sep=params->get_double("density_sep");
	double density_e=params->get_double("density_e");
	double density_Al=params->get_double("density_Al");
	double density_Cu=params->get_double("density_Cu");
	double density_binder=params->get_double("density_binder");
	
	double Cp_s_neg=params->get_double("Cp_s_neg");
	double Cp_s_pos=params->get_double("Cp_s_pos");
	double Cp_e=params->get_double("Cp_e");
	double Cp_sep=params->get_double("Cp_sep");
	double Cp_Al=params->get_double("Cp_Al");
	double Cp_Cu=params->get_double("Cp_Cu");
	double Cp_binder=params->get_double("Cp_binder");
	
	double h=params->get_double("h");
	params->leave_subsection();	
	
	std::vector<double> local_RVEdata(particle_num*5, 0);
	//initialize 
	this->reinitLinearSystem();
  hp::FEValues<dim> hp_fe_values (fe_collection, q_collection, update_values | update_quadrature_points  | update_JxW_values | update_gradients);  
  FEFaceValues<dim> activeMaterial_fe_face_values (*fe_system[active_material_fe], *common_face_quadrature, update_values | update_quadrature_points | update_JxW_values | update_normal_vectors | update_gradients);
  FEFaceValues<dim> currentCollector_fe_face_values (*fe_system[current_collector_fe], *common_face_quadrature, update_values | update_quadrature_points | update_JxW_values | update_normal_vectors | update_gradients);
	
  FullMatrix<double> local_matrix;
  Vector<double>            local_rhs;
	local_matrix = 0; local_rhs = 0; 
  std::vector<types::global_dof_index> local_dof_indices;
	
  //loop over cells
  PETScWrappers::Vector localized_U(solution);
  PETScWrappers::Vector localized_Un(solution_prev);
	
	ResidualEq->dt=current_dt;
  typename hp::DoFHandler<dim>::active_cell_iterator cell = hpFEM<dim>::dof_handler.begin_active(), endc=hpFEM<dim>::dof_handler.end();
  for (;cell!=endc; ++cell){
		if (cell->subdomain_id() == this_mpi_process){	
			const Point<dim> cell_center = cell->center();	
			hp_fe_values.reinit (cell);
			const unsigned int dofs_per_cell = cell->get_fe().dofs_per_cell;
    	const FEValues<dim> &fe_values = hp_fe_values.get_present_fe_values();
    	unsigned int n_q_points= fe_values.n_quadrature_points;
			
    	local_matrix.reinit (dofs_per_cell, dofs_per_cell);
    	local_rhs.reinit (dofs_per_cell);
    	local_dof_indices.resize (dofs_per_cell);
			cell->get_dof_indices (local_dof_indices);
			
    	//AD variables
    	Table<1, Sacado::Fad::DFad<double> > ULocal(dofs_per_cell);
			Table<1, double > ULocalConv(dofs_per_cell);
			Table<1, double > U0Local(dofs_per_cell);
			
	  	for (unsigned int i=0; i<dofs_per_cell; ++i){
				if (std::abs(localized_U(local_dof_indices[i]))<1.0e-16) ULocal[i]=0.0;
				else{ULocal[i]=localized_U(local_dof_indices[i]);}
				ULocal[i].diff (i, dofs_per_cell);
				ULocalConv[i]= localized_Un(local_dof_indices[i]);
	  	}
		
    	Table<1, Sacado::Fad::DFad<double> > R(dofs_per_cell); 
			for(unsigned int i=0;i<dofs_per_cell;i++) R[i]=0.0;
			/*
			*evaluate primary fields
			*/
 	 		deformationMap<Sacado::Fad::DFad<double>, dim> defMap(n_q_points); //displacement fields
			//deformationMap<Sacado::Fad::DFad<double>, dim> defMap_mesh(n_q_points); //displacement fields of mesh for electrolyte

			dealii::Table<2,Sacado::Fad::DFad<double> > velocity(n_q_points, dim);
			dealii::Table<3,Sacado::Fad::DFad<double> > velocity_grad(n_q_points, dim,dim);
			dealii::Table<1,double>  c_li_conv(n_q_points),c_li_plus_conv(n_q_points),  T_conv(n_q_points), phi_s_conv(n_q_points), phi_e_conv(n_q_points);
  		dealii::Table<1,Sacado::Fad::DFad<double> > Pressure(n_q_points), c_li(n_q_points),c_li_plus(n_q_points),T(n_q_points),  phi_s(n_q_points), phi_e(n_q_points);
  		dealii::Table<2,Sacado::Fad::DFad<double> >  c_li_grad(n_q_points, dim), c_li_plus_grad(n_q_points, dim),T_grad(n_q_points, dim),  phi_s_grad(n_q_points, dim),phi_e_grad(n_q_points, dim);
			
			dealii::Table<2,Sacado::Fad::DFad<double> > u_mesh(n_q_points, dim); 
			dealii::Table<3,Sacado::Fad::DFad<double> > u_mesh_grad(n_q_points,dim, dim); 
			
			dealii::Table<2,Sacado::Fad::DFad<double> > Fc(n_q_points, dim), Ft(n_q_points, dim);
			dealii::Table<1,Sacado::Fad::DFad<double>> D_li(n_q_points), D_li_plus(n_q_points),sigma_e(n_q_points);
			dealii::Table<1,Sacado::Fad::DFad<double> > Q_ohm(n_q_points);
			double sigma_s, lambda, density, Cp, omega_t, c_li_max, youngsModulus, nu;

			if(cell->material_id()==electrolyte_id) {
				getDeformationMap<Sacado::Fad::DFad<double>, dim>(fe_values, u_m_dof, ULocal, defMap);
			}
			else {
				getDeformationMap<Sacado::Fad::DFad<double>, dim>(fe_values, u_dof, ULocal, defMap);
			}
			
			evaluateVectorFunction<Sacado::Fad::DFad<double>,dim>(fe_values, v_dof, ULocal, velocity);
			evaluateVectorFunctionGradient<Sacado::Fad::DFad<double>,dim>(fe_values, v_dof, ULocal, velocity_grad, defMap);
			evaluateScalarFunction<Sacado::Fad::DFad<double>,dim>(fe_values, p_dof, ULocal, Pressure);
			
			evaluateScalarFunction<Sacado::Fad::DFad<double>,dim>(fe_values, c_li_dof, ULocal, c_li);
			evaluateScalarFunction<double,dim>(fe_values, c_li_dof, ULocalConv, c_li_conv);
			evaluateScalarFunctionGradient<Sacado::Fad::DFad<double>,dim>(fe_values, c_li_dof, ULocal, c_li_grad, defMap);
			
			evaluateScalarFunction<double,dim>(fe_values, phi_s_dof, ULocalConv, phi_s_conv);	
			evaluateScalarFunction<Sacado::Fad::DFad<double>,dim>(fe_values, phi_s_dof, ULocal, phi_s);
			evaluateScalarFunctionGradient<Sacado::Fad::DFad<double>,dim>(fe_values, phi_s_dof, ULocal, phi_s_grad, defMap);
			
			evaluateScalarFunction<double,dim>(fe_values, c_li_plus_dof, ULocalConv, c_li_plus_conv);	
			evaluateScalarFunction<Sacado::Fad::DFad<double>,dim>(fe_values, c_li_plus_dof, ULocal, c_li_plus);
			evaluateScalarFunctionGradient<Sacado::Fad::DFad<double>,dim>(fe_values, c_li_plus_dof, ULocal, c_li_plus_grad, defMap);
			
			evaluateScalarFunction<double,dim>(fe_values, phi_e_dof, ULocalConv, phi_e_conv);		
			evaluateScalarFunction<Sacado::Fad::DFad<double>,dim>(fe_values, phi_e_dof, ULocal, phi_e);
			evaluateScalarFunctionGradient<Sacado::Fad::DFad<double>,dim>(fe_values, phi_e_dof, ULocal, phi_e_grad, defMap);
			
			evaluateScalarFunction<double,dim>(fe_values, T_dof, ULocalConv, T_conv);
			evaluateScalarFunction<Sacado::Fad::DFad<double>,dim>(fe_values, T_dof, ULocal, T);
			evaluateScalarFunctionGradient<Sacado::Fad::DFad<double>,dim>(fe_values, T_dof, ULocal, T_grad, defMap);
			
			evaluateVectorFunction<Sacado::Fad::DFad<double>,dim>(fe_values, u_m_dof, ULocal, u_mesh);
			
			evaluateVectorFunctionGradient<Sacado::Fad::DFad<double>,dim>(fe_values, u_m_dof, ULocal, u_mesh_grad,defMap);
			
			//initialization
			for(unsigned int q=0; q<n_q_points;q++) Q_ohm[q]=0.0;
			
			if(cell->material_id()==active_material_id )	{			
		 		if (cell_center[1]<=electrode_Y1){
					for(unsigned int q=0; q<n_q_points;q++) D_li[q]=D_li_neg;
					sigma_s=sigma_s_neg;
					c_li_max=c_li_max_neg;
					lambda=lambda_s_neg;
					density=density_s_neg;
					Cp=Cp_s_neg;
					omega_t=omega_t_s_neg;
					youngsModulus=youngsModulus_s_neg;
					nu=nu_s_neg;
				}
				else if(cell_center[1]>=electrode_Y2){
					for(unsigned int q=0; q<n_q_points;q++) D_li[q]=D_li_pos;
					sigma_s=sigma_s_pos;
					c_li_max=c_li_max_pos;
					lambda=lambda_s_pos;
					density=density_s_pos;
					Cp=Cp_s_neg;
					omega_t=omega_t_s_pos;
					youngsModulus=youngsModulus_s_pos;
					nu=nu_s_pos;		
				}
			}
			else if(cell->material_id()==electrolyte_id )	{
				for(unsigned int q=0; q<n_q_points;q++){
					sigma_e[q]=electricChemoFormula->formula_conductivity_e(T[q],c_li_plus[q], 1);
					D_li_plus[q]=electricChemoFormula->formula_diffusivity_e(T[q],c_li_plus[q], 1);
				}
				lambda=lambda_e;
				density=density_e;
				Cp=Cp_e;
				omega_t=0;
				ResidualEq->viscosity=viscosity;
				ResidualEq->density=density_e;
				 
				youngsModulus=youngsModulus_mesh;
				nu=nu_mesh;
			}
			else if(cell->material_id()==current_collector_id )	{
				if (cell_center[1]<=electrode_Y1){
					sigma_s=sigma_Al;
					lambda=lambda_Al;
					density=density_Al;
					Cp=Cp_Al;
					omega_t=omega_t_Al;
					youngsModulus=youngsModulus_Al;
					nu=nu_Al;
				}
				else if (cell_center[1]>=electrode_Y2){
					sigma_s=sigma_Cu;
					lambda=lambda_Cu;
					density=density_Cu;
					Cp=Cp_Cu;
					omega_t=omega_t_Cu;
					youngsModulus=youngsModulus_Cu;
					nu=nu_Cu;
				}
				else {
					sigma_s=sigma_binder;
					lambda=lambda_binder;
					density=density_binder;
					Cp=Cp_binder;
					omega_t=omega_t_binder;
					youngsModulus=youngsModulus_binder;
					nu=nu_binder;
				}
			}
			else if(cell->material_id()==solid_id ){
				lambda=lambda_sep;
				density=density_sep;
				Cp=Cp_sep;
				omega_t=omega_t_sep;
				youngsModulus=youngsModulus_sep;
				nu=nu_sep;
			}

     //==================================================================================================================
			//elasticity
			ResidualEq->setLameParametersByYoungsModulusPoissonRatio(youngsModulus, nu);		
			if (cell->material_id()==active_material_id and cell_center[1]<=electrode_Y1+3){
				for(unsigned int q=0; q<n_q_points;q++){
				  Sacado::Fad::DFad<double>  omega_c_s=electricChemoFormula->solid_particle_expansion(c_li[q]/c_li_max, 1);
				  if (cell_center[1]<=electrode_Y1+3) omega_c_s-=0.0939;
					for(unsigned int i=0; i<dim;i++){
					  Fc[q][i]=std::pow(1+omega_c_s,1.0/3.0);
					  Ft[q][i]=std::pow(omega_t*(T[q]-T_ini)+1,1.0/3.0);
					}
				}
			}
			else {
				for(unsigned int q=0; q<n_q_points;q++){
					for(unsigned int i=0; i<dim;i++){
						Fc[q][i]=1;
						if (cell->material_id()==electrolyte_id ){Ft[q][i]=1;}
						else {Ft[q][i]=std::pow(omega_t*(T[q]-T_ini)+1,1.0/3.0);}
					}
				}
			}
			dealii::Table<3, Sacado::Fad::DFad<double> > P_stress(n_q_points,dim,dim);
			dealii::Table<3, Sacado::Fad::DFad<double> > Fe(n_q_points,dim,dim);
			dealii::Table<3, Sacado::Fad::DFad<double> > E(n_q_points,dim,dim);
			
     	for (unsigned int q=0; q<n_q_points; ++q){   
  			for (unsigned int i=0; i<dim; ++i){
    			for (unsigned int j=0; j<dim; ++j){
  					Fe[q][i][j]=defMap.F[q][i][j]/Fc[q][i]/Ft[q][i];;
  				}
				}
			}			
			bool infinitesimal_strain_indicator=false;
			int elasticity_dof=u_dof;
			if (cell->material_id()==electrolyte_id ) { bool infinitesimal_strain_indicator=true; elasticity_dof=u_m_dof;}
			
			ResidualEq->evaluateStrain(Fe, E, defMap, infinitesimal_strain_indicator);
			ResidualEq->evaluateSaint_Venant_KirchhoffStress(P_stress,Fe, E);
			ResidualEq->residualForMechanics(fe_values, elasticity_dof, R, P_stress);			
			
			//ElectroChemo
			if (cell->material_id()==active_material_id or cell->material_id()==current_collector_id ) {
				dealii::Table<1,Sacado::Fad::DFad<double> > RHS_phi_s(n_q_points);
				dealii::Table<2,Sacado::Fad::DFad<double> > j_li(n_q_points, dim), i_phi_s(n_q_points, dim);
      	for(unsigned int q=0; q<n_q_points;q++){
					for(unsigned int i=0; i<dim;i++){
						if(cell->material_id()==active_material_id ) j_li[q][i]=-D_li[q]*c_li_grad[q][i];
						i_phi_s[q][i]=-sigma_s*phi_s_grad[q][i];	
						Q_ohm[q]+=-i_phi_s[q][i]*phi_s_grad[q][i];
					}		
					RHS_phi_s[q]=0;
      	}

			  if (cell->material_id()==active_material_id ) ResidualEq->residualForDiffusionEq(fe_values, c_li_dof, R, defMap, c_li, c_li_conv, j_li);
				ResidualEq->residualForPoissonEq(fe_values, phi_s_dof, R,defMap, i_phi_s, RHS_phi_s);
		
			} 
			else if(cell->material_id()==electrolyte_id ){
				dealii::Table<1,Sacado::Fad::DFad<double> > RHS_phi_e(n_q_points), divVelocity(n_q_points);
				dealii::Table<2,Sacado::Fad::DFad<double> > j_li_plus(n_q_points, dim), i_phi_e(n_q_points, dim);
      	for(unsigned int q=0; q<n_q_points;q++){
					divVelocity[q]=0;
					for(unsigned int i=0; i<dim;i++){
						divVelocity[q]+=velocity_grad[q][i][i];
						i_phi_e[q][i]=-sigma_e[q]*phi_e_grad[q][i]+2*Rr*T[q]/F*sigma_e[q]*(1-t_0)/c_li_plus[q]*c_li_plus_grad[q][i];  
						j_li_plus[q][i]=-D_li_plus[q]*c_li_plus_grad[q][i]+t_0/F*i_phi_e[q][i]+c_li_plus[q]*velocity[q][i];
						Q_ohm[q]+=-i_phi_e[q][i]*phi_e_grad[q][i];
					}	
					RHS_phi_e[q]=0;
      	}
					
				ResidualEq->residualForDiffusionEq(fe_values, c_li_plus_dof, R, defMap, c_li_plus, c_li_plus_conv, j_li_plus);
				ResidualEq->residualForPoissonEq(fe_values, phi_e_dof, R, defMap, i_phi_e, RHS_phi_e);
				ResidualEq->residualForStokesEq(fe_values, v_dof, R,  defMap, velocity_grad, Pressure);
				ResidualEq->residualForContinuityEq(fe_values, p_dof, R, defMap, divVelocity);				
			}		
			
			//thermal
			dealii::Table<1,double>  LHS_T_conv(n_q_points);
			dealii::Table<1,Sacado::Fad::DFad<double> > LHS_T(n_q_points);
			dealii::Table<2,Sacado::Fad::DFad<double> > flux_T(n_q_points, dim);
			for(unsigned int q=0; q<n_q_points;q++){
			  //Q_ohm[q]=0;
				LHS_T[q]=density*Cp*T[q];
				LHS_T_conv[q]=density*Cp*T_conv[q];
				for(unsigned int i=0; i<dim;i++){
					flux_T[q][i]=-lambda*T_grad[q][i];
				}
			}
			ResidualEq->residualForDiff_ReacEq(fe_values, T_dof, R, defMap, LHS_T, LHS_T_conv, flux_T,Q_ohm);
			
			//heat dissipation and out_current
			for (unsigned int faceID=0; faceID<2*dim; faceID++){
				if(cell->face(faceID)->boundary_id()==dim*2 ){
					currentCollector_fe_face_values.reinit (cell, faceID);
					double current;
					if(cell_center[1]<=electrode_Y1) current=current_IpA;
					else current=-current_IpA;
					ResidualEq->residualForNeummanBC(fe_values, currentCollector_fe_face_values, phi_s_dof, R, current);
				}
				if(cell->face(faceID)->boundary_id()==dim or  cell->face(faceID)->boundary_id()==dim*2 ){
					currentCollector_fe_face_values.reinit (cell, faceID);
					const unsigned int n_face_quadrature_points = currentCollector_fe_face_values.n_quadrature_points;
					dealii::Table<1,Sacado::Fad::DFad<double> > T_face(n_face_quadrature_points),heat_transfer(n_face_quadrature_points);
					evaluateScalarFunction<Sacado::Fad::DFad<double>,dim>(fe_values, currentCollector_fe_face_values, T_dof, ULocal, T_face);
					for(unsigned int q=0;q<n_face_quadrature_points;q++) heat_transfer[q]=h*(T_face[q]-T_ini);
					ResidualEq->residualForNeummanBC(fe_values, currentCollector_fe_face_values, T_dof, R, heat_transfer);
				}
			}	
			
		 //==================================================================================================================
			//interface at activeMaterial/electrolyte
			if(cell->material_id()==active_material_id  ){
				for (unsigned int f=0; f<GeometryInfo<dim>::faces_per_cell; ++f){
					if (cell->at_boundary(f) == false){
						if(cell->neighbor(f)->material_id()==electrolyte_id  and cell->neighbor(f)->has_children() == false){
						  assemble_interface_activeMaterial_electrolyte_term(cell, f, ULocal, ULocalConv, localized_U, R);                                          
						}
						else if (cell->neighbor(f)->has_children() == true){
							for (unsigned int sf=0; sf<cell->face(f)->n_children(); ++sf){
								if (cell->neighbor_child_on_subface(f, sf)->material_id()==electrolyte_id ){
								  assemble_interface_activeMaterial_electrolyte_term(cell, f, ULocal, ULocalConv, localized_U, R);                                                                
									break;
								}
							}
						}
					}
				}
			}
			//interface at electrolyte/activeMaterial+separator
			if(cell->material_id()==electrolyte_id ){
	
				for (unsigned int f=0; f<GeometryInfo<dim>::faces_per_cell; ++f){
				  bool interface_electrode=false;
					if (cell->at_boundary(f) == false){
						if(cell->neighbor(f)->material_id()==active_material_id  and cell->neighbor(f)->has_children() == false){
						  assemble_interface_electrolyte_activeMaterial_term(cell, f, ULocal, ULocalConv, localized_U, R);   
							interface_electrode=true;                                       
						}
						else if (cell->neighbor(f)->has_children() == true){
							for (unsigned int sf=0; sf<cell->face(f)->n_children(); ++sf){
								if (cell->neighbor_child_on_subface(f, sf)->material_id()==active_material_id ){
								  assemble_interface_electrolyte_activeMaterial_term(cell, f, ULocal, ULocalConv, localized_U, R);  
									interface_electrode=true;                                                              
									break;
								}
							}
						}
					}	
				}
			}
 		 
			//interface at solid/electrolyte
			if(cell->material_id()==solid_id ){
				for (unsigned int f=0; f<GeometryInfo<dim>::faces_per_cell; ++f){
					if (cell->at_boundary(f) == false){
						if(cell->neighbor(f)->material_id()==electrolyte_id and cell->neighbor(f)->has_children() == false){
						  assemble_interface_solid_electrolyte_term(cell, f, ULocal, ULocalConv, localized_U, R);                                                                      
						}
						else if (cell->neighbor(f)->has_children() == true){
							for (unsigned int sf=0; sf<cell->face(f)->n_children(); ++sf){
								if (cell->neighbor_child_on_subface(f, sf)->material_id()==electrolyte_id ){
									assemble_interface_solid_electrolyte_term(cell, f, ULocal, ULocalConv, localized_U, R);                                                                
									break;
								}
							}
						}
					}
				}
			}
			//interface at currentCollector(binder)/electrolyte
			if(cell->material_id()==current_collector_id ){
				for (unsigned int f=0; f<GeometryInfo<dim>::faces_per_cell; ++f){
					if (cell->at_boundary(f) == false){
						if(cell->neighbor(f)->material_id()==electrolyte_id and cell->neighbor(f)->has_children() == false){
						  assemble_interface_currentCollector_electrolyte_term(cell, f, ULocal, ULocalConv, localized_U, R);                                                                      
						}
						else if (cell->neighbor(f)->has_children() == true){
							for (unsigned int sf=0; sf<cell->face(f)->n_children(); ++sf){
								if (cell->neighbor_child_on_subface(f, sf)->material_id()==electrolyte_id ){
									assemble_interface_currentCollector_electrolyte_term(cell, f, ULocal, ULocalConv, localized_U, R);                                                                
									break;
								}
							}
						}
					}
				}
			}
			
			{//scalling
			  ResidualEq->scalling(fe_values,T_dof,R,1e-3);
			  ResidualEq->scalling(fe_values,phi_s_dof,R,1e-3);
			}
			//==================================================================================================================
		  //volume integration
			double cell_frac=ResidualEq->volumeIntegration(fe_values, defMap.detF);
      double C_li_ave=ResidualEq->volumeIntegration(fe_values, c_li);
			double C_li_plus_ave=ResidualEq->volumeIntegration(fe_values, c_li_plus);
			for(unsigned int i=0;i<particle_num;i++){
				double bar1=i*7.5+currentCollector_Y1;
				double bar2=(i-8)*7.5+electrode_Y2;
				if(i<8 and cell_center[1]>=bar1 and cell_center[1]<bar1+7.5) {
					if(cell->material_id()==active_material_id ){
						local_RVEdata[i]+=cell_frac;
						local_RVEdata[i+3*particle_num]+=C_li_ave;
						 break;
					}
					if(cell->material_id()==current_collector_id ){
						local_RVEdata[i+particle_num]+=cell_frac;
						 break;
					}
					
					if(cell->material_id()==electrolyte_id){
						local_RVEdata[i+2*particle_num]+=cell_frac;
						local_RVEdata[i+4*particle_num]+=C_li_plus_ave;
						 break;
					}
				}
				
				else if(i>=8 and i<14 and cell_center[1]>=bar2 and cell_center[1]<=bar2+7.5) {
					if(cell->material_id()==active_material_id ){
						local_RVEdata[i]+=cell_frac;
						local_RVEdata[i+3*particle_num]+=C_li_ave;
						 break;
					}
					if(cell->material_id()==current_collector_id ){
						local_RVEdata[i+particle_num]+=cell_frac;
						 break;
					}
					if(cell->material_id()==electrolyte_id ){
						local_RVEdata[i+2*particle_num]+=cell_frac;
						local_RVEdata[i+4*particle_num]+=C_li_plus_ave;
						 break;
					}		
				}
				else if(cell_center[1]>electrode_Y1 and cell_center[1]<electrode_Y2){
					if(cell->material_id()==electrolyte_id ){
						local_RVEdata[5*particle_num-1]+=C_li_plus_ave;
						local_RVEdata[3*particle_num-1]+=cell_frac;
						 break;
					}
					if(cell->material_id()==solid_id ){
						local_RVEdata[particle_num-1]+=cell_frac;
  					break;
					}
				}						
			}	
			//==================================================================================================================
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
	MPI_Reduce(&local_RVEdata[0], &RVEdata[0], 5*particle_num, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
}

template class initBoundValProbs<1>;
template class initBoundValProbs<2>;
template class initBoundValProbs<3>;