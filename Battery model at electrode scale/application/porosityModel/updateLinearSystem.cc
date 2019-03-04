#include "initBoundValProbs.h"

template <int dim>
void initBoundValProbs<dim>::updateLinearSystem()
{	
	//initialize 
	params->enter_subsection("Initial condition");
  double c_li_max_neg=params->get_double("c_li_max_neg");
	double c_li_max_pos=params->get_double("c_li_max_pos");
	double c_li_100_neg=params->get_double("c_li_100_neg");
  double c_li_100_pos=params->get_double("c_li_100_pos");
	double c_li_plus_ini=params->get_double("c_li_plus_ini");
	
	double T_ini=params->get_double("T_0");
  params->leave_subsection();	
	
	params->enter_subsection("ElectroChemo");
	double F=params->get_double("F");
	double Rr=params->get_double("Rr");
	double sigma_neg=params->get_double("sigma_neg");
	double sigma_pos=params->get_double("sigma_pos");
	double t_0=params->get_double("t_0");
  double D_li_neg=params->get_double("D_li_neg");
	double D_li_pos=params->get_double("D_li_pos");
	
	double eps_s_0_neg=params->get_double("eps_s_0_neg");
	double eps_s_0_pos=params->get_double("eps_s_0_pos");
	double eps_s_0_sep=params->get_double("eps_s_0_sep");
	double eps_b_0_neg=params->get_double("eps_b_0_neg");
	double eps_b_0_pos=params->get_double("eps_b_0_pos");
	double eps_b_0_sep=params->get_double("eps_b_0_sep");
	
	double R_s_0_neg=params->get_double("R_s_0_neg");
	double R_s_0_pos=params->get_double("R_s_0_pos");
	double R_s_0_sep=params->get_double("R_s_0_sep");
	params->leave_subsection();	
	
	params->enter_subsection("Elasticity");
	double youngsModulus_neg=params->get_double("youngsModulus_neg");
	double youngsModulus_pos=params->get_double("youngsModulus_pos");
	double youngsModulus_sep=params->get_double("youngsModulus_sep");
	
	double nu_sep=params->get_double("nu_sep");
	double nu_neg=params->get_double("nu_neg");
	double nu_pos=params->get_double("nu_pos");
	
	double kappa_sep=params->get_double("kappa_neg");
	double kappa_neg=params->get_double("kappa_pos");
	double kappa_pos=params->get_double("kappa_s");
	double kappa_s=params->get_double("kappa_sep");
	
	double Pl=params->get_double("pl");
	double Pb=params->get_double("pb");
	
	double omega_neg=params->get_double("omega_neg");
	double omega_pos=params->get_double("omega_pos");
	double omega_sep=params->get_double("omega_sep");
	double omega_s=params->get_double("omega_s");
	params->leave_subsection();	
	
	params->enter_subsection("Thermal");
	double lambda_neg=params->get_double("lambda_neg");
	double lambda_pos=params->get_double("lambda_pos");
	double lambda_sep=params->get_double("lambda_sep");
	
	double density_neg=params->get_double("density_neg");
	double density_pos=params->get_double("density_pos");
	double density_sep=params->get_double("density_sep");
	
	double Cp_s_neg=params->get_double("Cp_s_neg");
	double Cp_s_pos=params->get_double("Cp_s_pos");
	double Cp_sep=params->get_double("Cp_sep");
	
	double h=params->get_double("h");
	params->leave_subsection();	
	
	//set rhs to zero
	this->reinitLinearSystem();
  hp::FEValues<dim> hp_fe_values (fe_collection, q_collection, update_values | update_quadrature_points  | update_JxW_values | update_gradients);  
  FEFaceValues<dim> electrode_fe_face_values (*fe_system[electrode_fe], *common_face_quadrature, update_values | update_quadrature_points | update_JxW_values | update_normal_vectors | update_gradients);
  FEFaceValues<dim> separator_fe_face_values (*fe_system[separator_fe], *common_face_quadrature, update_values | update_quadrature_points | update_JxW_values | update_normal_vectors | update_gradients);
		
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
			getDeformationMap<Sacado::Fad::DFad<double>, dim>(fe_values, u_dof, ULocal, defMap);
	 	 	
			deformationMap<double, dim> defMapConv(n_q_points); 
	 	 	getDeformationMap<double, dim>(fe_values, u_dof, ULocalConv, defMapConv);
			
	 	 	dealii::Table<1,double> c_li_conv(n_q_points), c_li_plus_conv(n_q_points),T_conv(n_q_points), phi_s_conv(n_q_points), phi_e_conv(n_q_points);
	  	dealii::Table<1,Sacado::Fad::DFad<double> > c_li(n_q_points), c_li_plus(n_q_points), T(n_q_points), phi_s(n_q_points), phi_e(n_q_points);
	  	dealii::Table<2,Sacado::Fad::DFad<double> > c_li_plus_grad(n_q_points, dim), T_grad(n_q_points, dim),phi_s_grad(n_q_points, dim),phi_e_grad(n_q_points, dim);
			

			evaluateScalarFunction<double,dim>(fe_values, c_li_dof, ULocalConv, c_li_conv);
			evaluateScalarFunction<double,dim>(fe_values, c_li_plus_dof, ULocalConv, c_li_plus_conv);
			evaluateScalarFunction<double,dim>(fe_values, T_dof, ULocalConv, T_conv);
			evaluateScalarFunction<double,dim>(fe_values, phi_s_dof, ULocalConv, phi_s_conv);
			evaluateScalarFunction<double,dim>(fe_values, phi_e_dof, ULocalConv, phi_e_conv);
		
			evaluateScalarFunction<Sacado::Fad::DFad<double>,dim>(fe_values, c_li_dof, ULocal, c_li);
			evaluateScalarFunction<Sacado::Fad::DFad<double>,dim>(fe_values, c_li_plus_dof, ULocal, c_li_plus);
			evaluateScalarFunction<Sacado::Fad::DFad<double>,dim>(fe_values, T_dof, ULocal, T);
			evaluateScalarFunction<Sacado::Fad::DFad<double>,dim>(fe_values, phi_s_dof, ULocal, phi_s);
			evaluateScalarFunction<Sacado::Fad::DFad<double>,dim>(fe_values, phi_e_dof, ULocal, phi_e);
		
			evaluateScalarFunctionGradient<Sacado::Fad::DFad<double>,dim>(fe_values, c_li_plus_dof, ULocal, c_li_plus_grad, defMap);
			evaluateScalarFunctionGradient<Sacado::Fad::DFad<double>,dim>(fe_values, T_dof, ULocal, T_grad, defMap);
			evaluateScalarFunctionGradient<Sacado::Fad::DFad<double>,dim>(fe_values, phi_s_dof, ULocal, phi_s_grad, defMap);
			evaluateScalarFunctionGradient<Sacado::Fad::DFad<double>,dim>(fe_values, phi_e_dof, ULocal, phi_e_grad, defMap);
			//assign paramters
	    dealii::Table<1,Sacado::Fad::DFad<double> > sigma_e(n_q_points), D_li(n_q_points),D_li_plus(n_q_points);
			int domain =0;			
			double sigma, c_li_max, youngsModulus, nu, omega, lambda, density, Cp, kappa, eps_s_0, eps_b_0, R_s_0;

			if(cell->material_id()==electrode_id )	{			
		 		if (cell_center[1]<=electrode_Y1){
					domain=-1;
					for(unsigned int q=0; q<n_q_points;q++) D_li[q]=D_li_neg;
					sigma=sigma_neg;
					c_li_max=c_li_max_neg;
					lambda=lambda_neg;
					density=density_neg;
					Cp=Cp_s_neg;
					omega=omega_neg;
					youngsModulus=youngsModulus_neg;
					nu=nu_neg;
					kappa=kappa_neg;
					
					eps_s_0=eps_s_0_neg;
					eps_b_0=eps_b_0_neg;
					R_s_0=R_s_0_neg;
				}
				else if(cell_center[1]>=electrode_Y2){
					domain=1;
					for(unsigned int q=0; q<n_q_points;q++) D_li[q]=D_li_pos;
					sigma=sigma_pos;
					c_li_max=c_li_max_pos;
					lambda=lambda_pos;
					density=density_pos;
					Cp=Cp_s_neg;
					omega=omega_pos;
					youngsModulus=youngsModulus_pos;
					nu=nu_pos;		
					kappa=kappa_pos;
					
					eps_s_0=eps_s_0_pos;
					eps_b_0=eps_b_0_pos;
					R_s_0=R_s_0_pos;
				}
			}
			else if(cell->material_id()==separator_id )	{
				domain=0;
				lambda=lambda_sep;
				density=density_sep;
				Cp=Cp_sep;
				omega=omega_sep;
				youngsModulus=youngsModulus_sep;
				nu=nu_sep;
				kappa=kappa_sep;
				
				eps_s_0=eps_s_0_sep;
				eps_b_0=eps_b_0_sep;
				R_s_0=R_s_0_sep;
			}
			for(unsigned int q=0; q<n_q_points;q++){
				sigma_e[q]=electricChemoFormula->formula_conductivity_e(T[q],c_li_plus[q], 1);
				//sigma_e[q]=1e8;
				D_li_plus[q]=electricChemoFormula->formula_diffusivity_e(T[q],c_li_plus[q], 1);
			}
     //==================================================================================================================
			//elasticity
			ResidualEq->setLameParametersByYoungsModulusPoissonRatio(youngsModulus, nu);		
			//omega_c_electrode
			dealii::Table<3, Sacado::Fad::DFad<double> > P_stress(n_q_points,dim,dim);
			dealii::Table<3, Sacado::Fad::DFad<double> > Fe(n_q_points,dim,dim);
			dealii::Table<3, Sacado::Fad::DFad<double> > E(n_q_points,dim,dim);
			
	    dealii::Table<1,Sacado::Fad::DFad<double> > eps_l(n_q_points),eps_s(n_q_points), R_e(n_q_points);
			dealii::Table<1,double> eps_l_conv(n_q_points),eps_s_conv(n_q_points);
			
			Fe=defMap.F;
			for(unsigned int q=0; q<n_q_points;q++){
				Sacado::Fad::DFad<double>  beta, beta_s, beta_t, beta_t_s;
				double  beta_conv, beta_s_conv, beta_t_conv, beta_t_s_conv;
				
				if (cell->material_id()==electrode_id){
					beta_s=electricChemoFormula->solid_particle_expansion(c_li[q]/c_li_max, 1);
				 	beta=electricChemoFormula->electrode_expansion(c_li[q]/c_li_max, 1);
					
					beta_s_conv=electricChemoFormula->solid_particle_expansion(c_li_conv[q]/c_li_max, 1).val();
				 	beta_conv=electricChemoFormula->electrode_expansion(c_li_conv[q]/c_li_max, 1).val();
				  if (cell_center[1]<=electrode_Y1+3) {
						beta_s-=0.0939; beta_s_conv-=0.0939;
						beta-=0.019; beta_conv-=0.019;
					}
				}
				else {beta =0; beta=0; beta_s=0; beta_s_conv=0;}
				
				beta_t=(T[q]-T_ini)*omega;
				beta_t_s=(T[q]-T_ini)*omega_s;
				
				beta_t_conv=(T_conv[q]-T_ini)*omega;
				beta_t_s_conv=(T_conv[q]-T_ini)*omega_s;
				
				Sacado::Fad::DFad<double> factor=((kappa*(defMap.detF[q]/(1+beta)/(1+beta_t)-1)-Pl-Pb)/kappa_s+1)*(1+beta_s)*(1+beta_t_s);
				eps_s[q]=factor/defMap.detF[q]*eps_s_0;
				eps_l[q]=1-eps_s[q]-eps_b_0/defMap.detF[q];
				R_e[q]=std::pow(factor,1/3)*R_s_0;
				//volume fraction at previous time step
				double factor_conv=((kappa*(defMapConv.detF[q]/(1+beta_conv)/(1+beta_t_conv)-1)-Pl-Pb)/kappa_s+1)*(1+beta_s_conv)*(1+beta_t_s_conv);
				eps_s_conv[q]=factor_conv/defMapConv.detF[q]*eps_s_0;
				eps_l_conv[q]=1-eps_s_conv[q]-eps_b_0/defMapConv.detF[q];
				
				for(unsigned int i=0; i<dim;i++){
					Fe[q][i][i]=defMap.F[q][i][i]/std::pow(1+beta,1.0/3.0)/std::pow(beta_t+1,1.0/3.0);					 
				}
			}
			bool infinitesimal_strain_indicator=true;
			
			ResidualEq->evaluateStrain(Fe, E, defMap, infinitesimal_strain_indicator);
			ResidualEq->evaluateSaint_Venant_KirchhoffStress(P_stress,Fe, E);
			ResidualEq->residualForMechanics(fe_values, u_dof, R, P_stress);			
			//==================================================================================================================
			//ElectroChemo
	  	dealii::Table<1,Sacado::Fad::DFad<double> > jn(n_q_points);	
			dealii::Table<1,Sacado::Fad::DFad<double> > LHS_c_li(n_q_points), LHS_c_li_plus(n_q_points), RHS_phi_s(n_q_points),RHS_phi_e(n_q_points);
			dealii::Table<1,Sacado::Fad::DFad<double> > reaction_c_li(n_q_points), reaction_c_li_plus(n_q_points);
			dealii::Table<1,double> LHS_c_li_conv(n_q_points), LHS_c_li_plus_conv(n_q_points);
			dealii::Table<2,Sacado::Fad::DFad<double> > j_li(n_q_points, dim), j_li_plus(n_q_points, dim), i_phi_s(n_q_points, dim),i_phi_e(n_q_points, dim);

      for(unsigned int q=0; q<n_q_points;q++){
				if (cell->material_id()==electrode_id) jn[q]=3/R_e[q]*electricChemoFormula->formula_jn(T[q], c_li[q], c_li_plus[q], phi_s[q], phi_e[q], domain);
				else jn[q]=0;
				
				LHS_c_li[q]=eps_s[q]*c_li[q];
				LHS_c_li_conv[q]=eps_s_conv[q]*c_li_conv[q];
				reaction_c_li[q]=-eps_s[q]*jn[q];
				RHS_phi_s[q]=-F*jn[q];
				
				LHS_c_li_plus[q]=eps_l[q]*c_li_plus[q];
				LHS_c_li_plus_conv[q]=eps_l_conv[q]*c_li_plus_conv[q];
				reaction_c_li_plus[q]=(1-t_0)*eps_s[q]*jn[q];
				RHS_phi_e[q]=F*jn[q];

				for(unsigned int i=0; i<dim;i++){
					if(cell->material_id()==electrode_id ) j_li[q][i]=0;
					j_li_plus[q][i]=-std::pow(eps_l[q],0.5)*D_li_plus[q]*c_li_plus_grad[q][i];
					i_phi_s[q][i]=-sigma*eps_s[q]*phi_s_grad[q][i];
					i_phi_e[q][i]=sigma_e[q]*std::pow(eps_l[q],1.5)*(-phi_e_grad[q][i]-2*Rr*T[q]/F*(1-t_0)/c_li_plus[q]*c_li_plus_grad[q][i]);
					//i_phi_e[q][i]=sigma_e[q]*std::pow(eps_l[q],1.5)*(-phi_e_grad[q][i]);  
				}		
			}
			if (cell->material_id()==electrode_id ){
			 	ResidualEq->residualForDiff_ReacEq(fe_values, c_li_dof, R, defMap, LHS_c_li, LHS_c_li_conv, j_li,reaction_c_li);
				ResidualEq->residualForPoissonEq(fe_values, phi_s_dof, R,defMap, i_phi_s, RHS_phi_s);
			} 
			ResidualEq->residualForDiff_ReacEq(fe_values, c_li_plus_dof, R, defMap, LHS_c_li_plus, LHS_c_li_plus_conv, j_li_plus,reaction_c_li_plus);
			ResidualEq->residualForPoissonEq(fe_values, phi_e_dof, R, defMap, i_phi_e, RHS_phi_e);					
			
			//thermal
			dealii::Table<1,double>  LHS_T_conv(n_q_points);
			dealii::Table<1,Sacado::Fad::DFad<double> > LHS_T(n_q_points), Q_rxn(n_q_points),Q_rev(n_q_points),Q_ohm(n_q_points),Q(n_q_points);
			dealii::Table<2,Sacado::Fad::DFad<double> > flux_T(n_q_points, dim);
			for(unsigned int q=0; q<n_q_points;q++){
				Sacado::Fad::DFad<double> usc=electricChemoFormula->formula_Usc(c_li[q]/c_li_max,domain);
			  Q_rxn[q]=F*eps_s[q]*jn[q]*(phi_s[q]-phi_e[q]-usc);
				Q_rev[q]=F*jn[q]*T[q]*electricChemoFormula->formula_dUdt(c_li[q]/c_li_max);
				Q_ohm[q]=0;
				for(unsigned int i=0; i<dim;i++){
					flux_T[q][i]=-lambda*T_grad[q][i];
					Q_ohm[q]+=-i_phi_s[q][i]*phi_s_grad[q][i]-i_phi_e[q][i]*phi_e_grad[q][i];
				}
				Q[q]= Q_rxn[q]+Q_rev[q]+Q_ohm[q];
				LHS_T[q]=density*Cp*T[q];
				LHS_T_conv[q]=density*Cp*T_conv[q];
				
				Q[q]=0;
			}
			ResidualEq->residualForDiff_ReacEq(fe_values, T_dof, R, defMap, LHS_T, LHS_T_conv, flux_T,Q);
			
			//heat dissipation and out_current
			for (unsigned int faceID=0; faceID<2*dim; faceID++){
				if(cell->face(faceID)->at_boundary()){
					FEFaceValues<dim>* fe_face_values;
					if(domain==1 or domain==-1) fe_face_values=&electrode_fe_face_values;
					else if(domain==0) fe_face_values=&separator_fe_face_values;
					fe_face_values->reinit (cell, faceID);
					const unsigned int n_face_quadrature_points = fe_face_values->n_quadrature_points;
					dealii::Table<1,Sacado::Fad::DFad<double> > T_face(n_face_quadrature_points),heat_transfer(n_face_quadrature_points);
					evaluateScalarFunction<Sacado::Fad::DFad<double>,dim>(fe_values, *fe_face_values, T_dof, ULocal, T_face);
					for(unsigned int q=0;q<n_face_quadrature_points;q++) heat_transfer[q]=h*(T_face[q]-T_ini);
					ResidualEq->residualForNeummanBC(fe_values, *fe_face_values, T_dof, R, heat_transfer);
					if(cell->face(faceID)->boundary_id()==dim*2 ){
						double current;
						if(cell_center[1]<=electrode_Y1) current=current_IpA;
						else current=-current_IpA;
						ResidualEq->residualForNeummanBC(fe_values, *fe_face_values, phi_s_dof, R, current);
					}
				}
			}	
			
			{//scalling
				
			  ResidualEq->scalling(fe_values,c_li_plus_dof,R,1e-3);
				ResidualEq->scalling(fe_values,phi_e_dof,R,1e-3);
				ResidualEq->scalling(fe_values,c_li_dof,R,1e-5);
				ResidualEq->scalling(fe_values,phi_s_dof,R,1e-6);
			  ResidualEq->scalling(fe_values,T_dof,R,1e-5);
				
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
}

template class initBoundValProbs<1>;
template class initBoundValProbs<2>;
template class initBoundValProbs<3>;