#include "initBoundValProbs.h"

template <int dim>
void initBoundValProbs<dim>::updateLinearSystem()
{
	params->enter_subsection("Mechanics");
	double youngsModulus=params->get_double("youngsModulus");
	double poissonRatio=params->get_double("poissonRatio");
	
	double saturation_matID_Cortex=params->get_double("saturation_matID_Cortex");
	double saturation_matID_Subcortex=params->get_double("saturation_matID_Subcortex");
	std::string GROWTH=params->get("GROWTH");
	params->leave_subsection();	
	
	params->enter_subsection("Concentration");
	std::string advection_type=params->get("advection_type");
	double mobility_c1=params->get_double("mobility_c1");
	double mobility_c2=params->get_double("mobility_c2");
	double reac_10=params->get_double("reac_10");
	double reac_20=params->get_double("reac_20");
	double inward_flux=params->get_double("inward_flux");
	params->leave_subsection();	
	//initialize 
	this->reinitLinearSystem();
	
	std::vector<double> local_dataStack(3, 0);
	
	//set current_dt for Residual class
	ResidualEq->dt=current_dt;
  hp::FEValues<dim> hp_fe_values (fe_collection, q_collection, update_values | update_quadrature_points  | update_JxW_values | update_gradients);  
  FEFaceValues<dim> Cortex_fe_face_values (*fe_system[Cortex_fe], *common_face_quadrature, update_values | update_quadrature_points | update_JxW_values | update_normal_vectors | update_gradients);
  FEFaceValues<dim> Subcortex_fe_face_values (*fe_system[Subcortex_fe], *common_face_quadrature, update_values | update_quadrature_points | update_JxW_values | update_normal_vectors | update_gradients);

  hp::FEValues<dim> hp_fe_values_add (fe_collection_add, q_collection_add, update_values | update_quadrature_points  | update_JxW_values | update_gradients);  

  FullMatrix<double> local_matrix;
  Vector<double>            local_rhs;
	local_matrix = 0; local_rhs = 0; 
  std::vector<types::global_dof_index> local_dof_indices;
	std::vector<types::global_dof_index> local_dof_indices_add;
	
  //loop over cells
  PETScWrappers::Vector localized_U(solution);
  PETScWrappers::Vector localized_Un(solution_prev);
	PETScWrappers::Vector localized_U0(solution_0);
	PETScWrappers::Vector localized_U_add(additional_data);
	
	ResidualEq->dt=current_dt;
	ResidualEq->setLameParametersByYoungsModulusPoissonRatio(youngsModulus, poissonRatio);
		
  typename hp::DoFHandler<dim>::active_cell_iterator cell = hpFEM<dim>::dof_handler.begin_active(), endc=hpFEM<dim>::dof_handler.end();
  typename hp::DoFHandler<dim>::active_cell_iterator cell_add = dof_handler_add->begin_active();
	
  for (;cell!=endc; ++cell, ++cell_add){
		/*nothing for the third domain, doing this will eliminate redundant opearation and speed up the code,
		* need to comment out if we have something for the third domain
		*/
		if(cell->material_id()==Ventricle_id) continue;
		
		if (cell->subdomain_id() == this_mpi_process){
		  const Point<dim> cell_center = cell->center();	
			hp_fe_values.reinit (cell);
			hp_fe_values_add.reinit (cell_add);
			
			const unsigned int dofs_per_cell = cell->get_fe().dofs_per_cell;
    	const FEValues<dim> &fe_values = hp_fe_values.get_present_fe_values();
    	unsigned int n_q_points= fe_values.n_quadrature_points;
			
			const unsigned int dofs_per_cell_add = cell_add->get_fe().dofs_per_cell;
    	const FEValues<dim> &fe_values_add = hp_fe_values_add.get_present_fe_values();
			
    	local_matrix.reinit (dofs_per_cell, dofs_per_cell);
    	local_rhs.reinit (dofs_per_cell);
    	local_dof_indices.resize (dofs_per_cell);
			cell->get_dof_indices (local_dof_indices);
			
    	local_dof_indices_add.resize (dofs_per_cell_add);
			cell_add->get_dof_indices (local_dof_indices_add);
			
    	//AD variables
    	Table<1, Sacado::Fad::DFad<double> > ULocal(dofs_per_cell);
			Table<1, double > ULocalConv(dofs_per_cell);
			Table<1, double > U0Local(dofs_per_cell);
			
			Table<1, double > ULocal_add(dofs_per_cell_add);
			
	  	for (unsigned int i=0; i<dofs_per_cell; ++i){
				if (std::abs(localized_U(local_dof_indices[i]))<1.0e-16) ULocal[i]=0.0;
				else{ULocal[i]=localized_U(local_dof_indices[i]);}
				ULocal[i].diff (i, dofs_per_cell);
				ULocalConv[i]= localized_Un(local_dof_indices[i]);
				U0Local[i]= localized_U0(local_dof_indices[i]);
	  	}
			
	  	for (unsigned int i=0; i<dofs_per_cell_add; ++i){
				if (std::abs(localized_U_add(local_dof_indices_add[i]))<1.0e-16) ULocal_add[i]=0.0;
				else{ULocal_add[i]=localized_U_add(local_dof_indices_add[i]);}
	  	}
    	Table<1, Sacado::Fad::DFad<double> > R(dofs_per_cell); 
			for(unsigned int i=0;i<dofs_per_cell;i++) R[i]=0.0;
			/*
			*evaluate primary fields
			*/
			deformationMap<Sacado::Fad::DFad<double>, dim> defMap(n_q_points); 
			getDeformationMap<Sacado::Fad::DFad<double>, dim>(fe_values, u_dof, ULocal, defMap);
			
			dealii::Table<1,double>  c_1_conv(n_q_points), c_2_conv(n_q_points),c_1_0(n_q_points),  c_2_0(n_q_points);
			dealii::Table<1,Sacado::Fad::DFad<double> > c_1(n_q_points), c_2(n_q_points);
			dealii::Table<2,Sacado::Fad::DFad<double> >  c_1_grad(n_q_points, dim), c_2_grad(n_q_points, dim);
			dealii::Table<2,double > advection_direction(n_q_points, dim);
			dealii::Table<1,double>  c_1_inverse(n_q_points);
			//evaluateVectorFunction<double,dim>(fe_values_add, 0, ULocal_add, advection_direction);
			
			evaluateScalarFunction<double,dim>(fe_values, c1_dof, ULocalConv, c_1_conv);
			evaluateScalarFunction<double,dim>(fe_values, c1_dof, U0Local, c_1_0);
			evaluateScalarFunction<Sacado::Fad::DFad<double>,dim>(fe_values, c1_dof, ULocal, c_1);	
			evaluateScalarFunctionGradient<Sacado::Fad::DFad<double>,dim>(fe_values, c1_dof, ULocal, c_1_grad,defMap);
			evaluateScalarFunction<double,dim>(fe_values_add, dim, ULocal_add, c_1_inverse);
			evaluateScalarFunction<double,dim>(fe_values, c2_dof, ULocalConv, c_2_conv);
			evaluateScalarFunction<double,dim>(fe_values, c2_dof, U0Local, c_2_0);
			evaluateScalarFunction<Sacado::Fad::DFad<double>,dim>(fe_values, c2_dof, ULocal, c_2);	
			evaluateScalarFunctionGradient<Sacado::Fad::DFad<double>,dim>(fe_values, c2_dof, ULocal, c_2_grad,defMap);
		
			/*
			*evaluate general diffusion and reaction term
			*/
			dealii::Table<1,Sacado::Fad::DFad<double> > c_1_reaction(n_q_points), c_2_reaction(n_q_points);
			dealii::Table<2,Sacado::Fad::DFad<double> > c_1_flux(n_q_points, dim),c_2_flux(n_q_points, dim);
			
			//if(cell->material_id()==Cortex_id) reac_10=0;			
			//chemo formula
			
			if(std::strcmp(advection_type.c_str(),"fromFile")==0){ evaluateVectorFunction<double,dim>(fe_values_add, 0, ULocal_add, advection_direction);}
			else if(std::strcmp(advection_type.c_str(),"radial")==0){ 
				for(unsigned int q=0; q<n_q_points;q++){
					const Point<dim> posR = fe_values.quadrature_point(q);
					for(unsigned int i=0; i<dim;i++){
						advection_direction[q][i]=posR[i]/std::sqrt(posR.square());
					}
				}
			}
			else if(std::strcmp(advection_type.c_str(),"noAdvection")==0){
				for(unsigned int q=0; q<n_q_points;q++){
					for(unsigned int i=0; i<dim;i++){
						advection_direction[q][i]=0;
					}
				}
			}
    	for(unsigned int q=0; q<n_q_points;q++){
				const Point<dim> posR = fe_values.quadrature_point(q);
				for(unsigned int i=0; i<dim;i++){
					c_1_flux[q][i]=-mobility_c1*c_1_grad[q][i];
					c_2_flux[q][i]=-mobility_c2*c_2_grad[q][i];
				}		
				c_1_reaction[q]=reac_10;
				c_2_reaction[q]=reac_20;
				
    	}			
			if(std::abs(inward_flux)>1.0e-12 and cell->material_id()==Subcortex_id){
				for (unsigned int f=0; f<GeometryInfo<dim>::faces_per_cell; ++f){
					bool at_interface=false;
					if(cell->face(f)->at_boundary()==false ){
						bool at_interface=false;
			    	if(cell->neighbor(f)->material_id()==Ventricle_id and cell->neighbor(f)->has_children() == false){
			      	at_interface=true;
			    	}
			    	// check for the cell that has children cells
			    	else if (cell->neighbor(f)->has_children() == true){
			      	for (unsigned int sf=0; sf<cell->face(f)->n_children(); ++sf){
								if (cell->neighbor_child_on_subface(f, sf)->material_id()==Ventricle_id){
				  				at_interface=true;                                                              
				  				break;
								}
			      	}
			    	}
			  	}
			
			  	if(at_interface){
			    	Subcortex_fe_face_values.reinit (cell, f);
			    	ResidualEq->residualForNeummanBC(fe_values, Subcortex_fe_face_values, c1_dof, R, inward_flux);
			  	}
				}
			}							
			//mechanics
			dealii::Table<3, Sacado::Fad::DFad<double> > P(n_q_points,dim,dim);
			dealii::Table<3, Sacado::Fad::DFad<double> > Fe(n_q_points,dim,dim);
			dealii::Table<1,Sacado::Fad::DFad<double> > fac(n_q_points);
			double sat;			

		  if (cell->material_id()==Cortex_id ){
				sat=saturation_matID_Cortex;
		  }
		  else if (cell->material_id()==Subcortex_id )
			{
				sat=saturation_matID_Subcortex;
		  }
			else{}
			for(unsigned int q=0; q<n_q_points;q++){
				for (unsigned int i=0; i<dim; ++i){
					for (unsigned int j=0; j<dim; ++j){
						Fe[q][i][j]=0.0;
					}
				}
				if (c_1_0[q] > 0){
					if(std::strcmp(GROWTH.c_str(),"Uniform")==0 ){
						fac[q] = 1.0; //Uniform growth
					} 
					else if(std::strcmp(GROWTH.c_str(),"Isotropic")==0 ){
						fac[q]=std::pow((c_1[q]/c_1_0[q]), 1.0/3.0); // Isotropic Growth
					} 
					else if(std::strcmp(GROWTH.c_str(),"Tangential")==0 ){
						fac[q]=std::pow((c_1[q]/c_1_0[q]), 1.0/2.0); //Tangential
					} 
					else{pcout << "Growth type not supported\n"; exit(1);}
				}
				else {fac[q] = 1.0;}
				
		    if (fac[q] <= 1.0e-15){
					printf("*************Non positive growth factor*************. Value %12.4e\n", fac[q].val());
		    }
		    
      		    if (fac[q] < sat){ fac[q] = 1.0; }
		    else{ fac[q] /= sat; }
				
		    if(std::strcmp(GROWTH.c_str(),"Uniform")==0 ){Fe=defMap.F;}
				else if(std::strcmp(GROWTH.c_str(),"Isotropic")==0 ){
					for (unsigned int i=0; i<dim; ++i){
						for (unsigned int j=0; j<dim; ++j){
					  	Fe[q][i][j]=defMap.F[q][i][j]/fac[q];
						}
					}
				}
				else if(std::strcmp(GROWTH.c_str(),"Tangential")==0 ){
		      Table<1, Sacado::Fad::DFad<double> > FeR(dim);
		      for (unsigned int i=0; i<dim; ++i){
						FeR[i] = 0.0;
			  		for (unsigned int j=0; j<dim; ++j){
			    		FeR[i]+=defMap.F[q][i][j]*advection_direction[q][j];
		        }
		      }
		      for (unsigned int i=0; i<dim; ++i){
						for (unsigned int j=0; j<dim; ++j){
						  //Fe[q][i][j]=(2.0 - fac[q])*defMap.F[q][i][j] + (fac[q] - 1.0)*FeR[i]*advection_direction[q][j];
						  Fe[q][i][j]=defMap.F[q][i][j]/fac[q] + (1.0 - 1.0/fac[q])*FeR[i]*advection_direction[q][j]; 
		        }
		      }
				}		
			}
			//youngsModulus, poissonRatio are not changing so move this to the top, only need to do it one time
			//ResidualEq->setLameParametersByYoungsModulusPoissonRatio(youngsModulus, poissonRatio);
			ResidualEq->evaluateNeoHookeanStress(P, Fe);
			
			ResidualEq->residualForMechanics(fe_values, u_dof, R, P);
			ResidualEq->residualForDiff_ReacEq(fe_values, c1_dof, R, defMap, c_1, c_1_conv, c_1_flux, c_1_reaction);
			ResidualEq->residualForDiff_ReacEq(fe_values, c2_dof, R, defMap,  c_2, c_2_conv, c_2_flux, c_2_reaction);
			
		
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
			if(cell->material_id()==Cortex_id) local_dataStack[1]+=ResidualEq->volumeIntegration(fe_values, defMap.detF);
			else if(cell->material_id()==Subcortex_id) local_dataStack[2]+=ResidualEq->volumeIntegration(fe_values, defMap.detF);
			//calculate surface area of Cortex
			if(cell->material_id()==Cortex_id ){
	      for (unsigned int f=0; f<GeometryInfo<dim>::faces_per_cell; ++f){
					if (cell->at_boundary(f) == true){
						Cortex_fe_face_values.reinit(cell, f);
						const Point<dim> X = Cortex_fe_face_values.quadrature_point(0);
						if(X[0]>0){
							const unsigned int n_face_quadrature_points = Cortex_fe_face_values.n_quadrature_points;
							deformationMap<Sacado::Fad::DFad<double>, dim> defMap_surface(n_face_quadrature_points); 
							getDeformationMap<Sacado::Fad::DFad<double>, dim>(fe_values, Cortex_fe_face_values, u_dof, ULocal, defMap_surface);
							dealii::Table<1,double > tem;
							local_dataStack[0]+=ResidualEq->surfaceIntegration(Cortex_fe_face_values, defMap_surface.detF);
						}
					}
				}
			}
		}
	}
	MPI_Reduce(&local_dataStack[0], &dataStack[0], local_dataStack.size(), MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
	this->LinearSystemCompressAdd();
}

template class initBoundValProbs<1>;
template class initBoundValProbs<2>;
template class initBoundValProbs<3>;
