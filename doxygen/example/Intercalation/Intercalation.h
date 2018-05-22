/**
 * @page Intercalation Example 1 : Intercalation with finte strain mechanics
 * \section Introduction
 * We solve two diffusion reaction equations:
 * \f[
	\frac{\partial C_\text{1}}{\partial t}+\nabla\cdot\boldsymbol{j}_1=r_1 \\
\frac{\partial C_\text{2}}{\partial t}+\nabla\cdot\boldsymbol{j}_2=r_2
\f]
* where \f$\boldsymbol{j}_1 \f$ and  \f$\boldsymbol{j}_2 \f$ are flux terms:
 * \f[
	 \boldsymbol{j}_1^{diff}=-M_1\nabla C_\text{1}; \quad  \boldsymbol{j}_2^{diff}=-M_2\nabla C_\text{2}\\
\f]
*\f$r_1\f$ and \f$r_2\f$ are reaction terms, in this example they are
 * \f[
	 r_1= reac_{10}; \quad  r_1= reac_{20}
\f]

 * Besides chemistry, we also solve elasticity problem at finite strain:
 * \f[
	\nabla\cdot\boldsymbol{T} = \boldsymbol{0}\\
\boldsymbol{T}= \frac{1}{\det{\boldsymbol{F}^{\text{e}}}}\frac{\partial W}{\partial \boldsymbol{F}^{\text{e}}}\boldsymbol{F}^{\text{e}}
\f]
 * To make it more interesting, we have mechanical deformation due to species intecalation.
 * \f[
\boldsymbol{F}=\boldsymbol{F}^{\text{e}}\boldsymbol{F}^{\text{g}}\\
\boldsymbol{F}^{\text{g}}=\left(\frac{C_\text{1}}{C_\text{10}}\right)^{\frac{1}{3}}\mathbb{1}
\f]
* For deomstraion we use a very simple mesh with three different domains. We will solve the abolve equations on domain 1 and 2.
\htmlonly <style>div.image img[src="domain.png"]{width:500px;}</style> \endhtmlonly 
 *@image html domain.png 
 *
 * \subsection define1 Define primary variables over different domains
 * In deal.ii to solve equations in different domains, \c Fe_Nothing is used and different \c Fe_System need to be defined. By using dealMutiphysics,
 * it can be easily taken care of by two user defined vector: \c <b>primary_variables</b> and \c <b>FE_support</b>.
*\code{.cpp}
		std::vector<std::vector<std::string> > primary_variables(3);		
		primary_variables[0].push_back("u"); primary_variables[0].push_back("component_is_vector");
	  	primary_variables[1].push_back("c1"); primary_variables[1].push_back("component_is_scalar");
	  	primary_variables[2].push_back("c2"); primary_variables[2].push_back("component_is_scalar");
		
		int number_domain=3;
		int mech_degree=1;
		int diff_degree=1;
		std::vector<std::vector<int> > FE_support(number_domain);// store order of polynomial basis functions, 0 means FE_Nothing	
		for(unsigned int i=0;i<2;i++){
			FE_support[i].push_back(mech_degree);
			FE_support[i].push_back(diff_degree);
			FE_support[i].push_back(diff_degree);
		}
		FE_support[2].push_back(0);
		FE_support[2].push_back(0);
		FE_support[2].push_back(0);
 *\endcode
 * Then user can use function 
*\code{.cpp}
 *  hpFEM<dim>::setup_FeSystem(fe_system, fe_collection, q_collection, primary_variables_dof,primary_variables,FE_support,volume_quadrature);
* \endcode
 * to setup \c fe_system, \c fe_collection, \c q_collection with quadarture points : \c volume_quadrature.\n
 * If the mesh is generated outside deal.ii (e.g. cubit), material ID can be pre-defined. To setup the \c Fe_System for each element, we can simily
 * use
*\code{.cpp}
 *  hpFEM<dim>::set_active_fe_indices (FE_support, hpFEM<dim>::dof_handler);
* \endcode 
* \subsection outPut Output and restart
* We can use class \c FEMdata to easily write output and have capability of re-start the code. After initializing the class we can set the output name by
*\code{.cpp}
 *  FEMdata<dim,vectorType>::set_output_name(primary_variables);
* \endcode 
* And write vtk file using 
*\code{.cpp}
  std::string output_path = output_directory+"output-0.vtk";
  FEMdata<dim,vectorType>::write_vtk(solution_0, output_path);
* \endcode 
* To use restart we first need to create a snapshot and resume the vector from ths snapshot
*\code{.cpp}
		std::string snapshot_path = snapshot_directory+"snapshot-"+std::to_string(current_increment)+".txt";
		FEMdata<dim,vectorType>::create_vector_snapshot(solution_old, snapshot_path);
		FEMdata<dim,vectorType>::resume_vector_from_snapshot(solution_new, snapshot_path);
* \endcode 
* \subsection assemble Assemble rsesidual functions of two diffuction reactions equations and elasticity at finite strain 
* In FEM modeling, we need to provide system_matrix (Jacobin matrix)\f$\frac{\partial \boldsymbol{R}}{\partial \boldsymbol{x}} \f$, and right hand side 
* vector \f$-\boldsymbol{R}\f$. They are usually achieved by assemble \c local_matrix and \c local_rhs over elements. We first need to overload abstract
* function \c updateLinearSystem. In this exmaple we need to solve two diffusion-reaction equations and one finite strain mechanics problem.
*\subsubsection DRq diffusion-reaction equations
* For diffusion-reaction problems we need to provide the flux and reaction term, and in this example we also have advection and few extra terms for stablization.
* All the term will be combined into \c flux and \c reaction terms.
 * \f[
	 \boldsymbol{j}_1^{diff}=-M_1\nabla C_\text{1}; \quad  \boldsymbol{j}_2^{diff}=-M_2\nabla C_\text{2}\\
\f]
The following code basically evaluates the terms shown above.
*\code{.cpp}
for(unsigned int q=0; q<n_q_points;q++){
	reac_11[q]=0; reac_12[q]=0; c_1_reac[q]=0; 
	reac_21[q]=0; reac_22[q]=0; c_2_reac[q]=0; 
	velDotGradSpatc_1[q]=0;
	velDotGradSpatc_2[q]=0;
	const Point<dim> posR = fe_values.quadrature_point(q);
	for(unsigned int i=0; i<dim;i++){
		c_1_flux[q][i]=-mobility_c1*c_1_grad[q][i];
		c_2_flux[q][i]=-mobility_c2*c_2_grad[q][i];
	}		
	c_1_reaction[q]=reac_10;
	c_2_reaction[q]=reac_20;
}
*\endcode
* Also in above code we need values for primary variables and their gradients (at current configuration), we can evaluate them using \c evaluation \c functions:
*\code{.cpp}
			evaluateScalarFunction(fe_values, primary_variables_dof[1], ULocalConv, c_1_conv);
			evaluateScalarFunction(fe_values, primary_variables_dof[1], U0Local, c_1_0);
			evaluateScalarFunction(fe_values, primary_variables_dof[1], ULocal, c_1);	
			evaluateScalarFunctionGradient(fe_values, primary_variables_dof[1], ULocal, c_1_grad,defMap);
			
			evaluateScalarFunction(fe_values, primary_variables_dof[2], ULocalConv, c_2_conv);
			evaluateScalarFunction(fe_values, primary_variables_dof[2], U0Local, c_2_0);
			evaluateScalarFunction(fe_values, primary_variables_dof[2], ULocal, c_2);	
			evaluateScalarFunctionGradient(fe_values, primary_variables_dof[2], ULocal, c_2_grad,defMap);
*\endcode
* The weak form of diffution-reaction equations can be written as
\f[
\mathscr{R}=\int_{\Omega_{\text{e}}}w\frac{\partial C}{\partial t}\text{d}v-\int_{\Omega_{\text{e}}} \nabla w \boldsymbol{j} \text{d}v+\int_{s}w\boldsymbol{j}\cdot\boldsymbol{n} \text{d}s=\int_{\Omega_{\text{e}}}w r \text{d}v
\f]
*After we have all fluxes and reaction terms we can simply call the corresponding residual function to evaluate the residual for \f$C_\text{1}\f$ and \f$C_\text{2}\f$ . 
*\code{.cpp}
			Residual<vectorType,dim>::residualForDiff_ReacEq(fe_values, primary_variables_dof[1], R, defMap, c_1, c_1_conv, c_1_flux, c_1_reaction);
			Residual<vectorType,dim>::residualForDiff_ReacEq(fe_values, primary_variables_dof[2], R, defMap, c_2, c_2_conv, c_2_flux, c_2_reaction);
*\endcode
*Please note the above functions will not do the the boundary integration (in this example we have trivial neumman boundary condition). To apply neumman boundary conditions
* please check \c Residual class which offer a easy way to apply a variety of tpye of neumman boundary conditions. 
* \subsubsection mechanics Finite strain mechanics using NeoHookean model
* For all mechanics problem we first need to define stress, strain and have a deformation map from which we obtain the strain
*\code{.cpp}
			dealii::Table<3, Sacado::Fad::DFad<double> > P(n_q_points,dim,dim);
			dealii::Table<3, Sacado::Fad::DFad<double> > Fe(n_q_points,dim,dim);
			deformationMap<Sacado::Fad::DFad<double>, dim> defMap(n_q_points); 
			getDeformationMap(fe_values, primary_variables_dof[0], ULocal, defMap);
*\endcode
* In this example, we consider no growth and isotropic growth:
*\code{.cpp}
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
					else{pcout << "Growth type not supported\n"; exit(1);}
				}
				else {fac[q] = 1.0;}
				
		    if (fac[q] <= 1.0e-15){
					printf("*************Non positive growth factor*************. Value %12.4e\n", fac[q].val());
		    }
				
		    if (fac[q] < sat){ fac[q] = 1.0; }
		    else{ fac[q] /= sat; }
				
				if(std::strcmp(GROWTH.c_str(),"Isotropic")==0 ){
					for (unsigned int i=0; i<dim; ++i){
						for (unsigned int j=0; j<dim; ++j){
					  	Fe[q][i][j]=defMap.F[q][i][j]/fac[q];
						}
					}
				}		
*\endcode
* After all this we can evaluate stress. For Neohookean model
 \f[
\boldsymbol{P}=\boldsymbol{F}^e\boldsymbol{S}\\
\boldsymbol{S}=0.5\lambda\det(\boldsymbol{C})\boldsymbol{C}^{-1}-(0.5\lambda+\mu)\boldsymbol{C}^{-1}+\mu\boldsymbol{1}
\f]
where \f$\boldsymbol{C}=\boldsymbol{F}^e(\boldsymbol{F}^e)^{T}\f$. If young's modulus and Possion ratio are provided, we need to set the Lame parameters first and evaluate the stress  
*\code{.cpp}
Residual<vectorType,dim>::setLameParametersByYoungsModulusPoissonRatio(youngsModulus, poissonRatio);
Residual<vectorType,dim>::evaluateNeoHookeanStress(P, Fe);
*\endcode
The weak form of elasticity problem is 
*\f[
\mathscr{R}_u=\int_{\Omega_{\text{e}}}\nabla \boldsymbol{w}\boldsymbol{T} \text{d}v- \int_{s}\boldsymbol{w} \boldsymbol{f} \cdot \boldsymbol{n} \text{d}s = \boldsymbol{0}
\f]
Again we just simply call the corresponding function to evalue the residual function for \f$\boldsymbol{u}\f$.
\code{.cpp}
Residual<vectorType,dim>::residualForMechanics(fe_values, primary_variables_dof[0], R, P);
*\endcode 
* 
* \subsection solve Solving nonlinear residual functions
* when we have assembled the system_matrix \f$\frac{\partial \boldsymbol{R}}{\partial \boldsymbol{x}} \f$, and right hand side 
* vector \f$-\boldsymbol{R}\f$, we can solve the Nonlinear residual functions using Newton family method. For dynamics problem we just need to specify the initial and total
* time, and then call \c nonlinearSolve
 *\code{.cpp}
  for (current_time=initial_time; current_time<=total_time; current_time+=dt){
    current_increment++;
  	t_solve = clock();
		solveClass< dim, matrixType, vectorType >::nonlinearSolve(solution);
		solution_prev=solution;
 }
*\endcode
* in the \c solveClass, it will take the system_matrix and right hand side 
* vector generated by \c updateLinearSystem we overloaded.\n \n
* Now let us discuss the parameter managements and how to choose parameters in solvers.
*\subsection parameter Parameterhandler
*Recall that 
 * \f[
	 \boldsymbol{j}_1^{diff}=-M_1\nabla C_\text{1}; \quad  \boldsymbol{j}_2^{diff}=-M_2\nabla C_\text{2}\\
\f]
 * \f[
	 r_1= reac_{10} \\
 r_1= reac_{20}
\f]
* In this example (and almost all real applications), there are many parameters. 
* DealShell code use \c Parameterhandler to manager all parameters, and \c Parameterhandler is also used in this example to manager the parameters. \c Parameterhandler allows
* us to modify parameters at runtime without re-compiling the code. To use \c Parameterhandler we first need to declare all parameters we will use. Please note
* parameters used in dealShell are already declared so we only need to declare parameters as we need.
In this example we nee
*\code{.cpp}
	params->declare_entry("dt","0",Patterns::Double() );
	params->declare_entry("totalTime","0",Patterns::Double() );
	
	params->declare_entry("mesh","1",Patterns::FileName() );
	params->declare_entry("output_directory","1",Patterns::DirectoryName() );
	params->declare_entry("snapshot_directory","1",Patterns::DirectoryName() );
	
	//declare paramters for mechanics
	params->enter_subsection("Mechanics");
	params->declare_entry("youngsModulus","0",Patterns::Double() );
	params->declare_entry("poissonRatio","0",Patterns::Double() );
	params->declare_entry("saturation_matID_0","0",Patterns::Double() );
	params->declare_entry("saturation_matID_1","0",Patterns::Double() );
	
	params->declare_entry("GROWTH","Isotropic",Patterns::Selection("Uniform|Isotropic") );
	params->leave_subsection();	
	
	//declare paramters for concentrations
	params->enter_subsection("Concentration");
	params->declare_entry("c1_ini","0",Patterns::Double() );
	params->declare_entry("c2_ini","0",Patterns::Double() );
	params->declare_entry("c1_ini_interface","0",Patterns::Double() );
	params->declare_entry("c2_ini_interface","0",Patterns::Double() );
	
	params->declare_entry("mobility_c1","0",Patterns::Double() );
	params->declare_entry("mobility_c2","0",Patterns::Double() );
	params->declare_entry("reac_10","0",Patterns::Double() );
	params->declare_entry("reac_20","0",Patterns::Double() );
	params->leave_subsection();	
*\endcode
* Then we need to create a \c .prm file to set values for all parameters. Also we need to set parameters for solvers in the file.
*\code{.cpp}
#parameters file

#set global parameters
set dt = 0.1
set totalTime = 1

set mesh = /Users/wzhenlin/GitLab/researchCode/brainMorph/mesh/testMesh.msh
set output_directory = output/
set snapshot_directory = snapshot/
# set equations <name />type(not useful now)

#set Mechanics = NeoHookean finite strain
#set Concentration = Diffusion reaction


subsection Mechanics
		set youngsModulus =  5.0e3
		set poissonRatio =  0.45
			
		set saturation_matID_0 = 10
		set saturation_matID_1 = 1
		
	 set GROWTH = Tangential
		#set GROWTH = Isotropic
end

subsection Concentration
		set c1_ini = 0.5
		set c2_ini = 0.5
		set c1_ini_interface = 1
		set c2_ini_interface = 1
						
		set mobility_c1 = 0.1
	 set mobility_c2 = 0.1
	 set reac_10 = 0
		set reac_20 = 0

end
						
#
# parameters reserved for deal.ii first level code:
#nonLinear_method : classicNewton
#solver_method (direct) : PETScsuperLU, PETScMUMPS
#solver_method (iterative) : PETScGMRES PETScBoomerAMG
#relative_norm_tolerance, absolute_norm_tolerance, max_iterations
#
subsection Nonlinear_solver
		set nonLinear_method = classicNewton
		set relative_norm_tolerance = 5.0e-12
		set absolute_norm_tolerance = 5.0e-12
		set max_iterations = 50
end
						
subsection Linear_solver
		set solver_method = PETScMUMPS
		set system_matrix_symmetricFlag = false # default is false
end
*\endcode
* Alteravely we can use GUI to set the parameters.
 * \section lib dealMutiphysics 
 * To use dealMutiphysics, we need to define a class derived from \c solveClass<int dim, class matrixType, class vectorType> and \c hpFEM<int dim>,
 * so that we can used pre-defined variables, functions, and overload the abstract virtual function \c updateLinearSystem() in parent class.
 * The last reason make the inherientance necessary. 
 * Also we usually include the following three classes
*\code{.cpp}
		ParameterHandler* params;	
		Residual<Sacado::Fad::DFad<double>,dim>* ResidualEq;
		FEMdata<dim, PETScWrappers::MPI::Vector>* FEMdata_out;
	\endcode
*\section results Results
\htmlonly <style>div.image img[src="c1.png"]{width:500px;}</style> \endhtmlonly 
\image html c1.png concentraion c1 at deformed configuration
\htmlonly <style>div.image img[src="u.png"]{width:500px;}</style> \endhtmlonly 
\image html u.png displacment magnitude
*\section com Complete code
 * The complete implementaion can be found at GitLab [https://gitlab.com/compPhysCode/brainMorph]. 
 */
