/**
 * @page growth Example 4 : Growth model
 * \section Introduction
* In this example, we first want to solve one diffusion equation:
 * \f[
	\frac{\partial C}{\partial t}+\nabla\cdot\boldsymbol{j}=0
\f]
* where \f$\boldsymbol{j}_1=-M\nabla C\f$ is the flux.
The boundary condition is
\f[
C=1 \text{ on }\Gamma_1; \quad \nabla \mu\cdot\boldsymbol{n}=0 \text{ on }\Gamma \backslash \Gamma_1
\f]
 * Besides chemistry, we also solve elasticity problem at finite strain:
 * \f[
	\nabla\cdot\boldsymbol{T} = \boldsymbol{0}\\
\boldsymbol{T}= \frac{1}{\det{\boldsymbol{F}^{\text{e}}}}\frac{\partial W}{\partial \boldsymbol{F}^{\text{e}}}\boldsymbol{F}^{\text{e}}
\f]
 * To make it more interesting, we have mechanical deformation induced by species intecalation.
 * \f[
\boldsymbol{F}=\boldsymbol{F}^{\text{e}}\boldsymbol{F}^{\text{g}}\\
\boldsymbol{F}^{\text{g}}=\left(\frac{C}{C_\text{0}}\right)^{\frac{1}{3}}\mathbb{1}
\f]
In this example, we have two domains/materials. The diffusion equation is solved over the whole domains, while we restrict the mechanics into the one of the domain.
We want to model that the species transports from domain 1 into domain 2, and casues expansion of domain 2. By this example, we demonstrate how to setup 
multiple domains, using <B>FE_Nothing</B> to exclude primary varialbe from certain doamin, and applying DOF constrains on the interface.

*\section imple Implementation
We first define the one scalar variable, and one vector variable for displacment:
*\code{.cpp}
std::vector<std::vector<std::string> > primary_variables(2);		
	  primary_variables[0].push_back("c1"); primary_variables[0].push_back("component_is_scalar");
	  primary_variables[1].push_back("u"); primary_variables[1].push_back("component_is_vector");
\endcode
We define two domains and basis order for each primal variables:
*\code{.cpp}
		int number_domain=2;
		int basis_order=1;
		std::vector<std::vector<int> > FE_support(number_domain);// store order of polynomial basis functions, 0 means FE_Nothing	
		FE_support[0].push_back(basis_order);
		FE_support[0].push_back(0);
		FE_support[1].push_back(basis_order);
		FE_support[1].push_back(basis_order);
\endcode
In domain 1, we set the order of polynomial basis to be zero for the second variable, i.e. the displacment, which will impose FE_Nothing to this field.
 
Before launching the <a href="../html/classinit_bound_val_probs.html">initBoundValProbs</a>, we need to initialize the <B>ParameterHandler</B> and declare all paramters we may use:
 *\code{.cpp}
ParameterHandler params;
		params.enter_subsection("parameters");
		params.declare_entry("youngsModulus","0",Patterns::Double() );
		params.declare_entry("poissonRatio","0",Patterns::Double() );
		params.declare_entry("c_ini","0",Patterns::Double() );
		params.declare_entry("M","0",Patterns::Double() );
		params.leave_subsection();	
\endcode
Now we just need to have class inherited from <a href="../html/classinit_bound_val_probs.html">initBoundValProbs</a> class, and overload several functions as discussed below.
function:
 *\code{.cpp}
template <int dim>
class growth: public initBoundValProbs<dim>
{
	public:
		growth(std::vector<std::vector<std::string> > _primary_variables, std::vector<std::vector<int> > _FE_support, ParameterHandler& _params);
		//this is a overloaded function 
		void get_residual(const typename hp::DoFHandler<dim>::active_cell_iterator &cell, const FEValues<dim>& fe_values, Table<1, Sacado::Fad::DFad<double> >& R, Table<1, Sacado::Fad::DFad<double>>& ULocal, Table<1, double >& ULocalConv);
		void setMultDomain();
		void setup_constraints();
		ParameterHandler* params;		
};
\endcode
In the first overloaded <B>get_residual</B> function, we define the residual for our problem. 
We use the pre-defined functions for diffusion equation and mechanics given the flux and stress P.
 *\code{.cpp}
this->ResidualEq.residualForDiffusionEq(fe_values, c_dof, R, c, c_conv, j_c);
 this->ResidualEq.residualForMechanics(fe_values, u_dof, R, P);	
\endcode
we can calculate the flux as before, here we only discuss how to evaluate the deformation gradient tensor with isotropic growth. 
First we need the total deformation graident tensor
*\code{.cpp}
	deformationMap<Sacado::Fad::DFad<double>, dim> defMap(n_q_points); 
	getDeformationMap<Sacado::Fad::DFad<double>, dim>(fe_values, u_dof, ULocal, defMap);
\endcode
Then the deformation gradient tensor induced by elasticity can be computed via the inverse of multiplicative decomposition, and the stress can be 
evaluated after chosing certain material model:
*\code{.cpp}
	dealii::Table<3, Sacado::Fad::DFad<double> > P(n_q_points,dim,dim), Fe(n_q_points,dim,dim);
	
	for(unsigned int q=0; q<n_q_points;q++){
		for (unsigned int i=0; i<dim; ++i){
			for (unsigned int j=0; j<dim; ++j){
	  		Fe[q][i][j]=defMap.F[q][i][j]/std::pow((c[q]/c_ini), 1.0/3.0); //Isotropic growth
			}
		}
	}
	this->ResidualEq.evaluateNeoHookeanStress(P, Fe);// NeoHookean model, Saint_Venant_Kirchhoff is also available 
\endcode
Setting multiple domains/materials is quite easy. We first need to overload the <B>setMultDomain</B> function, and pop out the cell iterator in <B>Deal.ii</B>
to set the material id for each cell/element. Remember to assign Fe_system to corresponding cells by calling the <B>set_active_fe_indices</B> afterwards.
*\code{.cpp}
template <int dim>
void growth<dim>::setMultDomain()
{
	this->pcout<<"setMultDomain"<<std::endl;
	
  for (typename Triangulation<dim>::active_cell_iterator cell = this->dof_handler.begin_active(); cell != this->dof_handler.end(); ++cell){
    Point<dim> cell_center = cell->center();
		if(cell_center[2]<0.5) cell->set_material_id(0);
		else cell->set_material_id(1);
	}
	//assign Fe_system to corresponding cells
	this->set_active_fe_indices (this->FE_support, this->dof_handler);
	
}
\endcode
To apply Dirchlet boundary condition and contraints in general, we need to overload the <B>setup_constraints</B> function.
We can apply the Dirchlet boundary condition easily by <B>constraints</B>. Similarly 
applying linear contraints on certain DOF can also be done by <B> constraints</B>, we only need to find the DOF where we want to apply the constraints:
*\code{.cpp}
template <int dim>
void growth<dim>::setup_constraints()
{
	this->pcout<<"setup_constraints"<<std::endl;
	hpFEM<dim>::constraints.clear ();
	DoFTools::make_hanging_node_constraints (this->dof_handler, hpFEM<dim>::constraints);
	
	int totalDOF=this->totalDOF(this->primary_variables);
  std::vector<bool> c_component (totalDOF, false); c_component[0]=true; 
	//apply constraints on boundary
	VectorTools:: interpolate_boundary_values (hpFEM<dim>::dof_handler, dim, ZeroFunction<dim> (totalDOF),hpFEM<dim>::constraints, c_component);
	//apply constraints on interface (domain 1 side)
	
	std::vector<types::global_dof_index> local_face_dof_indices_1 (this->fe_system[1]->dofs_per_face);
  typename hp::DoFHandler<dim>::active_cell_iterator cell = this->dof_handler.begin_active(), endc=this->dof_handler.end();
  for (;cell!=endc; ++cell){
		if(cell->material_id()==1 ){
    	for (unsigned int f=0; f<GeometryInfo<dim>::faces_per_cell; ++f){
				if (cell->at_boundary(f) == false){
					if(cell->neighbor(f)->material_id()==0 and cell->neighbor(f)->has_children() == false){
		  			cell->face(f)->get_dof_indices (local_face_dof_indices_1, 1);
		  			for (unsigned int i=0; i<local_face_dof_indices_1.size(); ++i){
					  	const unsigned int ck = this->fe_system[1]->face_system_to_component_index(i).first;
					  	if(ck>0) hpFEM<dim>::constraints.add_line (local_face_dof_indices_1[i]);//add constrain line for all u dofs
		  		 	}
					}
				}
			}
		}
	}
	
	hpFEM<dim>::constraints.close ();
}
\endcode



The last thing we need to define is the initial condition, we can simpily overload the  <a href="../html/class_initial_conditions.html#aa10cfdd7350c3810a8deab707f397657">vector_value()</a> function
of the <a href="../html/class_initial_conditions.html"> InitialConditions </a> class,
 *\code{.cpp}
void InitialConditions<dim>::vector_value (const Point<dim>   &p, Vector<double>   &values) const{
  Assert (values.size() == 2, ExcDimensionMismatch (values.size(), 2));
	if(p[2]==0) values(0)= 1;
  else values(0)= 0.5;
	values(1)=0;
}
\endcode

*\section results Results
\htmlonly <style>div.image img[src="E4.png"]{width:700px;}</style> \endhtmlonly 
\image html E4.png

*The results are generated using paramters shown below.
* The complete implementaion can be found at  <a href="https://github.com/mechanoChem/mechanoChemFEM/tree/example/Example4_growth">Github</a>. 
* 
*\code{.cpp}
#parameters file

subsection Problem
set print_parameter = true

set dt = 1
set totalTime = 5
set current_increment = 0
set off_output_index=0
set current_time = 0
set resuming_from_snapshot = false

#set mesh = /Users/wzhenlin/GitLab/researchCode/brainMorph/mesh/testMesh.msh
#set mesh = /home/wzhenlin/workspace/brainMorph/mesh/STA21_hex.msh
set output_directory = output/
set snapshot_directory = snapshot/

#FEM
set volume_quadrature = 3 
set face_quadrature = 2 

end

subsection Geometry
set X_0 = 0
set Y_0 = 0
set Z_0 = 0
set X_end = 5 
set Y_end = 5
set Z_end = 2 #no need to 2D

set element_div_x=5
set element_div_y=5
set element_div_z=4 #no need to 2D
end

subsection parameters
set c_ini =0.5
set youngsModulus =  5.0e3
set poissonRatio =  0.45
set M=1 
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
		set relative_norm_tolerance = 1.0e-12
		set absolute_norm_tolerance = 1.0e-12
		set max_iterations = 10
end
						
subsection Linear_solver
		set solver_method = PETScsuperLU
		set system_matrix_symmetricFlag = false # default is false
end

	\endcode
 */
