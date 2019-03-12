/**
 * @page CahnHilliard Example 2 : CahnHilliard with single species
 * \section Introduction
 * The Cahn-Hilliard equation is:
 * \f[
	\frac{\partial C}{\partial t}=\nabla\cdot(M\nabla\mu)\\
\mu=\frac{\partial g}{\partial C}-\nabla\cdot k\nabla C
\f]
where \f$g\f$ is a non-convex, ``homogeneous'' free energy density function, whose form has been chosen
\f[
g(C)=\omega(C-C_\alpha)^2(C-C_\beta)^2
\f]
The boundary condiiton is
* \f[
\nabla\mu\cdot\boldsymbol{n}=0; \nabla C\cdot\boldsymbol{n}=0 \text{ on }\Gamma
\f]
The double-well non-convex free energy density function, \f$g(C)\f$, drives segregation of the system into two distinct types.

Attention is called to the well-known fourth-order nature of this partial differential equation in the concentration \f$C\f$. The polynomial basis can only achieve
C0 continuity across the element. To overcome this difficuity, we split the equation into two equations:
 * \f[
	\frac{\partial C}{\partial t}+\nabla\cdot(-M\nabla\mu)=0\\
k\nabla^2 C=\frac{\partial g}{\partial C}-\mu
\f]
The first equation is diffusion equation, and the second one the Possion equation. 
*\section imple Implementation
We first define the two scalar primary variables:
*\code{.cpp}
std::vector<std::vector<std::string> > primary_variables(2);		
	  primary_variables[0].push_back("c1"); primary_variables[0].push_back("component_is_scalar");
	  primary_variables[1].push_back("mu"); primary_variables[1].push_back("component_is_scalar");
\endcode
and we solve both species in one domain. We define the domain and basis order for each primal variables:
*\code{.cpp}
		int number_domain=1;
		int diff_degree=1;
		std::vector<std::vector<int> > FE_support(number_domain);// store order of polynomial basis functions, 0 means FE_Nothing	
		FE_support[0].push_back(diff_degree);
		FE_support[0].push_back(diff_degree);
\endcode
Before launching the <a href="../html/classinit_bound_val_probs.html">initBoundValProbs</a>, we need to initialize the <B>ParameterHandler</B> and declare all paramters we may use:
 *\code{.cpp}
ParameterHandler params;
params.enter_subsection("Concentration");	
params.declare_entry("omega","0",Patterns::Double() );
params.declare_entry("c_alpha","0",Patterns::Double() );
//... more parameters 
params.leave_subsection();	
\endcode
Now we just need to have class inherited from <a href="../html/classinit_bound_val_probs.html">initBoundValProbs</a> class, and overload the <a href="../html/classinit_bound_val_probs.html#ac8f2c3e2a1040c70b709900dc3dfdaea">get_residual()</a>
function:
 *\code{.cpp}
template <int dim>
class CahnHilliard: public initBoundValProbs<dim>
{
	public:
		CahnHilliard(std::vector<std::vector<std::string> > _primary_variables, std::vector<std::vector<int> > _FE_support, ParameterHandler& _params);
		//this is a overloaded function 
		void get_residual(const typename hp::DoFHandler<dim>::active_cell_iterator &cell, const FEValues<dim>& fe_values, Table<1, Sacado::Fad::DFad<double> >& R, Table<1, Sacado::Fad::DFad<double>>& ULocal, Table<1, double >& ULocalConv);
		ParameterHandler* params;		
};
\endcode
In the overloaded <B>get_residual</B> function, we define the residual for our problem. As our equations are one standrad diffusion equation and Possion equation, we can simiply use 
the pre-defined model:
 *\code{.cpp}
	this->ResidualEq.residualForDiffusionEq(fe_values, c_dof, R, c_1, c_1_conv, j_c_1);
	this->ResidualEq.residualForPoissonEq(fe_values, mu_dof, R, kappa_c_1_grad, rhs_mu);
\endcode
Though before we call these two functions, we need the flux and reactions terms, we need to first evaluate the values of the primary fields and their spatial gradients.
We also need to evaluate the value of primary fields at previous time step for the Backward Euler time scheme:
 *\code{.cpp}
	dealii::Table<1,double>  c_1_conv(n_q_points);
	dealii::Table<1,Sacado::Fad::DFad<double> > c_1(n_q_points), mu(n_q_points);
	dealii::Table<2,Sacado::Fad::DFad<double> >  c_1_grad(n_q_points, dim), mu_grad(n_q_points, dim);
	
	evaluateScalarFunction<double,dim>(fe_values, c_dof, ULocalConv, c_1_conv);
	evaluateScalarFunction<Sacado::Fad::DFad<double>,dim>(fe_values, c_dof, ULocal, c_1);	
	evaluateScalarFunctionGradient<Sacado::Fad::DFad<double>,dim>(fe_values, c_dof, ULocal, c_1_grad);
	
	evaluateScalarFunction<Sacado::Fad::DFad<double>,dim>(fe_values, mu_dof, ULocal, mu);	
	evaluateScalarFunctionGradient<Sacado::Fad::DFad<double>,dim>(fe_values, mu_dof, ULocal, mu_grad);
	
	//evaluate diffusion and reaction term
	dealii::Table<1,Sacado::Fad::DFad<double> > rhs_mu(n_q_points);
	dealii::Table<2,Sacado::Fad::DFad<double> > j_c_1(n_q_points, dim), kappa_c_1_grad(n_q_points, dim);
	
	j_c_1=table_scaling<Sacado::Fad::DFad<double>, dim>(mu_grad,-M);//-D_1*c_1_grad
	kappa_c_1_grad=table_scaling<Sacado::Fad::DFad<double>, dim>(c_1_grad,kappa);
	
	for(unsigned int q=0; q<n_q_points;q++) rhs_mu[q]=2*omega*(c_1[q]-c_alpha)*(c_1[q]-c_beta)*(2*c_1[q]-c_alpha-c_beta)-mu[q];
\endcode

The last thing we need to define is the initial condition, we can simpily overload the  <a href="../html/class_initial_conditions.html#aa10cfdd7350c3810a8deab707f397657">vector_value()</a> function
of the <a href="../html/class_initial_conditions.html"> InitialConditions </a> class,
 *\code{.cpp}
void InitialConditions<dim>::vector_value (const Point<dim>   &p, Vector<double>   &values) const{
  Assert (values.size() == 2, ExcDimensionMismatch (values.size(), 2));
  values(1) = 0;    
 values(0)= 0.5 + 0.04*(static_cast <double> (rand())/(static_cast <double>(RAND_MAX))-0.5);
}
\endcode
*\section results Results
\htmlonly <style>div.image img[src="E2.png"]{width:400px;}</style> \endhtmlonly 
\image html E2.png

*The results are generated using paramters shown below.
* The complete implementaion can be found at  <a href="https://github.com/mechanoChem/mechanoChemFEM/tree/example/Example2_CahnHilliard">Github</a>. 
* 
*\code{.cpp}
#parameters file

subsection Problem
set print_parameter = true

set dt = 5
set totalTime = 250
set current_increment = 0
set off_output_index=0
set current_time = 0
set resuming_from_snapshot = false

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
set X_end = 1 
set Y_end = 1
set Z_end = 2.0 #no need to 2D

set element_div_x=50
set element_div_y=50
set element_div_z=5 #no need to 2D
end

subsection Concentration
set omega = 0.25
set c_alpha = 0.2
set c_beta = 0.8
set kappa = 0.002
set M=0.1
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
