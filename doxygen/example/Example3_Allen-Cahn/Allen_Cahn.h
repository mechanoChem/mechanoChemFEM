/**
 * @page Allen_Cahn Example 3 : Allen-Cahn equation
 * \section Introduction
 * The Allen Cahn equation is:
 * \f[
	\frac{\partial C}{\partial t}=-M\mu\\
\mu=\frac{\partial g}{\partial C}-\nabla\cdot k\nabla C
\f]
where \f$g\f$ is a non-convex, ``homogeneous'' free energy density function, whose form has been chosen
\f[
g(C)=\omega(C-C_\alpha)^2(C-C_\beta)^2
\f]
The boundary condition is
\f[
\nabla \mu\cdot\boldsymbol{n}=0 \text{ on }\Gamma 
\f]
The double-well non-convex free energy density function, \f$g(C)\f$, drives segregation of the system into two distinct types.

Before implementing the code, we rewrite the equation in the form of standard diffusion-reaction equations:
\f[
\frac{\partial C}{\partial t}+\nabla\cdot\boldsymbol{j}=r\\
\boldsymbol{j}=-Mk\nabla C\\
r=-M\frac{\partial g}{\partial C}=-2M\omega(C-C_\alpha)(C-C_\beta)(2C-C_\alpha-C_\beta)
\f]

*\section imple Implementation
We first define the single scalar primary variable:
*\code{.cpp}
		std::vector<std::vector<std::string> > primary_variables(1);		
	  primary_variables[0].push_back("mu"); primary_variables[0].push_back("component_is_scalar");
\endcode
and setup the order of basis function for it:
*\code{.cpp}
int number_domain=1;
int diff_degree=1;
std::vector<std::vector<int> > FE_support(number_domain);// store order of polynomial basis functions, 0 means FE_Nothing	
FE_support[0].push_back(diff_degree);
\endcode
Before launching the <a href="../html/classinit_bound_val_probs.html">initBoundValProbs</a>, we need to initialize the <B>ParameterHandler</B> and declare all paramters we may use:
 *\code{.cpp}
ParameterHandler params;
params.enter_subsection("Parameters");	
params.declare_entry("omega","0",Patterns::Double() );
params.declare_entry("c_alpha","0",Patterns::Double() );
//... more parameters 
params.leave_subsection();	
\endcode
Now we just need to have class inherited from <a href="../html/classinit_bound_val_probs.html">initBoundValProbs</a> class, and overload the <a href="../html/classinit_bound_val_probs.html#ac8f2c3e2a1040c70b709900dc3dfdaea">get_residual()</a>
function:
 *\code{.cpp}
template <int dim>
class AllenCahn: public initBoundValProbs<dim>
{
	public:
		AllenCahn(std::vector<std::vector<std::string> > _primary_variables, std::vector<std::vector<int> > _FE_support, ParameterHandler& _params);
		//this is a overloaded function 
		void get_residual(const typename hp::DoFHandler<dim>::active_cell_iterator &cell, const FEValues<dim>& fe_values, Table<1, Sacado::Fad::DFad<double> >& R, Table<1, Sacado::Fad::DFad<double>>& ULocal, Table<1, double >& ULocalConv);
		ParameterHandler* params;		
};
\endcode
In the overloaded <B>get_residual</B> function, we define the residual for our problem. As our equation is standrad diffusion-reaction, we can simiply use 
the pre-defined model:
 *\code{.cpp}
this->ResidualEq.residualForDiff_ReacEq(fe_values, 0, R, mu, mu_conv, j_mu, rhs_mu);
\endcode
Before we call these two functions, we need the flux and reactions term, we need to first evaluate the values of the primary field and the spatial gradient.
We also need to evaluate the value of primary field at previous time step for the Backward Euler time scheme:
 *\code{.cpp}
	dealii::Table<1,double>  mu_conv(n_q_points);
	dealii::Table<1,Sacado::Fad::DFad<double> >  mu(n_q_points);
	dealii::Table<2,Sacado::Fad::DFad<double> >  mu_grad(n_q_points, dim);
	
	evaluateScalarFunction<double,dim>(fe_values, 0, ULocalConv, mu_conv);
	evaluateScalarFunction<Sacado::Fad::DFad<double>,dim>(fe_values, 0, ULocal, mu);	
	evaluateScalarFunctionGradient<Sacado::Fad::DFad<double>,dim>(fe_values, 0, ULocal, mu_grad);
	
	//evaluate diffusion and reaction term
	dealii::Table<1,Sacado::Fad::DFad<double> > rhs_mu(n_q_points);
	dealii::Table<2,Sacado::Fad::DFad<double> > j_mu(n_q_points, dim);
	
	j_mu=table_scaling<Sacado::Fad::DFad<double>, dim>(mu_grad,-kappa*M);
	
	for(unsigned int q=0; q<n_q_points;q++){
		 rhs_mu[q]=-M*(2*omega*(mu[q]-c_alpha)*(mu[q]-c_beta)*(2*mu[q]-c_alpha-c_beta));
	 }
\endcode
The last thing we need to define is the initial condition, we can simpily overload the  <a href="../html/class_initial_conditions.html#a369cea7ba74f8cd0a6ca12e0c164ff74">value()</a> function
of the <a href="../html/class_initial_conditions.html"> InitialConditions </a> class to define the initial condition for single scalar field. 
 *\code{.cpp}
double InitialConditions<dim>::value(const Point<dim>   &p, const unsigned int 	component) const{
  return 0.5 + 0.2*(static_cast <double> (rand())/(static_cast <double>(RAND_MAX))-0.5); 
}
\endcode

*\section results Results
\htmlonly <style>div.image img[src="E3.png"]{width:600px;}</style> \endhtmlonly 
\image html E3.png

*The results are generated using paramters shown below.
* The complete implementaion can be found at  <a href="https://github.com/mechanoChem/mechanoChemFEM/tree/example/Example3_Allen-Cahn">Github</a>. 
* 
*\code{.cpp}
#parameters file

subsection Problem
set print_parameter = true

set dt = 1
set totalTime = 50
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

set element_div_x=100
set element_div_y=100
set element_div_z=5 #no need to 2D
end

subsection Parameters
set omega = 0.8
set c_alpha = 0.2
set c_beta = 0.8
set kappa = 0.001
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
