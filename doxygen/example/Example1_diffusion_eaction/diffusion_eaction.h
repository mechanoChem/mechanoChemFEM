/**
 * @page diffusion_reaction Example 1 : Turing patterns: diffusion-reaction systems
 * \section Introduction
 * We solve two diffusion reaction equations:
 * \f[
	\frac{\partial C_\text{1}}{\partial t}+\nabla\cdot\boldsymbol{j}_1=r_1 \\
\frac{\partial C_\text{2}}{\partial t}+\nabla\cdot\boldsymbol{j}_2=r_2
\f]
* where \f$\boldsymbol{j}_1 \f$ and  \f$\boldsymbol{j}_2 \f$ are flux terms:
 * \f[
	 \boldsymbol{j}_1=-M_1\nabla C_\text{1}; \quad  \boldsymbol{j}_2=-M_2\nabla C_\text{2}\\
\f]
*\f$r_1\f$ and \f$r_2\f$ are reaction terms:
 * \f[
	 r_1= R_{10}+R_{11}C_1+R_{13}C_1^2C_2; \quad  r_1= R_{20}+R_{21}C_1^2C_2
\f]
The boundary condiiton is
* \f[
\boldsymbol{j}_1\cdot\boldsymbol{n}=j_n \text{ on }\Gamma_2;\quad  \quad 
\boldsymbol{j}_1\cdot\boldsymbol{n}=0 \text{ on }\Gamma \backslash \Gamma_2;  \quad  \quad 
 \boldsymbol{j}_2\cdot\boldsymbol{n}=0 \text{ on }\Gamma 
\f]
The coupled diffusion-reaction equations for two species follow Schnakenberg kinetics.
For an activator-inhibitor species pair, these equations use auto-inhibition with cross-activation of a short range species, and auto-activation with cross-inhibition of a long range species to form so-called Turing patterns.

*\section imple Implementation
We first define the two scalar primary variables:
*\code{.cpp}
std::vector<std::vector<std::string> > primary_variables(2);		
	  primary_variables[0].push_back("c1"); primary_variables[0].push_back("component_is_scalar");
	  primary_variables[1].push_back("c2"); primary_variables[1].push_back("component_is_scalar");
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
params.declare_entry("D_1","0",Patterns::Double() );
params.declare_entry("D_2","0",Patterns::Double() );
//... more parameters 
params.leave_subsection();	
\endcode
Now we just need to have class inherited from <a href="../html/classinit_bound_val_probs.html">initBoundValProbs</a> class, and overload the <a href="../html/classinit_bound_val_probs.html#ac8f2c3e2a1040c70b709900dc3dfdaea">get_residual()</a>
function:
 *\code{.cpp}
template <int dim>
class diffusion_reaction: public initBoundValProbs<dim>
{
	public:
		diffusion_reaction(std::vector<std::vector<std::string> > _primary_variables, std::vector<std::vector<int> > _FE_support, ParameterHandler& _params);
		//this is a overloaded function 
		void get_residual(const typename hp::DoFHandler<dim>::active_cell_iterator &cell, const FEValues<dim>& fe_values, Table<1, Sacado::Fad::DFad<double> >& R, Table<1, Sacado::Fad::DFad<double>>& ULocal, Table<1, double >& ULocalConv);
		ParameterHandler* params;		
};
\endcode
In the overloaded <B>get_residual</B> function, we define the residual for our problem. As our equations are standrad diffusion-reaction, we can simiply use 
the pre-defined model:
 *\code{.cpp}
this->ResidualEq.residualForDiff_ReacEq(fe_values, c_1_dof, R, c_1, c_1_conv, j_c_1, reaction_1);
this->ResidualEq.residualForDiff_ReacEq(fe_values, c_2_dof, R, c_2, c_2_conv, j_c_2, reaction_2);
\endcode
Though before we call these two functions, we need the flux and reactions terms, we need to first evaluate the values of the primary fields and their spatial gradients.
We also need to evaluate the value of primary fields at previous time step for the Backward Euler time scheme:
 *\code{.cpp}
	dealii::Table<1,double>  c_1_conv(n_q_points), c_2_conv(n_q_points);
	dealii::Table<1,Sacado::Fad::DFad<double> > c_1(n_q_points), c_2(n_q_points);
	dealii::Table<2,Sacado::Fad::DFad<double> >  c_1_grad(n_q_points, dim), c_2_grad(n_q_points, dim);
	
	evaluateScalarFunction<double,dim>(fe_values, c_1_dof, ULocalConv, c_1_conv);
	evaluateScalarFunction<Sacado::Fad::DFad<double>,dim>(fe_values, c_1_dof, ULocal, c_1);	
	evaluateScalarFunctionGradient<Sacado::Fad::DFad<double>,dim>(fe_values, c_1_dof, ULocal, c_1_grad);
	
	evaluateScalarFunction<double,dim>(fe_values, c_2_dof, ULocalConv, c_2_conv);
	evaluateScalarFunction<Sacado::Fad::DFad<double>,dim>(fe_values, c_2_dof, ULocal, c_2);	
	evaluateScalarFunctionGradient<Sacado::Fad::DFad<double>,dim>(fe_values, c_2_dof, ULocal, c_2_grad);

	
	//evaluate diffusion and reaction term
	dealii::Table<1,Sacado::Fad::DFad<double> > reaction_1(n_q_points), reaction_2(n_q_points);
	dealii::Table<2,Sacado::Fad::DFad<double> > j_c_1(n_q_points, dim),j_c_2(n_q_points, dim);
	
	j_c_1=table_scaling<Sacado::Fad::DFad<double>, dim>(c_1_grad,-D_1);//-D_1*c_1_grad
	j_c_2=table_scaling<Sacado::Fad::DFad<double>, dim>(c_2_grad,-D_2);//-D_2*c_2_grad
	
	for(unsigned int q=0; q<n_q_points;q++){
		reaction_1[q]=R_10+R_11*c_1[q]+R_12*c_2[q]+R_13*c_1[q]*c_1[q]*c_2[q];
		reaction_2[q]=R_20+R_21*c_1[q]+R_22*c_2[q]+R_23*c_1[q]*c_1[q]*c_2[q];
	}
\endcode
Besides the residual for the PDEs, we have the boundary conditions on one surface:
 *\code{.cpp}
	for (unsigned int faceID=0; faceID<2*dim; faceID++){
		if(cell->face(faceID)->boundary_id()==dim*2 ){
		  FEFaceValues<dim> fe_face_values(fe_values.get_fe(), *(this->common_face_quadrature), update_values | update_quadrature_points | update_JxW_values);
			fe_face_values.reinit(cell,faceID);
			this->ResidualEq.residualForNeummanBC(fe_values, fe_face_values, c_1_dof, R, jn);
		}
	}
\endcode
The last thing we need to define is the initial condition, we can simpily overload the  <a href="../html/class_initial_conditions.html#aa10cfdd7350c3810a8deab707f397657">vector_value()</a> function
of the <a href="../html/class_initial_conditions.html"> InitialConditions </a> class,
 *\code{.cpp}
void InitialConditions<dim>::vector_value (const Point<dim>   &p, Vector<double>   &values) const{
  Assert (values.size() == 2, ExcDimensionMismatch (values.size(), 2));
  values(1) = 0;    
  values(0)= 0.5 + 0.1*static_cast <double> (rand())/(static_cast <double>(RAND_MAX/2.0))/2;
}
\endcode
*\section results Results
The right plot shows the patterns of the Schnakenberg kinetics.
\htmlonly <style>div.image img[src="E1.png"]{width:400px;}</style> \endhtmlonly 
\image html E1.png

*The results are generated using paramters shown below.
* The complete implementaion can be found at  <a href="https://github.com/mechanoChem/mechanoChemFEM/tree/example/Example1_diffusion_eaction">Github</a>.
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
set X_end = 10
set Y_end = 10
set Z_end = 2.0 #no need to 2D

set element_div_x=50
set element_div_y=50
set element_div_z=5 #no need to 2D
end

subsection Concentration

set D_1 = 0.1
set D_2 = 4.0
set R_10 = 0.1
set R_11 = -1
set R_13 = 1
set R_20 = 0.9
set R_21 = 0
set R_23 = -1
set jn=-0.01
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
		set max_iterations = 20
end
						
subsection Linear_solver
		set solver_method = PETScsuperLU
		set system_matrix_symmetricFlag = false # default is false
end
	\endcode
 */