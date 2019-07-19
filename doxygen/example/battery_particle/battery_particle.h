/**
 * @page battery_particle Battery model at particle scale
 * \section Introduction
A battery cell consists of porous, positive and negative electrodes, a separator and current collectors(see Figure 1). 
The electrode consists of active-material particles and carbon-binders, and the porespace is filled with the electrolyte.
The two electrodes are isolated by a separator which is usuallymade of polyolefin for Li-ion batteries.
The porous separator is perfused with electrolyte, allowingionic transport between the electrodes.
Metallic current collectors are located at either end of the battery.

This example present a coupled continuum formulation for the electrostatic,  chemical,  thermal, mechanical and fluid physics in battery materials.
The treatment is at the particle scale, at whic hthe active particles held together by binders, the porous separator, current collectors and the perfusing electrolyte are explicitly modeled.

This code is used to generate results for paper:
<B>A multi-physics battery model with particle scale resolution ofporosity evolution driven by intercalation strain and electrolyte flow</B>,
Z. Wang, K. Garikipati, Journal of the Electrochemical Society, Vol. 165: A2421-A438, 2018, doi:10.1149/2.0141811jes 
*\htmlonly <style>div.image img[src="domain_binder.png"]{width:700px;}</style> \endhtmlonly 
*\image html domain_binder.png A schematic showing each component of a battery cell. The entire battery cell is a rolled or folded structure of a layer.

\subsection section1 The multi-physics particle scale model
We lay down the governing equations for primary variables in three dimensions.
As is the convention in continuum mechanics, the solids (current collector,active-material particle and polymeric separator particles) are described in a Lagrangian setting,
while the fluid (electrolyte) flow is described in an Eulerian setting. Electrolytic flow around thedeforming solid components alters the fluid domain.
The fluid mesh therefore must be remapped asthe computation proceeds. The Arbitrary Lagrangian-Eulerian (ALE) framework is adopted herefor this purpose.
\subsubsection subsub1 Electro-chemo-thermal equations

In the active material we have conservation of lithium:
\f[
\frac{\partial C_\text{Li}}{\partial t}+\nabla\cdot\boldsymbol{j}=0\\
\boldsymbol{j}=-D\nabla C_\text{Li}
\f]
where \f$D\f$ is diffusivity.
Conservation of of Li\f$^+\f$  cations in the electrolyte is written similarly:
\f[
\frac{\partial C_{\text{Li}^+}}{\partial t}+\nabla\cdot\boldsymbol{j}_+=0
\f]
In general, the flux of each species can depend on the concentration gradients of the other species. This is a special case of the Onsager reciprocity relations with the dependence on chemical potential gradients reduced to concentration gradient dependence. Following Newman and Thomas-Alyea\f$^{\text{J. Newman and K. E. Thomas-Alyea,John Wiley and Sons, Inc, (2004).}}\f$, we have 
\f[
\boldsymbol{j}_+=-D_{+}\nabla C_{\text{Li}^+} +\frac{t_+}{F}\boldsymbol{i}_\text{E}+C_{\text{Li}^+}\boldsymbol{v}
\f]
for a binary solution. Here, \f$D_{+}\f$ is the diffusivity, \f$t_+\f$ is the transference number of the cation, \f$\boldsymbol{v}\f$ is the electrolyte velocity, \f$F\f$ is the Faraday constant, and \f$\boldsymbol{i}_\text{E}\f$ is the total current in the electrolyte phase which will be derived in what follows.

The total current in the current collector and active material is governed by Ohm's law:
\f[
\boldsymbol{i}_\text{S}=-\sigma_\text{S} \nabla\phi_\text{S}
\f]
where \f$\phi_\text{S}\f$ is the electric potential, and the conductivity \f$\sigma_\text{S}\f$ differs between the current collector and active material. In the electrolyte we have
\f[
\boldsymbol{i}_\text{E}=-\sigma_\text{E}\nabla\phi_\text{E}-\gamma_\text{D}\nabla\ln C_{\text{Li}^+}
\f]
where \f$\sigma_\text{E}\f$ is the electrolyte's conductivity, and \f$\gamma_\text{D}\f$ is the diffusion conductivity evaluated as:
\f[
\gamma_\text{D}=\frac{2R\theta \sigma_\text{E}}{F}\left(1-t_+\right)\left(1+\frac{d \ln f}{d\ln  C_{\text{Li}^+}}\right),
\f]
with \f$f\f$ being the mean molar activity coefficient of the electrolyte, which is assumed to be constant. Thus the relation simplifies to
\f[
\gamma_\text{D}=\frac{2 R\theta \sigma_\text{E}}{F}\left(1-t_+\right).
\f]
Under the electroneutrality approximation, we have
\f[
\nabla\cdot (- \sigma_\text{s}\nabla\phi_\text{S})=0
\f]
in the active material, and
\f[
\nabla\cdot\left(-\sigma_\text{E}\nabla\phi_\text{E}-\gamma_\text{D}\nabla\ln C_{\text{Li}^+}\right)=0
\f]
in the electrolyte. These two equations properly describe the electric potentials with the electroneutrality approximation, such that the double layer effect is neglected.

\subsection section2 The standard thermal equations
Heat generation and transport are governed by the heat equation, which is derived from the first law of thermodynamics. Since the velocity of the electrolyte is low, with a Peclet number Pe \f$\sim 1.0\times 10^{-7}\f$, we neglect heat flux by advection. For the temperature \f$\theta\f$, we have the standard form of the heat equation:
\f[
\rho C_p\frac{\text{d}\theta}{\text{d}t}+\nabla \cdot\boldsymbol{q}=0
\f]
where \f$\rho\f$ is  the mass density of the electrode, \f$C_p\f$ is the specific heat and \f$\boldsymbol{q}\f$ is the heat flux. The heat flux can be expressed as
\f[
\boldsymbol{q}=\phi F\sum_i z_i\boldsymbol{j}_i-\lambda\nabla\theta+\boldsymbol{q}^\text{D}
\f]
where \f$\lambda\f$ is the thermal conductivity. The first term on right is associated with Joule heating and the last term \f$\boldsymbol{q}^\text{D}\f$ is the Dufour effect, which is ignored in this work. Again, with the electroneutrality approximation we have $\nabla\cdot\sum z_i \boldsymbol{j}_i=0$, we have
\f[
\rho C_p\frac{\text{d}\theta}{\text{d}t}-\lambda \nabla^2\theta+\nabla\phi_\text{S}\cdot\boldsymbol{i}_\text{S}=0
\f]
in the electrode, and 
\f[
\rho C_p\frac{\text{d}\theta}{\text{d}t}-\lambda \nabla^2\theta+\nabla\phi_\text{E}\cdot\boldsymbol{i}_\text{E}=0
\f]
in the electrolyte. 

\subsection section3 Finite strain mechanics
Lithium intercalation and de-intercalation induce expansion and contraction, respectively, of the active material. Additionally, the active material, binder, porous separator and current collector undergo thermal expansion. 
The kinematics of finite strain leads to the following decomposition in the active material:
\f[
\boldsymbol{F}=\boldsymbol{F}^\text{e}\boldsymbol{F}^\text{c}\boldsymbol{F}^{\theta}.
\f]
For other solid sub-domains (binder and polymeric separator particles) we have
\f[
\boldsymbol{F}=\boldsymbol{F}^\text{e}\boldsymbol{F}^{\theta}.
\f]
Here, \f$\boldsymbol{F} = \boldsymbol{I} + \partial\boldsymbol{u}/\partial\boldsymbol{X}\f$, is the total deformation gradient tensor in each of the relevant solid regions. Its multiplicative components \f$\boldsymbol{F}^\text{e}\f$, \f$\boldsymbol{F}^\text{c}\f$ and \f$\boldsymbol{F}^{\theta}\f$, are, respectively, the elastic, chemical (induced by lithium intercalation) and thermal components. In the absence of a body force the strong form of the mechanics problem in the current configuration is 
\f[
\nabla\cdot\boldsymbol{T} = \boldsymbol{0},
\f]
\f[
\text{for}\quad \boldsymbol{T}= \frac{1}{\det{\boldsymbol{F}^\text{e}}}\frac{\partial W}{\partial \boldsymbol{F}^\text{e}}\boldsymbol{F}^{\text{e}^\text{T}},
\f]
where \f$\boldsymbol{T}\f$  is the Cauchy stress tensor and \f$W\f$ is the strain energy density function.
The chemical and thermal expansion components of the multiplicative decomposition of \f$\boldsymbol{F}\f$ are modelled as isotropic, e.g.:
\f[
F^{\text{c}}_{iJ}=(1+\beta^{\text{c}})^{1/3}\delta_{iJ} \\
F^{\theta}_{iJ}=(1+\beta^{\theta})^{1/3}\delta_{iJ}
\f]
The lithiation swelling response function is parameterized by the lithium concentration. We write
\f[
\beta^{\theta}(\theta)=\Omega_{\theta} (\theta-\theta_0)
\f]
for the thermal expansion functions.

\subsection section4 Incompressible fluid model
We model the electrolyte as an incompressible, creeping flow. In this regime, with inertia and body forces being neglected, the Stokes equations are: 
\f[
-2\eta\nabla\cdot\varepsilon(\boldsymbol{v})+\nabla p=0
\f]
\f[
\nabla\cdot\boldsymbol{v}=0
\f]
where \f$\varepsilon(\boldsymbol{v})=\frac{1}{2}(\nabla\boldsymbol{v}+\nabla\boldsymbol{v}^\text{T})\f$ is the strain rate tensor, and \f$\eta\f$ is the dynamic viscosity.

\subsection section5 The Arbitrary Lagrangian-Eulerian framework
In the fluid (electrolyte) phase, whose description is Eulerian, the mesh must be mapped to evolve with the  fluid sub-domain as it flows around the deforming solid components. If the mesh displacement is \f$\boldsymbol{u}_\text{m}\f$, the spatial gradient operator acting on a variable \f$\tau\f$ in the Eulerian setting is
\f[
\nabla \tau=\frac{\partial\tau}{\partial\boldsymbol{X}}\boldsymbol{F}_\text{m}^{-1}
\f]
where \f$\boldsymbol{F}_\text{m}= \boldsymbol{I} +\partial \boldsymbol{u}_\text{m}/\partial \boldsymbol{X}\f$, is the  deformation gradient tensor of the mesh. In the solid phase, \f$\boldsymbol{u}_\text{m}\f$ coincides with the displacement of material points, \f$\boldsymbol{u}\f$, and \f$\boldsymbol{F}_\text{m}\f$ is replaced by the deformation gradient tensor of the solid, \f$\boldsymbol{F}\f$. In the fluid phase the description of fluid mesh deformation could be arbitrary, but its choice proves to be critical in solid-fluid interaction problems. The mesh displacement in the fluid phase can be solved by Poisson equations, arbitrary elasticity equations or bi-harmonic equations.

The large lithium intercalation strain may lead to a dramatic volume change of the active material. The large deformation of active material may cause extreme distortions of the mesh in the fluid. For this reason, we apply  adaptive mesh rezoning schemes. The idea is to introduce a constraint condition over two adjacent element nodes \f$i\f$ and \f$j\f$ such that the relative difference of net displacements is less than the element length \f$h^e\f$:
\f[
|\boldsymbol{u}_m^i-\boldsymbol{u}_m^j|\le\alpha h^e,
\f]
which is equivalent to
\f[
    |\nabla\boldsymbol{u}_m|\le\alpha
\f]
where \f$\alpha\in[0,1)\f$ is the tolerance for element distortion. This constraint is applied element wise as a penalty term:
\f[
    \nabla^2(1+\tau_m\boldsymbol{u}_m)=\boldsymbol{0}
\f]
where \f$\tau_m\f$ is the weight function that imposes spatially varying stiffening effects for the mesh. The key idea for choosing \f$\tau_m\f$ is to enforce this value to be high for smaller elements and small for large elements. In this work we choose
\f[
   \tau_m=\left(\frac{\det \boldsymbol{F}_m^0}{\det \boldsymbol{F}_m} \right)^\delta
\f]
where \f$\boldsymbol{F}_m^0\f$ is the initial deformation gradient tensor of the mesh and \f$\delta\f$ is a constant value.

\subsection subsection3 Boundary and interface condition
The multi-physics character of this problem extends to the restriction of certain physics--and therefore the corresponding partial differential equations--to specific sub-domains. The conventional application of boundary conditions translates to interface conditions where the sub-domains meet.
During discharging or charging, the following chemical reactions occur at the interface between the active material and electrolyte:
\f[
\text{Li}\to \text{Li}^++q^-, \quad \text{-ve electrode during discharge/+ve electrode during charge,}\nonumber\\
 \text{Li}^++q^- \to\text{Li},\quad \text{+ve electrode during discharge/-ve electrode during charge}.\nonumber
\f]
The reaction rate is given by the Butler-Volmer model
\f[
j_n=j_0\left(\text{exp}\left(\frac{\alpha_aF}{R\theta}(\phi_\text{S}-\phi_\text{E}-U)\right)-\text{exp}\left(-\frac{\alpha_aF}{R\theta}(\phi_\text{S}-\phi_\text{E}-U)\right)\right) \label{eq:BVeuqation}\\
j_0=k_0(C_{\text{Li}^+})^{\alpha_a} \frac{(C_\text{Li}^{\text{max}}-C_\text{Li})^{\alpha_a}}{(C_1^{\text{max}})^{\alpha_a}}  \frac{(C_\text{Li})^{\alpha_c}}{{(C_1^\text{max}})^{\alpha_c}}
\f]
where \f$\alpha_a, \alpha_c\f$ are transfer coefficients and \f$k_0\f$ is a kinetic rate constant.
The interface condition for lithium and lithium ions:
\f[
\boldsymbol{j}\cdot \boldsymbol{n}_\text{s}=j_n \quad \text{on } \Gamma_\text{s-e},\\
\boldsymbol{j}_+\cdot \boldsymbol{n}_\text{e}=-j_n\ \quad \text{on } \Gamma_\text{s-e},
\f]
The interface condition for the electric potential:
\f[
 \boldsymbol{i}_\text{S}\cdot\boldsymbol{n}_\text{s}=Fj_n  \quad \text{on } \Gamma_\text{s-e}\\
 \boldsymbol{i}_\text{E}\cdot\boldsymbol{n}_\text{e}=-Fj_n \ \quad \text{on } \Gamma_\text{s-e}
\f]
At the interface, the chemical reactions also generate heat:
\f[
Q_{\Gamma_\text{s-e}}=q_{\text{rxn}}+q_{\text{rev}}
\f]
 where
\f[
q_{\text{rxn}}=Fj_n(\phi_\text{S}-\phi_\text{E}-U) ;\quad \text{irreversible entropic heat,} \\
q_{\text{rev}}=Fj_n\theta\frac{\partial U}{\partial \theta};\quad\text{reversible entropic heat.}
\f]
The electrolytic fluid is bounded within the solid materials: active material, binder, polymeric separator particles, and current collector. No-slip conditions give  
\f[
\boldsymbol{v}=\dot{\boldsymbol{u}} \quad \text{on } \Gamma_\text{s-e}, \Gamma_\text{e-c}, \Gamma_\text{e-p}
\f]
Traction continuity on all solid-fluid interfaces gives:
\f[
\boldsymbol{T}\boldsymbol{n}_\text{s}=-(2\eta\cdot\varepsilon(\boldsymbol{v})-p\boldsymbol{1})\boldsymbol{n}_\text{e} \quad \text{on } \Gamma_\text{s-e}
\\
\boldsymbol{T}\boldsymbol{n}_\text{c}=-(2\eta\cdot\varepsilon(\boldsymbol{v})-p\boldsymbol{1})\boldsymbol{n}_\text{e} \quad \text{on } \Gamma_\text{e-c}
\\
\boldsymbol{T}\boldsymbol{n}_\text{p}=-(2\eta\cdot\varepsilon(\boldsymbol{v})-p\boldsymbol{1})\boldsymbol{n}_\text{e} \quad \text{on } \Gamma_\text{e-p}
\f]
Compatibility at the solid-fliud interfaces equates the ALE mesh deformation, \f$\boldsymbol{u}_\text{m}\f$ in the electrolyte sub-domain to the solid deformations:
\f[
    \boldsymbol{u}_\text{m}=\boldsymbol{u} \quad \text{on } \Gamma_\text{s-e}, \Gamma_\text{c-e}, \Gamma_\text{p-e}.
\f]
\htmlonly <style>div.image img[src="boundaryCondition.png"]{width:500px;}</style> \endhtmlonly 
\image html boundaryCondition.png A schematic showing boundary and interface conditions.

\section Implementation Implementation: level 0 developer
The code is very similar to \ref battery_electrodeScale. Here we mainly discuss how to apply interface condition and apply constraints during solve step.

<B>Applying boundary and interface conditions</B><br>
During updating linear system, we apply these boundary and interface conditions discussed above. For interface coditions, we first need to find these interface.
For example to find interface bewteen active_material and electrolyte we can loop the element face and check the domain id of elements on both sides.
\code{.cpp}
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
\endcode
Once we find the interface, the corresponding function for assemble interfacial residual is called. Now let's take a look of one of these functions.
The following code iitialize the finite elemnt system for the elements on both sides and their faces.
\code{.cpp}
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
\endcode
Then we can evaluate all fields at surface quadrature points using \link evaluateScalarFunction<T, dim> evaluate functions \endlink
\code{.cpp}
	evaluateVectorFunctionGradient<Sacado::Fad::DFad<double>,dim>(electrolyte_fe_values, electrolyte_fe_face_values, v_dof, ULocal_electrolyte, velocity_grad,defMap_mesh_surface);
	evaluateScalarFunction<Sacado::Fad::DFad<double>,dim>(electrolyte_fe_values, electrolyte_fe_face_values, p_dof, ULocal_electrolyte, Pressure);
	evaluateScalarFunction<Sacado::Fad::DFad<double>,dim>(electrolyte_fe_values, electrolyte_fe_face_values, c_li_plus_dof, ULocal_electrolyte, c_li_plus);
	evaluateScalarFunction<Sacado::Fad::DFad<double>,dim>(electrolyte_fe_values, electrolyte_fe_face_values, phi_e_dof, ULocal_electrolyte, phi_e);
	
	evaluateScalarFunction<Sacado::Fad::DFad<double>,dim>(activeMaterial_fe_values, activeMaterial_fe_face_values, T_dof, ULocal, T);
	evaluateScalarFunction<Sacado::Fad::DFad<double>,dim>(activeMaterial_fe_values, activeMaterial_fe_face_values, c_li_dof, ULocal, c_li);
	evaluateScalarFunction<Sacado::Fad::DFad<double>,dim>(activeMaterial_fe_values, activeMaterial_fe_face_values, phi_s_dof, ULocal, phi_s);
\endcode
These interface conditions are Neumman conditions. And the residual is assembled by calling the residual functions as we discuss before. 

<B>Applying constraints during solve step</B><br>
The interface condition \f$\boldsymbol{u}_\text{m}=\boldsymbol{u}\f$ is a master-slave interface, meanning \f$\boldsymbol{u}_\text{m}\f$ does not affect \f$\boldsymbol{u}\f$ but just take its value as a Dirichlet boundary condition.
So we cannot use the constraint materix to apply this constraint. Here we apply this constrain manually during the solve step. First we need to find the pair of DOFs for \f$\boldsymbol{u}\f$ and \f$\boldsymbol{u}_\text{m}\f$ at same node. 
And save it as constrain_index.
\code{.cpp}
				    if (ck1>=0 && ck1<dim ){
							int dof_per_node=local_face_dof_indices2.size()/4;
							hpFEM<dim>::constraints.add_line (local_face_dof_indices1[i]);//u_mesh														
							for (unsigned int j=(node_u_m_applied/2)*dof_per_node; j<local_face_dof_indices2.size(); ++j){
								const unsigned int ck3 = fe_system[current_collector_fe]->face_system_to_component_index(j).first;
								if(ck3==ck1){
									constrain_index[0].push_back(local_face_dof_indices1[i]);
									constrain_index[1].push_back(local_face_dof_indices2[j]);
									node_u_m_applied++;
									break;
								}
							}						
						}	
\endcode
This constraint is equivalent to \f$d\boldsymbol{u}_\text{m}=d\boldsymbol{u}\f$ which can be applied during the nonlinear solve step. We only need to overload function \link solveClass< dim, matrixType, vectorType >::apply_dU_constrain apply_dU_constrain \endlink.
\code{.cpp}
template <int dim>
void initBoundValProbs<dim>::apply_dU_constrain(PETScWrappers::MPI::Vector& dU)
{
	PETScWrappers::Vector localized_dU(this->dU);	
  for(unsigned int i=0;i<constrain_index[0].size();i++){
		localized_dU(constrain_index[0][i])=localized_dU(constrain_index[1][i]);
	}
	dU = localized_dU;
}
\endcode
Also use another \link solveClass< dim, matrixType, vectorType >::nonlinearSolve(vectorType &U, vectorType &dU) nonlinear solver function \endlink which take a extra argument <B>vectorType& dU</B>.

*\section results Results
We present few results here by using the following Parameter file:
\code{.cpp}
#parameters file
#1. Thermal modeling of cylindrical lithium ion battery during discharge cycle, Dong Hyup Jeon ⇑,1, Seung Man Baek,2011,Energy Conversion and Management LiC6
# 2. Numerical study on lithium titanate battery thermal response under adiabatic condition,Qiujuan Sun a, Qingsong Wanga 2015, Energy Conversion and Management
#3.white Theoretical Analysis of Stresses in a Lithium Ion Cell    
#4 A Computational Model of the Mechanical Behavior within Reconstructed LixCoO2 Li-ion Battery Cathode Particles, Veruska Malavé, J.R. Berger, EA
#5. A pseudo three-dimensional electrochemicalethermal model of a prismatic LiFePO4 battery during discharge process
#6. solid diffusion Single-Particle Model for a Lithium-Ion Cell: Thermal Behavior, Meng Guo,a Godfrey Sikha,b,* and Ralph E. Whitea
#global parameters


#declare problem setting
subsection Problem
set print_parameter = true
set dt = 0.5
set totalTime = 370
set step_load = true

set first_domain_id = 1
set active_material_id=1
set electrolyte_id=2
set current_collector_id=3
set binder_id=3 
set solid_id=4

set active_material_fe=0
set electrolyte_fe=1
set current_collector_fe=2
set binder_fe=2//not necessary
set solid_fe=3


#directory
set mesh = /home/wzhenlin/workspace/batteryCode/mesh/2D/2D_R10_labeled.inp
set output_directory = output/
set snapshot_directory = snapshot/ 
#FEM
set volume_quadrature = 4 
set face_quadrature = 3 
#applied current
set IpA = -350 
end

# some useful geometry information beforehand
subsection Geometry
set X_0 = -13
set Y_0 = 0
set Z_0 = 0
set X_end = 13 
set Y_end = 153
set Z_end = 0

set electrode_Y1 = 76
set electrode_Y2 = 99
set currentCollector_Y1 = 16
set currentCollector_Y2 = 144
set particle_R = 10
set particle_number = 15
end

# initial condition
subsection Initial condition
set c_li_max_neg = 28.7e-3
set c_li_max_pos = 37.5e-3
set c_li_100_neg = 0.915
set c_li_100_pos = 0.022
set c_li_0_neg = 0.02
set c_li_0_pos = 0.98
set c_li_plus_ini = 1.0e-3
set T_0 = 298
end

# parameter for elastiticy equations
subsection Elasticity
set youngsModulus_Al = 70e3 #cite5
set youngsModulus_Cu = 120e3 #cite5	
set youngsModulus_binder = 2.3e3
set youngsModulus_s_neg = 12e3 #cite 3 in porosity paper
set youngsModulus_s_pos = 370e3 #cite4
set youngsModulus_sep = 0.5e3 #cite 3 in porosity paper

set nu_Al = 0.35
set nu_Cu = 0.34
set nu_binder = 0.3 
set nu_sep = 0.35 #cite 3 in porosity paper
set nu_s_neg = 0.3 #cite 3 in porosity paper
set nu_s_pos = 0.2 #cite4 

set omega_c = 3.5 #not used but 3.5 is from linear coff in porosity paper
#following from cite 3 in porosity paper
set omega_t_s_neg = 8.62e-6 
set omega_t_s_pos = 4.06e-6
set omega_t_sep = 13.32e-5
set omega_t_Al = 23.6e-6
set omega_t_Cu = 17e-6
set omega_t_binder = 190e-6 
end

#parameter for fluid equations
subsection Fluid
set viscosity = 1.0 #original data is 1.0e-9, but scaled to 1.0
set youngsModulus_mesh = 1.0e3
set nu_mesh = 0.3
end

# parameter for electro-chemo equations
subsection ElectroChemo #cite6
set sigma_s_neg = 1.5e8
set sigma_s_pos = 0.5e8
set sigma_Al = 3.5e9
set sigma_Cu = 6e9
set sigma_binder =1e8

set t_0 = 0.2
set D_li_neg = 5e-1 #cite 40 in porosity paper in code use expression cite 3 
set D_li_pos = 1.0e-1 #cite3

# parameter for Butler-Volmer equations at ElectricChemo class
set F = 96485.3329
set Rr = 8.3144598
set alpha_a = 0.5
set alpha_c = 0.5
set k_neg = 0.8 #to match porosity paper
set k_pos = 0.8 #to match porosity paper
set c_max_neg = 28.7e-3
set c_max_pos = 37.5e-3

end

# parameter for thermal equations
subsection Thermal
set lambda_s_neg = 5e6 #cite1
set lambda_s_pos = 1.85e6  #cite1
set lambda_e = 0.33e6  #cite2
set lambda_sep = 1e6  #cite1
set lambda_Al = 160e6  #cite5
set lambda_Cu = 400e6  #cite5
set lambda_binder = 0.19e6 #PVDF

set density_s_neg = 5.03e-15  #cite1
set density_s_pos = 2.29e-15  #cite1
set density_sep = 1.2e15  #cite1
set density_e = 1e-15  #cite2
set density_Al = 2.7e-15  #cite5
set density_Cu = 8.9e-15  #cite5
set density_binder = 1.78e-15  #cite5

set Cp_s_neg = 0.7e12  #cite1
set Cp_s_pos = 1.17e12  #cite1
set Cp_e = 1978.16e12  #cite2
set Cp_sep = 700e12  #cite1
set Cp_Al = 903e12  #cite5
set Cp_Cu = 385e12  #cite5
set Cp_binder = 700e12 

set h = 5 #Heat transfer coefficient
end
			
					
#==============================================================================
# parameters reserved for deal.ii first level code:
#nonLinear_method : classicNewton
#solver_method (direct) : PETScsuperLU, PETScMUMPS
#solver_method (iterative) : PETScGMRES PETScBoomerAMG
#relative_norm_tolerance, absolute_norm_tolerance, max_iterations
#
subsection Nonlinear_solver
		set nonLinear_method = classicNewton
		set relative_norm_tolerance = 1.0e-16
		set absolute_norm_tolerance = 1.0e-7
		set max_iterations = 60
end
						
subsection Linear_solver
		set solver_method = PETScsuperLU
		set system_matrix_symmetricFlag = false # default is false
end
\endcode
\htmlonly <style>div.image img[src="c_li.png"]{width:500px;}</style> \endhtmlonly 
\image html c_li.png Lithium concentration under different far-field boundary conditions and intercalation strain as-sumptions

\htmlonly <style>div.image img[src="c_li_plus.png"]{width:500px;}</style> \endhtmlonly 
\image html c_li_plus.png Lithium ion concentration under different far-field boundary conditions and intercalation strain as-sumptions
		
\htmlonly <style>div.image img[src="V.png"]{width:500px;}</style> \endhtmlonly 
\image html V.png Velocity of electrolyte under different far-field boundary conditions and intercalation strain as-sumptions
 * The complete implementaion can be found at  <a href="https://github.com/mechanoChem/mechanoChemFEM/tree/example/Battery%20model%20at%20particle%20scale">Github</a>. 
 */