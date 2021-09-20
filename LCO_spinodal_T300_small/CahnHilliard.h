#include <deal.II/numerics/error_estimator.h>
#include "mechanoChemFEM.h"
#include "DNN.h"
template <int dim>
class CahnHilliard: public mechanoChemFEM<dim>
{
 public:
  CahnHilliard();
  //this is a overloaded function 
  void get_residual(const typename hp::DoFHandler<dim>::active_cell_iterator &cell, const FEValues<dim>& fe_values, Table<1, Sacado::Fad::DFad<double> >& R, Table<1, Sacado::Fad::DFad<double>>& ULocal, Table<1, double >& ULocalConv);
  void make_grid();
  void refine_grid_once();
  void solve_ibvp();

  DNN freeEn;
  int iter_count=0;
  int max_refine_level=5;
  int min_refine_level=2;

  ParameterHandler* params;		
};

template <int dim>
CahnHilliard<dim>::CahnHilliard()
{
  //This let you use one params to get all parameters pre-defined in the mechanoChemFEM
  params=this->params_mechanoChemFEM;
  params->enter_subsection("Concentration");
  params->declare_entry("c_ini","0",Patterns::Double() );
  params->declare_entry("kappa","0",Patterns::Double() );
  params->declare_entry("M","0",Patterns::Double() );
  params->leave_subsection();		
  params->enter_subsection("Geometry");
  params->declare_entry("R","1",Patterns::Double() );
  params->declare_entry("n_refine","1",Patterns::Integer() );
  params->leave_subsection();		
	
  //Declear the parameters before load it
  this->load_parameters("../parameters.prm");
	
  //define main fields from parameter file.
  this->define_primary_fields();
  //Set up the ibvp.
  this->init_ibvp();

  freeEn.reinit(2,"",true);

}

template <int dim>
void CahnHilliard<dim>::get_residual(const typename hp::DoFHandler<dim>::active_cell_iterator &cell, const FEValues<dim>& fe_values, Table<1, Sacado::Fad::DFad<double> >& R, Table<1, Sacado::Fad::DFad<double>>& ULocal, Table<1, double >& ULocalConv)
{
  //evaluate primary fields
  params->enter_subsection("Concentration");
  //double M=params->get_double("M");
  Sacado::Fad::DFad<double> M;
  double kappa=params->get_double("kappa");
  params->leave_subsection();
  double jn = 0.;

  unsigned int n_q_points= fe_values.n_quadrature_points;
  int c_dof=0, mu_dof=1;
		
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
	
  //j_c_1=table_scaling<Sacado::Fad::DFad<double>, dim>(mu_grad,-M);//-D_1*c_1_grad
  double kB = 8.617333262145e-5; // eV/K
  double T = 300; // K
  for(unsigned int q=0; q<n_q_points; q++){
    M = 1.e-2*std::exp(-274.0*(1.05-c_1[q])*(0.47-c_1[q])*(1.-c_1[q]))*c_1[q]/(kB*T); //micron^2/(eV*s*unit_cell)
    for (int i=0; i<dim; ++i){
      j_c_1[q][i] = -M*mu_grad[q][i];
      kappa_c_1_grad[q][i] = kappa*c_1_grad[q][i];
    }
  }
  //kappa_c_1_grad=table_scaling<Sacado::Fad::DFad<double>, dim>(c_1_grad,kappa);

  //DNN free energy evaluation
  std::vector<Sacado::Fad::DFad<double> > features(1), y;
  std::vector<std::vector<Sacado::Fad::DFad<double> > > dy_dx;

  for(unsigned int q=0; q<n_q_points; q++){
    features[0]=c_1[q];
    freeEn.eval(features,y,dy_dx);
    rhs_mu[q]=dy_dx[0][0]-mu[q];
  }
	
  //call residual functions
  this->ResidualEq.residualForDiffusionEq(fe_values, c_dof, R, c_1, c_1_conv, j_c_1);
  this->ResidualEq.residualForPoissonEq(fe_values, mu_dof, R, kappa_c_1_grad, rhs_mu);

  for (unsigned int faceID=0; faceID<2*dim; faceID++){
    if(cell->face(faceID)->boundary_id()==0 ){
      FEFaceValues<dim> fe_face_values(fe_values.get_fe(), *(this->common_face_quadrature), update_values | update_quadrature_points | update_JxW_values);
      fe_face_values.reinit(cell,faceID);
      this->ResidualEq.residualForNeummanBC(fe_values, fe_face_values, c_dof, R, jn);
    }
  }
	
}


template <int dim>
void CahnHilliard<dim>::make_grid()
{

  params=this->params_mechanoChemFEM;	
  params->enter_subsection("Geometry");
  double R=params->get_double("R");
  int n_refine=params->get_integer("n_refine");
  params->leave_subsection();	

  GridGenerator::hyper_ball (this->triangulation,Point<dim>(),R);
  static const SphericalManifold<dim> boundary;
  this->triangulation.set_all_manifold_ids_on_boundary(0);
  this->triangulation.set_manifold (0, boundary);
  this->triangulation.refine_global (n_refine);
}

/**************************************
 * Locally refine the mesh once
 *************************************/
template <int dim>
void CahnHilliard<dim>::refine_grid_once(){

  PETScWrappers::Vector localized_U(this->solution);
  Vector<float> estimated_error_per_cell (hpFEM<dim>::triangulation.n_active_cells());
  KellyErrorEstimator<dim>::estimate (hpFEM<dim>::dof_handler,
				      QGauss<dim-1>(3),
				      typename FunctionMap<dim>::type(),
				      localized_U,
				      estimated_error_per_cell);
  //GridRefinement::refine_and_coarsen_optimize (hpFEM<dim>::triangulation, estimated_error_per_cell,2);
  GridRefinement::refine_and_coarsen_fixed_fraction (hpFEM<dim>::triangulation, estimated_error_per_cell, 0.18, 0.12);
  //GridRefinement::refine_and_coarsen_fixed_fraction (hpFEM<dim>::triangulation, estimated_error_per_cell, 0.25, 0.12);
  if (hpFEM<dim>::triangulation.n_levels() > max_refine_level){
    for (typename Triangulation<dim>::active_cell_iterator cell = hpFEM<dim>::triangulation.begin_active(max_refine_level); cell != hpFEM<dim>::triangulation.end(); ++cell){
      if (cell->is_locally_owned()){
	cell->clear_refine_flag ();
      }
    }
  }


  SolutionTransfer<dim,PETScWrappers::Vector,hp::DoFHandler<dim> > solution_trans(hpFEM<dim>::dof_handler);
  hpFEM<dim>::triangulation.prepare_coarsening_and_refinement();

  solution_trans.prepare_for_coarsening_and_refinement(localized_U);
  hpFEM<dim>::triangulation.execute_coarsening_and_refinement();

  this->mark_boundary();
  this->setup_linear_system();
  this->apply_boundary_condition();

  PETScWrappers::Vector localized_Unew(this->solution);
  solution_trans.interpolate(localized_U,localized_Unew);
  this->solution = localized_Unew;
  this->constraints_mechanoChemFEM->distribute (this->solution);

  this->pcout << "   Number of active cells:       " << hpFEM<dim>::triangulation.n_active_cells() << std::endl;
  this->pcout << "   Number of degrees of freedom: " << hpFEM<dim>::dof_handler.n_dofs() << std::endl; 

}


/**************************************
 * Implement adaptive time stepping
 *************************************/
template <int dim>
void CahnHilliard<dim>::solve_ibvp()
{
  int converged;
  while(true){
    converged=this->nonlinearSolve(this->solution);
    //reset iter_count to 0
    if (converged > -1){
      break;
    }
    else{
      iter_count=0;
      this->pcout<<"Not converged, reduce dt."<<std::endl;
      this->current_dt /= std::pow(10.,0.25);
    }
  }
  if (converged < 4){
    iter_count++;
  }
  else{
    iter_count = 0;
  }
  // Check if the current_time is a multiple of the doubled timestep
  //bool multiple = (std::fmod(this->current_time,2.*this->current_dt) < 1.e-8);
  //if((iter_count>5) && multiple && (this->current_dt < 1.)) {
  if((iter_count>5)) {
    this->pcout<<"Increasing timestep \n";
    iter_count=0;
    this->current_dt *= std::pow(10.,0.125);//double dt
  }

  // Adaptive mesh refinement
  //refine_grid_once();

  this->solution_prev=this->solution;
}


template <int dim>
void InitialConditions<dim>::vector_value (const Point<dim>   &p, Vector<double>   &values) const{
  Assert (values.size() == 2, ExcDimensionMismatch (values.size(), 2));
  values(1) = 0;    
  values(0)= 0.53 + 0.01*(static_cast <double> (rand())/(static_cast <double>(RAND_MAX))-0.5);
}
template class InitialConditions<1>;
template class InitialConditions<2>;
template class InitialConditions<3>;
