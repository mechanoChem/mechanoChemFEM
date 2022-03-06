/*************************************
 * Allen-Cahn, Cahn-Hilliard, mechanics w/ eigenstrain
 *************************************/

#include <deal.II/numerics/error_estimator.h>
#include "extraFunctions.h"
#include "mechanoChemFEM.h"
#include "DNN.h"

#define PI std::acos(-1)

/**************************************
 * Class declaration
 *************************************/
template <int dim>
class precipitate_PF: public mechanoChemFEM<dim>
{

 public:
  precipitate_PF(std::vector<std::vector<std::string> > _primary_variables,
		 std::vector<std::vector<int> > _FE_support,
		 ParameterHandler& _params);
  void setup_model();
  void get_residual(const typename hp::DoFHandler<dim>::active_cell_iterator &cell,
		    const FEValues<dim>& fe_values,
		    Table<1, Sacado::Fad::DFad<double> >& R,
		    Table<1, Sacado::Fad::DFad<double> >& ULocal,
		    Table<1, double >& ULocalConv);
  void setup_constraints();
  void refine_grid_once();
  void refine_grid();
  void solve_ibvp();
  void update_post_TS();
  ParameterHandler* params;		

  //Chemistry parameters
  double M, L, c_avg, kappa1, kappa2;
  DNN freeEn;

  //Mechanics parameters
  DNN psi_Ni, psi_Ni3Al;

  //Mesh refinement
  int n_init_global_refine, n_init_local_refine, max_refine_level;

  //Other
  int iter_count=0;
  unsigned int n_q_points;
  int c_dof, mu_dof, eta1_dof, eta2_dof, eta3_dof, u_dof;
  dealii::Table<2,Sacado::Fad::DFad<double> > eye;

};

/**************************************
 * Class constructor
 *************************************/
template <int dim>
precipitate_PF<dim>::precipitate_PF(std::vector<std::vector<std::string> > _primary_variables,
				    std::vector<std::vector<int> > _FE_support,
				    ParameterHandler& _params)
:
mechanoChemFEM<dim>(_primary_variables, _FE_support, _params),
  params(&_params),
  eye(dim,dim)
{

  setup_model();
  refine_grid();

}

/**************************************
 * Define parameters
 *************************************/
template <int dim>
void precipitate_PF<dim>::setup_model()
{
  freeEn.reinit(2,"cfe_",false);
  psi_Ni.reinit(2,"Ni_",true);
  psi_Ni3Al.reinit(2,"Ni3Al_",true);

  params->enter_subsection("parameters");

  //Chemistry parameters
  M = params->get_double("mobility");
  L = params->get_double("kinetic_coeff");
  kappa1 = std::pow(params->get_double("sqrt_kappa1"),2.);
  kappa2 = std::pow(params->get_double("sqrt_kappa2"),2.);
  c_avg=params->get_double("c_avg");
  params->leave_subsection();

  //Mesh refinement parameters
  params->enter_subsection("mesh_refinement");
  n_init_global_refine = params->get_integer("n_init_global_refine");
  n_init_local_refine = params->get_integer("n_init_local_refine");
  max_refine_level = params->get_integer("max_refine_level");
  params->leave_subsection();

  for (unsigned int i=0; i<dim; ++i){
    eye[i][i] = 1.;
  }

  c_dof=0; mu_dof = 1;
  eta1_dof=2; eta2_dof=3; eta3_dof=4;
  u_dof=5;

}

template<typename T>
T a(T eta0){
  
  double temp = 326.85, //C
    a_g = 0.354590 + 5.74055e-6*temp - 1.00986e-10*temp*temp,
    a_gp = 0.356276 + 6.16199e-6*temp - 11.32198e-10*temp*temp,
    c_g = 0.12,
    c_gp = 0.23;

  return (a_gp - a_g)/(c_gp - c_g)*(eta0 - c_g) + a_g;

}

double a_c(){
  
  double temp = 326.85, //C
    a_g = 0.354590 + 5.74055e-6*temp - 1.00986e-10*temp*temp,
    a_gp = 0.356276 + 6.16199e-6*temp - 11.32198e-10*temp*temp,
    c_g = 0.12,
    c_gp = 0.23;

  return (a_gp - a_g)/(c_gp - c_g);

}

template<typename T>
void chem_features(T c,T eta1,T eta2,T eta3,std::vector<T> &features){

  features.resize(4);
  features[0] = c;
  features[1] = 16.*eta1*eta2*eta3;
  features[2] = 4.*(eta1*eta1 + eta2*eta2 + eta3*eta3);
  features[3] = 64.*(eta1*eta1*eta2*eta2 + eta2*eta2*eta3*eta3 + eta3*eta3*eta1*eta1);
}

template<typename T>
void chem_features_der(T c,T eta1,T eta2,T eta3,std::vector<std::vector<T> > &dhdeta){

  dhdeta.resize(4,std::vector<T>(4,0.));
  dhdeta[0][0] = 1.;
  dhdeta[1][1] = 16.*eta2*eta3;
  dhdeta[1][2] = 16.*eta1*eta3;
  dhdeta[1][3] = 16.*eta1*eta2;
  dhdeta[2][1] = 8.*eta1;
  dhdeta[2][2] = 8.*eta2;
  dhdeta[2][3] = 8.*eta3;
  dhdeta[3][1] = 64.*2.*eta1*(eta2*eta2+eta3*eta3);
  dhdeta[3][2] = 64.*2.*eta2*(eta1*eta1+eta3*eta3);
  dhdeta[3][3] = 64.*2.*eta3*(eta1*eta1+eta2*eta2);
}

template<typename T>
void strain_features(const std::vector<T> &e,std::vector<T> &features){

  features.resize(6);
  features[0] = e[0];
  features[1] = std::sqrt(0.5)*(e[1]*e[1] + e[2]*e[2]);
  features[2] = std::sqrt(1./3.)*(e[3]*e[3] + e[4]*e[4] + e[5]*e[5]);
  features[3] = 0.5*(e[2]*e[2]*e[2] - 3.*e[2]*e[1]*e[1]);
  features[4] = e[2]*(2.*e[3]*e[3] - e[4]*e[4] - e[5]*e[5])/2. - std::sqrt(3.)*e[1]*(e[4]*e[4] - e[5]*e[5])/2.;
  features[5] = std::sqrt(6.)*e[3]*e[4]*e[5];
}

template<typename T>
void strain_features_der(const std::vector<T> &e,std::vector<std::vector<T> > &dhde){

  dhde.resize(6,std::vector<T>(6,0.));
  dhde[0][0] = 1.;
  dhde[1][1] = 2.*std::sqrt(0.5)*e[1];
  dhde[1][2] = 2.*std::sqrt(0.5)*e[2];
  dhde[2][3] = 2.*std::sqrt(1./3.)*e[3];
  dhde[2][4] = 2.*std::sqrt(1./3.)*e[4];
  dhde[2][5] = 2.*std::sqrt(1./3.)*e[5];
  dhde[3][1] = -3.*e[2]*e[1];
  dhde[3][2] = 1.5*(e[2]*e[2] - e[1]*e[1]);
  dhde[4][1] = -std::sqrt(3.)/2.*(e[4]*e[4] - e[5]*e[5]);
  dhde[4][2] = 0.5*(2.*e[3]*e[3] - e[4]*e[4] - e[5]*e[5]);
  dhde[4][3] = 2.*e[2]*e[3];
  dhde[4][4] = -e[2]*e[4] - std::sqrt(3.)*e[1]*e[4];
  dhde[4][5] = -e[2]*e[5] + std::sqrt(3.)*e[1]*e[5];
  dhde[5][3] = std::sqrt(6.)*e[4]*e[5];
  dhde[5][4] = std::sqrt(6.)*e[3]*e[5];
  dhde[5][5] = std::sqrt(6.)*e[3]*e[4];
}

template<typename T>
void strain_params(Table<2,T> Ee,std::vector<T> &e){

  e.resize(6);
  if (Ee.size(0) == 3){
    e[0] = (Ee[0][0] + Ee[1][1] + Ee[2][2])/std::sqrt(3.);
    e[1] = (Ee[0][0] - Ee[1][1])/std::sqrt(2.);
    e[2] = (2.*Ee[2][2] - Ee[0][0] - Ee[1][1])/std::sqrt(6.);
    e[3] = std::sqrt(2.)*Ee[1][2];
    e[4] = std::sqrt(2.)*Ee[0][2];
    e[5] = std::sqrt(2.)*Ee[0][1];
  }
  else if (Ee.size(0) == 2){
    e[0] = (Ee[0][0] + Ee[1][1])/std::sqrt(3.);
    e[1] = (Ee[0][0] - Ee[1][1])/std::sqrt(2.);
    e[2] = -(Ee[0][0] + Ee[1][1])/std::sqrt(6.);
    e[3] = 0.;
    e[4] = 0.;
    e[5] = std::sqrt(2.)*Ee[0][1];
  }
}

template<typename T,int dim>
void strain_params_der(std::vector<Table<2,T> > &dedE,std::vector<T> e){

  dedE.resize(6,Table<2,T>(dim,dim));
  if (dim == 3){
    dedE[0][0][0] = dedE[0][1][1] = dedE[0][2][2] = 1./std::sqrt(3.);
    dedE[1][0][0] = 1./std::sqrt(2.);
    dedE[1][1][1] = -1./std::sqrt(2.);
    dedE[2][0][0] = dedE[2][1][1] = -1./std::sqrt(6.);
    dedE[2][2][2] = 2./std::sqrt(6.);
    dedE[3][1][2] = dedE[3][2][1] = std::sqrt(2.);
    dedE[4][0][2] = dedE[4][2][0] = std::sqrt(2.);
    dedE[5][0][1] = dedE[5][1][0] = std::sqrt(2.);
  }
  else if (dim == 2){
    dedE[0][0][0] = dedE[0][1][1] = 1./std::sqrt(3.);
    dedE[1][0][0] = 1./std::sqrt(2.);
    dedE[1][1][1] = -1./std::sqrt(2.);
    dedE[2][0][0] = dedE[2][1][1] = -1./std::sqrt(6.);
    dedE[5][0][1] = dedE[5][1][0] = std::sqrt(2.);
  }
}

template<typename T>
void SVK_eval(std::vector<T> h,std::vector<T> &p_Ni,std::vector<std::vector<T> > &dpNi_dh){

  double C11 = 222., C12 = 152.5, C44 = 112.; 
  //0.5*(C11 + 2C12)*h1^2 + 0.5*(C11 - C12)*sqrt(2)*h2 + C44*sqrt(3)*h3
  p_Ni.resize(1);
  dpNi_dh.resize(6,std::vector<T>(1,0.));
  p_Ni[0] = 0.5*(C11 + 2.*C12)*h[0]*h[0] + 0.5*(C11 - C12)*std::sqrt(2.)*h[1] + C44*std::sqrt(3.)*h[2];
  dpNi_dh[0][0] = (C11 + 2.*C12)*h[0];
  dpNi_dh[1][0] = 0.5*(C11 - C12)*std::sqrt(2.);
  dpNi_dh[2][0] = C44*std::sqrt(3.);

}


/**************************************
 * Define the physics
 *************************************/
template <int dim>
void precipitate_PF<dim>::get_residual(const typename hp::DoFHandler<dim>::active_cell_iterator &cell,
				       const FEValues<dim>& fe_values, 
				       Table<1, Sacado::Fad::DFad<double> >& R, 
				       Table<1, Sacado::Fad::DFad<double> >& ULocal, 
				       Table<1, double >& ULocalConv)
{
  n_q_points = fe_values.n_quadrature_points;

  /**************************************
   * Define variables for fields and derivatives
   *************************************/
  //Cahn-Hilliard
  dealii::Table<1,double>  c_conv(n_q_points);
  dealii::Table<1,Sacado::Fad::DFad<double> > c(n_q_points), mu(n_q_points), rhs_mu(n_q_points);
  dealii::Table<2,Sacado::Fad::DFad<double> > c_grad(n_q_points, dim), j_c(n_q_points, dim),
    mu_grad(n_q_points, dim), kappa1_cgrad(n_q_points, dim);

  //Allen-Cahn
  dealii::Table<1,double>  eta1_conv(n_q_points), eta2_conv(n_q_points), eta3_conv(n_q_points);
  dealii::Table<1,Sacado::Fad::DFad<double> > eta1(n_q_points), eta2(n_q_points), eta3(n_q_points);
  dealii::Table<2,Sacado::Fad::DFad<double> > eta1_grad(n_q_points, dim), j_eta1(n_q_points, dim),
    eta2_grad(n_q_points, dim), j_eta2(n_q_points, dim),
    eta3_grad(n_q_points, dim), j_eta3(n_q_points, dim);
  dealii::Table<1,Sacado::Fad::DFad<double> > reacAC1(n_q_points);
  dealii::Table<1,Sacado::Fad::DFad<double> > reacAC2(n_q_points);
  dealii::Table<1,Sacado::Fad::DFad<double> > reacAC3(n_q_points);

  evaluateScalarFunction<double,dim>(fe_values, c_dof, ULocalConv, c_conv);
  evaluateScalarFunction<Sacado::Fad::DFad<double>,dim>(fe_values, c_dof, ULocal, c);	
  evaluateScalarFunctionGradient<Sacado::Fad::DFad<double>,dim>(fe_values, c_dof, ULocal, c_grad);//at current configuration 

  evaluateScalarFunction<Sacado::Fad::DFad<double>,dim>(fe_values, mu_dof, ULocal, mu);	
  evaluateScalarFunctionGradient<Sacado::Fad::DFad<double>,dim>(fe_values, mu_dof, ULocal, mu_grad);//at current configuration 

  evaluateScalarFunction<double,dim>(fe_values, eta1_dof, ULocalConv, eta1_conv);
  evaluateScalarFunction<Sacado::Fad::DFad<double>,dim>(fe_values, eta1_dof, ULocal, eta1);	
  evaluateScalarFunctionGradient<Sacado::Fad::DFad<double>,dim>(fe_values, eta1_dof, ULocal, eta1_grad);//at current configuration 

  evaluateScalarFunction<double,dim>(fe_values, eta2_dof, ULocalConv, eta2_conv);
  evaluateScalarFunction<Sacado::Fad::DFad<double>,dim>(fe_values, eta2_dof, ULocal, eta2);	
  evaluateScalarFunctionGradient<Sacado::Fad::DFad<double>,dim>(fe_values, eta2_dof, ULocal, eta2_grad);//at current configuration 

  evaluateScalarFunction<double,dim>(fe_values, eta3_dof, ULocalConv, eta3_conv);
  evaluateScalarFunction<Sacado::Fad::DFad<double>,dim>(fe_values, eta3_dof, ULocal, eta3);	
  evaluateScalarFunctionGradient<Sacado::Fad::DFad<double>,dim>(fe_values, eta3_dof, ULocal, eta3_grad);//at current configuration 
  
  //mechanics
  deformationMap<Sacado::Fad::DFad<double>, dim> defMap(n_q_points); 
  getDeformationMap<Sacado::Fad::DFad<double>, dim>(fe_values, u_dof, ULocal, defMap);
  dealii::Table<3, Sacado::Fad::DFad<double> > P(n_q_points,dim,dim);
  
  //Loop over quadrature points
  for (unsigned int q=0; q<n_q_points; ++q){

  /**************************************
   * Define mechanics
   *************************************/
    
    Table<2, Sacado::Fad::DFad<double> > F (dim, dim);
    Table<2, Sacado::Fad::DFad<double> > Ee (dim, dim), S (dim, dim);
    Sacado::Fad::DFad<double> lambda = a(c[q])/a(c_avg), tmp1;
    //Sacado::Fad::DFad<double> lambda = 1., tmp1;
    
    F = copy_from_q(defMap.F,q);
    tmp1 = 1./(lambda*lambda);
    Ee = 0.5*(tmp1*(transpose(F)*F) - eye);

    std::vector<Sacado::Fad::DFad<double> > e, h, p_Ni, p_Ni3Al;
    std::vector<std::vector<Sacado::Fad::DFad<double> > > dpNi_dh, dpNi3Al_dh, dh_de;
    std::vector<Table<2,Sacado::Fad::DFad<double> > > de_dEe;

    strain_params(Ee,e);
    strain_params_der<Sacado::Fad::DFad<double>,dim>(de_dEe,e);
    strain_features(e,h);
    strain_features_der(e,dh_de);
    psi_Ni.eval(h,p_Ni,dpNi_dh); // 0.5*(C11 + 2C12)e1^2 + 0.5*(C11 - C12)(e2^2 + e3^2) + C44(e4^2 + e5^2 + e6^2) = 0.5*(C11 + 2C12)*h1^2 + 0.5*(C11 - C12)*sqrt(2)*h2 + C44*sqrt(3)*h3
    psi_Ni3Al.eval(h,p_Ni3Al,dpNi3Al_dh);
    //SVK_eval(h,p_Ni,dpNi_dh);
    //SVK_eval(h,p_Ni3Al,dpNi3Al_dh);
    for (unsigned int k=0; k<dim; ++k){
      for (unsigned int l=0; l<dim; ++l){
	S[k][l] = 0.;
	for (unsigned int i=0; i<6; ++i){
	  for (unsigned int j=0; j<6; ++j){
	    S[k][l] += ((1. - 4.*c[q])*dpNi_dh[i][0] + (4.*c[q])*dpNi3Al_dh[i][0])*dh_de[i][j]*de_dEe[j][k][l];
	    //S[k][l] += dpNi3Al_dh[i][0]*dh_de[i][j]*de_dEe[j][k][l];
	  }
	}
      }
    }
    copy_to_q(P,tmp1*(F*S),q); //P = Fe*S, but we solve P{F^\lambda}^{-1}, where Fe = (1/lambda)*F and {F^\lambda}^{-1} = (1/lambda)\identity
    
  /**************************************
   * Define chemistry
   *************************************/
     
    Sacado::Fad::DFad<double> lambda_c = a_c()/a(c_avg);
    Table<2, Sacado::Fad::DFad<double> > Ee_c;
    tmp1 = -std::pow(lambda,-3)*lambda_c;
    Ee_c = tmp1*(transpose(F)*F);
    
    std::vector<Sacado::Fad::DFad<double> > features(4), y;
    std::vector<std::vector<Sacado::Fad::DFad<double> > > dy_dx, dhdeta;

    chem_features(c[q],eta1[q],eta2[q],eta3[q],features);
    freeEn.eval(features,y,dy_dx);
    
    Sacado::Fad::DFad<double> f_c=0., psi_c,
      f_eta1=0., f_eta2=0., f_eta3=0.;
    
    psi_c = 4.*(p_Ni3Al[0] - p_Ni[0]);
    psi_c += double_contract(S,Ee_c);

    chem_features_der(c[q],eta1[q],eta2[q],eta3[q],dhdeta);
    for (unsigned int i=0; i<4; ++i){
      f_c += 14.2*dy_dx[i][0]*dhdeta[i][0];  //~14.2 GPa/(eV/atom)
      f_eta1 += 14.2*dy_dx[i][0]*dhdeta[i][1];
      f_eta2 += 14.2*dy_dx[i][0]*dhdeta[i][2];
      f_eta3 += 14.2*dy_dx[i][0]*dhdeta[i][3];
    }
    f_c += psi_c;

    //Allen-Cahn
    copy_to_q(j_eta1,-L*kappa2*copy_from_q(eta1_grad,q),q);
    copy_to_q(j_eta2,-L*kappa2*copy_from_q(eta2_grad,q),q);
    copy_to_q(j_eta3,-L*kappa2*copy_from_q(eta3_grad,q),q);
    reacAC1[q] = -L*f_eta1;
    reacAC2[q] = -L*f_eta2;
    reacAC3[q] = -L*f_eta3;

    //Cahn-Hilliard
    rhs_mu[q] = f_c - mu[q];    
    copy_to_q(j_c,-M*copy_from_q(mu_grad,q),q);
    copy_to_q(kappa1_cgrad,kappa1*copy_from_q(c_grad,q),q);
    
  }


  /**************************************
   * Call residual functions
   *************************************/
  this->ResidualEq.residualForDiffusionEq(fe_values, c_dof, R, c, c_conv, j_c); //CahnHilliard (part 1)
  this->ResidualEq.residualForPoissonEq(fe_values, mu_dof, R, kappa1_cgrad, rhs_mu); //CahnHilliard (part 2)
  this->ResidualEq.residualForDiff_ReacEq(fe_values, eta1_dof, R, eta1, eta1_conv, j_eta1, reacAC1); //Allen-Cahn
  this->ResidualEq.residualForDiff_ReacEq(fe_values, eta2_dof, R, eta2, eta2_conv, j_eta2, reacAC2); //Allen-Cahn
  this->ResidualEq.residualForDiff_ReacEq(fe_values, eta3_dof, R, eta3, eta3_conv, j_eta3, reacAC3); //Allen-Cahn
  this->ResidualEq.residualForMechanics(fe_values, u_dof, R, P); //Finite strain elasticity (need to correct this)
}

/**************************************
 * Define Dirichlet boundary conditions
 *************************************/
template <int dim>
void precipitate_PF<dim>::setup_constraints(){

  this->pcout<<"setup_constraints"<<std::endl;
  hpFEM<dim>::constraints.clear ();
  DoFTools::make_hanging_node_constraints (this->dof_handler, hpFEM<dim>::constraints);
  
  int totalDOF=this->totalDOF(this->primary_variables);
  std::vector<bool> elas (totalDOF, false);
  for (unsigned int i=0; i<dim; ++i){
    elas[totalDOF-1-i]=true;
  } 
  //apply constraints on boundary
  VectorTools:: interpolate_boundary_values (hpFEM<dim>::dof_handler, 1, ZeroFunction<dim> (totalDOF),hpFEM<dim>::constraints, elas);
  VectorTools:: interpolate_boundary_values (hpFEM<dim>::dof_handler, 2, ZeroFunction<dim> (totalDOF),hpFEM<dim>::constraints, elas);
  VectorTools:: interpolate_boundary_values (hpFEM<dim>::dof_handler, 3, ZeroFunction<dim> (totalDOF),hpFEM<dim>::constraints, elas);
  VectorTools:: interpolate_boundary_values (hpFEM<dim>::dof_handler, 4, ZeroFunction<dim> (totalDOF),hpFEM<dim>::constraints, elas);
  if (dim == 3){
       VectorTools:: interpolate_boundary_values (hpFEM<dim>::dof_handler, 5, ZeroFunction<dim> (totalDOF),hpFEM<dim>::constraints, elas);
       VectorTools:: interpolate_boundary_values (hpFEM<dim>::dof_handler, 6, ZeroFunction<dim> (totalDOF),hpFEM<dim>::constraints, elas);
  }
  
  hpFEM<dim>::constraints.close ();

}

/**************************************
 * Locally refine the mesh once
 *************************************/
template <int dim>
void precipitate_PF<dim>::refine_grid_once(){

  /*
  if(this->current_increment < 40){
    max_refine_level = 3;
  }
  else{
    max_refine_level = 5;
  }
  */

  dealii::Vector<double> localized_U(this->solution);
  Vector<float> estimated_error_per_cell (hpFEM<dim>::triangulation.n_active_cells());
  KellyErrorEstimator<dim>::estimate (hpFEM<dim>::dof_handler,
				      QGauss<dim-1>(3),
				      typename FunctionMap<dim>::type(),
				      localized_U,
				      estimated_error_per_cell);
  GridRefinement::refine_and_coarsen_fixed_fraction (hpFEM<dim>::triangulation, estimated_error_per_cell, 0.25, 0.12);
  if (hpFEM<dim>::triangulation.n_levels() > max_refine_level){
    for (typename Triangulation<dim>::active_cell_iterator cell = hpFEM<dim>::triangulation.begin_active(max_refine_level); cell != hpFEM<dim>::triangulation.end(); ++cell){
      if (cell->is_locally_owned()){
	cell->clear_refine_flag ();
      }
    }
  }

  SolutionTransfer<dim,dealii::Vector<double>,hp::DoFHandler<dim> > solution_trans(hpFEM<dim>::dof_handler);
  hpFEM<dim>::triangulation.prepare_coarsening_and_refinement();

  solution_trans.prepare_for_coarsening_and_refinement(localized_U);
  hpFEM<dim>::triangulation.execute_coarsening_and_refinement();

  this->mark_boundary();
  this->setup_linear_system();
  this->setup_constraints();

  dealii::Vector<double> localized_Unew(this->solution);
  solution_trans.interpolate(localized_U,localized_Unew);
  this->solution = localized_Unew;
  hpFEM<dim>::constraints.distribute (this->solution);

  this->pcout << "   Number of active cells:       " << hpFEM<dim>::triangulation.n_active_cells() << std::endl;
  this->pcout << "   Number of degrees of freedom: " << hpFEM<dim>::dof_handler.n_dofs() << std::endl; 

}

/**************************************
 * Initial refinement of the mesh
 *************************************/
template <int dim>
void precipitate_PF<dim>::refine_grid(){

  hpFEM<dim>::triangulation.refine_global(n_init_global_refine);

  this->setup_linear_system();
  this->apply_initial_condition(); // Random ICs - don't apply after refinements
  //hpFEM<dim>::triangulation.refine_global(n_init_global_refine);
  /*
  for(unsigned int i=0; i<n_init_local_refine; ++i){
    this->setup_linear_system();
    this->apply_initial_condition();
    
    refine_grid_once();
  }
  / */
}

/**************************************
 * Implement adaptive time stepping
 *************************************/
template <int dim>
void precipitate_PF<dim>::solve_ibvp()
{
  int reduce_time=0;
  int converge_iter;
  while(reduce_time<8) {
    converge_iter=this->nonlinearSolve(this->solution);
    //reset iter_count to 0
    if (converge_iter>0) break;
    else{
      iter_count=0;
      this->pcout<<"not converge, reduce dt by half"<<std::endl;
      this->current_dt *= 0.5;
    }
  }
  if(converge_iter<4){
    iter_count++;
  }
  else{
    iter_count = 0;
  }
  if(iter_count>3) {
    iter_count=0;
    this->current_dt *= 2.;//sqrt(2.);//double dt
  }

  //*
  if(this->current_increment < 6 || this->current_increment%3 == 0){ 
    refine_grid_once();
  }
  // */
  /*
  if(this->current_increment > 4 && (this->current_increment < 10 || this->current_increment%3 == 0)){ 
    refine_grid_once();
  }
  */

  this->solution_prev=this->solution;
}

/**************************************
 * Adaptive mesh refinement
 *************************************/
template <int dim>
void precipitate_PF<dim>::update_post_TS(){

  //if(this->current_increment > 10){
  //  this->current_dt *= 1.3;
  //}
  /*
  if(this->current_increment < 6 || this->current_increment%3 == 0){ 
    refine_grid_once();
  }
  */
}

template <int dim>
void InitialConditions<dim>::vector_value (const Point<dim>   &p, Vector<double>   &values) const{
  Assert (values.size() == totalDOF, ExcDimensionMismatch (values.size(), totalDOF));
  /**************************************
   * Read in necessary parameters
   *************************************/
  params->enter_subsection("parameters");
  double c_avg=params->get_double("c_avg");
  params->leave_subsection();

  //*
  {
    values(0) = c_avg  + 0.03*(0.5 - (double)(rand() % 100 )/100.0); //composition
    values(1) = 0.; //mu
    values(2) = 0.02*(0.5 - (double)(rand() % 100 )/100.0); //order parameter
    values(3) = 0.02*(0.5 - (double)(rand() % 100 )/100.0); //order parameter
    values(4) = 0.02*(0.5 - (double)(rand() % 100 )/100.0); //order parameter
  }
  // */

  /**************************************
   * Define initial conditions for composition, order parameter (spherical precipitates)
   *************************************/
  /*
  // Initial precipitate radius
  double r = 1.;
  double r1 = r - 0.1, r2 = r + 0.1;

  // Define number and location of initial precipitate seeds
  unsigned int n_points = 1;
  std::vector<double> x_points = {-0.85583912};
  std::vector<double> y_points = {-3.99068818};
  std::vector<double> z_points = {0.};
  std::vector<int> var = {1};

  // Initial precipitate and matrix compositions for a given average composition
  double c_p = 0.226;
  double c_matrix = 0.07;//c_avg*domainVol - c_p*n_points*(4./3.*PI*std::pow(r2,3) - 4.*PI/(r2-r1)*(1./4.*std::pow(r2,4) - 1./3.*r1*std::pow(r2,3) + 1./12.*std::pow(r1,4)));
  //c_matrix /= domainVol - n_points*(4./3.*PI*std::pow(r2,3) - 4.*PI/(r2-r1)*(1./4.*std::pow(r2,4) - 1./3.*r1*std::pow(r2,3) + 1./12.*std::pow(r1,4)));

  double eta1_p, eta2_p, eta3_p;

  // Assign values
  double tmp_c, tmp_eta1, tmp_eta2, tmp_eta3;
  values(0) = c_matrix; //composition
  values(1) = 0.;
  values(2) = 0.; //order parameter
  values(3) = 0.; //order parameter
  values(4) = 0.; //order parameter
  for(unsigned int i=0; i<n_points; ++i){

    if(var[i] == 1){
      eta1_p = 0.22;
      eta2_p = 0.22;
      eta3_p = 0.22;
    }
    if(var[i] == 2){
      eta1_p = -0.22;
      eta2_p = 0.22;
      eta3_p = -0.22;
    }
    if(var[i] == 3){
      eta1_p = -0.22;
      eta2_p = -0.22;
      eta3_p = 0.22;
    }
    if(var[i] == 4){
      eta1_p = 0.22;
      eta2_p = -0.22;
      eta3_p = -0.22;
    }

    Point<dim> center(x_points[i],y_points[i]);//,z_points[i]);
    double rho = p.distance(center);
    if(rho < r1){ //Completely inside precipitate
      values(0) = c_p; //c
      values(2) = eta1_p; //eta1
      values(3) = eta2_p; //eta2
      values(4) = eta3_p; //eta3
      break;
    }
    else if(rho < r2){ //Within a (linearly) diffuse interface
      tmp_c = (c_matrix - c_p)/(r2 - r1)*(rho - r1) + c_p;
      tmp_eta1 = (0.0 - eta1_p)/(r2 - r1)*(rho - r1) + eta1_p;
      tmp_eta2 = (0.0 - eta2_p)/(r2 - r1)*(rho - r1) + eta2_p;
      tmp_eta3 = (0.0 - eta3_p)/(r2 - r1)*(rho - r1) + eta3_p;

      values(0) = tmp_c ;//c
      values(2) = tmp_eta1; //eta1
      values(3) = tmp_eta2; //eta2
      values(4) = tmp_eta3; //eta3
    }
  }
  */

}

template class InitialConditions<1>;
template class InitialConditions<2>;
template class InitialConditions<3>;
