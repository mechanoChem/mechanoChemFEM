#include <deal.II/numerics/error_estimator.h>
#include <deal.II/grid/grid_in.h>
#include <deal.II/base/mpi.h>
#include <random>
#include "mechanoChemFEM.h"
#include "DNN.h"
//#include "rand_grid.h"
#include "rand_seed.h"

using namespace dealii;

//Initial conditions
template <int dim>
class InitCond: public Function<dim>{
public:
  InitCond (int _totalDOF,
	    std::vector<std::vector<std::string> >& _primary_variables,
	    std::vector<unsigned int >& _primary_variables_dof,
	    ParameterHandler& _params);
  int totalDOF;
  std::vector<std::vector<std::string> > primary_variables;
  std::vector<unsigned int > primary_variables_dof;
  ParameterHandler* params;
  virtual double value(const Point<dim>   &p, const unsigned int 	component=0) const;
  virtual void vector_value (const Point<dim>   &p, Vector<double>   &values) const;
  //RandGrid rand_grid, rand_grid2;
  RandSeed rand_seed;
  std::vector<std::vector<double> > invQ;

};

template <int dim>
InitCond<dim>::InitCond (int _totalDOF,
			 std::vector<std::vector<std::string> >& _primary_variables,
			 std::vector<unsigned int >& _primary_variables_dof,
			 ParameterHandler& _params)
  :
  Function<dim>(_totalDOF),
  primary_variables(_primary_variables),
  primary_variables_dof(_primary_variables_dof),
  totalDOF(_totalDOF),
  params(&_params),
  rand_seed(20),
  invQ(32,std::vector<double>(7))
{
  std::ifstream file;
  file.open("../invQ_LCO.txt");
  for(unsigned int i=0; i<32; ++i){
    for(unsigned int j=0; j<7; ++j){
      file >> invQ[i][j];
    }
  }
  file.close();
  
}

template <int dim>
double InitCond<dim>::value(const Point<dim>   &p, const unsigned int 	component) const{
  return 0;   
}

template <int dim>
void InitCond<dim>::vector_value (const Point<dim>   &p, Vector<double>   &values) const{
  Assert (values.size() == 8, ExcDimensionMismatch (values.size(), 8));

  // A couple approaches:
  /* Random, coarse field: define a random field of numbers on a coarse grid (one-time, separately).
     Simple to convert from position to indices on the grid. Can interpolate.
     Random circles: similarly, define random set of circle centers. For each position,
     check if lies within radius of any center. Also need to previously define random
     index for each circle to define which eta is nonzero. 
  */
  
  //values(0)= 0.5 + 0.0*(static_cast <double> (rand())/(static_cast <double>(RAND_MAX))-0.5); // c (eta0)
  values(0)= 0.53 + 0.0*(static_cast <double> (rand())/(static_cast <double>(RAND_MAX))-0.5); // c (eta0)
  //values(0)= 0.525 + rand_grid2.eval(p[0],p[1],0); // c (eta0)
  values(1) = 0.0; // mu_c
  /*
  values(2) = 0.0*(static_cast <double> (rand())/(static_cast <double>(RAND_MAX))-0.5); //eta1
  values(3) = 0.0*(static_cast <double> (rand())/(static_cast <double>(RAND_MAX))-0.5); //eta2
  values(4) = 0.0*(static_cast <double> (rand())/(static_cast <double>(RAND_MAX))-0.5); //eta3
  values(5) = 0.0*(static_cast <double> (rand())/(static_cast <double>(RAND_MAX))-0.5); //eta4
  values(6) = 0.0*(static_cast <double> (rand())/(static_cast <double>(RAND_MAX))-0.5); //eta5
  values(7) = 0.0*(static_cast <double> (rand())/(static_cast <double>(RAND_MAX))-0.5); //eta6
  */

  /*
  if (sqrt(p[0]*p[0] + p[1]*p[1]) < 0.2){
    values(0) = 0.5;
    values(2) = 0.4;
  }
  else{
    values(0) = 0.53;
    values(2) = 0.;
  }
  */
  //int ind = (std::rand() % 6) + 2;
  //int ind = 2;//(std::rand() % 2) + 2;

  //values(ind) += 0.9*(static_cast <double> (rand())/(static_cast <double>(RAND_MAX))-0.5);
  //values(ind) += 0.4*(static_cast <double> (rand())/(static_cast <double>(RAND_MAX)));
  //values(ind) = 0.15/(exp(-sqrt(2/2.e-5)*(0.075 - sqrt(p[0]*p[0] + p[1]*p[1]))) + 1.);
  //values(2) = 0.9*(static_cast <double> (rand())/(static_cast <double>(RAND_MAX)) - 0.5);//
  values(2) = rand_seed.eval(p[0],p[1],0);
  values(3) = rand_seed.eval(p[0],p[1],1);
  values(4) = rand_seed.eval(p[0],p[1],2);
  values(5) = rand_seed.eval(p[0],p[1],3);
  values(6) = rand_seed.eval(p[0],p[1],4);
  values(7) = rand_seed.eval(p[0],p[1],5);


  std::vector<int> conv = {0,2,3,4,5,6,7};
  // Now, make sure it lies in a physically realistic range
  double check;
  for(unsigned int i=0; i<invQ.size(); ++i){
    check = 0.;
    for(unsigned int j=0; j<invQ[0].size(); ++j){
      check += invQ[i][j]*values(conv[j]);
    }
    if ((check > 1) || (check < 0)){
      std::cout << "Trouble!\n";
    }
  }

  /*
    if (p[1] > 0){
    values(ind) += 0.3;//0.*(static_cast <double> (rand())/(static_cast <double>(RAND_MAX))-0.5);
    }
    else{
    values(ind) += -0.3;//0.*(static_cast <double> (rand())/(static_cast <double>(RAND_MAX))-0.5);
    }
  */
}

template <int dim>
class CahnHilliard: public mechanoChemFEM<dim>
{
public:
  CahnHilliard();
  //this is a overloaded function 
  void get_residual(const typename hp::DoFHandler<dim>::active_cell_iterator &cell,
		    const FEValues<dim>& fe_values,
		    Table<1, Sacado::Fad::DFad<double> >& R,
		    Table<1, Sacado::Fad::DFad<double>>& ULocal,
		    Table<1, double >& ULocalConv);
  void make_grid();
  void apply_initial_condition();  
  void refine_grid_once();
  void ini_updateLinearSystem();
  void solve_ibvp();
  bool checkBounds();

  DNN freeEn;
  int iter_count=0;
  int max_refine_level=4;

  ParameterHandler* params;

  int c_dof, mu_dof, eta1_dof;

  std::vector<std::vector<std::vector<double> > > current_noise;
  FEMdata<dim, Vector<double> > FEMdata_refine;

  std::vector<std::vector<double> > invQ;
  //bool outOfBounds;

  std::vector<double> gradEnergy_eta;
  double gradEnergy_0;
  double c_avg;
  double area;
};

template <int dim>
CahnHilliard<dim>::CahnHilliard()
  :
  invQ(32,std::vector<double>(7)),
  gradEnergy_eta(6)
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

  freeEn.reinit(3,"T300_7D_",true);

  std::ifstream file;
  file.open("../invQ_LCO.txt");
  for(unsigned int i=0; i<32; ++i){
    for(unsigned int j=0; j<7; ++j){
      file >> invQ[i][j];
    }
  }
  file.close();

  
  c_dof=0; mu_dof = 1; eta1_dof=2;
}

template<typename T>
void chem_features(T c,T e1,T e2,T e3,T e4,T e5,T e6,std::vector<T> &features){

  features.resize(12);
  features[0] = c;
  features[1] = 2./3.*(e1*e1 + e2*e2 + e3*e3 +
		       e4*e4 + e5*e5 + e6*e6);
  features[2] = 8./3.*(pow(e1,4) + pow(e2,4) + pow(e3,4) +
		       pow(e4,4) + pow(e5,4) + pow(e6,4));
  features[3] = 4./3.*((e1*e1 + e2*e2)*
		       (e3*e3 + e4*e4 + e5*e5 + e6*e6) +
		       (e3*e3 + e6*e6)*(e4*e4 + e5*e5));
  features[4] = 16./3.*(e1*e1*e2*e2 + e3*e3*e6*e6 + e4*e4*e5*e5);
  features[5] = 32./3.*(pow(e1,6) + pow(e2,6) + pow(e3,6) +
			pow(e4,6) + pow(e5,6) + pow(e6,6));
  features[6] = 8./3.*((pow(e1,4) + pow(e2,4))*
		       (e3*e3 + e4*e4 + e5*e5 + e6*e6) +
		       (pow(e3,4) + pow(e6,4))*(e4*e4 + e5*e5) +
		       (e1*e1 + e2*e2)*
		       (pow(e3,4) + pow(e4,4) + pow(e5,4) + pow(e6,4)) +
		       (e3*e3 + e6*e6)*(pow(e4,4) + pow(e5,4)));
  features[7] = 16./3.*(e1*e1*e2*e2*(e3*e3 + e4*e4 + e5*e5 + e6*e6) +
			e3*e3*e6*e6*(e1*e1 + e2*e2 + e4*e4 + e5*e5) +
			e4*e4*e5*e5*(e1*e1 + e2*e2 + e3*e3 + e6*e6));
  features[8] = 32./3.*(e1*e1*pow(e2,4) + e3*e3*pow(e6,4) + e4*e4*pow(e5,4) +
			pow(e1,4)*e2*e2 + pow(e3,4)*e6*e6 + pow(e4,4)*e5*e5);
  features[9] = 8.*(e1*e1 + e2*e2)*(e3*e3 + e6*e6)*(e4*e4 + e5*e5);
  features[10] = 64./5.*(e1*e2*(e3*e3 - e6*e6)*(e4*e4 - e5*e5) +
			 e3*e6*(e1*e1 - e2*e2)*(e4*e4 - e5*e5) +
			 e4*e5*(e1*e1 - e2*e2)*(e3*e3 - e6*e6));
  features[11] = 64.*sqrt(5.)*e1*e2*e3*e4*e5*e6;
  
}

template<typename T>
void chem_features_der(T c,T e1,T e2,T e3,T e4,T e5,T e6,std::vector<std::vector<T> > &dhde){

  dhde.resize(12,std::vector<T>(7,0.));

  dhde[0][0] = 1.;

  dhde[1][1] = 4./3.*e1;
  dhde[1][2] = 4./3.*e2;
  dhde[1][3] = 4./3.*e3;
  dhde[1][4] = 4./3.*e4;
  dhde[1][5] = 4./3.*e5;
  dhde[1][6] = 4./3.*e6;
  
  dhde[2][1] = 32./3.*pow(e1,3);
  dhde[2][2] = 32./3.*pow(e2,3);
  dhde[2][3] = 32./3.*pow(e3,3);
  dhde[2][4] = 32./3.*pow(e4,3);
  dhde[2][5] = 32./3.*pow(e5,3);
  dhde[2][6] = 32./3.*pow(e6,3);

  dhde[3][1] = 8./3.*e1*(e3*e3 + e4*e4 + e5*e5 + e6*e6);
  dhde[3][2] = 8./3.*e2*(e3*e3 + e4*e4 + e5*e5 + e6*e6);
  dhde[3][3] = 8./3.*e3*(e1*e1 + e2*e2 + e4*e4 + e5*e5);
  dhde[3][4] = 8./3.*e4*(e1*e1 + e2*e2 + e3*e3 + e6*e6);
  dhde[3][5] = 8./3.*e5*(e1*e1 + e2*e2 + e3*e3 + e6*e6);
  dhde[3][6] = 8./3.*e6*(e1*e1 + e2*e2 + e4*e4 + e5*e5);
  
  dhde[4][1] = 32./3.*e1*e2*e2;
  dhde[4][2] = 32./3.*e1*e1*e2;
  dhde[4][3] = 32./3.*e3*e6*e6;
  dhde[4][4] = 32./3.*e4*e5*e5;
  dhde[4][5] = 32./3.*e4*e4*e5;
  dhde[4][6] = 32./3.*e3*e3*e6;

  dhde[5][1] = 64.*pow(e1,5);
  dhde[5][2] = 64.*pow(e2,5);
  dhde[5][3] = 64.*pow(e3,5);
  dhde[5][4] = 64.*pow(e4,5);
  dhde[5][5] = 64.*pow(e5,5);
  dhde[5][6] = 64.*pow(e6,5);

  dhde[6][1] = 8./3.*(4.*pow(e1,3)*(e3*e3 + e4*e4 + e5*e5 + e6*e6) +
		      2.*e1*(pow(e3,4) + pow(e4,4) + pow(e5,4) + pow(e6,4)));
  dhde[6][2] = 8./3.*(4.*pow(e2,3)*(e3*e3 + e4*e4 + e5*e5 + e6*e6) +
		      2.*e2*(pow(e3,4) + pow(e4,4) + pow(e5,4) + pow(e6,4)));
  dhde[6][3] = 8./3.*(4.*pow(e3,3)*(e1*e1 + e2*e2 + e4*e4 + e5*e5) +
		      2.*e3*(pow(e1,4) + pow(e2,4) + pow(e4,4) + pow(e5,4)));
  dhde[6][4] = 8./3.*(4.*pow(e4,3)*(e1*e1 + e2*e2 + e3*e3 + e6*e6) +
		      2.*e4*(pow(e1,4) + pow(e2,4) + pow(e3,4) + pow(e6,4)));
  dhde[6][5] = 8./3.*(4.*pow(e5,3)*(e1*e1 + e2*e2 + e3*e3 + e6*e6) +
		      2.*e5*(pow(e1,4) + pow(e2,4) + pow(e3,4) + pow(e6,4)));
  dhde[6][6] = 8./3.*(4.*pow(e6,3)*(e1*e1 + e2*e2 + e4*e4 + e5*e5) +
		      2.*e6*(pow(e1,4) + pow(e2,4) + pow(e4,4) + pow(e5,4)));
		      
  dhde[7][1] = 32./3.*e1*(e2*e2*(e3*e3 + e4*e4 + e5*e5 + e6*e6) +
			  e3*e3*e6*e6 + e4*e4*e5*e5);
  dhde[7][2] = 32./3.*e2*(e1*e1*(e3*e3 + e4*e4 + e5*e5 + e6*e6) +
			  e3*e3*e6*e6 + e4*e4*e5*e5);
  dhde[7][3] = 32./3.*e3*(e6*e6*(e1*e1 + e2*e2 + e4*e4 + e5*e5) +
			  e1*e1*e2*e2 + e4*e4*e5*e5);
  dhde[7][4] = 32./3.*e4*(e5*e5*(e1*e1 + e2*e2 + e3*e3 + e6*e6) +
			  e1*e1*e2*e2 + e3*e3*e6*e6);
  dhde[7][5] = 32./3.*e5*(e4*e4*(e1*e1 + e2*e2 + e3*e3 + e6*e6) +
			  e1*e1*e2*e2 + e3*e3*e6*e6);
  dhde[7][6] = 32./3.*e6*(e3*e3*(e1*e1 + e2*e2 + e4*e4 + e5*e5) +
			  e1*e1*e2*e2 + e4*e4*e5*e5);
  
  dhde[8][1] = 32./3.*(2.*e1*pow(e2,4) + 4.*pow(e1,3)*e2*e2);
  dhde[8][2] = 32./3.*(2.*e2*pow(e1,4) + 4.*pow(e2,3)*e1*e1);
  dhde[8][3] = 32./3.*(2.*e3*pow(e6,4) + 4.*pow(e3,3)*e6*e6);
  dhde[8][4] = 32./3.*(2.*e4*pow(e5,4) + 4.*pow(e4,3)*e5*e5);
  dhde[8][5] = 32./3.*(2.*e5*pow(e4,4) + 4.*pow(e5,3)*e4*e4);
  dhde[8][6] = 32./3.*(2.*e6*pow(e3,4) + 4.*pow(e6,3)*e3*e3);

  dhde[9][1] = 16.*e1*(e3*e3 + e6*e6)*(e4*e4 + e5*e5);
  dhde[9][2] = 16.*e2*(e3*e3 + e6*e6)*(e4*e4 + e5*e5);
  dhde[9][3] = 16.*e3*(e1*e1 + e2*e2)*(e4*e4 + e5*e5);
  dhde[9][4] = 16.*e4*(e1*e1 + e2*e2)*(e3*e3 + e6*e6);
  dhde[9][5] = 16.*e5*(e1*e1 + e2*e2)*(e3*e3 + e6*e6);
  dhde[9][6] = 16.*e6*(e1*e1 + e2*e2)*(e4*e4 + e5*e5);

  dhde[10][1] = 64./5.*(e2*(e3*e3 - e6*e6)*(e4*e4 - e5*e5) +
			2.*e1*(e3*e6*(e4*e4 - e5*e5) +
			       e4*e5*(e3*e3 - e6*e6)));
  dhde[10][2] = 64./5.*(e1*(e3*e3 - e6*e6)*(e4*e4 - e5*e5) -
			2.*e2*(e3*e6*(e4*e4 - e5*e5) +
			       e4*e5*(e3*e3 - e6*e6)));
  dhde[10][3] = 64./5.*(e6*(e1*e1 - e2*e2)*(e4*e4 - e5*e5) +
			2.*e3*(e1*e2*(e4*e4 - e5*e5) +
			       e4*e5*(e1*e1 - e2*e2)));
  dhde[10][4] = 64./5.*(e5*(e1*e1 - e2*e2)*(e3*e3 - e6*e6) +
			2.*e4*(e3*e6*(e1*e1 - e2*e2) +
			       e1*e2*(e3*e3 - e6*e6)));
  dhde[10][5] = 64./5.*(e4*(e1*e1 - e2*e2)*(e3*e3 - e6*e6) -
			2.*e5*(e3*e6*(e1*e1 - e2*e2) +
			       e1*e2*(e3*e3 - e6*e6)));
  dhde[10][6] = 64./5.*(e3*(e1*e1 - e2*e2)*(e4*e4 - e5*e5) -
			2.*e6*(e1*e2*(e4*e4 - e5*e5) +
			       e4*e5*(e1*e1 - e2*e2)));
  
  dhde[11][1] = 64.*sqrt(5.)*e2*e3*e4*e5*e6;
  dhde[11][2] = 64.*sqrt(5.)*e1*e3*e4*e5*e6;
  dhde[11][3] = 64.*sqrt(5.)*e1*e2*e4*e5*e6;
  dhde[11][4] = 64.*sqrt(5.)*e1*e2*e3*e5*e6;
  dhde[11][5] = 64.*sqrt(5.)*e1*e2*e3*e4*e6;
  dhde[11][6] = 64.*sqrt(5.)*e1*e2*e3*e4*e5;
  
}

template <int dim>
void CahnHilliard<dim>::get_residual(const typename hp::DoFHandler<dim>::active_cell_iterator &cell,
				     const FEValues<dim>& fe_values,
				     Table<1, Sacado::Fad::DFad<double> >& R,
				     Table<1, Sacado::Fad::DFad<double>>& ULocal,
				     Table<1, double >& ULocalConv)
{
  //evaluate primary fields
  params->enter_subsection("Concentration");
  //double M=params->get_double("M");
  Sacado::Fad::DFad<double> M, L;
  double kappa1=params->get_double("kappa"); //Cahn-Hilliard
  double kappa2=(1./256.)*params->get_double("kappa"); //Allen-Cahn
  params->leave_subsection();
  double jn = 0.;
  
  unsigned int n_q_points= fe_values.n_quadrature_points;

  //Cahn-Hilliard
  dealii::Table<1,double>  c_conv(n_q_points);
  dealii::Table<1,Sacado::Fad::DFad<double> > c(n_q_points), mu(n_q_points);
  dealii::Table<2,Sacado::Fad::DFad<double> >  c_grad(n_q_points, dim), mu_grad(n_q_points, dim);

  evaluateScalarFunction<double,dim>(fe_values, c_dof, ULocalConv, c_conv);
  evaluateScalarFunction<Sacado::Fad::DFad<double>,dim>(fe_values, c_dof, ULocal, c);	
  evaluateScalarFunctionGradient<Sacado::Fad::DFad<double>,dim>(fe_values, c_dof, ULocal, c_grad);
	
  evaluateScalarFunction<Sacado::Fad::DFad<double>,dim>(fe_values, mu_dof, ULocal, mu);	
  evaluateScalarFunctionGradient<Sacado::Fad::DFad<double>,dim>(fe_values, mu_dof, ULocal, mu_grad);
  
  //Allen-Cahn
  std::vector<dealii::Table<1,double> > eta_conv(6,dealii::Table<1,double>(n_q_points));
  std::vector<dealii::Table<1,Sacado::Fad::DFad<double> > > eta(6,dealii::Table<1,Sacado::Fad::DFad<double> >(n_q_points));
  std::vector<dealii::Table<2,Sacado::Fad::DFad<double> > > eta_grad(6,dealii::Table<2,Sacado::Fad::DFad<double> >(n_q_points, dim)),
    j_eta(6,dealii::Table<2,Sacado::Fad::DFad<double> >(n_q_points, dim));
  std::vector<dealii::Table<1,Sacado::Fad::DFad<double> > > reacAC(6,dealii::Table<1,Sacado::Fad::DFad<double> >(n_q_points));

  for (unsigned int i=0; i<6; i++){
    evaluateScalarFunction<double,dim>(fe_values, eta1_dof+i, ULocalConv, eta_conv[i]);
    evaluateScalarFunction<Sacado::Fad::DFad<double>,dim>(fe_values, eta1_dof+i, ULocal, eta[i]);	
    evaluateScalarFunctionGradient<Sacado::Fad::DFad<double>,dim>(fe_values, eta1_dof+i, ULocal, eta_grad[i]);//at current configuration 
  }
  
  //evaluate diffusion and reaction term
  dealii::Table<1,Sacado::Fad::DFad<double> > rhs_mu(n_q_points);
  dealii::Table<2,Sacado::Fad::DFad<double> > j_c(n_q_points, dim), kappa_c_grad(n_q_points, dim);

  //DNN free energy evaluation
  std::vector<Sacado::Fad::DFad<double> > features(12), f, dfdeta(7);
  std::vector<std::vector<Sacado::Fad::DFad<double> > > dfdh, dhdeta;
  
  //j_c=table_scaling<Sacado::Fad::DFad<double>, dim>(mu_grad,-M);//-D_1*c_grad
  double kB = 8.617333262145e-5; // eV/K
  double T = 300; // K
  for(unsigned int q=0; q<n_q_points; q++){

    /*
    std::vector<double> values = {c[q].val(),eta[0][q].val(),eta[1][q].val(),eta[2][q].val(),eta[3][q].val(),eta[4][q].val(),eta[5][q].val()};
    // Now, make sure it lies in a physically realistic range
    double check;
    for(unsigned int i=0; i<invQ.size(); ++i){
      check = 0.;
      for(unsigned int j=0; j<invQ[0].size(); ++j){
	check += invQ[i][j]*values[j];
      }
      if ((check > 1) || (check < 0)){
	outOfBounds = true;
	//this->pcout << "Out of bounds...\n";
	//break;
      }
    }
    */

    chem_features(c[q],eta[0][q],eta[1][q],eta[2][q],eta[3][q],eta[4][q],eta[5][q],features);
    /*
    Sacado::Fad::DFad<double> zero, eta1_eq, deta1_dc;
    zero = 0.;
    if (c[q] < 0.52){
      eta1_eq = 0.4;
      deta1_dc = 0.;
    }
    else if (c[q] > 0.54){
      eta1_eq = 0.;
      deta1_dc = 0.;
    }
    else{
      eta1_eq = 0.4*(3.*pow((0.54-c[q])/0.02,2) - 2.*pow((0.54-c[q])/0.02,3));
      deta1_dc = -0.4/0.02*(6.*(0.54-c[q])/0.02 - 6.*pow((0.54-c[q])/0.02,2));
    }
    */
    //chem_features(c[q],eta1_eq,zero,zero,zero,zero,zero,features);
    //chem_features(c[q],zero,zero,zero,zero,zero,zero,features);
    freeEn.eval(features,f,dfdh);

    // Get the derivative of the free energy w.r.t. each input
    // (use chain rule because of symmetry input functions)
    chem_features_der(c[q],eta[0][q],eta[1][q],eta[2][q],eta[3][q],eta[4][q],eta[5][q],dhdeta);
    //chem_features_der(c[q],eta1_eq,zero,zero,zero,zero,zero,dhdeta);
    //chem_features_der(c[q],zero,zero,zero,zero,zero,zero,dhdeta);
    for (unsigned int i=0; i<dfdeta.size(); ++i){
      dfdeta[i] = 0;
      for (unsigned int j=0; j<dhdeta.size(); ++j){
	dfdeta[i] += dfdh[j][0]*dhdeta[j][i];
      }
    }

    // Cahn-Hilliard
    M = 1.e-2*std::exp(-274.0*(1.05-c[q])*(0.47-c[q])*(1.-c[q]))*c[q]/(kB*T); //micron^2/(eV*s*unit_cell)
    L = M*1.e3;
    //L = 1.e5;
    for (int j=0; j<6; ++j){
      //L = std::min(L,100.*exp(-4.*log(1000.)*(eta[j][q] - 0.5)*(eta[j][q] + 0.5))-100.);
      //L = std::min(L,10.*exp(-4.*log(100.)*(eta[j][q] - 0.5)*(eta[j][q] + 0.5))-10.);
      //L = std::min(L,-4.*1000.*(eta[j][q] - 0.5)*(eta[j][q] + 0.5));
      //L = std::max(L,0.);
      //L = std::min(L,1.e3*exp(-4.*log(100.)*(eta[j][q] - 0.5)*(eta[j][q] + 0.5)));
    }
    //M = 0.;//L*1.e-3; //We'll improve this later; just a check for now

    for (int i=0; i<dim; ++i){
      j_c[q][i] = -M*mu_grad[q][i];
      kappa_c_grad[q][i] = kappa1*c_grad[q][i];
    }
    rhs_mu[q] = dfdeta[0] - mu[q];
    //rhs_mu[q] = (dfdeta[0] + dfdeta[1]*deta1_dc) - mu[q];

    // Allen-Cahn
    int cell_label = cell->user_index();
    for (int j=0; j<6; ++j){
      //L = 1.e3*exp(-4.*log(100.)*(eta[j][q] - 0.5)*(eta[j][q] + 0.5));
      for (int i=0; i<dim; ++i){
	j_eta[j][q][i] = -L*kappa2*eta_grad[j][q][i];
      }
      reacAC[j][q] = -L*dfdeta[j+1];// + sqrt(2.*kB*T*L)*current_noise[cell_label][j][q];
      //reacAC[j][q] = -L*dfdeta[j+1] + sqrt(2.*L)*current_noise[cell_label][j][q];
      //reacAC[j][q] = -L*(dfdeta[j+1] - 0.*current_noise[cell_label][j][q]);
    }

    //integrate interfacial energy
    for (int i=0; i<dim; ++i){
      gradEnergy_0 += 0.5*kappa1*c_grad[q][i].val()*c_grad[q][i].val()*fe_values.JxW(q);
      for (int j=0; j<6; ++j){
	gradEnergy_eta[j] += 0.5*kappa2*eta_grad[j][q][i].val()*eta_grad[j][q][i].val()*fe_values.JxW(q);
      }
    }
    c_avg += c[q].val()*fe_values.JxW(q);
    area += fe_values.JxW(q);
      
  }
	
  //call residual functions
  this->ResidualEq.residualForDiffusionEq(fe_values, c_dof, R, c, c_conv, j_c); //CahnHilliard (part 1)
  this->ResidualEq.residualForPoissonEq(fe_values, mu_dof, R, kappa_c_grad, rhs_mu); //CahnHilliard (part 2)
  for (int j=0; j<6; ++j){
    this->ResidualEq.residualForDiff_ReacEq(fe_values, eta1_dof+j, R, eta[j], eta_conv[j], j_eta[j], reacAC[j]); //Allen-Cahn
  }
  
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
  /*
    GridGenerator::hyper_ball (this->triangulation,Point<dim>(),R);
    static const SphericalManifold<dim> boundary;
    this->triangulation.set_all_manifold_ids_on_boundary(0);
    this->triangulation.set_manifold (0, boundary);
  */
  std::ifstream fin("initial_mesh.msh");
  GridIn<dim> grid_in;
  grid_in.attach_triangulation(this->triangulation);
  grid_in.read_msh(fin);
  this->triangulation.refine_global (n_refine);

  if (this->resuming_from_snapshot){
    Vector<double> est_err_in;
    
    for (int i=1; i<=this->current_increment+this->off_output_index; ++i){
      std::string refine_path = this->snapshot_directory+"refinement-"+std::to_string(i)+".dat";
      std::ifstream in(refine_path);
      est_err_in.block_read(in);
      //GridRefinement::refine_and_coarsen_fixed_fraction (hpFEM<dim>::triangulation, est_err_in, 0.18, 0.12);
      GridRefinement::refine_and_coarsen_fixed_fraction (hpFEM<dim>::triangulation, est_err_in, 0.14, 0.12);
      if (hpFEM<dim>::triangulation.n_levels() > max_refine_level){
	for (typename Triangulation<dim>::active_cell_iterator cell = hpFEM<dim>::triangulation.begin_active(max_refine_level); cell != hpFEM<dim>::triangulation.end(); ++cell){
	  if (cell->is_locally_owned()){
	    cell->clear_refine_flag ();
	  }
	}
      }
      hpFEM<dim>::triangulation.prepare_coarsening_and_refinement();
      hpFEM<dim>::triangulation.execute_coarsening_and_refinement();
    }
  }

}

template <int dim>
void CahnHilliard<dim>::apply_initial_condition()
{ 
  this->pcout << "applying initial condition\n";
  int totalDOF=this->totalDOF(this->primary_variables);
  InitCond<dim> init_cond(totalDOF, this->primary_variables, this->primary_variables_dof, *params);
  //init_cond.rand_grid.reinit(8,6,-0.4,0.4,-0.5,0.5);
  //init_cond.rand_grid.reinit(11,6,-0.3,0.3,-0.5,0.5);
  //init_cond.rand_grid.reinit(16,6,-0.4,0.4,-0.5,0.5);
  init_cond.rand_seed.reinit(30,6,0.4,0.04,0.5);
  //init_cond.rand_grid.reinit(21,6,-0.4,0.4,-0.5,0.5);
  //init_cond.rand_grid2.reinit(21,1,-0.03,0.03,-0.5,0.5);
  VectorTools::interpolate(this->dof_handler, init_cond, this->solution_prev); 

  this->solution_prev.compress(VectorOperation::insert);
  this->solution=this->solution_prev;
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
  //GridRefinement::refine_and_coarsen_fixed_fraction (hpFEM<dim>::triangulation, estimated_error_per_cell, 0.10, 0.15);

  // Save estimated_error_per_cell to replicate the mesh refinement when restarting
  if(this->save_snapshot){
    
    Vector<double> est_err_out(estimated_error_per_cell);
    std::string snapshot_path = this->snapshot_directory+"refinement-"+std::to_string(this->current_increment+this->off_output_index)+".dat";
    //FEMdata_refine.create_vector_snapshot(est_err_out, snapshot_path);
    //const char *snappath = snapshot_path.c_str(); 
    std::ofstream out(snapshot_path);
    est_err_out.block_write(out);
  }
  
  GridRefinement::refine_and_coarsen_fixed_fraction (hpFEM<dim>::triangulation, estimated_error_per_cell, 0.14, 0.12);
  //GridRefinement::refine_and_coarsen_fixed_fraction (hpFEM<dim>::triangulation, estimated_error_per_cell, 0.10, 0.12);
  //GridRefinement::refine_and_coarsen_fixed_fraction (hpFEM<dim>::triangulation, estimated_error_per_cell, 0.18, 0.12);
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

template <int dim>
void CahnHilliard<dim>::ini_updateLinearSystem()
{
  //this->pcout << "Clearing grad energy\n";
  for (int j=0; j<6; ++j){
    gradEnergy_eta[j] = 0.;
  }
  gradEnergy_0 = 0.;
  c_avg = 0.;
  area = 0.;
}


/**************************************
 * Implement adaptive time stepping
 *************************************/
template <int dim>
void CahnHilliard<dim>::solve_ibvp()
{

  // Set up Langevin noise for the current time
  /*
    std::random_device rd;
    std::mt19937 gen(rd());
    std::normal_distribution<double> normal_dist(0,1);

    current_noise.resize(0);
    int n_q_points = this->volume_quadrature->size();
    int cell_label = 0;
    typename hp::DoFHandler<dim>::active_cell_iterator cell = hpFEM<dim>::dof_handler.begin_active(), endc=hpFEM<dim>::dof_handler.end();
    for (;cell!=endc; ++cell){
    if (cell->is_locally_owned()){
    std::vector<std::vector<double> > tmp(6,std::vector<double>(n_q_points));
    for (int j=0; j<6; ++j){
    for (int i=0; i<n_q_points; ++i){
    tmp[j][i] = normal_dist(gen);
    }
    }
    current_noise.push_back(tmp);
    cell->set_user_index(cell_label);
    cell_label++;
    }
    }
  */    
  
  int converged;
  while(true){
    //outOfBounds = false;
    converged=this->nonlinearSolve(this->solution);
    //int allOutOfBounds = Utilities::MPI::max(int(outOfBounds),this->mpi_communicator);
    //int allOutOfBounds = Utilities::MPI::max(int(checkBounds()),this->mpi_communicator);
    //bool allOutOfBounds = checkBounds();
    //this->pcout << "Out of bounds: " << allOutOfBounds << std::endl;
    //reset iter_count to 0
    if ((converged > -1)){// && (allOutOfBounds==0)){
      break;
    }
    else{
      iter_count=0;
      this->pcout<<"Not converged or out of bounds, reduce dt."<<std::endl;
      this->current_dt /= std::pow(10.,0.25);
      this->solution=this->solution_prev;
    }
  }
  if (converged < 4){
    iter_count++;
  }
  else{
    iter_count = 0;
  }
  // Check if the current_time is a multiple of the doubled timestep
  bool multiple = (std::fmod(this->current_time,2.*this->current_dt) < 1.e-8);
  if(iter_count>4) {
    this->pcout<<"Increasing timestep \n";
    iter_count=0;
    this->current_dt *= std::pow(10.,0.125);//double dt
  }
  std::vector<double> totalGradEnergy_eta(6,0);
  for (int j=0; j<6; ++j){
    totalGradEnergy_eta[j] = Utilities::MPI::sum(gradEnergy_eta[j],this->mpi_communicator);
  }
  double totalGradEnergy_0 = Utilities::MPI::sum(gradEnergy_0,this->mpi_communicator);
  double total_c_avg = Utilities::MPI::sum(c_avg,this->mpi_communicator);
  double total_area = Utilities::MPI::sum(area,this->mpi_communicator);
  this->pcout << "Grad energy 0: " << totalGradEnergy_0 << "\narea: " << total_area << "\nc_avg: " << total_c_avg/total_area << "\n";
  for (int j=0; j<6; ++j){
    this->pcout << "Grad energy eta " << j << ": " << totalGradEnergy_eta[j] << ", ";
  }
  this->pcout << "\n";
  
  // Adaptive mesh refinement
  refine_grid_once();

  this->solution_prev=this->solution;
}



/**************************************
 * Check if solution is in bounds
 *************************************/
template <int dim>
bool CahnHilliard<dim>::checkBounds()
{

  bool outBounds = false;
  dealii::Vector<double> loc_U(this->solution);
  /*
  //this->pcout << "Size of loc_U: " << loc_U.size() << std::endl;
  for (unsigned int k=0; k<loc_U.size(); k+=8){  //The number of dof is hardcoded right now
    std::vector<double> values = {loc_U[k],loc_U[k+2],loc_U[k+3],loc_U[k+4],loc_U[k+5],loc_U[k+6],loc_U[k+7]};
    // Now, make sure it lies in a physically realistic range
    double check;
    for(unsigned int i=0; i<invQ.size(); ++i){
      check = 0.;
      for(unsigned int j=0; j<invQ[0].size(); ++j){
	//this->pcout << "Value: " << values[j] << std::endl;
	check += invQ[i][j]*values[j];
      }
      if ((check > 1) || (check < 0)){
	outBounds = true;
	break;
      }
    }
    if (outBounds){
      break;
    }
  }
  // */
  //*
  std::vector<types::global_dof_index> local_dof_indices;
  typename hp::DoFHandler<dim>::active_cell_iterator cell = hpFEM<dim>::dof_handler.begin_active(), endc=hpFEM<dim>::dof_handler.end();
  for (;cell!=endc; ++cell){
    if (cell->is_locally_owned()){
      //if (cell->subdomain_id() == this->this_mpi_process){
      const unsigned int dofs_per_cell = cell->get_fe().dofs_per_cell;
      Table<1, double > ULocal(dofs_per_cell);
      local_dof_indices.resize (dofs_per_cell);
      cell->get_dof_indices (local_dof_indices);
      for (unsigned int i=0; i<dofs_per_cell; ++i){
	//std::cout << "Here1a" << loc_U.size() << " " << local_dof_indices[i] << " " << i << "\n";
	ULocal[i]=loc_U(local_dof_indices[i]);
      }
      //this->pcout << "Here2\n";
      for (unsigned int k=0; k<dofs_per_cell; k+=8){
	std::vector<double> values = {ULocal[k],ULocal[k+2],ULocal[k+3],ULocal[k+4],ULocal[k+5],ULocal[k+6],ULocal[k+7]};
	//this->pcout << "Here3\n";
	// Now, make sure it lies in a physically realistic range
	double check;
	for(unsigned int i=0; i<invQ.size(); ++i){
	  check = 0.;
	  for(unsigned int j=0; j<invQ[0].size(); ++j){
	    //this->pcout << "Value: " << values[j] << std::endl;
	    check += invQ[i][j]*values[j];
	  }
	  if ((check > 1) || (check < 0)){
	    outBounds = true;
	    break;
	  }
	}
	//this->pcout << "Here4\n";
	if (outBounds){
	  break;
	}
      }
    }
    if (outBounds){
      break;
    }
  }
  // */
  
  
  return outBounds;

}
