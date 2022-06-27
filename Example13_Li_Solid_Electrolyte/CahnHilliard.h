/*
zhenlin wang 2019
*CahnHilliard
*/
#include "mechanoChemFEM.h"
#include <math.h>       /* exp */
template <int dim>
class CahnHilliard: public mechanoChemFEM<dim>
{
	public:
		CahnHilliard();
		//this is a overloaded function 
		void get_residual(const typename hp::DoFHandler<dim>::active_cell_iterator &cell, const FEValues<dim>& fe_values, Table<1, Sacado::Fad::DFad<double> >& R, Table<1, Sacado::Fad::DFad<double>>& ULocal, Table<1, double >& ULocalConv);
		ParameterHandler* params;		
		void apply_initial_condition();
		ConstraintMatrix* constraints;

	  int c_dof=0, mu_dof=1, phi_dof=2, zeta_dof=3;
};
template <int dim>
CahnHilliard<dim>::CahnHilliard()
{
	constraints=this->constraints_mechanoChemFEM;
	//This let you use one params to get all parameters pre-defined in the mechanoChemFEM
	params=this->params_mechanoChemFEM;
	params->enter_subsection("Concentration");
	params->declare_entry("c_ini","0",Patterns::Double() );

	params->declare_entry("omega","0",Patterns::Double() );
	params->declare_entry("c_alpha","0",Patterns::Double() );
	params->declare_entry("c_beta","0",Patterns::Double() );
	params->declare_entry("kappa","0",Patterns::Double() );
	params->declare_entry("M","0",Patterns::Double() );
	params->declare_entry("phi_0","0",Patterns::Double() );
	params->declare_entry("zeta_0","0",Patterns::Double() );
	params->declare_entry("flux_0","0",Patterns::Double() );
	params->declare_entry("m_ratio_0","1.0",Patterns::Double() );
	params->leave_subsection();		
	
	//Declear the parameters before load it
	this->load_parameters("parameters.prm");
	
	//define main fields from parameter file.
	this->define_primary_fields();
	//Set up the ibvp.
	this->init_ibvp();
}


template <int dim>
void CahnHilliard<dim>::apply_initial_condition()
{
  std::cout << "applying initial condition (new)\n";

  // 0, 1, 2, 3 for different dof

  typename hp::DoFHandler<dim>::active_cell_iterator cell = this->dof_handler.begin_active(), endc=this->dof_handler.end();
  for (;cell!=endc; ++cell){
    if (cell->subdomain_id() == this->this_mpi_process){
      Point<dim> center=cell->center();
      hp::FEValues<dim> hp_fe_values (this->fe_collection, this->q_collection, update_values | update_quadrature_points);
      hp_fe_values.reinit (cell);
      const FEValues<dim> &fe_values = hp_fe_values.get_present_fe_values();
      const unsigned int dofs_per_cell = cell->get_fe().dofs_per_cell;
      std::vector<unsigned int> local_dof_indices (dofs_per_cell);
      cell->get_dof_indices (local_dof_indices);

      int vertex_id = 0;
      for (unsigned int i=0; i<dofs_per_cell; ++i) {
        int ck = fe_values.get_fe().system_to_component_index(i).first;
        //if (ck==c_dof) this->solution_prev(local_dof_indices[i]) = exp(-10.0 * cell->vertex(vertex_id)[0]) * exp(-10.0 * (1.0-cell->vertex(vertex_id)[1]));
        if (ck==c_dof) this->solution_prev(local_dof_indices[i]) = exp(-5.0 * cell->vertex(vertex_id)[0]) * exp(-5.0 * (1.0-cell->vertex(vertex_id)[1]));
        if (ck==c_dof and this->solution_prev(local_dof_indices[i]) <0.01) this->solution_prev(local_dof_indices[i]) = 0.01;
        if (ck==mu_dof) this->solution_prev(local_dof_indices[i]) = 0.0;
        if (ck==phi_dof) this->solution_prev(local_dof_indices[i]) = exp(-30.0 * cell->vertex(vertex_id)[0]); //phi
        if (ck==zeta_dof) this->solution_prev(local_dof_indices[i]) = exp(-30.0 * (1.0-cell->vertex(vertex_id)[1])); // zeta
        if (ck==zeta_dof) vertex_id += 1;
      }//dofs_per_cell
    } // this_mpi_process
  }//cell

	
	this->solution_prev.compress(VectorOperation::insert);
	this->solution=this->solution_prev;

	constraints->clear ();
	DoFTools::make_hanging_node_constraints (this->dof_handler, *constraints);

  {
    hp::FEValues<dim> hp_fe_values (this->fe_collection, this->q_collection, update_values | update_quadrature_points);
    typename hp::DoFHandler<dim>::active_cell_iterator cell = this->dof_handler.begin_active(), endc = this->dof_handler.end();
    for (; cell != endc; ++cell) {
      if (cell->subdomain_id() == this->this_mpi_process) {
        //int cell_id = cell->active_cell_index();
        Point<dim> center=cell->center();
        hp_fe_values.reinit (cell);
        const FEValues<dim> &fe_values = hp_fe_values.get_present_fe_values();

        const unsigned int dofs_per_cell = cell->get_fe().dofs_per_cell;
        std::vector<unsigned int> local_dof_indices(dofs_per_cell);
        cell->get_dof_indices(local_dof_indices);
        //-------------- fix all displacement
        for (unsigned int i=0; i<dofs_per_cell; ++i) {
          const unsigned int ck = fe_values.get_fe().system_to_component_index(i).first;
          //std::cout << cell_id << " ck " << ck << std::endl;
          if (ck==phi_dof or ck==zeta_dof)
          {
            auto globalDOF = local_dof_indices[i];
            constraints->add_line(globalDOF);
            constraints->set_inhomogeneity(globalDOF, 0.0);
          }

        }

      }
    }
  }

	constraints->close ();
  std::cout << " end of constraint " << std::endl;
}

template <int dim>
void CahnHilliard<dim>::get_residual(const typename hp::DoFHandler<dim>::active_cell_iterator &cell, const FEValues<dim>& fe_values, Table<1, Sacado::Fad::DFad<double> >& R, Table<1, Sacado::Fad::DFad<double>>& ULocal, Table<1, double >& ULocalConv)
{
	//evaluate primary fields
	params->enter_subsection("Concentration");
	double M=params->get_double("M");
	double omega=params->get_double("omega");
	double c_alpha=params->get_double("c_alpha");
	double c_beta=params->get_double("c_beta");
	double kappa=params->get_double("kappa");
	double phi_0=params->get_double("phi_0");
	double zeta_0=params->get_double("zeta_0");
	double flux_0=params->get_double("flux_0");
	double M_ratio=params->get_double("m_ratio_0");
	params->leave_subsection();
	unsigned int n_q_points= fe_values.n_quadrature_points;
		
	dealii::Table<1,double>  c_1_conv(n_q_points);
	dealii::Table<1,Sacado::Fad::DFad<double> > c_1(n_q_points), mu(n_q_points), phi(n_q_points), zeta(n_q_points);
	dealii::Table<2,Sacado::Fad::DFad<double> >  c_1_grad(n_q_points, dim), mu_grad(n_q_points, dim), phi_grad(n_q_points, dim), zeta_grad(n_q_points, dim);
	dealii::Table<2,double>  c_1_grad_conv(n_q_points, dim);
	
	evaluateScalarFunction<double,dim>(fe_values, c_dof, ULocalConv, c_1_conv);
	evaluateScalarFunctionGradient<double,dim>(fe_values, c_dof, ULocalConv, c_1_grad_conv);

	evaluateScalarFunction<Sacado::Fad::DFad<double>,dim>(fe_values, c_dof, ULocal, c_1);	
	evaluateScalarFunctionGradient<Sacado::Fad::DFad<double>,dim>(fe_values, c_dof, ULocal, c_1_grad);
	
	evaluateScalarFunction<Sacado::Fad::DFad<double>,dim>(fe_values, mu_dof, ULocal, mu);	
	evaluateScalarFunctionGradient<Sacado::Fad::DFad<double>,dim>(fe_values, mu_dof, ULocal, mu_grad);

	evaluateScalarFunction<Sacado::Fad::DFad<double>,dim>(fe_values, phi_dof, ULocal, phi);	
	evaluateScalarFunctionGradient<Sacado::Fad::DFad<double>,dim>(fe_values, phi_dof, ULocal, phi_grad);

	evaluateScalarFunction<Sacado::Fad::DFad<double>,dim>(fe_values, zeta_dof, ULocal, zeta);	
	evaluateScalarFunctionGradient<Sacado::Fad::DFad<double>,dim>(fe_values, zeta_dof, ULocal, zeta_grad);
	
	//evaluate diffusion and reaction term
	dealii::Table<1,Sacado::Fad::DFad<double> > rhs_mu(n_q_points);
	dealii::Table<2,Sacado::Fad::DFad<double> > j_c_1(n_q_points, dim), kappa_c_1_grad(n_q_points, dim);
	
	kappa_c_1_grad=table_scaling<dim, Sacado::Fad::DFad<double>, Sacado::Fad::DFad<double> >(c_1_grad,kappa);


  //scalar M
  //j_c_1=table_scaling<dim, Sacado::Fad::DFad<double>, Sacado::Fad::DFad<double> >(mu_grad,-M);//-D_1*c_1_grad
  //for(unsigned int q=0; q<n_q_points;q++) 
  //{
    //for(unsigned int i0=0; i0<dim;i0++) 
    //{
      //j_c_1[q][i0] = mu_grad[q][i0] * (-M) * c_1[q] * (1.0-c_1[q]);
    //}
  //}

  // tensor M based on old direction
  double N_grad_c[dim] = {0.0};
  double T_grad_c[dim] = {0.0};
  double M_tensor[dim][dim] = {0.0};

  //j_c_1=table_scaling<dim, Sacado::Fad::DFad<double>, Sacado::Fad::DFad<double> >(mu_grad,-M);//-D_1*c_1_grad
  for(unsigned int q=0; q<n_q_points;q++) 
  {
    for(unsigned int i0=0; i0<dim;i0++) 
    {
      for(unsigned int j0=0; j0<dim;j0++) 
      {
        M_tensor[i0][j0] = 0.0;
      }
    }
    if (c_1_conv[q]>0.1 and c_1_conv[q]<0.9)
    {
      double len = 0.0;
      for(unsigned int i0=0; i0<dim;i0++) 
      {
        len += c_1_grad_conv[q][i0]*c_1_grad_conv[q][i0];
      }
      len = sqrt(len);
      for(unsigned int i0=0; i0<dim;i0++) 
      {
        N_grad_c[i0] = c_1_grad_conv[q][i0] / len; 
      }
      //std::cout << " c_1_grad_conv: dc_dx=" << c_1_grad_conv[q][0] << "; dc_dy=" << c_1_grad_conv[q][1] << " len=" << len << " x " << N_grad_c[0] << " y " << N_grad_c[1]<< std::endl;
      T_grad_c[0] = N_grad_c[1];
      T_grad_c[1] = -N_grad_c[0];
      for(unsigned int i0=0; i0<dim;i0++) 
      {
        for(unsigned int j0=0; j0<dim;j0++) 
        {
          M_tensor[i0][j0] = M * N_grad_c[i0] * N_grad_c[j0] + M * M_ratio * T_grad_c[i0] * T_grad_c[j0];
        }
      }
      //std::cout << " M_tensor: 00=" << M_tensor[0][0] << "; 01=" << M_tensor[0][1] << "; 10=" << M_tensor[1][0] << "; 11=" << M_tensor[1][1] << std::endl;
    }
    else
    {
      for(unsigned int i0=0; i0<dim;i0++) 
      {
        M_tensor[i0][i0] = M;
      }
    }

    // tensor M
    for(unsigned int i0=0; i0<dim;i0++) 
    {
      j_c_1[q][i0] = 0.0;
      for(unsigned int j0=0; j0<dim;j0++) 
      {
        j_c_1[q][i0] += (-M_tensor[i0][j0]) * mu_grad[q][j0] * c_1[q] * (1.0-c_1[q]);
      }
    }
  }

  //// tensor M based on new direction
  //Sacado::Fad::DFad<double> N_grad_c[dim] = {0.0};
  //Sacado::Fad::DFad<double> T_grad_c[dim] = {0.0};
  //Sacado::Fad::DFad<double> M_tensor[dim][dim] = {0.0};
  //double M_ratio = 1.00001;  // interface diffusion is faster than the bulk diffusion

  ////j_c_1=table_scaling<dim, Sacado::Fad::DFad<double>, Sacado::Fad::DFad<double> >(mu_grad,-M);//-D_1*c_1_grad
  //for(unsigned int q=0; q<n_q_points;q++) 
  //{
    //for(unsigned int i0=0; i0<dim;i0++) 
    //{
      //for(unsigned int j0=0; j0<dim;j0++) 
      //{
        //M_tensor[i0][j0] = 0.0;
      //}
    //}

    //if (c_1[q]>0.1 and c_1[q]<0.9)
    //{
      //Sacado::Fad::DFad<double> len = 0.0;
      //for(unsigned int i0=0; i0<dim;i0++) 
      //{
        //len += c_1_grad[q][i0]*c_1_grad[q][i0];
      //}
      //len = sqrt(len);
      //for(unsigned int i0=0; i0<dim;i0++) 
      //{
        //N_grad_c[i0] = c_1_grad[q][i0] / len; 
      //}
      ////std::cout << " c_1_grad_conv: dc_dx=" << c_1_grad_conv[q][0] << "; dc_dy=" << c_1_grad_conv[q][1] << " len=" << len << " x " << N_grad_c[0] << " y " << N_grad_c[1]<< std::endl;
      //T_grad_c[0] = N_grad_c[1];
      //T_grad_c[1] = -N_grad_c[0];
      //for(unsigned int i0=0; i0<dim;i0++) 
      //{
        //for(unsigned int j0=0; j0<dim;j0++) 
        //{
          //M_tensor[i0][j0] = M * N_grad_c[i0] * N_grad_c[j0] + M * M_ratio * T_grad_c[i0] * T_grad_c[j0];
        //}
      //}
      ////std::cout << " M_tensor: 00=" << M_tensor[0][0] << "; 01=" << M_tensor[0][1] << "; 10=" << M_tensor[1][0] << "; 11=" << M_tensor[1][1] << std::endl;
    //}
    //else
    //{
      //for(unsigned int i0=0; i0<dim;i0++) 
      //{
        //M_tensor[i0][i0] = M;
      //}
    //}
    //// tensor M
    //for(unsigned int i0=0; i0<dim;i0++) 
    //{
      //j_c_1[q][i0] = 0.0;
      //for(unsigned int j0=0; j0<dim;j0++) 
      //{
        //j_c_1[q][i0] += (-M_tensor[i0][j0]) * mu_grad[q][j0] * c_1[q] * (1.0-c_1[q]);
      //}
    //}
  //}
	
	for(unsigned int q=0; q<n_q_points;q++) 
    rhs_mu[q]=
      //(1.0-phi[q])*(1.0-zeta[q])*2.0*omega*(c_1[q]-c_alpha)*(c_1[q]-c_beta)*(2*c_1[q]-c_alpha-c_beta) 
      //omega=A, B=3.5*A
      (1.0-phi[q])*(1.0-zeta[q])*(omega*(log(c_1[q]/(1.0-c_1[q]))) + 2.5 *omega*(1.0-2.0*c_1[q]))
      //+ phi_0 * c_1[q] * (phi_grad[q][0]*phi_grad[q][0]+phi_grad[q][1]*phi_grad[q][1]) 
      //+ zeta_0 * c_1[q]* (zeta_grad[q][0]*zeta_grad[q][0]+zeta_grad[q][1]*zeta_grad[q][1]) 
      + phi_0 * (phi_grad[q][0]*phi_grad[q][0]+phi_grad[q][1]*phi_grad[q][1])  // free energy linearly depends on c
      + zeta_0 * (zeta_grad[q][0]*zeta_grad[q][0]+zeta_grad[q][1]*zeta_grad[q][1]) 
      -mu[q];
	
	//call residual functions
	this->ResidualEq.residualForDiffusionEq(fe_values, c_dof, R, c_1, c_1_conv, j_c_1);
	this->ResidualEq.residualForPoissonEq(fe_values, mu_dof, R, kappa_c_1_grad, rhs_mu);

  //if (cell->active_cell_index() == 0) std::cout << " kappa " << kappa << " M " << M  << " omega "<< omega << std::endl;


	//apply_Neumann_boundary_condition();
  //
  //std::vector<types::global_dof_index> local_face_dof_indices(fe_values.get_fe().dofs_per_face);
  //std::vector<types::global_dof_index> local_dof_indices;
  //local_dof_indices.resize(dofs_per_cell);
  //cell->get_dof_indices(local_dof_indices);
	//BC
	for (unsigned int faceID=0; faceID<2*dim; faceID++){
		if(cell->face(faceID)->at_boundary()==true ){
			//if(cell->face(faceID)->center()[0] <= 0.001 and cell->face(faceID)->center()[1] >= 0.97)
		  if(cell->face(faceID)->center()[0] <= 0.001)
      {
        //std::cout << " c_1[q] " << c_1[0].val() << " " << c_1[1].val() << " " << c_1[2].val() << " " << c_1[3].val() << std::endl;
        if (c_1[0] > 0.35 and c_1[0] < 0.65)
        {
		      FEFaceValues<dim> fe_face_values(fe_values.get_fe(), *(this->common_face_quadrature), update_values | update_quadrature_points | update_JxW_values);
		  	  fe_face_values.reinit(cell,faceID);
		  	  this->ResidualEq.residualForNeummanBC(fe_values, fe_face_values, c_dof , R, flux_0);
        }
		  }
    }
	}

  const unsigned int dofs_per_cell = cell->get_fe().dofs_per_cell;
	for (unsigned int i0=0; i0<dofs_per_cell; ++i0){
    if (R[i0] != R[i0])
    {
      std::cout << " R[i0] " << R[i0] << std::endl;
      exit(0);
    }
  }
	
}

template <int dim>
void InitialConditions<dim>::vector_value (const Point<dim>   &p, Vector<double>   &values) const{
  Assert (values.size() == 2, ExcDimensionMismatch (values.size(), 2));
  values(1) = 0;    
 values(0)= 0.5 + 0.04*(static_cast <double> (rand())/(static_cast <double>(RAND_MAX))-0.5);
}
template class InitialConditions<1>;
template class InitialConditions<2>;
template class InitialConditions<3>;
