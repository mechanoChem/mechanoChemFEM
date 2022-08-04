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
		void output_results();
		ConstraintMatrix* constraints;

	  int c_dof=0, mu_dof=1, u_dof=2;
    double total_phi_e = 0.0;
    double total_phi_i = 0.0;
    double total_phi_b = 0.0;

	  Vector<double> cell_phi_e; 
	  Vector<double> cell_phi_i; 
	  Vector<double> cell_phi_b; 

	  Vector<double> cell_sig11; 
	  Vector<double> cell_sig22; 
	  Vector<double> cell_sig12; 
};
template <int dim>
CahnHilliard<dim>::CahnHilliard()
{
	constraints=this->constraints_mechanoChemFEM;
	//This let you use one params to get all parameters pre-defined in the mechanoChemFEM
	params=this->params_mechanoChemFEM;
	params->enter_subsection("Parameters");
	params->declare_entry("L_0","0",Patterns::Double() );
	params->declare_entry("T_0","0",Patterns::Double() );
	params->declare_entry("E_0","0",Patterns::Double() );
	params->declare_entry("nu_0","0",Patterns::Double() );
	params->declare_entry("h_0","0",Patterns::Double() );
	params->declare_entry("f_0","0",Patterns::Double() );
	params->declare_entry("M_0","0",Patterns::Double() );
	params->declare_entry("d_1","0",Patterns::Double() );
	params->declare_entry("d_2","0",Patterns::Double() );
  params->declare_entry("h_bar_y","0",Patterns::Double() );

	params->declare_entry("lambda_bar","0",Patterns::Double() );
	params->declare_entry("gamma_bar","0",Patterns::Double() );
	params->declare_entry("gamma_E","0",Patterns::Double() );
	params->declare_entry("D_1","0",Patterns::Double() );
	params->declare_entry("D_2","0",Patterns::Double() );
	params->declare_entry("D_3","0",Patterns::Double() );
	params->declare_entry("D_4","0",Patterns::Double() );
	params->declare_entry("D_5","0",Patterns::Double() );

  params->declare_entry("enable_coupling","0",Patterns::Integer());

	params->leave_subsection();		
	
	//Declear the parameters before load it
	this->load_parameters("parameters.prm");
	
	//define main fields from parameter file.
	this->define_primary_fields();
	//Set up the ibvp.
	this->init_ibvp();

  cell_phi_e.reinit(this->triangulation.n_active_cells());
  cell_phi_i.reinit(this->triangulation.n_active_cells());
  cell_phi_b.reinit(this->triangulation.n_active_cells());
  cell_sig11.reinit(this->triangulation.n_active_cells());
  cell_sig22.reinit(this->triangulation.n_active_cells());
  cell_sig12.reinit(this->triangulation.n_active_cells());
}


template <int dim>
void CahnHilliard<dim>::apply_initial_condition()
{
  std::cout << "applying initial condition (new)\n";

  // 0, 1, 2, 3 for different dof

  typename hp::DoFHandler<dim>::active_cell_iterator cell = this->dof_handler.begin_active(), endc=this->dof_handler.end();
  for (;cell!=endc; ++cell){
    if (cell->subdomain_id() == this->this_mpi_process)
    {
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
        if (ck==c_dof) this->solution_prev(local_dof_indices[i]) = 0.35 + 0.0001*(static_cast <double> (rand())/(static_cast <double>(RAND_MAX))-0.5);
        //if (ck==c_dof) this->solution_prev(local_dof_indices[i]) = exp(-5.0 * cell->vertex(vertex_id)[0]) * exp(-5.0 * (1.0-cell->vertex(vertex_id)[1]));
        //if (ck==c_dof and this->solution_prev(local_dof_indices[i]) <0.01) this->solution_prev(local_dof_indices[i]) = 0.01;
        if (ck==mu_dof) this->solution_prev(local_dof_indices[i]) = 0.0;
        if (ck==u_dof or ck==u_dof+1 or ck==u_dof+2) this->solution_prev(local_dof_indices[i]) = 0.0;
      }//dofs_per_cell
    } // this_mpi_process
  }//cell

	
	this->solution_prev.compress(VectorOperation::insert);
	this->solution=this->solution_prev;

	constraints->clear ();
	DoFTools::make_hanging_node_constraints (this->dof_handler, *constraints);

  // TODO : add proper elasticity constraint here
  {
	  params->enter_subsection("Geometry");
	  double x_min=params->get_double("x_min");
	  double x_max=params->get_double("x_max");
	  params->leave_subsection();

    hp::FEValues<dim> hp_fe_values (this->fe_collection, this->q_collection, update_values | update_quadrature_points);
    typename hp::DoFHandler<dim>::active_cell_iterator cell = this->dof_handler.begin_active(), endc = this->dof_handler.end();
    for (; cell != endc; ++cell) {
      if (cell->subdomain_id() == this->this_mpi_process) 
      {
        //int cell_id = cell->active_cell_index();
        Point<dim> center=cell->center();
        hp_fe_values.reinit (cell);
        const FEValues<dim> &fe_values = hp_fe_values.get_present_fe_values();

        const unsigned int dofs_per_cell = cell->get_fe().dofs_per_cell;
        std::vector<unsigned int> local_dof_indices(dofs_per_cell);
        cell->get_dof_indices(local_dof_indices);

        int _vertex_id = -1;
        for (unsigned int i=0; i<dofs_per_cell; ++i) {
          const unsigned int ck = fe_values.get_fe().system_to_component_index(i).first;
          if (ck == 0) _vertex_id += 1;
          //std::cout << cell_id << " ck " << ck << std::endl;
          if ( (ck==u_dof or ck==u_dof+1) and (std::abs(cell->vertex(_vertex_id)[0] - x_min) < 1e-3))
          //if (ck==u_dof or ck==u_dof+1)
          {
            //std::cout << " i (fix x y) " << i << " center " << cell->vertex(_vertex_id) << std::endl;
            auto globalDOF = local_dof_indices[i];
            constraints->add_line(globalDOF);
            constraints->set_inhomogeneity(globalDOF, 0.0);
          }

          // if (ck==c_dof or ck==mu_dof) 
          // {// fix mu and c
          //   //std::cout << " i (fix x y) " << i << " center " << cell->vertex(_vertex_id) << std::endl;
          //   auto globalDOF = local_dof_indices[i];
          //   constraints->add_line(globalDOF);
          //   constraints->set_inhomogeneity(globalDOF, 0.0);
          // }
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
  int cell_id = cell->active_cell_index();
  this->cell_phi_e[cell_id] = 0.0;
  this->cell_phi_i[cell_id] = 0.0;
  this->cell_phi_b[cell_id] = 0.0;

  //if (cell_id == 0) 
  //{
    //this->total_phi_e = 0.0;
    //this->total_phi_i = 0.0;
    //this->total_phi_b = 0.0;
  //}
	//evaluate primary fields
	params->enter_subsection("Parameters");

	double L_0=params->get_double("L_0");
	double T_0=params->get_double("T_0");
	double E_0=params->get_double("E_0");
	double nu_0=params->get_double("nu_0");
	double h_0=params->get_double("h_0");
	double f_0=params->get_double("f_0");
	double M_0=params->get_double("M_0");
	double d_1=params->get_double("d_1");
	double d_2=params->get_double("d_2");

  double h_bar_y=params->get_double("h_bar_y");


	double lambda_bar=params->get_double("lambda_bar");
	double gamma_bar=params->get_double("gamma_bar");
	double gamma_E=params->get_double("gamma_E");

	double D_1=params->get_double("D_1");
	double D_2=params->get_double("D_2");
	double D_3=params->get_double("D_3");
	double D_4=params->get_double("D_4");
	double D_5=params->get_double("D_5");

	int enable_coupling =params->get_double("enable_coupling");

	params->leave_subsection();


	params->enter_subsection("Geometry");
	double x_min=params->get_double("x_min");
	double x_max=params->get_double("x_max");
	double y_min=params->get_double("y_min");
	double y_max=params->get_double("y_max");
	params->leave_subsection();

	unsigned int n_q_points= fe_values.n_quadrature_points;
  const unsigned int dofs_per_cell = cell->get_fe().dofs_per_cell;

	////mechanics
	//deformationMap<Sacado::Fad::DFad<double>, dim> defMap(n_q_points); 
	//getDeformationMap<Sacado::Fad::DFad<double>, dim>(fe_values, u_dof, ULocal, defMap);
		
	dealii::Table<1,Sacado::Fad::DFad<double> > c_1(n_q_points), mu(n_q_points);
	dealii::Table<2,Sacado::Fad::DFad<double> > c_1_grad(n_q_points, dim), mu_grad(n_q_points, dim);
	dealii::Table<1,double>  c_1_conv(n_q_points);
	dealii::Table<2,double>  c_1_grad_conv(n_q_points, dim);
	
	evaluateScalarFunction<double,dim>(fe_values, c_dof, ULocalConv, c_1_conv);
	evaluateScalarFunctionGradient<double,dim>(fe_values, c_dof, ULocalConv, c_1_grad_conv);

	evaluateScalarFunction<Sacado::Fad::DFad<double>,dim>(fe_values, c_dof, ULocal, c_1);	
	evaluateScalarFunctionGradient<Sacado::Fad::DFad<double>,dim>(fe_values, c_dof, ULocal, c_1_grad);
	
	evaluateScalarFunction<Sacado::Fad::DFad<double>,dim>(fe_values, mu_dof, ULocal, mu);	
	evaluateScalarFunctionGradient<Sacado::Fad::DFad<double>,dim>(fe_values, mu_dof, ULocal, mu_grad);
	
	//evaluate diffusion and reaction term
	dealii::Table<1,Sacado::Fad::DFad<double> > rhs_mu(n_q_points);
	dealii::Table<2,Sacado::Fad::DFad<double> > j_c_1(n_q_points, dim), kappa_c_1_grad(n_q_points, dim);

  //scalar M
  j_c_1=table_scaling<dim, Sacado::Fad::DFad<double>, Sacado::Fad::DFad<double> >(mu_grad,-M_0); 

	//call residual functions
	this->ResidualEq.residualForDiffusionEq(fe_values, c_dof, R, c_1, c_1_conv, j_c_1); // D_1=1 is implicitly used here.
  //std::cout << "done with rho" << std::endl;


  // mechanics
  Sacado::Fad::DFad<double>  phi_e[n_q_points];
  Sacado::Fad::DFad<double>  phi_b[n_q_points];
  Sacado::Fad::DFad<double>  phi_i[n_q_points];

  for (unsigned int q = 0; q < n_q_points; ++q) {
    phi_e[q] = 0.0;
    phi_b[q] = 0.0;
    phi_i[q] = 0.0;

    Sacado::Fad::DFad<double> rho = c_1[q];
    //Sacado::Fad::DFad<double> g_rho_fcn = std::pow(1.0/(1.0 + std::exp(-10.0*(rho - 0.5))), 3.0);
    Sacado::Fad::DFad<double> g_rho_fcn = 1.0/(1.0 + std::exp(-10.0*(rho - 0.5))); // without power 3
    if (not enable_coupling) 
    {
      g_rho_fcn = 1.0;
      D_3 = 0.0;
    }

    //std::cout << "done with g_rho " << g_rho_fcn<< std::endl;

    Sacado::Fad::DFad<double> grad_u[dim][dim] = {0.0};
    Sacado::Fad::DFad<double> strain[dim][dim] = {0.0};
    Sacado::Fad::DFad<double> stress[dim][dim] = {0.0};

    for (unsigned int c = 0; c < dim; c++) {
      for (unsigned int d = 0; d < dim; d++) {
        strain[c][d] = 0.0;
      }
    }
  //std::cout << "done with strain 1" << std::endl;
    for (unsigned int i = 0; i < dofs_per_cell; ++i) {
      //std::cout << " ULocal[i] " << i << " " << ULocal[i] << std::endl;
      const unsigned int c = fe_values.get_fe().system_to_component_index(i).first;

      if ( c >= u_dof and c <u_dof+dim){
        for (unsigned int d = 0; d < dim; d++) {
          grad_u[c-u_dof][d] += ULocal[i] * fe_values.shape_grad_component(i, q, c)[d];  // B(3x3x8)@q=B(c,d,i)@q
          //std::cout << " grad_u " << grad_u[0][0] << " "<< fe_values.shape_grad_component(i, q, c)[d] << std::endl;
        }
      }
    }

    for (unsigned int c = 0; c < dim; c++)
      for (unsigned int d = 0; d < dim; d++) strain[c][d] = 0.5 * (grad_u[c][d] + grad_u[d][c]);  
    //std::cout << "done with strain " << strain[0][0] << std::endl;

    // E_0 is replaced with 1.0
    //double mu = E_0/2.0/(1.0+nu_0);
    //double lambda = (E_0*nu_0)/(1.0+nu_0)/(1.0-2.0*nu_0);
    double mu = 1.0/2.0/(1.0+nu_0);
    double lambda = (1.0*nu_0)/(1.0+nu_0)/(1.0-2.0*nu_0);

    // compute stress
    Sacado::Fad::DFad<double> tr_e = 0.0;
    for (unsigned int c = 0; c < dim; c++) tr_e += strain[c][c];

    for (unsigned int c = 0; c < dim; c++) {
      for (unsigned int d = 0; d < dim; d++) {
        stress[c][d] = g_rho_fcn * 2.0 * mu * strain[c][d]; // check if the g_rho_fcn is correct.
        if (c == d) stress[c][d] += g_rho_fcn * lambda * tr_e;
      }
    }
    this->cell_sig11[cell_id] = stress[0][0].val();
    this->cell_sig22[cell_id] = stress[1][1].val();
    this->cell_sig12[cell_id] = stress[0][1].val();
  //std::cout << "done with stress " << stress[0][0] << std::endl;

    for (unsigned int c = 0; c < dim; c++) {
      for (unsigned int d = 0; d < dim; d++) {
        phi_e[q] += stress[c][d] * strain[c][d];
      }
    }

  //std::cout << "done with phi_e" << std::endl;
    for (unsigned i = 0; i < dofs_per_cell; ++i) {
      const unsigned int c = fe_values.get_fe().system_to_component_index(i).first;
      if ( c >= u_dof and c < u_dof + dim){
        for (unsigned int d = 0; d < dim; d++)
        {
          R[i] -= D_4 * stress[c-u_dof][d] * fe_values.shape_grad_component(i, q, c)[d] * fe_values.JxW(q);
          //std::cout << "done with R[i] " << i << " " << R[i] << std::endl;
        }
      }
    }

  } // q for mechanics
  //std::cout << "done with mechanics" << std::endl;

	//kappa_c_1_grad=table_scaling<dim, Sacado::Fad::DFad<double>, Sacado::Fad::DFad<double> >(c_1_grad,kappa);
	kappa_c_1_grad=table_scaling<dim, Sacado::Fad::DFad<double>, Sacado::Fad::DFad<double> >(c_1_grad,L_0*L_0); // dimensionless parameters
  for(unsigned int q=0; q<n_q_points;q++) 
  {
    rhs_mu[q]= 
      - D_3 * ( 1484.13*std::exp(10*c_1[q])/(148.413+std::exp(10*c_1[q]))/(148.413+std::exp(10*c_1[q]))*phi_e[q] ) // elastic part
      + D_2 * (4*c_1[q]*c_1[q]*c_1[q] - 6*c_1[q]*c_1[q] + 2*c_1[q] - 57.5646*std::pow(10, -50.0*c_1[q]) + 5.75646*std::pow(10,-49)*std::pow(10,50.0*c_1[q]))
      -mu[q];

    //if (rhs_mu[q] != rhs_mu[q]) std::cout << " rhs_mu[q] " << q << " " << rhs_mu[q] << " c_1 " << c_1[q] << " phi_e[q] " << phi_e[q] << std::endl;
  }
	this->ResidualEq.residualForPoissonEq(fe_values, mu_dof, R, kappa_c_1_grad, rhs_mu);

  for(unsigned int q=0; q<n_q_points;q++) 
  {
    phi_b[q] = c_1[q]*c_1[q]*(1.0-c_1[q])*(1.0-c_1[q]) + 0.5*(std::pow(10, -50.0*c_1[q]) + std::pow(10, 50*(c_1[q]-1)));
    phi_i[q] = lambda_bar * 0.00625 * 0.00625 * 0.5 * (c_1_grad[q][0]*c_1_grad[q][0] + c_1_grad[q][1]*c_1_grad[q][1]);
    this->cell_phi_e[cell_id] += gamma_bar * gamma_E * phi_e[q].val() * fe_values.JxW(q);
    this->cell_phi_b[cell_id] += phi_b[q].val() * fe_values.JxW(q);
    this->cell_phi_i[cell_id] += phi_i[q].val() * fe_values.JxW(q);
  }

  //std::cout << "done with mu" << std::endl;
  // discretize the mechanical loading to "# of load_steps".
  int load_steps = 1;

	//apply Neumann boundary condition
  for (unsigned int faceID=0; faceID<2*dim; faceID++){
    if(cell->face(faceID)->at_boundary()==true ){
      if( std::abs(cell->face(faceID)->center()[0] - x_max) <= 1e-4 )
      {
        if( cell->face(faceID)->center()[1] >= (y_max - d_1 - d_2) and cell->face(faceID)->center()[1] <= (y_max - d_1) )
        {
          FEFaceValues<dim> fe_face_values(fe_values.get_fe(), *(this->common_face_quadrature), update_values | update_quadrature_points | update_JxW_values);
          fe_face_values.reinit(cell,faceID);

          unsigned int dofs_per_cell= fe_values.dofs_per_cell;
          unsigned int n_face_q_points=fe_face_values.n_quadrature_points;
          unsigned int ck;
          for (unsigned int i=0; i<dofs_per_cell; ++i) {
            ck = fe_values.get_fe().system_to_component_index(i).first;
            if(ck>=u_dof && ck<u_dof + dim){
              for (unsigned int q=0; q<n_face_q_points; ++q){
                if (ck == u_dof + 1)
                {
                  // check if the y traction is done correctly or not
                  R[i] += fe_face_values.shape_value(i, q)* h_bar_y * 1.0 / load_steps * std::min(this->current_increment, load_steps) * fe_face_values.JxW(q); // traction in the y direction.
                } // y direction traction
              } // loop of q
            } // ck condition
          }	 // loop of dofs_per_cell
        } // at y d_1
      } // at x_max
    } // at boundary
  }

  //if (cell_id == 0) std::cout << " load factor " <<  std::min(this->current_increment, load_steps)*1.0/load_steps << std::endl;

  //std::cout << "done with BCs" << std::endl;

	for (unsigned int i0=0; i0<dofs_per_cell; ++i0){
    int cell_id = cell->active_cell_index();
    Point<dim> center=cell->center();
    if (cell_id == -1)     std::cout << " R[i0] " << i0 << " "  << R[i0] << std::endl;
    //std::cout << i0 << " R[i0] " << R[i0] << " cell_id " << cell_id << " center " << center  << std::endl;
    if (R[i0] != R[i0])
    {
      std::cout << " R[i0] " << R[i0] << " cell_id " << cell_id << " center " << center  << std::endl;
      exit(0);
    }
  }

  ////std::cout << cell_id << std::endl;
  //if (cell_id == 12799) 
    //std::cout 
      //<< " total_phi_e " << this->total_phi_e 
      //<< " total_phi_b " << this->total_phi_b
      //<< " total_phi_i " << this->total_phi_i 
      //<< std::endl;
	
}

template <int dim>
void CahnHilliard<dim>::output_results()
{
	Vector<float> subdomain_id(this->triangulation.n_active_cells()); 
  typename hp::DoFHandler<dim>::active_cell_iterator elem = this->dof_handler.begin_active(), endc = this->dof_handler.end();             
  unsigned int _j = 0;                                                                                                      
  for (;elem!=endc; ++elem){                                                                                                
		    subdomain_id(_j++) = elem->subdomain_id();
	}
  //std::cout << " after elem id " << std::endl;
	Vector<double> _cell_phi_e(this->triangulation.n_active_cells()); 
	Vector<double> _cell_phi_i(this->triangulation.n_active_cells()); 
	Vector<double> _cell_phi_b(this->triangulation.n_active_cells()); 

	Vector<double> _cell_sig11(this->triangulation.n_active_cells()); 
	Vector<double> _cell_sig22(this->triangulation.n_active_cells()); 
	Vector<double> _cell_sig12(this->triangulation.n_active_cells()); 

  Utilities::MPI::sum(cell_phi_e, MPI_COMM_WORLD, _cell_phi_e);
  Utilities::MPI::sum(cell_phi_i, MPI_COMM_WORLD, _cell_phi_i);
  Utilities::MPI::sum(cell_phi_b, MPI_COMM_WORLD, _cell_phi_b);
  Utilities::MPI::sum(cell_sig11, MPI_COMM_WORLD, _cell_sig11);
  Utilities::MPI::sum(cell_sig22, MPI_COMM_WORLD, _cell_sig22);
  Utilities::MPI::sum(cell_sig12, MPI_COMM_WORLD, _cell_sig12);

  total_phi_e = _cell_phi_e.l1_norm();
  total_phi_i = _cell_phi_i.l1_norm();
  total_phi_b = _cell_phi_b.l1_norm();

	//write vtk and snapshot for solution
	if(this->save_output){ 
		std::string output_path = this->output_directory+"output-"+std::to_string(this->current_increment)+".vtk";
		this->FEMdata_out.clear_data_vectors();
    this->FEMdata_out.data_out.add_data_vector(subdomain_id, "sub_id");

	  this->FEMdata_out.data_out.add_data_vector(_cell_phi_e, "phi_e");
	  this->FEMdata_out.data_out.add_data_vector(_cell_phi_i, "phi_i");
	  this->FEMdata_out.data_out.add_data_vector(_cell_phi_b, "phi_b");
	  this->FEMdata_out.data_out.add_data_vector(_cell_sig11, "sig11");
	  this->FEMdata_out.data_out.add_data_vector(_cell_sig22, "sig22");
	  this->FEMdata_out.data_out.add_data_vector(_cell_sig12, "sig12");

		if(this->current_increment%this->skip_output==0) this->FEMdata_out.write_vtk(this->solution_prev, output_path);	

	}
	if(this->save_snapshot){
		std::string snapshot_path = this->snapshot_directory+"snapshot-"+std::to_string(this->current_increment+this->off_output_index)+".dat";
    this->pcout << " save to " << snapshot_path << std::endl;
		this->FEMdata_out.create_vector_snapshot(this->solution, snapshot_path);
	}

  this->pcout 
    << " total_phi_e " << this->total_phi_e 
    << " total_phi_b " << this->total_phi_b
    << " total_phi_i " << this->total_phi_i 
    << std::endl;
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
