/*
zhenlin wang 2019
*coupled diffusion reaction
*/
#include "mechanoChemFEM.h"
template <int dim>
class growth: public mechanoChemFEM<dim>
{
	public:
		growth();
		//this is a overloaded function 
		void setMultDomain();
		void apply_boundary_condition();
		void get_residual(const typename hp::DoFHandler<dim>::active_cell_iterator &cell, const FEValues<dim>& fe_values, Table<1, Sacado::Fad::DFad<double> >& R, Table<1, Sacado::Fad::DFad<double>>& ULocal, Table<1, double >& ULocalConv);
		ParameterHandler* params;		
		
		ConstraintMatrix* constraints;
};
template <int dim>
growth<dim>::growth()
{
	//pass the pointer to "constraints" in mechanoChemFEM
	constraints=this->constraints_mechanoChemFEM;
	//This let you use one params to get all parameters pre-defined in the mechanoChemFEM
	params=this->params_mechanoChemFEM;
	params->enter_subsection("parameters");
	params->declare_entry("youngsModulus","0",Patterns::Double() );
	params->declare_entry("poissonRatio","0",Patterns::Double() );
	params->declare_entry("c_ini","0",Patterns::Double() );
	params->declare_entry("M","0",Patterns::Double() );
	params->leave_subsection();	
	
	//Declear the parameters before load it
	this->load_parameters("../parameters.prm");
	
	//define main fields from parameter file.
	this->define_primary_fields();
	//Set up the ibvp.
	this->init_ibvp();
}

template <int dim>
void growth<dim>::get_residual(const typename hp::DoFHandler<dim>::active_cell_iterator &cell, const FEValues<dim>& fe_values, Table<1, Sacado::Fad::DFad<double> >& R, Table<1, Sacado::Fad::DFad<double>>& ULocal, Table<1, double >& ULocalConv)
{
//evaluate primary fields
	params->enter_subsection("parameters");
	double M=params->get_double("M");
	double c_ini=params->get_double("c_ini");
	
	double youngsModulus=params->get_double("youngsModulus");
	double poissonRatio=params->get_double("poissonRatio");
	params->leave_subsection();
	
	unsigned int n_q_points= fe_values.n_quadrature_points;
	int c_dof=0, u_dof=1;
	  
	//mechanics
	deformationMap<Sacado::Fad::DFad<double>, dim> defMap(n_q_points); 
	getDeformationMap<Sacado::Fad::DFad<double>, dim>(fe_values, u_dof, ULocal, defMap);
	
	//chemo
	dealii::Table<1,double>  c_conv(n_q_points);
	dealii::Table<1,Sacado::Fad::DFad<double> > c(n_q_points);
	dealii::Table<2,Sacado::Fad::DFad<double> >  c_grad(n_q_points, dim), j_c(n_q_points, dim);
	
	evaluateScalarFunction<double,dim>(fe_values, c_dof, ULocalConv, c_conv);
	evaluateScalarFunction<Sacado::Fad::DFad<double>,dim>(fe_values, c_dof, ULocal, c);	
	evaluateScalarFunctionGradient<Sacado::Fad::DFad<double>,dim>(fe_values, c_dof, ULocal, c_grad);//at current configuration 

	j_c=table_scaling<Sacado::Fad::DFad<double>, dim>(c_grad,-M);//-D_1*c_1_grad	
	
	dealii::Table<3, Sacado::Fad::DFad<double> > P(n_q_points,dim,dim), Fe(n_q_points,dim,dim);
	
	for(unsigned int q=0; q<n_q_points;q++){
		for (unsigned int i=0; i<dim; ++i){
			for (unsigned int j=0; j<dim; ++j){
	  		Fe[q][i][j]=defMap.F[q][i][j]/std::pow((c[q]/c_ini), 1.0/3.0); //Isotropic growth
			}
		}
	}
	
	//call residual functions
	this->ResidualEq.residualForDiffusionEq(fe_values, c_dof, R, c, c_conv, j_c);
	
  this->ResidualEq.setLameParametersByYoungsModulusPoissonRatio(youngsModulus, poissonRatio);
	this->ResidualEq.evaluateNeoHookeanStress(P, Fe);// NeoHookean model, Saint_Venant_Kirchhoff is also available 
  this->ResidualEq.residualForMechanics(fe_values, u_dof, R, P);	
	
}


	
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

//set Dirichlet BC
template <int dim>
void growth<dim>::apply_boundary_condition()
{
	this->pcout<<"setup_constraints"<<std::endl;
	constraints->clear ();
	DoFTools::make_hanging_node_constraints (this->dof_handler, *constraints);
	
	int totalDOF=this->totalDOF(this->primary_variables);
  std::vector<bool> c_component (totalDOF, false); c_component[0]=true; 
	//apply constraints on boundary
	VectorTools:: interpolate_boundary_values (this->dof_handler, dim, ZeroFunction<dim> (totalDOF),*constraints, c_component);
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
					  	if(ck>0) constraints->add_line (local_face_dof_indices_1[i]);//add constrain line for all u dofs
		  		 	}
					}
				}
			}
		}
	}
	
	constraints->close ();
}

template <int dim>
void InitialConditions<dim>::vector_value (const Point<dim>   &p, Vector<double>   &values) const{
	Assert (values.size() == totalDOF, ExcDimensionMismatch (values.size(), totalDOF));
	if(p[2]==0) values(0)= 1;
  else values(0)= 0.5;
	values(1)=0;
}

template class InitialConditions<1>;
template class InitialConditions<2>;
template class InitialConditions<3>;