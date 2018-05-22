#ifndef hpFEM_h
#define hpFEM_h


//minimum include file
#include <deal.II/lac/constraint_matrix.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/hp/dof_handler.h>
#include <deal.II/grid/tria.h>

#include <deal.II/fe/fe_dgq.h>
#include <deal.II/fe/fe_nothing.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/fe_values.h>

#include <deal.II/hp/dof_handler.h>
#include <deal.II/hp/q_collection.h>
#include <deal.II/hp/fe_values.h>
#include <deal.II/hp/fe_collection.h>

using namespace dealii;
template<int dim>
/**
*This class is intended to be base/abstract class. It contains three essential
*objects in deal.ii:  \c Triangulation<dim>, \c hp::DoFHandler<dim> and \c ConstraintMatrix. Those three will be shared in other classes: \code{.cpp}solveClass< dim, matrixType, vectorType > and FEMdata< dim, vectorType >\endcode
*Class hpFEM< dim > also provide functions to
*/
class hpFEM
{
	
  public:
		/**
		* abstract class, do nothing specifically except initializing dof_handler with triangulation.
		*/
    hpFEM();
    ~hpFEM();
		void setup_FeSystem(std::vector<std::shared_ptr<FESystem<dim>> >& fe_system, hp::FECollection<dim>& fe_collection, hp::QCollection<dim>& q_collection, std::vector<unsigned int >& primary_variables_dof, 
									 std::vector<std::vector<std::string> >& primary_variables, std::vector<std::vector<int> >& FE_support, const QGauss<dim>& volume_quadrature);
		void set_active_fe_indices(std::vector<std::vector<int> >& FE_support, hp::DoFHandler<dim>& local_handler, int domain=0);
		int totalDOF(std::vector<std::vector<std::string> >& primary_variables);
		Triangulation<dim>    triangulation;
		hp::DoFHandler<dim>   dof_handler;
		ConstraintMatrix constraints;

};

#endif
