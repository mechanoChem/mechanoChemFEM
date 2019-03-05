/*
zhenlin wang 2019
*/
#ifndef initial_conds_
#define initial_conds_
#include <cstdlib>
#include <deal.II/numerics/vector_tools.h>
#include <deal.II/base/point.h>
#include <deal.II/base/function.h>
#include <deal.II/base/parameter_handler.h>
//
using namespace dealii;

//Initial conditions
template <int dim>
class InitialConditions: public Function<dim>{
public:
  InitialConditions (int _totalDOF,ParameterHandler& _params);
	
	int totalDOF;
	ParameterHandler* params;
  virtual void vector_value (const Point<dim>   &p, Vector<double>   &values) const;
	virtual double value(const Point<dim>   &p, const unsigned int 	component=0) const;
};

#endif
