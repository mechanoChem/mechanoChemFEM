/*
zhenlin wang 2019
*/
#ifndef initial_conds_
#define initial_conds_
#include <cstdlib>
#include <deal.II/numerics/vector_tools.h>
#include <deal.II/base/point.h>
#include <deal.II/base/function.h>
//
using namespace dealii;

//Initial conditions
template <int dim>
class InitialConditions: public Function<dim>{
public:
  InitialConditions (int _totalDOF);
	
	int totalDOF;
  virtual void vector_value (const Point<dim>   &p, Vector<double>   &values) const;
};

#endif
