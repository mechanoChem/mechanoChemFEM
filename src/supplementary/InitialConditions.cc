/**
 * @defgroup EvaluationFunctions Evaluation Functions
 */
#include"../../include/supplementary/InitialConditions.h"
template <int dim>
InitialConditions<dim>::InitialConditions (int _totalDOF): Function<dim>(_totalDOF), totalDOF(_totalDOF) {}

template <int dim>
void InitialConditions<dim>::vector_value (const Point<dim>   &p, Vector<double>   &values) const{
  values = 0;    
}
template <int dim>
double InitialConditions<dim>::value(const Point<dim>   &p, const unsigned int 	component) const{
  return 0;   
}
template class InitialConditions<1>;
template class InitialConditions<2>;
template class InitialConditions<3>;