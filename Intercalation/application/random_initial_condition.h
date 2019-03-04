/* Iniital Conditions for the coupled transport-mechanics code
   Written by Krishna Garikipati
*/
#ifndef initial_conds_
#define initial_conds_
#include <cstdlib>
//
using namespace dealii;

//Initial conditions
template <int dim>
class InitialConditions: public Function<dim>{
public:
  InitialConditions (int _totalDOF): Function<dim>(_totalDOF), totalDOF(_totalDOF) {}
	
	int totalDOF;
  void vector_value (const Point<dim>   &p, Vector<double>   &values) const{
    Assert (values.size() == totalDOF, ExcDimensionMismatch (values.size(), totalDOF));

    values(0) = -1.0 + static_cast <double> (rand())/(static_cast <double>(RAND_MAX/2.0));//Randomized over [-1.0,1.0]
    values(1) = 0;    

  }
};

#endif