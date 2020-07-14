/*
zhenlin wang 2019
*/

#include"../../include/mechanoChemFEM.h"
template <int dim>
std::vector<double> mechanoChemFEM<dim>::get_solution()
{	
	
	Vector<double> localized_U(this->solution);
	int size=localized_U.size();
	std::vector<double> solution_std(size);
	for (unsigned int i = 0; i < size; ++i)
	  solution_std[i] = localized_U[i];
	return solution_std;
}
template class mechanoChemFEM<1>;
template class mechanoChemFEM<2>;
template class mechanoChemFEM<3>;