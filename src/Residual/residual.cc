#include"../../include/Residual.h"

template <class T, int dim>
Residual<T, dim>::Residual(){};

template <class T, int dim>
Residual<T, dim>::~Residual(){};

template class Residual<Sacado::Fad::DFad<double>, 1>;
template class Residual<Sacado::Fad::DFad<double>, 2>;
template class Residual<Sacado::Fad::DFad<double>, 3>;
template class Residual<double, 1>;
template class Residual<double, 2>;
template class Residual<double, 3>;