#ifndef dataStruct_H_
#define dataStruct_H_
#include <deal.II/base/table.h>
#include <Sacado.hpp>
template <class T, int dim>
  struct deformationMap{
  deformationMap(unsigned int _n_q_points): n_q_points(_n_q_points), F(_n_q_points, dim, dim),  invF(_n_q_points, dim, dim), detF(_n_q_points){}
    dealii::Table<3, T> F, invF;
    dealii::Table<1, T> detF;
		unsigned int n_q_points;
 };

template <class T, int dim>
  struct deformationMapwithGrad{
  deformationMapwithGrad(unsigned int n_q_points): F(n_q_points, dim, dim),  invF(n_q_points, dim, dim), gradF(n_q_points, dim, dim, dim), detF(n_q_points){}
    dealii::Table<3, T> F, invF;
    dealii::Table<4, T> gradF;
    dealii::Table<1, T> detF;
 };
	

#endif	