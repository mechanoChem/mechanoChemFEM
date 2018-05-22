/*
 * functionEvaluations.h
 *
 *  Created on: Sept. 1, 2016, zhenlin wang
 *  (Intend to combine IGA and deal.ii Fe_value, when using own table class)
 * Currently only for deal.ii use
 */

#ifndef FUNCTIONEVALUATIONS_H_
#define FUNCTIONEVALUATIONS_H_
#include <deal.II/fe/fe_values.h>
#include <deal.II/base/table.h>
#include <Sacado.hpp>
#include "supplementaryFunctions.h"
#include "dataStruct.h"

using namespace dealii;

template <class T, int dim>
void evaluateScalarFunction(const FEValues<dim>& fe_values, unsigned int DOF, Table<1, T>& ULocal, Table<1, T>& U);
template <class T, int dim>
void evaluateScalarFunction(const FEValues<dim>& fe_values,const FEFaceValues<dim>& fe_face_values, unsigned int DOF, Table<1, T>& ULocal, Table<1, T>& U);

template <class T, int dim>
void evaluateScalarFunctionGradient(const FEValues<dim>& fe_values, unsigned int DOF, Table<1, T>& ULocal, Table<2, T>& gradU);
template <class T, int dim>
void evaluateScalarFunctionGradient(const FEValues<dim>& fe_values, unsigned int DOF, Table<1, T>& ULocal, Table<2, T>& gradU, deformationMap<T, dim>& defMap);

template <class T, int dim>
void evaluateScalarFunctionGradient(const FEValues<dim>& fe_values,const FEFaceValues<dim>& fe_face_values, unsigned int DOF, Table<1, T>& ULocal, Table<2, T>& gradU);
template <class T, int dim>
void evaluateScalarFunctionGradient(const FEValues<dim>& fe_values,const FEFaceValues<dim>& fe_face_values, unsigned int DOF, Table<1, T>& ULocal, Table<2, T>& gradU, deformationMap<T, dim>& defMap);

template <class T, int dim>
void evaluateVectorFunction(const FEValues<dim>& fe_values, unsigned int DOF, Table<1, T>& ULocal, Table<2, T>& U);

template <class T, int dim>
void evaluateVectorFunctionGradient(const FEValues<dim>& fe_values, unsigned int DOF, Table<1, T>& ULocal, Table<3, T>& gradU);
template <class T, int dim>
void evaluateVectorFunctionGradient(const FEValues<dim>& fe_values, unsigned int DOF, Table<1, T>& ULocal, Table<3, T>& gradU, deformationMap<T, dim>& defMap);


template <class T, int dim>
void evaluateVectorFunctionGradient(const FEValues<dim>& fe_values, const FEFaceValues<dim>& fe_face_values, unsigned int DOF, Table<1, T>& ULocal, Table<3, T>& gradU);
template <class T, int dim>
void evaluateVectorFunctionGradient(const FEValues<dim>& fe_values, const FEFaceValues<dim>& fe_face_values, unsigned int DOF, Table<1, T>& ULocal, Table<3, T>& gradU, deformationMap<T, dim>& defMap);


template <class T, int dim>
void getDeformationMap(const FEValues<dim>& fe_values, unsigned int DOF, Table<1, T>& ULocal, deformationMap<T, dim>& defMap, unsigned int iteration=1);
template <class T, int dim>
void getDeformationMap(const FEValues<dim>& fe_values, const FEFaceValues<dim>& fe_face_values, unsigned int DOF, Table<1, T>& ULocal, deformationMap<T, dim>& defMap,unsigned int iteration=1);

template <class T, int dim>
void evaluateSpatialgradients(Table<2, T>& gradU, Table<2, T> gradU_spat, deformationMap<T, dim>& defMap);

template <class T, int dim>
void evaluateSpatialgradients(Table<3, T>& gradU, Table<3, T> gradU_spat, deformationMap<T, dim>& defMap);

#endif
