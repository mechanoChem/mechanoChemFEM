/*
 * supplementaryFunctions.h
 *
 *  Created on: May 11, 2011
 */

#ifndef SUPPLEMENTARYFUNCTIONS_H_
#define SUPPLEMENTARYFUNCTIONS_H_
#include <deal.II/lac/block_vector.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/block_sparse_matrix.h>
#include <deal.II/lac/sparse_direct.h>
#include <deal.II/lac/precondition.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/solver_bicgstab.h>
#include <deal.II/lac/solver_gmres.h>
#include <deal.II/lac/constraint_matrix.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/vector_tools.h>
#include <deal.II/numerics/matrix_tools.h>
#include <iostream>
#include <string>

using namespace dealii;

template <class T, int dim>
inline T determinantOfMinor(unsigned int theRowHeightY, unsigned int theColumnWidthX, Table<2, T>& matrix){
  unsigned int x1 = theColumnWidthX == 0 ? 1 : 0;  /* always either 0 or 1 */
  unsigned int x2 = theColumnWidthX == 2 ? 1 : 2;  /* always either 1 or 2 */
  unsigned int y1 = theRowHeightY   == 0 ? 1 : 0;  /* always either 0 or 1 */
  unsigned int y2 = theRowHeightY   == 2 ? 1 : 2;  /* always either 1 or 2 */
  return matrix[y1][x1]*matrix[y2][x2] - matrix[y1][x2]*matrix[y2][x1];
}


template <class T, int dim>
inline void getInverse(Table<2, T>& matrix, Table<2, T>& invMatrix, T& det){
	if (dim==1){
		det=matrix[0][0];
		invMatrix[0][0]=1.0/det;
	}
	else if(dim==2){
		det=matrix[0][0]*matrix[1][1]-matrix[0][1]*matrix[1][0];
		invMatrix[0][0]=matrix[1][1]/det;
		invMatrix[1][0]=-matrix[1][0]/det;
		invMatrix[0][1]=-matrix[0][1]/det;
		invMatrix[1][1]=matrix[0][0]/det;
	}
	else if(dim==3){
		det=  matrix[0][0]*determinantOfMinor<T, dim>(0, 0, matrix) - matrix[0][1]*determinantOfMinor<T, dim>(0, 1, matrix) +  matrix[0][2]*determinantOfMinor<T, dim>(0, 2, matrix);
		for (int y=0;  y< dim;  y++){
			for (int x=0; x< dim;  x++){
				invMatrix[y][x] = determinantOfMinor<T, dim>(x, y, matrix)/det;
				if( ((x + y) % 2)==1){invMatrix[y][x]*=-1;}
			}
		}
	}
	else throw "dim>3";
	if (std::abs(det)< 1.0e-15){
		printf("**************Near zero determinant in Matrix inversion***********************\n"); throw "Near zero determinant in Matrix inversion";
	}
}

template <class T, int dim>
inline Table<1, T> table_scaling(Table<1, T> matrix, double a){
	unsigned int n_q_points= matrix.size(0);
	Table<1,T> value(n_q_points);
	for(unsigned int q=0; q<n_q_points;q++){
		value[q]=a*matrix[q];	
	}
	return value;
}

template <class T, int dim>
inline Table<2, T> table_scaling(Table<2, T> matrix, double a){
	unsigned int n_q_points= matrix.size(0);
	Table<2, T> value(n_q_points,dim);
	for(unsigned int q=0; q<n_q_points;q++){
		for(unsigned int i=0; i<dim;i++){
			value[q][i]=a*matrix[q][i];
		}	
	}
	return value;
}

template <class T, int dim>
inline Table<1, T> table_scaling(Table<1, T> matrix, T a){
	unsigned int n_q_points= matrix.size(0);
	Table<1,T> value(n_q_points);
	for(unsigned int q=0; q<n_q_points;q++){
		value[q]=a*matrix[q];	
	}
	return value;
}

template <class T, int dim>
inline Table<2, T> table_scaling(Table<2, T> matrix, T a){
	unsigned int n_q_points= matrix.size(0);
	Table<2, T> value(n_q_points,dim);
	for(unsigned int q=0; q<n_q_points;q++){
		for(unsigned int i=0; i<dim;i++){
			value[q][i]=a*matrix[q][i];
		}	
	}
	return value;
}



template <int dim>
void print_mesh_info(const Triangulation<dim> &tria);

template <int dim>
void output_mesh(const Triangulation<dim> &tria, std::string path);

void move_file (const std::string &old_name, const std::string &new_name);

#endif /* SUPPLEMENTARYFUNCTIONS_H_ */

