To run the example, use the following command line: 

mpirun -np 4 ./main -pc_type lu -pc_factor_mat_solver_package superlu_dist

One can change the total number of CPUs(=4, currently) to a different number.
