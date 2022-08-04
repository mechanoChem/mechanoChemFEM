# phase field method for topology optimization
A Cahn-Hilliard formulation coupled with elasticity equation for TO.

For parallel run, use the following command:

mpirun -np 4 main -pc_type lu -pc_factor_mat_solver_type superlu_dist
