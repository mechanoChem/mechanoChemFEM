<B>mechanoChemFEM</B><br>
=======================================================================
mechanoChemFEM is a libarary for modeling of mechano-chemical problems using the finite element method. It consists of [classes](https://htmlpreview.github.io/?https://raw.githubusercontent.com/mechanoChem/mechanoChemFEM/master/doxygen/html/annotated.html) and [functions](https://htmlpreview.github.io/?https://raw.githubusercontent.com/mechanoChem/mechanoChemFEM/master/doxygen/html/modules.html) based on [Deal.ii](https://www.dealii.org/), and examples.

<B>Master branch</B> contains the current stable version of the code and documentaion.<br>
<B>example branch</B> contains all example code.<br>
<B>Other branches</B>  contain the older versions.


<B>List of contributors:</B><br>

Zhenlin Wang (lead developer）<br>

Krishna Garikipati<br>

[Computational Physics Group, University of Michigan](http://umich.edu/~compphys/index.html)


[<B>Code documentation</B>](https://htmlpreview.github.io/?https://raw.githubusercontent.com/mechanoChem/mechanoChemFEM/master/doxygen/html/index.html)

<B>Overview</B><br>
=======================================================================
[deal.ii](http://www.dealii.org) is a robust finite element library.
The dealMultiphysics code is a light library for enhancement of using deal.ii for FEM modeling. It contains five modules:


	hpFEM<T, dim> : offer high level functions dedicated to mulitphysics modeling using deal.ii
	
	Residual< T, dim > : offer residual functions of a variety of equations and boundary conditions

	solveClass< dim, matrixType, vectorType > : offer nonlinear solvers based on deal.ii linear solver

	FEMdata< dim, vectorType > : high level functions for data output and restart(checkpoint)
	
	supplementary:a collection of extra data structures and many helper functions

	
ParameterHandler is used for parameter management. 


<B>Version information</B><br>
=======================================================================
This is version 0.3, the intial release of the code.


<B>License</B><br>
=======================================================================
GNU Lesser General Public License (LGPL). Please see the file LICENSE for details.



<B>Acknowledgements</B><br>
=======================================================================
This code has been developed under the support of the following: <br>

- NSF DMREF grant: DMR1436154 "DMREF: Integrated Computational Framework for Designing Dynamically Controlled Alloy-Oxide Heterostructures" <br>
- DOE BES, Division of Materials Sciences and Engineering: Award #DE-SC0008637 that funds the PRedictive Integrated Structural Materials Science (PRISMS) Center at University of Michigan <br>
- Toyota Research Institute, Award #849910, "Computational framework for data-driven, predictive, multi-scale and multi-physics modeling of battery materials" <br>

<B>Installation with cmake:</B><br>
=======================================================================
1. <B>Install pre-required libs</B><br>


  1) Install CMake [http://www.cmake.org/download/]<br>
  2) Install deal.II [www.dealii.org/download.html] with<br>
		 - Trilinos [https://trilinos.org/]<br>
		 - PetSc [https://www.mcs.anl.gov/petsc/download/index.html]<br>
     Deal.II OSX binaries include full packages of deal.ii with Trillions and other useful libs.


2. <B>Install mechanoChemFEM lib</B><br>


  1) Goes into “build” folder<br>
  2) Modify CMakeList.txt for path of pre-required libs: deal.ii (with Trilinos, Petsc) <br>
  3) cmake CMakeLists.txt <br>
  4) $ make install or do $ make release install<br>
  5) $ make run <br>
     - test(optional) of installation which will run “main” in build folder <br>


<B>Usage</B><br>
=======================================================================
Please see examples.


data:01/20/2018
