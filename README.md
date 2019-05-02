# 3dVPM
3D Unsteady Vortex Panel Method, forked from pranavladkat/3dVPM

This code has been mainly modified with the intent to simulate the flow over multibladed wind-turbines. The additional features with respect to the original code are:

* Several surfaces can now interact through the extensive use of vectors

* Deformation of the blades can be imposed, such as an oscillating motion around the beam, a flapping motion, etc.

* The original code implemented a simple Morino condition at the trailing edge. As a consequence the pressure at the trailing edge had different values at the upper and lower sections. 
An iterative Runge-Kutta condition has been implemented (based on [1]) and displays a great deal of improvement with respect to the experimental data.

Requirements

* petsc 3.10: the system of equations to solve the doublets on the surfaces is based on the library Petsc.

* Intel MKL BLAS (not mandatory) if you don't have it, you can use the alternative code put in comment, below the calling to the MKL functions.

[1] Wang, Youjiang and Abdel-Maksoud, Mousatafa and Song, Baowei, "A fast method to realize the pressure Kutta condition in boundary element method for lifting bodies", Ocean Engineering 130 (2017) 398-406