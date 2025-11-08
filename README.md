
Summary:

This repository contains a C++ implementation of a Brownian dynamics simulation examining the diffusion of rod-like nanoparticles embedded in polymeric gel networks.

Project Overview:

In many soft‐matter and nanoscience contexts, understanding how anisotropic nanoparticles (rods) move through polymeric gel media is important for applications such as drug delivery, nanocomposite design, filtration, and sensing. This simulation tool models rod‐like particles undergoing Brownian motion while interacting stericly with a surrounding gel environment represented as rigid bodies.
In this code, you’ll find a focus on a rod‐particle geometry, remedies for anisotropic diffusion, and simulation of steric hindrance by a fixed gel network.

Repository Structure:

- RNP_diffusion_AR3_LJ_steric.cpp – The main C++ source file implementing the simulation algorithm.

- ar3.sh – A shell script (e.g., compile + run) wrapper to build and launch the simulation on HPC.

- README.md – (this file) provides usage information, context, and guidance.

Features & Capabilities:

•	Simulation of rod‐like nanoparticles (rigid rods) undergoing Brownian motion in a viscous medium.

•	Steric interaction (via Lennard-Jones potential) modelling between rods and gel network obstacles.

•	Adjustable simulation parameters (rod length, diameter, fiber length, fiber density, time step, total time).

•	Output of trajectories or statistics of diffusion (e.g., mean squared displacement).

•	Shell script for ease of compilation and execution.


Prerequisites:

•	A C++ compiler (e.g., g++, clang++) supporting C++11 or later.

•	A Unix‐style shell environment to run ar3.sh.

•	Optionally, plotting tools (e.g., Python + matplotlib, GNUplot) to analyze output trajectory/statistics.


  Expected Outcomes:
  
•	All important results are saved in txt files. 

•	Parameters.txt saves all variables, constants, and parameters of the system.

•	Nanorods coordianates are saved in these files: save_qx, save_qy, save_qz, save_qx_total, save_qy_total, save_qz_total.






Citation and further details:
- Rokhforouz, Mohammad-Reza, et al. "Brownian dynamics simulation of the diffusion of rod-like nanoparticles in polymeric gels." Soft Matter 21.27 (2025): 5529-5541.








