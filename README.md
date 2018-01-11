# CylindricalNanowire

# This simulator is still under contruction !!

1. Eigen library is required.
http://eigen.tuxfamily.org/index.php?title=Main_Page

  Download Eigen, copy folder "Eigen" to /usr/local/include  

2. openmp and gsl library is resuqired.

  Install gsl library

> sudo apt-get install libgsl0ldbl

  Enable openmp and gsl, Qmake flags

> QMAKE_CXXFLAGS += -fopenmp -lgsl -lgslcblas

> LIBS += -fopenmp -lgsl -lgslcblas


Nanowire simulator using cylindrical coordinate

This simulator includes 2D Poisson solver, 2D Continuity equation, 1D Schrodinger Solver.

Couple the Poisson solver and Continuity equation = Classical DD model

Couple the Poisson solver and Schrodinger solver = Schrodinger-Poisson Solver

Couple the Poisson solver, Schrodinger solver and Continuity equation = Quantum DD model

This simulator using Robin type boundary to simulate the surface potential, we do not explicitly simulate the oxide.
> Reference : A 2-D3-D Schrödinger-Poisson Drift-Diffusion Numerical Simulation of Radially-Symmetric Nanowire MOSFETs
