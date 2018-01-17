# CylindricalNanowire

Nanowire simulator using cylindrical coordinate

# The simulator is still under contruction !!

This simulator includes 2D Poisson solver(done), 2D Continuity equation(done), 1D Schrodinger Solver(done).

Couple the Poisson solver and Continuity equation = Classical DD model(done)

Couple the Poisson solver and Schrodinger solver = Schrodinger-Poisson Solver

Couple the Poisson solver, Schrodinger solver and Continuity equation = Quantum DD model

This simulator explicitly simulate the oxide (not using Robin type boundary approach).
So it could simulate the wave function penetration in Tox.

Btw, I have no idea how to handle 1/r when r=0 in cylindrical coordinate.
So the the Hamiltonian matrix is not include r=0.

> Reference : A 2-D3-D Schrödinger-Poisson Drift-Diffusion Numerical Simulation of Radially-Symmetric Nanowire MOSFETs

1. Eigen library is required.
http://eigen.tuxfamily.org/index.php?title=Main_Page

  Download Eigen, copy folder "Eigen" to /usr/local/include  

2. openmp and gsl library is resuqired.

  Install gsl library

> sudo apt-get install libgsl0ldbl

  Enable openmp and gsl, Qmake flags

> QMAKE_CXXFLAGS += -fopenmp -lgsl -lgslcblas

> LIBS += -fopenmp -lgsl -lgslcblas


