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
So the the Hamiltonian matrix does not include r=0.

Schrodinger solver & Poisson solver

> Reference 1 : A 2-D3-D Schrödinger-Poisson Drift-Diffusion Numerical Simulation of Radially-Symmetric Nanowire MOSFETs


Electron Continuity solver

> Reference 2 : Finite difference discretization of semiconductor drift-diffusion equations for nanowire solar cells.pdf

Th electron continuity solver from ref.1 does not work, because they didn't use Scharfetter-Gummel method.

Therefore electron continuity solver is come from ref.2.

Scharfetter-Gummel method

> Reference 3 : ScharfGum.pdf



Eigen library is required.
http://eigen.tuxfamily.org/index.php?title=Main_Page

  Download Eigen, copy folder "Eigen" to /usr/local/include  

openmp and gsl library is resuqired.

  Install gsl library

> sudo apt-get install libgsl0ldbl

  Enable openmp and gsl, Qmake flags

> QMAKE_CXXFLAGS += -fopenmp -lgsl -lgslcblas

> LIBS += -fopenmp -lgsl -lgslcblas


