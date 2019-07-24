High order layer potential evaluation library in R^2 

Fortran
=======

*_dr.f - driver files
*.make - make files

to compile a driver, please type 'make -f filename.make', e.g,

make -f hqsuppquad.make 


---

chunkmatc.f - 2d discretizer + geometry support routines

inter2dn.f - interaction kernel library, Laplace and Helmholtz equations in R^2

hqsuppquad.f - auxiliary quadratures for smooth, log, hilbert, and
               hypersingular kernels (up to 10th order)

cadavect.f - adaptive integration on a segment (smooth functions)


---

test8.f - Dirichlet solver for the Helmholtz equation in R^2, first kind
          integral equation. Put a unit charge inside/outside the circle,
          compare the direct potential with the potential produced by
          the solution charge density. Check for the resonances.

test8dn.f - Dirichlet solver for the Helmholtz equation in R^2, first kind
          integral equation. Put a unit charge inside/outside the circle,
          compare the direct potential with the potential produced by
          the solution charge density. Check for the resonances.
	  Check Dirichlet-Neumann map.

test9.f - Dirichlet solver for the Helmholtz equation in R^2, second kind
          integral equation. Put a unit charge inside/outside the circle,
          compare the direct potential with the potential produced by
          the solution charge density. Check for the resonances.

test9dn.f - Dirichlet solver for the Helmholtz equation in R^2, second kind
          integral equation. Put a unit charge inside/outside the circle,
          compare the direct potential with the potential produced by
          the solution charge density. Check for the resonances.
	  Check Dirichlet-Neumann map (requires hypersingular layer evaluation)


Matlab
======

matlab/

hqsuppquad*.m - auxiliary quadratures for smooth, log, hilbert, and
                hypersingular kernels (1, 2, 3, 4, 5, 6, 8, and 10th order)
                for all Gauss-Legendre support nodes

