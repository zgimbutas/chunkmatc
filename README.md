# High order layer potential evaluation library in R^2 

Author: Zydrunas Gimbutas

Date: August 12, 2025

Version 0.2

### License

```
Copyright (C) 2010-2012: Zydrunas Gimbutas

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met: 

1. Redistributions of source code must retain the above copyright notice, this
   list of conditions and the following disclaimer. 

2. Redistributions in binary form must reproduce the above copyright
   notice, this list of conditions and the following disclaimer in the
   documentation and/or other materials provided with the distribution.

3. Neither the name of the copyright holders nor the names of its
   contributors may be used to endorse or promote products derived
   from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR
ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
```

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

