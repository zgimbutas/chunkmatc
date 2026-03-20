cc SPDX-License-Identifier: BSD-3-Clause
cc Copyright (C) 2009-2026: Leslie Greengard and Zydrunas Gimbutas
cc See LICENSE file for full license text
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c       Complex GMRES routines (Arnoldi + Givens rotations)
c
c       Based on the netlib GMRES template:
c       "Templates for the Solution of Linear Systems: Building Blocks
c       for Iterative Methods", Barrett et al., SIAM, 1993.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
c
c
c
        subroutine cgmres
     $     (ier,n,a,multa,p1,p2,p3,p4,p5,p6,y,eps,numit,
     1     x,niter,errs,ngmrec,w)
        implicit real *8 (a-h,o-z)
        dimension errs(1)
        complex *16 a(n,n),x(1),y(1),w(1)
c
        external multa
c
c     This subroutine solves a complex linear system Ax=y by means
c     of restarted GMRES(m) algorithm with Arnoldi orthogonalization
c     and Givens rotations.
c
c     Input parameters:
c
c     n - the dimensionality of the linear system
c     y - the right hand side
c     a - the matrix of the system (or whatever other parameter)
c     eps - the required relative accuracy
c     multa - the user-defined matrix-vector multiplication subroutine
c
c     the calling sequence for multa must be
c
c     multa(a,p1,p2,p3,p4,p5,p6,x,y,n)               (1)
c
c       in (1), a is a matrix the system or whatever other parameter,
c       par1, par2 are whatever other parameters, x is an input vector,
c       y is a product Ax, and n is the dimensionality of a, x, and y.
c
c     numit - the maximum number of outer iterations (restarts) permitted
c     ngmrec - the number of inner iterations per restart cycle
c
c     w - must be at least n*(ngmrec+1) + (ngmrec+1)*(ngmrec+1)
c         + 3*(ngmrec+1) + n complex *16 elements long
c
c     Output parameters:
c
c     ier - error return code
c        ier=0 normal execution of the subroutine
c        ier=4 means that the maximum number of iterations
c           has been reached without achieving the required accuracy eps
c
c     x - the solution of the system
c     niter - the number of iterations performed
c     errs - the array of errors (relative) produced by the algorithm.
c
c
c       ... memory management
c
c       V: n x (ngmrec+1) Arnoldi basis
c
        iv=1
        liv=n*(ngmrec+1)
c
c       H: (ngmrec+1) x ngmrec upper Hessenberg matrix
c
        ih=iv+liv
        lih=(ngmrec+1)*ngmrec
c
c       cs: ngmrec Givens cosines
c
        ics=ih+lih
        lics=ngmrec
c
c       sn: ngmrec Givens sines
c
        isn=ics+lics
        lisn=ngmrec
c
c       s: (ngmrec+1) residual vector
c
        is=isn+lisn
        lis=ngmrec+1
c
c       r: n residual work vector
c
        ir=is+lis
        lir=n
c
        call cgmres1
     $     (ier,n,a,multa,p1,p2,p3,p4,p5,p6,y,eps,numit,x,niter,errs,
     1     ngmrec,w(iv),w(ih),w(ics),w(isn),w(is),w(ir))
c
        return
        end
c
c
c
c
c
        subroutine cgmres1
     $     (ier,n,a,multa,p1,p2,p3,p4,p5,p6,y,eps,numit,
     1     x,niter,errs,m,v,h,cs,sn,s,r)
        implicit real *8 (a-h,o-z)
c
        dimension errs(1)
        complex *16 a(n,n),x(n),y(n)
        complex *16 v(n,m+1),h(m+1,m),cs(m),sn(m),s(m+1),r(n)
        complex *16 d,temp
c
        external multa
c
        ier=0
        niter=0
c
c       ... compute norm of rhs
c
        bnrm2=0
        do 1000 i=1,n
        bnrm2=bnrm2+dble(dconjg(y(i))*y(i))
 1000   continue
        bnrm2=sqrt(bnrm2)
        if( bnrm2 .eq. 0.0d0 ) bnrm2=1.0d0
c
c       ... set x=0
c
        do 1100 i=1,n
        x(i)=0
 1100   continue
c
c       ... compute initial residual r = y - A*x = y
c
        call multa(a,p1,p2,p3,p4,p5,p6,x,r,n)
        do 1200 i=1,n
        r(i)=y(i)-r(i)
 1200   continue
c
        rnrm2=0
        do 1300 i=1,n
        rnrm2=rnrm2+dble(dconjg(r(i))*r(i))
 1300   continue
        rnrm2=sqrt(rnrm2)
c
        error=rnrm2/bnrm2
        if( error .lt. eps ) return
c
c
c       ============================================================
c       ... outer iteration (restarts)
c       ============================================================
c
        do 6000 iter=1,numit
c
c       ... compute residual r = y - A*x
c
        call multa(a,p1,p2,p3,p4,p5,p6,x,r,n)
        do 2000 i=1,n
        r(i)=y(i)-r(i)
 2000   continue
c
        rnrm2=0
        do 2100 i=1,n
        rnrm2=rnrm2+dble(dconjg(r(i))*r(i))
 2100   continue
        rnrm2=sqrt(rnrm2)
c
c       ... V(:,1) = r / norm(r)
c
        do 2200 i=1,n
        v(i,1)=r(i)/rnrm2
 2200   continue
c
c       ... s = norm(r) * e1
c
        s(1)=rnrm2
        do 2300 i=2,m+1
        s(i)=0
 2300   continue
c
c       ... zero out H
c
        do 2500 j=1,m
        do 2400 i=1,m+1
        h(i,j)=0
 2400   continue
 2500   continue
c
c
c       ============================================================
c       ... inner iteration (Arnoldi + Givens)
c       ============================================================
c
        do 5000 i=1,m
c
        niter=niter+1
c
c       ... w = A * V(:,i)
c
        call multa(a,p1,p2,p3,p4,p5,p6,v(1,i),v(1,i+1),n)
c
c       ... Gram-Schmidt orthogonalization
c
        do 3200 k=1,i
        d=0
        do 3000 j=1,n
        d=d+dconjg(v(j,k))*v(j,i+1)
 3000   continue
        h(k,i)=d
        do 3100 j=1,n
        v(j,i+1)=v(j,i+1)-d*v(j,k)
 3100   continue
 3200   continue
c
c       ... h(i+1,i) = norm(w)
c
        d=0
        do 3300 j=1,n
        d=d+dconjg(v(j,i+1))*v(j,i+1)
 3300   continue
        h(i+1,i)=sqrt(d)
c
c       ... V(:,i+1) = w / h(i+1,i)
c
        temp=1.0d0/h(i+1,i)
        do 3400 j=1,n
        v(j,i+1)=v(j,i+1)*temp
 3400   continue
c
c       ... apply previous Givens rotations to column i of H
c
        do 3600 k=1,i-1
        temp     = cs(k)*h(k,i) + sn(k)*h(k+1,i)
        h(k+1,i)=-dconjg(sn(k))*h(k,i) + cs(k)*h(k+1,i)
        h(k,i)  = temp
 3600   continue
c
c       ... compute the i-th Givens rotation
c
        call cgmres_rotmat(h(i,i),h(i+1,i),cs(i),sn(i))
c
c       ... apply the new rotation to s and H
c
        temp  = cs(i)*s(i)
        s(i+1)=-dconjg(sn(i))*s(i)
        s(i)  = temp
        h(i,i)= cs(i)*h(i,i) + sn(i)*h(i+1,i)
        h(i+1,i)=0
c
c       ... check convergence
c
        error=abs(s(i+1))/bnrm2
        errs(niter)=error
c
        write(*,1010) niter, error
 1010   format(' in cgmres, iter =',i4,', rel error =',e12.5)
c
        if( error .le. eps ) then
c
c       ... converged: back-solve and update x
c
        call cgmres_backsolve(h,s,i,m+1)
        do 4200 j=1,i
        do 4100 k=1,n
        x(k)=x(k)+s(j)*v(k,j)
 4100   continue
 4200   continue
c
        return
        endif
c
 5000   continue
c
c       ... end of inner iteration, update x with current best
c
        call cgmres_backsolve(h,s,m,m+1)
        do 5200 j=1,m
        do 5100 k=1,n
        x(k)=x(k)+s(j)*v(k,j)
 5100   continue
 5200   continue
c
c       ... recompute residual for convergence check
c
        call multa(a,p1,p2,p3,p4,p5,p6,x,r,n)
        do 5300 i=1,n
        r(i)=y(i)-r(i)
 5300   continue
c
        rnrm2=0
        do 5400 i=1,n
        rnrm2=rnrm2+dble(dconjg(r(i))*r(i))
 5400   continue
        rnrm2=sqrt(rnrm2)
c
        error=rnrm2/bnrm2
        if( error .le. eps ) return
c
 6000   continue
c
c       ... the maximum number of iterations has been reached
c
        ier=4
c
        return
        end
c
c
c
c
c
        subroutine cgmres_rotmat(a,b,c,s)
        implicit real *8 (a-h,o-z)
        complex *16 a,b,c,s,temp
c
c       Compute unitary Givens rotation parameters for a and b.
c       [c, s; -conjg(s), c] * [a; b] = [r; 0]
c       with c real >= 0 and |c|^2 + |s|^2 = 1.
c
        if( abs(b) .eq. 0.0d0 ) then
        c=1.0d0
        s=0.0d0
        elseif( abs(a) .eq. 0.0d0 ) then
        c=0.0d0
        s=1.0d0
        elseif( abs(a) .ge. abs(b) ) then
        temp=b/a
        d=sqrt(1.0d0+abs(temp)**2)
        c=1.0d0/d
        s=dconjg(temp)/d
        else
        temp=a/b
        d=sqrt(1.0d0+abs(temp)**2)
        c=abs(temp)/d
        s=(temp/abs(temp))/d
        endif
c
        return
        end
c
c
c
c
c
        subroutine cgmres_backsolve(h,s,k,ldh)
        implicit real *8 (a-h,o-z)
        complex *16 h(ldh,1),s(1)
c
c       Back-solve the upper triangular system H(1:k,1:k) * y = s(1:k)
c       Solution overwrites s(1:k).
c
        do 2000 i=k,1,-1
        s(i)=s(i)/h(i,i)
        do 1000 j=1,i-1
        s(j)=s(j)-h(j,i)*s(i)
 1000   continue
 2000   continue
c
        return
        end
c
c
c
c
c
