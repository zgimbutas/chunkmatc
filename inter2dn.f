cc Copyright (C) 2009-2011: Leslie Greengard and Zydrunas Gimbutas
cc Contact: greengard@cims.nyu.edu
cc 
cc This program is free software; you can redistribute it and/or modify 
cc it under the terms of the GNU General Public License as published by 
cc the Free Software Foundation; either version 2 of the License, or 
cc (at your option) any later version.  This program is distributed in 
cc the hope that it will be useful, but WITHOUT ANY WARRANTY; without 
cc even the implied warranty of MERCHANTABILITY or FITNESS FOR A 
cc PARTICULAR PURPOSE.  See the GNU General Public License for more 
cc details. You should have received a copy of the GNU General Public 
cc License along with this program; 
cc if not, see <http://www.gnu.org/licenses/>.
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c    $Date$
c    $Revision$
c
c       
c     This file contains interaction kernels for the Laplace and
c     Helmholtz equations in R^2.
c
c
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c       This is the end of the debugging routines and the beginning of
c       the interaction routines in R^2
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c       The calling sequence is:
c
c       subroutine interact(srcinfo,targinfo,cout,par1,par2,par3,par4)
c
c       srcinfo - geometry information at the source
c       targinfo - geometry information at the target
c       cout - the complex *16 value of interaction kernel
c       par1, par2, par3, par4 - extra parameters for future extensions
c
c       geometry information is encoded as a linear array containing
c          location, normal, and two tangent vectors, for a total of 6
c          real *8 elements, e.g. srcinfo 
c
c          src(1)=srcinfo(1)
c          src(2)=srcinfo(2)
c          srcnorm(1)=srcinfo(3)
c          srcnorm(2)=srcinfo(4)
c          srctang(1)=srcinfo(5)
c          srctang(2)=srcinfo(6)
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
c       ... Laplace kernels
c
c
        subroutine lfinter1(srcinfo,targinfo,cout,par1,par2,par3,par4)
        implicit real *8 (a-h,o-z)
        dimension srcinfo(*),targinfo(*)
        dimension src(2),targ(2)
        complex *16 cout
c
c       ... single layer potential, S_0
c
        src(1)=srcinfo(1)
        src(2)=srcinfo(2)
c
        targ(1)=targinfo(1)
        targ(2)=targinfo(2)
c
        dx=targ(1)-src(1)
        dy=targ(2)-src(2)
c
        r=sqrt(dx**2+dy**2)
c
ccc        cout=log(1/r)
        cout=-log(r)
c
        return
        end
c
c
c
c
c
        subroutine lfinter2(srcinfo,targinfo,cout,par1,par2,par3,par4)
        implicit real *8 (a-h,o-z)
        dimension srcinfo(*),targinfo(*)
        dimension src(2),targ(2)
        dimension srcnorm(2),targnorm(2)
        complex *16 cout
c
c       ... double layer potential, D_0
c
        src(1)=srcinfo(1)
        src(2)=srcinfo(2)
        srcnorm(1)=srcinfo(3)
        srcnorm(2)=srcinfo(4)
c
        targ(1)=targinfo(1)
        targ(2)=targinfo(2)
        targnorm(1)=targinfo(3)
        targnorm(2)=targinfo(4)
c
        dx=targ(1)-src(1)
        dy=targ(2)-src(2)
c
        r=sqrt(dx**2+dy**2)
        d=dx*srcnorm(1)+dy*srcnorm(2)
c
        cout=d/r**2
c
        return
        end
c
c
c
c
c
        subroutine lfinter3(srcinfo,targinfo,cout,par1,par2,par3,par4)
        implicit real *8 (a-h,o-z)
        dimension srcinfo(*),targinfo(*)
        dimension src(2),targ(2)
        dimension srcnorm(2),targnorm(2)
        complex *16 cout
c
c       ... derivative of single layer potential at the target
c       warning, this routine actually computes -S_0'
c
        src(1)=srcinfo(1)
        src(2)=srcinfo(2)
        srcnorm(1)=srcinfo(3)
        srcnorm(2)=srcinfo(4)
c
        targ(1)=targinfo(1)
        targ(2)=targinfo(2)
        targnorm(1)=targinfo(3)
        targnorm(2)=targinfo(4)
c
        dx=targ(1)-src(1)
        dy=targ(2)-src(2)
c
        r=sqrt(dx**2+dy**2)
        d=dx*targnorm(1)+dy*targnorm(2)
c
        cout=d/r**2
c
        return
        end
c
c
c
c
c
        subroutine lfinter4(srcinfo,targinfo,cout,par1,par2,par3,par4)
        implicit real *8 (a-h,o-z)
        dimension srcinfo(*),targinfo(*)
        dimension src(2),targ(2)
        dimension srcnorm(2),targnorm(2)
        complex *16 cout
c
c       ... derivative of double layer potential at the target
c       warning, this routine actually computes -D_0'
c
        src(1)=srcinfo(1)
        src(2)=srcinfo(2)
        srcnorm(1)=srcinfo(3)
        srcnorm(2)=srcinfo(4)
c
        targ(1)=targinfo(1)
        targ(2)=targinfo(2)
        targnorm(1)=targinfo(3)
        targnorm(2)=targinfo(4)
c
        dx=targ(1)-src(1)
        dy=targ(2)-src(2)
c
        r=sqrt(dx**2+dy**2)
c
        d=0
        d=d+(dx**2-dy**2)*targnorm(1)*srcnorm(1)
        d=d+(dy**2-dx**2)*targnorm(2)*srcnorm(2)
c
        d=d+2*(dx*dy)*targnorm(1)*srcnorm(2)
        d=d+2*(dy*dx)*targnorm(2)*srcnorm(1)
c
        cout=d/r**4
c
        return
        end
c
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
c       ... Helmholtz kernels
c
c
        subroutine hfinter1(srcinfo,targinfo,cout,zk,par2,par3,par4)
        implicit real *8 (a-h,o-z)
        complex *16 zk,ima
        dimension srcinfo(*),targinfo(*)
        dimension src(2),targ(2)
        complex *16 cout,z,h0,h1
        data ima/(0.0d0,1.0d0)/
c       
c       ... single layer potential, S_k
c
        src(1)=srcinfo(1)
        src(2)=srcinfo(2)
c
        targ(1)=targinfo(1)
        targ(2)=targinfo(2)
c
        dx=targ(1)-src(1)
        dy=targ(2)-src(2)
c
        r=sqrt(dx**2+dy**2)
c
        ifexpon=1
        z=zk*r
        call hank103(z,h0,h1,ifexpon)
c
        cout=h0
c
        return
        end
c
c
c
c
c
        subroutine hfinter2(srcinfo,targinfo,cout,zk,par2,par3,par4)
        implicit real *8 (a-h,o-z)
        complex *16 zk,ima
        dimension srcinfo(*),targinfo(*)
        dimension src(3),targ(3)
        dimension srcnorm(3),targnorm(3)
        complex *16 cout,z,h0,h1
        data ima/(0.0d0,1.0d0)/
c
c       ... double layer potential, D_k
c
        src(1)=srcinfo(1)
        src(2)=srcinfo(2)
        srcnorm(1)=srcinfo(3)
        srcnorm(2)=srcinfo(4)
c
        targ(1)=targinfo(1)
        targ(2)=targinfo(2)
        targnorm(1)=targinfo(3)
        targnorm(2)=targinfo(4)
c
        dx=targ(1)-src(1)
        dy=targ(2)-src(2)
c
        r=sqrt(dx**2+dy**2)
c
        ifexpon=1
        z=zk*r
        call hank103(z,h0,h1,ifexpon)
c
        d=dx*srcnorm(1)+dy*srcnorm(2)
        cout=d*zk*h1/r
c
        return
        end
c
c
c
c
c
        subroutine hfinter3(srcinfo,targinfo,cout,zk,par2,par3,par4)
        implicit real *8 (a-h,o-z)
        complex *16 zk,ima
        dimension srcinfo(*),targinfo(*)
        dimension src(3),targ(3)
        dimension srcnorm(3),targnorm(3)
        complex *16 cout,z,h0,h1
        data ima/(0.0d0,1.0d0)/
c
c       ... derivative of single layer potential at the target
c       warning, this routine actually computes -S_k'
c
        src(1)=srcinfo(1)
        src(2)=srcinfo(2)
        srcnorm(1)=srcinfo(3)
        srcnorm(2)=srcinfo(4)
c
        targ(1)=targinfo(1)
        targ(2)=targinfo(2)
        targnorm(1)=targinfo(3)
        targnorm(2)=targinfo(4)
c
        dx=targ(1)-src(1)
        dy=targ(2)-src(2)
c
        r=sqrt(dx**2+dy**2)
c
        ifexpon=1
        z=zk*r
        call hank103(z,h0,h1,ifexpon)
c
        d=dx*targnorm(1)+dy*targnorm(2)
        cout=d*zk*h1/r
c
        return
        end
c
c
c
c
c
        subroutine hfinter4(srcinfo,targinfo,cout,zk,par2,par3,par4)
        implicit real *8 (a-h,o-z)
        complex *16 zk,ima
        dimension srcinfo(*),targinfo(*)
        dimension src(3),targ(3)
        dimension srcnorm(3),targnorm(3)
        complex *16 cout,z,h0,h1,cd
        data ima/(0.0d0,1.0d0)/
c
c       ... derivative of double layer potential at the target
c       warning, this routine actually computes -D_k'
c
        src(1)=srcinfo(1)
        src(2)=srcinfo(2)
        srcnorm(1)=srcinfo(3)
        srcnorm(2)=srcinfo(4)
c
        targ(1)=targinfo(1)
        targ(2)=targinfo(2)
        targnorm(1)=targinfo(3)
        targnorm(2)=targinfo(4)
c
        dx=targ(1)-src(1)
        dy=targ(2)-src(2)
c
        r=sqrt(dx**2+dy**2)
c
        ifexpon=1
        z=zk*r
        call hank103(z,h0,h1,ifexpon)

c
c> simplify(diff(HankelH1(0,k*r),x,y));
c     2    2 1/2                 2    2 1/2                        2    2 1/2
c- ((x  + y )    HankelH1(0, k (x  + y )   ) k - 2 HankelH1(1, k (x  + y )   ))
c
c            /   2    2 3/2
c    x y k  /  (x  + y )
c          /
c
c> simplify(diff(HankelH1(0,k*r),x,x));
c    2   2    2 1/2                 2    2 1/2
c- (x  (x  + y )    HankelH1(0, k (x  + y )   ) k
c
c        2                 2    2 1/2     2                 2    2 1/2       /
c     - x  HankelH1(1, k (x  + y )   ) + y  HankelH1(1, k (x  + y )   )) k  /
c                                                                          /
c
c      2    2 3/2
c    (x  + y )
c

        cd=0
        cd=cd+((dx**2-dy**2)*h1-dx**2*r*h0*zk)*targnorm(1)*srcnorm(1)
        cd=cd+((dy**2-dx**2)*h1-dy**2*r*h0*zk)*targnorm(2)*srcnorm(2)
c
        cd=cd+(2*h1-r*h0*zk)*(dx*dy)*targnorm(1)*srcnorm(2)
        cd=cd+(2*h1-r*h0*zk)*(dy*dx)*targnorm(2)*srcnorm(1)
c
        cout=cd*zk/r**3
c
        return
        end
c
c
c
c
c
