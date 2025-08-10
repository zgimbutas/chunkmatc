cc Copyright (C) 2012: Zydrunas Gimbutas 
cc 
cc This software is being released under a modified FreeBSD license
cc (see COPYING in home directory). 
c
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c       This file contains a suite of routines for discretizing
c       integral equations in R^2
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
c       Geometry descriptor routines:
c       
c       fpatchpnt - flat segment in R^2
c       qpatchpnt - quadratic segment in R^2
c       spatchpnt - unit circle in R^2
c       epatchpnt - ellipse in R^2
c
c
c       Geometry refinement routines:
c
c       genrefineinfo - uniform refinement
c       genrefineinfo_dyadic - two-sided dyadic refinement
c       genrefineinfo_dyadic_url - uniform + two-sided dyadic refinement
c       genrefineinfo_dyadic_ul - uniform + left-sided dyadic refinement
c       genrefineinfo_dyadic_ur - uniform + right-sided dyadic refinement
c
c       rpatchpnt - refined geometry descriptor
c
c
c       Geometry processing routines:
c
c       chunkgeo - maps the standard simplex chunk at the location (t)
c           and returns the location of the point in R^2, the
c           derivatives with respect to parametrization, the length
c           element, the normal, and the tangent.
c
c       chunkallpnts - construct all discretization points, normals and tangents
c       chunkallwhts - construct all discretization weights
c
c       
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
c       Matrix generation routines:
c
c       chunkmatc - generate the (complex) matrix of interactions, on a
c           user-defined geometry, described by chunkpnt subroutine.
c
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
c       The calling sequence for the geometry descriptor routine is:
c
c       subroutine chunkpnt(ichunk,t,xy,dxydt,par1,par2,par3,par4)
c
c       Input parameters:
c
c       ichunk - the index of the chunk
c       t - parametrization parameter of the chunk on interval [0,1]
c       par1, par2, par3, par4 - extra parameters for future extensions
c
c       Output parameters: 
c
c       xy - the location of the point in R^2: real*8 xy(2)
c       dxydt - the derivatives of coordinate functions with
c               respect to (t) parameterization, real*8 dxydt(2)
c
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
c       The calling sequence for the interaction routine is:
c
c       subroutine interact(srcinfo,targinfo,cout,par1,par2,par3,par4)
c       
c       Input parameters:
c
c       srcinfo - geometry information at the source
c       targinfo - geometry information at the target
c       par1, par2, par3, par4 - extra parameters for future extensions
c
c       Output parameters: 
c
c       cout - the complex *16 value of interaction kernel
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
c
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c       the chunk geometry descriptor routines 
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
c
        subroutine fchunkpnt(ichunk,t,xy,dxydt,seginfo,par2,par3,par4)
        implicit real *8 (a-h,o-z)
        dimension xy(2),dxydt(2),seginfo(2,2,*)
c
c       This subroutine maps the standard interval (0,1)
c       into R^2
c
c       ... setup a flat segment in R^2
c
c       0 .. . .. 1
c
c
        x0=seginfo(1,1,ichunk)
        y0=seginfo(2,1,ichunk)
c
        x1=seginfo(1,2,ichunk)
        y1=seginfo(2,2,ichunk)
c
c
c       ... process the geometry, return the point location and the
c       derivatives with respect to t
c
        x=x0+t*(x1-x0)
        y=y0+t*(y1-y0)
c
        xy(1)=x
        xy(2)=y
c
        dxydt(1)=x1-x0
        dxydt(2)=y1-y0
c
        return
        end
c
c
c
c
c
        subroutine qchunkpnt(ichunk,t,xy,dxydt,seginfo,par2,par3,par4)
        implicit real *8 (a-h,o-z)
        dimension xy(2),dxydt(2),seginfo(2,3,*)
c
c       This subroutine maps the standard interval (0,1)
c       into R^2
c
c       ... setup a quadratic segment in R^2
c
c       0 .. A .. 1
c
c
        x0=seginfo(1,1,ichunk)
        y0=seginfo(2,1,ichunk)
c
        x1=seginfo(1,2,ichunk)
        y1=seginfo(2,2,ichunk)
c
        xa=seginfo(1,3,ichunk)
        ya=seginfo(2,3,ichunk)
c
c
c       ... process the geometry, return the point location and the
c       derivatives with respect to t
c
c
        xt=-3*x0+4*xa-x1
        yt=-3*y0+4*ya-y1
c
        xtt=x0-2*xa+x1
        ytt=y0-2*ya+y1
c
        x=x0+t*xt+2*t*t*xtt
        y=y0+t*yt+2*t*t*ytt
c
        xy(1)=x
        xy(2)=y
c
        dxydt(1)=xt+4*t*xtt
        dxydt(2)=yt+4*t*ytt
c
        return
        end
c
c
c
c
c
        subroutine schunkpnt(ichunk,t,xy,dxydt,par1,par2,par3,par4)
        implicit real *8 (a-h,o-z)
        dimension xy(2),dxydt(2)
c
c       This subroutine maps the standard interval (0,1)
c       into R^2
c
c       ... setup a circle in R^2
c
c
c
c       ... process the geometry, return the point location and the
c       derivatives with respect to t
c
c
        done=1
        pi=4*atan(done)
c
        x=cos(2*pi*t)
        y=sin(2*pi*t)
c
        xy(1)=x
        xy(2)=y
c
        dxydt(1)=2*pi*(-sin(2*pi*t))
        dxydt(2)=2*pi*(+cos(2*pi*t))
c
        return
c
        sx = 1.2d0
        sy = 1.3d0
c
c       ... ellipse
c
        xy(1)=xy(1)*sx
        xy(2)=xy(2)*sy
c
        dxydt(1)=dxydt(1)*sx
        dxydt(2)=dxydt(2)*sy
c
        return
        end
c
c
c
c
c
        subroutine echunkpnt(ichunk,t,xy,dxydt,par1,par2,par3,par4)
        implicit real *8 (a-h,o-z)
        dimension xy(2),dxydt(2)
c
c       This subroutine maps the standard interval (0,1)
c       into R^2
c
c       ... setup an ellipse in R^2
c
c
c
c       ... process the geometry, return the point location and the
c       derivatives with respect to t
c
c
        done=1
        pi=4*atan(done)
c
        x=cos(2*pi*t)
        y=sin(2*pi*t)
c
        xy(1)=x
        xy(2)=y
c
        dxydt(1)=2*pi*(-sin(2*pi*t))
        dxydt(2)=2*pi*(+cos(2*pi*t))
c
ccc        return
c
c       ... ellipse
c
        sx = 2.0d0
        sy = 2.0d0
c
        xy(1)=xy(1)*sx
        xy(2)=xy(2)*sy
c
        dxydt(1)=dxydt(1)*sx
        dxydt(2)=dxydt(2)*sy
c
        return
        end
c
c
c
c
c
        subroutine wchunkpnt(ichunk,t,xy,dxydt,par1,par2,par3,par4)
        implicit real *8 (a-h,o-z)
        dimension xy(2),dxydt(2)
c
c       This subroutine maps the standard interval (0,1)
c       into R^2. This is just an example, one segment only.
c
c       ... setup a flat segment in R^2
c
        x0=1
        y0=0
c
        x1=0
        y1=1
c
c       ... process the geometry, return the point location and the
c       derivatives with respect to t
c
        x=x0+t*(x1-x0)
        y=y0+t*(y1-y0)
c
        xy(1)=x
        xy(2)=y
c
        dxydt(1)=x1-x0
        dxydt(2)=y1-y0
c
        return
        end
c
c
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c       the chunk geometry refinement routines:
c       uniform, dyadic, left/right-sided dyadic
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
c
        subroutine genrefineinfo(noversamp,nchunks,
     $     nchunksout,ichunkinfo,refineinfo)
        implicit real *8 (a-h,o-z)
        dimension ichunkinfo(*),refineinfo(2,*)
c
c       This subroutine generates a refinement table, used by chunk
c       descriptor routines.
c
c       i--+--i
c
        nchunksout=nchunks*noversamp
c
        done=1
        dx=done/noversamp
        kk=0
c
        do 1400 i=1,nchunks
c
        do 1200 ix=1,noversamp
c
        kk=kk+1
c
c       ... uniform subdivision
c
        ichunkinfo(kk)=i
        refineinfo(1,kk)=dx
        refineinfo(2,kk)=dx*(ix-1)
c
 1200   continue
c
 1400   continue
c
c
ccc        call prinf('ichunkinfo=*',ichunkinfo,kk)
ccc        call prin2('refineinfo=*',refineinfo,2*kk)
ccc        pause
c
        return
        end
c
c
c
c
c
        subroutine genrefineinfo_dyadic(noversamp,nchunks,
     $     nchunksout,ichunkinfo,refineinfo)
        implicit real *8 (a-h,o-z)
        dimension ichunkinfo(*),refineinfo(2,*)
c
c       This subroutine generates a refinement table, used by chunk
c       descriptor routines.
c
c       i--+--i
c
        nchunksout=nchunks*noversamp
c
        done=1
        kk=0
c
        do 1400 i=1,nchunks
c
        sx=0
        dx=done/2.0d0**(noversamp/2)
c
        kk=kk+1
        ichunkinfo(kk)=i
        refineinfo(1,kk)=dx
        refineinfo(2,kk)=sx
c
        sx=sx+dx
        ntot=noversamp/2
c
        do 1200 ix=2,ntot
c
        kk=kk+1
c
c       ... dyadic subdivision
c
        ichunkinfo(kk)=i
        refineinfo(1,kk)=dx
        refineinfo(2,kk)=sx
c
ccc        write(*,*) ix, sx,dx, sx+dx
c
        sx=sx+dx
        dx=dx*2
c
 1200   continue
c
        dx=dx/2
c
        do 1300 ix=ntot+1,2*ntot-1
c
        kk=kk+1
c
c       ... dyadic subdivision
c
        ichunkinfo(kk)=i
        refineinfo(1,kk)=dx
        refineinfo(2,kk)=sx
c
ccc        write(*,*) ix, sx,dx, sx+dx
c
        sx=sx+dx
        dx=dx/2
c
 1300   continue
c
        dx=dx*2
c
        kk=kk+1
c
        ichunkinfo(kk)=i
        refineinfo(1,kk)=dx
        refineinfo(2,kk)=sx
c
ccc        write(*,*) sx,dx, sx+dx
c
 1400   continue
c
c
ccc        call prinf('ichunkinfo=*',ichunkinfo,kk)
ccc        call prin2('refineinfo=*',refineinfo,2*kk)
ccc        pause
c
        return
        end
c
c
c
c
c
        subroutine genrefineinfo_dyadic_ulr(noversamp,nmiddle,nchunks,
     $     nchunksout,ichunkinfo,refineinfo)
        implicit real *8 (a-h,o-z)
        dimension ichunkinfo(*),refineinfo(2,*)
c
c       This subroutine generates a refinement table, used by chunk
c       descriptor routines.
c
c       i--+--i
c
        nchunksout=nchunks*(noversamp+nmiddle)
c
        done=1
        kk=0
c
        do 1400 i=1,nchunks
c
        sx=0
        dx=done/2.0d0**(noversamp/4)
c
        kk=kk+1
        ichunkinfo(kk)=i
        refineinfo(1,kk)=dx
        refineinfo(2,kk)=sx
c
        kstart=kk
c
        sx=sx+dx
        ntot=noversamp/2
c
        do 1200 ix=2,ntot
c
        kk=kk+1
c
c       ... dyadic subdivision
c
        ichunkinfo(kk)=i
        refineinfo(1,kk)=dx
        refineinfo(2,kk)=sx
c
ccc        write(*,*) ix, sx,dx, sx+dx
c
        sx=sx+dx
        dx=dx*2
c
 1200   continue
c
        do 1250 ix=1,nmiddle
c
        kk=kk+1
c
c       ... uniform middle intervals
c
        ichunkinfo(kk)=i
        refineinfo(1,kk)=dx
        refineinfo(2,kk)=sx
c
        sx=sx+dx
c
 1250   continue
c
        dx=dx/2
c
        do 1300 ix=ntot+1,2*ntot-1
c
        kk=kk+1
c
c       ... dyadic subdivision
c
        ichunkinfo(kk)=i
        refineinfo(1,kk)=dx
        refineinfo(2,kk)=sx
c
ccc        write(*,*) ix, sx,dx, sx+dx
c
        sx=sx+dx
        dx=dx/2
c
 1300   continue
c
        dx=dx*2
c
        kk=kk+1
c
        ichunkinfo(kk)=i
        refineinfo(1,kk)=dx
        refineinfo(2,kk)=sx
c
ccc        write(*,*) sx,dx, sx+dx
c
        sx=sx+dx
c
ccc        write(*,*) sx
c        
        do 1350 ix=kstart,kstart+2*ntot+nmiddle-1
c
c       ... rescale all intervals
c
        refineinfo(1,ix)=refineinfo(1,ix)/sx
        refineinfo(2,ix)=refineinfo(2,ix)/sx
ccc        write(*,*) ix, refineinfo(1,ix),refineinfo(2,ix)
c
 1350   continue
c
 1400   continue
c
c
ccc        call prinf('ichunkinfo=*',ichunkinfo,kk)
ccc        call prin2('refineinfo=*',refineinfo,2*kk)
ccc        pause
c
        return
        end
c
c
c
c
c
        subroutine genrefineinfo_dyadic_ul(noversamp,nmiddle,nchunks,
     $     nchunksout,ichunkinfo,refineinfo)
        implicit real *8 (a-h,o-z)
        dimension ichunkinfo(*),refineinfo(2,*)
c
c       This subroutine generates a refinement table, used by chunk
c       descriptor routines.
c
c       i--+--i
c
        nchunksout=nchunks*(noversamp+nmiddle)
c
        done=1
        kk=0
c
        do 1400 i=1,nchunks
c
        sx=0
        dx=done/2.0d0**(noversamp)
c
        kk=kk+1
        ichunkinfo(kk)=i
        refineinfo(1,kk)=dx
        refineinfo(2,kk)=sx
c
        kstart=kk
c
        sx=sx+dx
        ntot=noversamp
c
        do 1200 ix=2,ntot
c
        kk=kk+1
c
c       ... dyadic subdivision
c
        ichunkinfo(kk)=i
        refineinfo(1,kk)=dx
        refineinfo(2,kk)=sx
c
ccc        write(*,*) ix, sx,dx, sx+dx
c
        sx=sx+dx
        dx=dx*2
c
 1200   continue
c
        do 1250 ix=1,nmiddle
c
        kk=kk+1
c
c       ... uniform middle intervals
c
        ichunkinfo(kk)=i
        refineinfo(1,kk)=dx
        refineinfo(2,kk)=sx
c
        sx=sx+dx
c
 1250   continue
c
c        
        do 1350 ix=kstart,kstart+noversamp+nmiddle-1
c
c       ... rescale all intervals
c
        refineinfo(1,ix)=refineinfo(1,ix)/sx
        refineinfo(2,ix)=refineinfo(2,ix)/sx
ccc        write(*,*) ix, refineinfo(1,ix),refineinfo(2,ix)
c
 1350   continue
c
 1400   continue
c
c
ccc        call prinf('ichunkinfo=*',ichunkinfo,kk)
ccc        call prin2('refineinfo=*',refineinfo,2*kk)
ccc        pause
c
        return
        end
c
c
c
c
c
        subroutine genrefineinfo_dyadic_ur(noversamp,nmiddle,nchunks,
     $     nchunksout,ichunkinfo,refineinfo)
        implicit real *8 (a-h,o-z)
        dimension ichunkinfo(*),refineinfo(2,*)
c
c       This subroutine generates a refinement table, used by chunk
c       descriptor routines.
c
c       i--+--i
c
        nchunksout=nchunks*(noversamp+nmiddle)
c
        done=1
        kk=0
c
        do 1400 i=1,nchunks
c
        sx=0
        dx=done/2.0d0**(noversamp)
c
c
        sx=0
        dx=1.0d0/(nmiddle+3)
c
        do 1250 ix=1,nmiddle
c
        kk=kk+1
        if(ix .eq. 1 ) kstart=kk
c
c       ... uniform middle intervals
c
        ichunkinfo(kk)=i
        refineinfo(1,kk)=dx
        refineinfo(2,kk)=sx
c
        sx=sx+dx
c
 1250   continue
c
        dx=dx/2
c
        do 1300 ix=nmiddle+1,nmiddle+noversamp-1
c
        kk=kk+1
c
c       ... dyadic subdivision
c
        ichunkinfo(kk)=i
        refineinfo(1,kk)=dx
        refineinfo(2,kk)=sx
c
ccc        write(*,*) ix, sx,dx, sx+dx
c
        sx=sx+dx
        dx=dx/2
c
 1300   continue
c
        dx=dx*2
c
        kk=kk+1
c
        ichunkinfo(kk)=i
        refineinfo(1,kk)=dx
        refineinfo(2,kk)=sx
c
ccc        write(*,*) sx,dx, sx+dx
c
        sx=sx+dx
c
ccc        write(*,*) sx
c        
        do 1350 ix=kstart,kstart+noversamp+nmiddle-1
c
c       ... rescale all intervals
c
        refineinfo(1,ix)=refineinfo(1,ix)/sx
        refineinfo(2,ix)=refineinfo(2,ix)/sx
ccc        write(*,*) ix, refineinfo(1,ix),refineinfo(2,ix)
c
 1350   continue
c
 1400   continue
c
c
ccc        call prinf('ichunkinfo=*',ichunkinfo,kk)
ccc        call prin2('refineinfo=*',refineinfo,2*kk)
ccc        pause
c
        return
        end
c
c
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
c
        subroutine rchunkpnt
     $     (ichunk,t,xy,dxydt,par1,chunkpnt,ichunkinfo,refineinfo)
        implicit real *8 (a-h,o-z)
        dimension xy(2),dxydt(2)        
        dimension ichunkinfo(*),refineinfo(2,*)
        external chunkpnt
c
c       This subroutine provides on-fly refinement routines 
c       for arbitrary chunks
c
c
c       ... refine the parametrization
c       
        jchunk=ichunkinfo(ichunk)
        tref=t*refineinfo(1,ichunk)+refineinfo(2,ichunk)
c
        call chunkpnt(jchunk,tref,xy,dxydt,par1,
     $     xpar2,xpar3,xpar4)
c
c       ... adjust all derivatives
c       
        dxydt(1)=dxydt(1)*refineinfo(1,ichunk)
        dxydt(2)=dxydt(2)*refineinfo(1,ichunk)
c
        return
        end
c
c
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c       the discretization routines in R^2
c       geometry processing
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c       
c       
        subroutine chunkgeo(chunkpnt,ichunk,t,par1,par2,par3,par4,
     $     xy,dxydt,ds,xynorm,xytang)
        implicit real *8 (a-h,o-z)
        external chunkpnt
        dimension xy(2),dxydt(2)
        dimension xynorm(2),xytang(2)
c
c       This subroutine maps the standard simplex chunk ichunk 
c       at the location (t) and returns the location of the point in R^2,
c       the derivatives with respect to parametrization, the length element,
c       the normal, and the tangent.
c
c       Input parameters:
c
c       chunkpnt: external: must be of the form
c                 chunkpnt(ichunk,t,xy,dxydt,par1,par2,par3,par4)
c       ichunk: integer: the index of the chunk
c       t: real *8: parametrization of the standard simplex chunk
c       par1,par2,par3,par4: extra parameters
c
c       Output parameters:
c       
c       xy: real*8(2): the location of the point in R^2
c       dxydt: real*8(2): the derivatives of coordinate functions with
c             respect to (t) parameterization
c       ds: real *8: the length element
c       xynorm: real*8(2): the normal vector
c       xytang: real*8(2): the tangent vector
c
c
c       ... retrieve a point on the chunk
c
        call chunkpnt(ichunk,t,xy,dxydt,par1,par2,par3,par4)
c
ccc        call prin2('xy=*',xy,2)
ccc        call prin2('dxydt=*',dxydt,2)
c
c       ... find the normal
c
        xynorm(1)=-dxydt(2)
        xynorm(2)=+dxydt(1)
c
c       ... find the length element
c
        d=sqrt(xynorm(1)**2+xynorm(2)**2)
c
        ds=d
c
ccc        call prin2('ds=*',ds,1)
c
        xynorm(1)=xynorm(1)/ds
        xynorm(2)=xynorm(2)/ds
c
c
c       ... construct the tangent vector
c       
        xytang(1)=+xynorm(2)
        xytang(2)=-xynorm(1)
c
ccc        call prin2('xynorm=*',xynorm,2)
ccc        call prin2('xytang=*',xytang,2)
c
        return
        end
c
c
c
c
c
        subroutine chunkinfo(xy,xynorm,xytang,xyinfo)
        implicit real *8 (a-h,o-z)
        dimension xy(2),xynorm(2),xytang(2)
        dimension xyinfo(6)
c
c       This subroutine compresses the geometry information into
c       the standard linear array.
c
c       Input parameters:
c       
c       xy: real*8(2): the location of the point in R^2
c       xynorm: real*8(2): the normal vector
c       xytang: real*8(2): the tangent vector
c
c       Output parameters:
c
c       xyinfo: real*8(6): the compressed geometry structure
c
c       ... compress geometry information
c
        xyinfo(1)=xy(1)
        xyinfo(2)=xy(2)
c
        xyinfo(3)=xynorm(1)
        xyinfo(4)=xynorm(2)
c
        xyinfo(5)=xytang(1)
        xyinfo(6)=xytang(2)
c
        return
        end
c
c
c
c
c
        subroutine chunkallpnts(nchunks,chunkpnt,par1,par2,par3,par4,
     $     npols,ts,ixys,xys,xynorms,xytangs,npts)
        implicit real *8 (a-h,o-z)
        external chunkpnt
        dimension ts(*)
        dimension xy(2),dxydt(2)
        dimension xynorm(2),xytang(2)
c
        dimension ixys(2,*),xys(2,*),xynorms(2,*),xytangs(2,*)
c
c       This subroutine return all points, normals and tangents from
c       geometry descriptor
c
c       Input parameters:
c
c       nchunks: integer: the number of chunks
c       chunkpnt: external: must be of the form
c                 chunkpnt(ichunk,t,xy,dxydt,par1,par2,par3,par4)
c       par1,par2,par3,par4: extra parameters
c       npols: integer: the total number of polynomials for each chunk
c       ts: real *8(*): local t-discretization points for each chunk
c
c       Output parameters:
c       
c       ixys: integer(2,npts): the index array for each chunk
c       xys: real*8(3,npts): the location of the points in R^2
c       xynorms: real*8(3,npts): the normal vectors
c       xytangs: real*8(3,npts): the first standard tangent vectors
c       npts: integer: the total number of points in discretization
c
c
        npts=0
        kk=1
c
        do 1400 ichunk=1,nchunks
        ixys(1,ichunk)=kk
        ixys(2,ichunk)=npols
        do 1200 i=1,npols
c
        t=ts(i)
c
        call chunkgeo(chunkpnt,ichunk,t,
     $     par1,par2,par3,par4,
     $     xy,dxydt,ds,xynorm,xytang)
c
        xys(1,kk)=xy(1)
        xys(2,kk)=xy(2)
c
        xynorms(1,kk)=xynorm(1)
        xynorms(2,kk)=xynorm(2)
c
        xytangs(1,kk)=xytang(1)
        xytangs(2,kk)=xytang(2)
c
        npts=npts+1
        kk=kk+1
c
 1200   continue
 1400   continue
c
c
        return
        end
c
c
c
c
c
        subroutine chunkallwhts(nchunks,chunkpnt,par1,par2,par3,par4,
     $     npols,ts,ws,whts,npts)
        implicit real *8 (a-h,o-z)
        external chunkpnt
        dimension ts(*),ws(*)
        dimension xy(2),dxydt(2)
        dimension xynorm(2),xytang(2)
c
c       This subroutine return the discretization weights
c
c       Input parameters:
c
c       nchunks: integer: the number of chunks
c       chunkpnt: external: must be of the form
c                 chunkpnt(ichunk,t,xy,dxydt,par1,par2,par3,par4)
c       par1,par2,par3,par4: extra parameters
c       npols: integer: the total number of polynomials for each chunk
c       ts: real *8(*): local t-discretization points for each chunk
c       ws: real *8(*): local integration weights for each chunk
c
c       Output parameters:
c       
c       whts: real*8(npts): the discretization weights
c       npts: integer: the total number of points in discretization
c
c
        dimension whts(*)
c
        npts=0
        kk=1
c
        do 1400 ichunk=1,nchunks
        do 1200 i=1,npols
c
        t=ts(i)
c
        call chunkgeo(chunkpnt,ichunk,t,
     $     par1,par2,par3,par4,
     $     xy,dxydt,ds,xynorm,xytang)
c
        whts(kk)=ds*ws(i)
c
        npts=npts+1
        kk=kk+1
c
 1200   continue
 1400   continue
c
c
        return
        end
c
c
c
c
c
        subroutine chunkallrads(nchunks,chunkpnt,par1,par2,par3,par4,
     $     npols,ts,centers,radii)
        implicit real *8 (a-h,o-z)
        external chunkpnt
        dimension ts(*)
        dimension xy(2),dxydt(2)
        dimension xynorm(2),xytang(2)
c
c       This subroutine return the discretization weights
c
c       Input parameters:
c
c       nchunks: integer: the number of chunks
c       chunkpnt: external: must be of the form
c                 chunkpnt(ichunk,t,xy,dxydt,par1,par2,par3,par4)
c       par1,par2,par3,par4: extra parameters
c       npols: integer: the total number of polynomials for each chunk
c       ts: real *8(*): local t-discretization points for each chunk
c       ws: real *8(*): local integration weights for each chunk
c
c       Output parameters:
c       
c       whts: real*8(npts): the discretization weights
c       npts: integer: the total number of points in discretization
c
c
        dimension centers(2,*),radii(*)
c
        npts=0
        kk=1
c
        do 1400 ichunk=1,nchunks
c
        t=0.5d0
c
        call chunkgeo(chunkpnt,ichunk,t,
     $     par1,par2,par3,par4,
     $     xy,dxydt,ds,xynorm,xytang)
c
        centers(1,ichunk)=xy(1)
        centers(2,ichunk)=xy(2)
c
        t=0.0d0
c
        call chunkgeo(chunkpnt,ichunk,t,
     $     par1,par2,par3,par4,
     $     xy,dxydt,ds,xynorm,xytang)
c
        r0 = sqrt((centers(1,ichunk)-xy(1))**2+
     $     (centers(2,ichunk)-xy(2))**2)
c
        t=1.0d0
c
        call chunkgeo(chunkpnt,ichunk,t,
     $     par1,par2,par3,par4,
     $     xy,dxydt,ds,xynorm,xytang)
c
        r1 = sqrt((centers(1,ichunk)-xy(1))**2+
     $     (centers(2,ichunk)-xy(2))**2)
c
        radii(ichunk) = max(r0,r1)
c        
 1400   continue
c
c
        return
        end
c
c
c
c
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c       the discretization routines in R^2
c       matrix generation
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c       
c
c       
        subroutine chunkmatc(nchunks,chunkpnt,par1,par2,par3,par4,
     $     norder,npols,ts,umatr,vmatr,
     $     ixys,xys,xynorms,xytangs,npts,
     $     interact,par5,par6,par7,par8,
     $     cmatr,w,lw,lused,ier)
        implicit real *8 (a-h,o-z)
c
c
c       This subroutine generates the (complex) matrix of interactions
c       on a user-defined geometry, described by chunkpnt subroutine.
c
c       Input parameters:
c
c       nchunks: integer: the number of chunks
c       chunkpnt: external: must be of the form
c                 chunkpnt(ichunk,t,xy,dxydt,par1,par2,par3,par4)
c       par1,par2,par3,par4: extra parameters
c       norder: integer: the order of interpolation
c       npols: integer: the total number of polynomials for each chunk
c       ts: real *8(*): local t-discretization points for each chunk
c       umatr: real *8(npols,npols): first interpolation matrix
c       vmatr: real *8(npols,npols): second interpolation matrix
c       interact: external: must be of the form
c            interact(srcinfo,targinfo,cout,par5,par6,par7,par8)
c
c       Output parameters:
c       
c       cmatr: complex*16(npts,npts): the interaction matrix
c
c       lused: integer: the total 
c       ier: integer: the error code
c
c       Work arrays:
c
c       w: real *8(*): work array, must be at least 2*10000 real *8 elements
c       lw: integer: the length of the work array
c
c
        external chunkpnt,interact
        dimension ts(*)
        dimension xy(2),dxydt(2)
        dimension xynorm(2),xytang(2)
c
        dimension ixys(2,*),xys(2,*),xynorms(2,*)
        dimension xytangs(2,*)
c
        dimension w(*)
        complex *16 cmatr(npts,npts)
c
        ier=0
c
c       ... allocate work arrays
c        
c       ... max 100 points per chunk
c
        nmax=100
c
        itmatr=1
        ltmatr=2*nmax*nmax
c
        lused7=ltmatr
c
        call chunkmatc0(nchunks,chunkpnt,par1,par2,par3,par4,
     $     norder,npols,ts,umatr,vmatr,
     $     ixys,xys,xynorms,xytangs,npts,
     $     interact,par5,par6,par7,par8,
     $     cmatr,w(itmatr),w(1+lused7),lw-lused7,lused8,ier)        
c
        lused=lused7+lused8
c
        return
        end
c
c
c
c
c
        subroutine chunkmatc0(nchunks,chunkpnt,par1,par2,par3,par4,
     $     norder,npols,ts,umatr,vmatr,
     $     ixys,xys,xynorms,xytangs,npts,
     $     interact,par5,par6,par7,par8,
     $     cmatr,tmatr,w,lw,lused,ier)
        implicit real *8 (a-h,o-z)
        external chunkpnt,interact
        dimension ts(*)
        dimension xy(2),dxydt(2)
        dimension xynorm(2),xytang(2)
c
        dimension ixys(2,*),xys(2,*),xynorms(2,*)
        dimension xytangs(2,*)
c
        dimension w(*)
        real *8, allocatable :: w_omp(:)
        complex *16 tmatr(*)
        complex *16 cmatr(npts,npts)
        complex *16 tmatr_omp(10000)
c
c
c       ... construct the off-diagonal blocks
c
ccc        do 1400 k=1,1
C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP$PRIVATE(i,j,ii,ipols,jj,jpols,tmatr_omp,lused,ier,w_omp) 
        do 1400 j=1,nchunks
c
        allocate(w_omp(2000000))
c
ccc        call prinf('od,j=*',j,1)
c
        do 1200 i=1,nchunks
        if ( i .eq. j ) goto 1200
c        
c       ... (i,j), i index - target, j index - source
c
        ii=ixys(1,i)
        jj=ixys(1,j)
        ipols=ixys(2,i)
        jpols=ixys(2,j)
c
cc        call prinf('i=*',i,1)
cc        call prinf('j=*',j,1)
c
        call chunkmatc_od(i,ipols,j,jpols,
     $     nchunks,chunkpnt,par1,par2,par3,par4,
     $     norder,npols,ts,umatr,vmatr,
     $     ixys,xys,xynorms,xytangs,npts,
     $     interact,par5,par6,par7,par8,
     $     tmatr_omp,w_omp,lw,lused,ier)
c
        call chunksubcpy(npts,cmatr,tmatr_omp,ii,ipols,jj,jpols)
c        
 1200   continue
        deallocate(w_omp)
 1400   continue        
C$OMP END PARALLEL DO
c
c
ccc        return
c
c       
c       ... construct the diagonal (self interaction) blocks
c
C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP$PRIVATE(i,ii,ipols,jj,jpols,tmatr_omp,lused,ier,w_omp) 
        do 2200 i=1,nchunks
c
        allocate(w_omp(2000000))
c
c       ... (i,i) 
c
        ii=ixys(1,i)
        ipols=ixys(2,i)
c
ccc        call prinf('dd,i=*',i,1)
c
        call chunkmatc_dd(i,ipols,i,ipols,
     $     nchunks,chunkpnt,par1,par2,par3,par4,
     $     norder,npols,ts,umatr,vmatr,
     $     ixys,xys,xynorms,xytangs,npts,
     $     interact,par5,par6,par7,par8,
     $     tmatr_omp,w_omp,lw,lused,ier)
c
        call chunksubcpy(npts,cmatr,tmatr_omp,ii,ipols,ii,ipols)
c              
        deallocate(w_omp)
 2200   continue
C$OMP END PARALLEL DO
c
        return
        end
c
c
c
c
c
        subroutine chunksubcpy(npts,cmatr,tmatr,ii,ipols,jj,jpols)
        implicit real *8 (a-h,o-z)
c
        complex *16 cmatr(npts,npts)
        complex *16 tmatr(ipols,jpols)
c
        do 1400 j=1,jpols
        do 1200 i=1,ipols
c
        cmatr(ii+i-1,jj+j-1)=tmatr(i,j)
 1200   continue
 1400   continue
c       
        return
        end
c
c
c
c
c
        subroutine chunkmatc_od(ichunk,ipols,jchunk,jpols,
     $     nchunks,chunkpnt,par1,par2,par3,par4,
     $     norder,npols,ts,umatr,vmatr,
     $     ixys,xys,xynorms,xytangs,npts,
     $     interact,par5,par6,par7,par8,
     $     tmatr,w,lw,lused,ier)
        implicit real *8 (a-h,o-z)
c
c       ... generate the off-diagonal block of interaction matrix
c
        external chunkpnt,interact
        dimension ts(*)
        dimension xy(2),dxydt(2)
        dimension xynorm(2),xytang(2)
c
        dimension ixys(2,*),xys(2,*),xynorms(2,*)
        dimension xytangs(2,*)
c
        dimension targinfo(6),info(20)
        dimension xpar1(20),xpar2(20)
c
        dimension w(*)
        complex *16 tmatr(ipols,jpols)
        complex *16 cvals(1000),coefs(1000)
c       
        dimension vert1(2),vert2(2),vert3(2)
        external chunkfun3
c
ccc        write(*,*) '.', ichunk, ipols, jchunk, jpols
c
        ii=ixys(1,ichunk)
        jj=ixys(1,jchunk)
c
ccc        write(*,*) '.', ii, jj
c
c       ... construct one off-diagonal block via collocation
c       
c       ... (i,j), i index - target, j index - source
c
        if( ipols .ne. npols ) then
        write(*,*) 'ipols .ne. npols'
        endif
        if( jpols .ne. npols ) then
        write(*,*) 'jpols .ne. npols'
        endif
c        
        do 1400 i=1,npols
c
c       ... on j-th chunk integrate all basis functions multiplied 
c       by interaction function at the target point ts(i)
c
c       ... first, initialize function to be integrated
c
        xpar1(1)=norder
        xpar1(2)=npols
        xpar1(3)=jchunk
c
        call chunkinfo(xys(1,ii+i-1),xynorms(1,ii+i-1),
     $     xytangs(1,ii+i-1),targinfo)
c
ccc        call prin2('=====================*',targinfo,0)
ccc        call prin2('inside _od, targinfo=*',targinfo,6)
        do j=1,6
        xpar2(j)=targinfo(j)
        enddo
c
c       ... then, call adaptive gaussian integration routine 
c       
ccc        m=6
ccc        eps=1d-6
c
        m=16
        eps=1d-13
c
        iquadtype=2
        maxrec=20
c
c
c
        a=0
        b=1
c
        if( iquadtype .eq. 1 ) then
c
        call cadavect(ier,a,b,chunkfun3,npols,
     $     chunkpnt,par1,par2,par3,par4,
     $     interact,par5,par6,par7,par8,xpar1,xpar2,
     $     m,eps,coefs,maxrec,numfunev,w,lw)
c
ccc        call prinf('numfunev=*',numfunev,1)
        if( ier .eq. 8 ) then
        write(*,*) 'maximum recursion depth of 200 has been reached'
        write(*,*) 'abort'
        stop
        endif

        endif
c
        if( iquadtype .eq. 2 ) then
c
        call cadachunk(ier,a,b,chunkfun3,npols,
     $     chunkpnt,par1,par2,par3,par4,
     $     interact,par5,par6,par7,par8,xpar1,xpar2,
     $     m,eps,coefs,maxrec,numfunev,w,lw,
     $     jchunk,targinfo,info)
c
ccc        call prinf('numfunev=*',numfunev,1)
ccc        call prinf('info=*',info,2)
        if( ier .eq. 8 ) then
        write(*,*) 'maximum recursion depth of 200 has been reached'
        write(*,*) 'abort'
        stop
        endif

        endif
c
ccc        call prinf('numfunev=*',numfunev,1)
ccc        call prinf('info=*',info,2)
ccc        call prin2('coefs=*',coefs,2*npols)
c
c
c       ... finally, convert the linear form of integral values to the
c       pointwise interation matrix, we will need umatr and vmatr for
c       this operation
c
        call chunkcoefs2cvals(npols,umatr,vmatr,coefs,cvals)
c
ccc        call prin2('cvals=*',cvals,2*npols)
c
        do 1200 j=1,npols
        tmatr(i,j)=cvals(j)
 1200   continue
c
 1400   continue        
c
c
        return
        end
c
c
c
c
c
        subroutine chunkmatc_dd(ichunk,ipols,jchunk,jpols,
     $     nchunks,chunkpnt,par1,par2,par3,par4,
     $     norder,npols,ts,umatr,vmatr,
     $     ixys,xys,xynorms,xytangs,npts,
     $     interact,par5,par6,par7,par8,
     $     tmatr,w,lw,lused,ier)
        implicit real *8 (a-h,o-z)
c
c       ... generate the diagonal block of interaction matrix,
c       self-interaction
c
        external chunkpnt,interact
        dimension ts(*)
        dimension xy(2),dxydt(2)
        dimension xynorm(2),xytang(2)
c
        dimension ixys(2,*),xys(2,*),xynorms(2,*)
        dimension xytangs(2,*)
c
        dimension targinfo(6),info(20)
        dimension xpar1(20),xpar2(20)
c
        dimension w(*)
        complex *16 tmatr(ipols,jpols)
        complex *16 cvals(1000),coefs(1000)
        complex *16 coefs1(1000),coefs2(1000),coefs3(1000)
c
        dimension xs(10000),ys(10000),ws(10000)
c       
        dimension vert1(2),vert2(2),vert3(2)
        dimension vert1a(2,3)
        external chunkfun3
c
c
ccc        write(*,*) '.', ichunk, ipols, jchunk, jpols
c
        ii=ixys(1,ichunk)
        jj=ixys(1,jchunk)
c
ccc        write(*,*) '.', ii, jj
c
c       ... construct one off-diagonal block via collocation
c       
c       ... (i,j), i index - target, j index - source
c
        if( ipols .ne. npols ) then
        write(*,*) 'ipols .ne. npols'
        endif
        if( jpols .ne. npols ) then
        write(*,*) 'jpols .ne. npols'
        endif
c        
        do 1400 i=1,npols
c
c       ... on j-th chunk integrate all basis functions multiplied 
c       by interaction function at the target point ts(i)
c
c       ... first, initialize function to be integrated
c
        xpar1(1)=norder
        xpar1(2)=npols
        xpar1(3)=jchunk
c
        call chunkinfo(xys(1,ii+i-1),xynorms(1,ii+i-1),
     $     xytangs(1,ii+i-1),targinfo)
c
ccc        call prin2('=====================*',targinfo,0)
ccc        call prin2('inside _dd, targinfo=*',targinfo,6)
        do j=1,6
        xpar2(j)=targinfo(j)
        enddo
c
cccc        write(*,*) i,ts(i)
c
c       ... then, call adaptive gaussian integration routine 
c
ccc        m=6
ccc        eps=1d-6
c
        m=16
        eps=1d-13
c
        iquadtype=2
        nrec=20
c
        if( norder .eq. 1 ) iquadtype=3
        if( norder .eq. 2 ) iquadtype=3
        if( norder .eq. 3 ) iquadtype=3
        if( norder .eq. 4 ) iquadtype=3
        if( norder .eq. 5 ) iquadtype=3
        if( norder .eq. 6 ) iquadtype=3
        if( norder .eq. 8 ) iquadtype=3
        if( norder .eq. 10 ) iquadtype=3
        if( norder .eq. 12 ) iquadtype=3
        if( norder .eq. 16 ) iquadtype=3
        if( norder .eq. 20 ) iquadtype=3
c
c
c
        if( iquadtype .eq. 1 ) then

        a=0
        b=ts(i)
c
        call cadavect(ier,a,b,chunkfun3,npols,
     $     chunkpnt,par1,par2,par3,par4,
     $     interact,par5,par6,par7,par8,xpar1,xpar2,
     $     m,eps,coefs1,maxrec,numfunev,w,lw)
c
        if( ier .eq. 8 ) then
        write(*,*) 'maximum recursion depth of 200 has been reached'
        write(*,*) 'abort'
        stop
        endif
c
        a=ts(i)
        b=1
c
        call cadavect(ier,a,b,chunkfun3,npols,
     $     chunkpnt,par1,par2,par3,par4,
     $     interact,par5,par6,par7,par8,xpar1,xpar2,
     $     m,eps,coefs2,maxrec,numfunev,w,lw)
c
        if( ier .eq. 8 ) then
        write(*,*) 'maximum recursion depth of 200 has been reached'
        write(*,*) 'abort'
        stop
        endif
c
        do j=1,npols
        coefs(j)=coefs1(j)+coefs2(j)
        enddo
c
        endif
c
c
        if( iquadtype .eq. 2 ) then

        a=0
        b=ts(i)
c
        call cadachunk(ier,a,b,chunkfun3,npols,
     $     chunkpnt,par1,par2,par3,par4,
     $     interact,par5,par6,par7,par8,xpar1,xpar2,
     $     m,eps,coefs1,maxrec,numfunev,w,lw,
     $     jchunk,targinfo,info)
c
        if( ier .eq. 8 ) then
        write(*,*) 'maximum recursion depth of 200 has been reached'
        write(*,*) 'abort'
        stop
        endif
c
        a=ts(i)
        b=1
c
        call cadachunk(ier,a,b,chunkfun3,npols,
     $     chunkpnt,par1,par2,par3,par4,
     $     interact,par5,par6,par7,par8,xpar1,xpar2,
     $     m,eps,coefs2,maxrec,numfunev,w,lw,
     $     jchunk,targinfo,info)
c
        if( ier .eq. 8 ) then
        write(*,*) 'maximum recursion depth of 200 has been reached'
        write(*,*) 'abort'
        stop
        endif
c
        do j=1,npols
        coefs(j)=coefs1(j)+coefs2(j)
        enddo
c
        endif
c
c
        if( iquadtype .eq. 3 ) then 
c
        inode=i
        t=ts(inode)
        call chunkgeo(chunkpnt,jchunk,t,
     $     par1,par2,par3,par4,
     $     xy,dxydt,ds,xynorm,xytang)
c
        call hqsuppquad(norder,inode,xs,ws,ns)
c       
ccc        call prin2('xs=*',xs,ns)
ccc        call prinf('ns=*',ns,1)
c
        do j=1,npols
        coefs(j)=0
        enddo
c
        do k=1,ns
        call chunkfun3(xs(k),npols,cvals,
     $     chunkpnt,par1,par2,par3,par4,
     $     interact,par5,par6,par7,par8,
     $     xpar1,xpar2)
        do j=1,npols
        coefs(j)=coefs(j)+ws(k)*cvals(j)
        enddo
        enddo
c
        endif
c
c
ccc        call prin2('coefs=*',coefs,2*npols)
c
c       ... finally, convert the linear form of integral values to the
c       pointwise interation matrix, we will need umatr and vmatr for
c       this operation
c
        call chunkcoefs2cvals(npols,umatr,vmatr,coefs,cvals)
c
ccc        call prin2('cvals=*',cvals,2*npols)
c
        do 1200 j=1,npols
        tmatr(i,j)=cvals(j)
 1200   continue
c
 1400   continue        
c
c
        return
        end
c
c
c
c
c
        subroutine chunkcoefs2cvals(npols,umatr,vmatr,coefs,cvals)
        implicit real *8 (a-h,o-z)
        dimension umatr(npols,npols),vmatr(npols,npols)
        complex *16 cvals(1000),coefs(1000),cd
c
        do 1400 i=1,npols
        cd=0
        do 1200 j=1,npols
        cd=cd+umatr(j,i)*coefs(j)
 1200   continue
        cvals(i)=cd
 1400   continue
c
        return
        end
c
c
c
c
c
        subroutine chunkfun3(t,npols,cvals,
     $     chunkpnt,par1,par2,par3,par4,
     $     interact,par5,par6,par7,par8,
     $     xpar1,xpar2)
        implicit real *8 (a-h,o-z)
c
c       ... provides functions to be integrated 
c      must be initialized via arrays xpar1 and xpar2
c
c
        external chunkpnt,interact
        complex *16 cvals(*),cout
        dimension pols(1000)
        dimension xy(2),dxydt(2)
        dimension xynorm(2),xytang(2)
        dimension srcinfo(6),targinfo(6)
        dimension xpar1(20),xpar2(20)
c
c       ... initialize the function evaluator
c
        norder=xpar1(1)
        jchunk=xpar1(3)
c
        do i=1,6
        targinfo(i)=xpar2(i)
        enddo
c
c
ccc        ntimes=ntimes+1
ccc        if( mod(ntimes,1000) .eq. 0 ) call prinf('ntimes=*',ntimes,1)
c
        call legepols(2*t-1,norder,pols)
c
ccc        call prin2('pols=*',pols,npols)
c
        call chunkgeo(chunkpnt,jchunk,t,
     $     par1,par2,par3,par4,
     $     xy,dxydt,ds,xynorm,xytang)
c
ccc        call prin2('xy=*',xy,3)
c
        call chunkinfo(xy,xynorm,xytang,srcinfo)
c
        call interact(srcinfo,targinfo,cout,par5,par6,par7,par8)
c
c        call prin2('srcinfo=*',srcinfo,6)
c        call prin2('ds=*',ds,1)
c        call prin2('targinfo=*',targinfo,6)
c        call prin2('cout=*',cout,1)
c
        do 1200 i=1,npols
        cvals(i)=cout*pols(i)*ds
ccc        cvals(i)=pols(i)*ds
 1200   continue
c
ccc        call prin2('cvals=*',cvals,2*npols)
c
        return
        end
c
c
c
c
c
