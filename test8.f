cc Copyright (C) 2012: Zydrunas Gimbutas 
cc 
cc This software is being released under a modified FreeBSD license
cc (see COPYING in home directory). 
c
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c       
c
c
        implicit real *8 (a-h,o-z)
        dimension xy(2),dxydt(2)
        dimension xynorm(2),xytang(2)
c    
        dimension tsout(100),wsout(100)
        dimension umatr(10000),vmatr(10000)
c
        dimension ichunkinfo(50000),refineinfo(2,50000)
c
        dimension ixys(2,50000),xys(2,50000),xynorms(2,50000)
        dimension xytangs(2,50000)
        dimension whts(50000)
c
        external fchunkpnt,qchunkpnt
        external schunkpnt,echunkpnt
        external rchunkpnt
c
        external lfinter1,lfinter2,lfinter3,lfinter4
        external hfinter1,hfinter2,hfinter3
c
        complex *16 rk
        complex *16 ima
c
        complex *16, allocatable :: cmatr(:,:)
        complex *16, allocatable :: rhs(:)
        complex *16, allocatable :: sol(:)
c
        complex *16, allocatable :: wgmres(:)
        complex *16, allocatable :: wcqr(:)
c
        complex *16 cd
c
        real *8, allocatable :: w(:)
c
        dimension source(2),target(2)
        complex *16 cpot,cpot0
c       
        external h2dmulta,h2dmultb
        dimension errs(10000)
c
        data ima/(0.0d0,1.0d0)/
c
c
c       SET ALL PARAMETERS
c        
        call prini(6,13)
c
        done=1
        pi=4*atan(done)
c
        lw=60 000 000
        allocate(w(lw))
c
cccc        rk=0.0d0
cccc        rk=1.0d0*pi
        rk=.1d0
c
        call prin2('rk=*',rk,2)
c       
c       ... retrieve the interpolation nodes
c
        itype=2
        norder=10
        
        call legeexps(itype,norder,tsout,umatr,vmatr,wsout)
        npols=norder
c
        do i=1,npols
        tsout(i)=(tsout(i)+1)/2
        wsout(i)=wsout(i)/2
        enddo
c
        call prinf('norder=*',norder,1)
        call prinf('npols=*',npols,1)
        call prin2('tsout=*',tsout,npols)
        call prin2('wsout=*',wsout,npols)
c
        d=0
        do 1100 i=1,npols
        d=d+wsout(i)
 1100   continue
c
        call prin2('sum of weights=*',d,1)
c
c
        itest=2
        call prinf('itest=*',itest,1)
c
c
        nchunks=1
c
c       ... refine the parametrization
c
        noversamp=8
        call genrefineinfo(noversamp,nchunks,
     $     nchunksout,ichunkinfo,refineinfo)
c
        nchunks=nchunksout
c
        call prinf('after oversampling, nchunks=*',nchunks,1)
ccc        call prinf('ichunkinfo=*',ichunkinfo,nchunks)
ccc        call prin2('refineinfo=*',refineinfo,2*nchunks)
c       
c
c       ... map the interpolation points on the chunk into R^2
c        
cccc        ichunk=1
c
        npts=npols*nchunks
        call prinf('npts=*',npts,1)
        call prinf('nchunks=*',nchunks,1)
        call prinf('npols=*',npols,1)
c
        allocate( cmatr(npts,npts) )
        allocate( rhs(npts) )
        allocate( sol(npts) )
c
c       ... plot all discretization nodes and normals
c
c       
        call prinf('============================*',i,0)
c
c
c
        do 1400 j=1,nchunks
        do 1200 i=1,npols
c
        ichunk=j
c
        t=tsout(i)
c
        call chunkgeo(rchunkpnt,ichunk,t,
     $     seginfo,schunkpnt,ichunkinfo,refineinfo,
     $     xy,dxydt,ds,xynorm,xytang)

c
cc        call prin2('xy=*',xy,2)
cc        call prin2('dxydt=*',dxydt,2)
c
cc        call prin2('ds=*',ds,1)
c
cc        call prin2('xynorm=*',xynorm,2)
cc        call prin2('xytang=*',xytang,2)
c
c
        write(16,*) xy(1),xy(2)
        write(17,*) xy(1),xy(2)
        write(17,*) 
     $     xy(1)+xynorm(1)/3,
     $     xy(2)+xynorm(2)/3
        write(17,*) 
        write(17,*) 
c
        write(117,*) 
     $     xynorm(1)/3,
     $     xynorm(2)/3
        write(117,*) 
        write(117,*) 
c
        write(18,*) xy(1),xy(2)
        write(18,*) 
     $     xy(1)+xytang(1)/3,
     $     xy(2)+xytang(2)/3
        write(18,*) 
        write(18,*) 
c
        write(118,*) 
     $     xytang(1)/3,
     $     xytang(2)/3
        write(118,*) 
        write(118,*) 
c
 1200   continue
 1400   continue
c
c
c       ... generate all discretization nodes and normals
c
c       
        call prinf('============================*',i,0)
c
        call chunkallpnts(nchunks,rchunkpnt,
     $     seginfo,schunkpnt,ichunkinfo,refineinfo,
     $     npols,tsout,ixys,xys,xynorms,
     $     xytangs,npts)
c
c
        call prinf('npts=*',npts,1)
        call prinf('nchunks=*',nchunks,1)
        call prinf('npols=*',npols,1)
c        call prinf('ixys=*',ixys,2*nchunks)
c        call prin2('xys=*',xys,2*npts)
c        call prin2('xynorms=*',xynorms,2*npts)
c
c
ccc     stop
c
c
c       ... call chunkmatc
c
c       
        call prinf('============================*',i,0)
c
        t1=second()
C$        t1=omp_get_wtime()
c
        call chunkallpnts(nchunks,rchunkpnt,
     $     seginfo,schunkpnt,ichunkinfo,refineinfo,
     $     npols,tsout,ixys,xys,xynorms,
     $     xytangs,npts)
        call chunkallwhts(nchunks,rchunkpnt,
     $     seginfo,schunkpnt,ichunkinfo,refineinfo,
     $     npols,tsout,wsout,whts,npts)
        call chunkmatc(nchunks,rchunkpnt,
     $     seginfo,schunkpnt,ichunkinfo,refineinfo,
     $     norder,npols,tsout,umatr,vmatr,
     $     ixys,xys,xynorms,xytangs,npts,
     $     hfinter1,rk,par6,par7,par8,
     $     cmatr,w,lw,lused,ier)
c
        t2=second()
C$        t2=omp_get_wtime()
c
ccc        call prin2('whts=*',whts,npts)
c
        cd=0
        do i=1,npts
        cd=cd+whts(i)
        enddo
c
        call prin2('sum, whts=*',cd,2)
        call prin2('whts-2 pi=*',cd-2*pi,2)
c
        call prinf('npts=*',npts,1)
ccc        call prin2('after chunkmatc, cmatr=*',cmatr,2*npts*npts)
        call prinf('after chunkmatc, ier=*',ier,1)
        call prin2('in chunkmatc, time=*',t2-t1,1)
c
c        do i=1,npts
c        rhs(i)=1
c        enddo
c
        ifexterior=0
        call prinf('ifexterior=*',ifexterior,1)
c
        if( ifexterior .eq. 0 ) then
c
        target(1)=0.2
        target(2)=-0.1
c
        source(1)=10
        source(2)=20
c
        endif
c
        if( ifexterior .eq. 1 ) then
c
        source(1)=0.2
        source(2)=-0.1
c
        target(1)=10
        target(2)=20
c
        endif
c
c
        call h2getrhs(rk,source,xys,rhs,npts)
c
        call prinf('npts=*',npts,1)
c
c       ... add the diagonal term
c
        call chunkdiag(cmatr,npts,0.0d0)
c
        ifprec = 1
c
        if( ifprec .eq. 1 ) then
        do i=1,npts
        do j=1,npts
        cmatr(i,j)=cmatr(i,j)*sqrt(whts(i))/sqrt(whts(j))
        enddo
        rhs(i)=rhs(i)*sqrt(whts(i))
        enddo
        endif
c
c
        ifsolve = 1
c
        if( ifsolve .eq. 1 ) then
c
        call prinf('entering cqrdecom, npts=*',npts,1)
c
        allocate( wcqr(2*npts**2+npts+10) )
c
        t1=second()
C$        t1=omp_get_wtime()
        call cqrdecom(cmatr,npts,wcqr,rcond)
        t2=second()
C$        t2=omp_get_wtime()
c
        call prin2('after cqrdecom, rcond=*',rcond,1)
        call prin2('in cqrdecom, time=*',t2-t1,1)
c
        t1=second()
C$        t1=omp_get_wtime()
        call cqrsolve(npts,wcqr,rhs,sol)
        t2=second()
C$        t2=omp_get_wtime()
ccc        call prin2('in cqrsolve, time=*',t2-t1,1)
c
        endif
c
c
        if( ifsolve .eq. 2 ) then
c
        call prinf('entering cgmres, npts=*',npts,1)
c
        t1=second()
C$        t1=omp_get_wtime()
        eps=1e-14
        numit=40
        ngmrec=40
c
        allocate( wgmres((ngmrec*2+4)*npts) )
c
        call cgmres(ier,npts,cmatr,
     $     h2dmulta,p1,p2,p3,p4,p5,p6,rhs,eps,numit,
     1     sol,niter,errs,ngmrec,wgmres)
        t2=second()
C$        t2=omp_get_wtime()
c
        call prinf('after cgmres, ier=*',ier,1)
        call prinf('after cgmres, niter=*',niter,1)
        call prin2('after cgmres, errs=*',errs,niter)
        call prin2('in cgmres, time=*',t2-t1,1)
c
        endif
c
c
        if( ifprec .eq. 1 ) then
        do i=1,npts
        sol(i)=sol(i)/sqrt(whts(i))
        rhs(i)=rhs(i)/sqrt(whts(i))
        enddo
        endif
c
c
ccc        call prin2('rhs=*',rhs,2*npts)
ccc        call prin2('sol=*',sol,2*npts)
c
        call prin2('rhs=*',rhs,2*min(npts,30))
        call prin2('sol=*',sol,2*min(npts,30))
c
        cd=0
        do i=1,npts
        cd=cd+sol(i)*whts(i)
        enddo
c
ccc        call prin2('capacitance, integral of the solution=*',cd,2)
ccc        call prin2('capacitance-1=*',cd-1,2)
c
        call h2soleva(rk,target,xys,sol,whts,npts,cpot)
        call prin2('cpot=*',cpot,2)
c
        call h2direva(rk,source,target,cpot0)
        call prin2('directly, cpot0=*',cpot0,2)
c
        call prin2('error=*',(cpot-cpot0)/cpot0,2)
c
        stop
        end
c
c
c
c
c
        subroutine h2getrhs(rk,source,xys,rhs,npts)
        implicit real *8 (a-h,o-z)
        dimension source(2),xys(2,1)
        complex *16 rhs(npts),rk,ima,z,h0,h1
        data ima/(0.0d0,1.0d0)/
c
        do 1200 i=1,npts
        dx=xys(1,i)-source(1)
        dy=xys(2,i)-source(2)
        r=sqrt(dx**2+dy**2)
c
        ifexpon=1
        z=rk*r
        call hank103(z,h0,h1,ifexpon)
c
        rhs(i)=h0
 1200   continue
c        
        return
        end
c
c
c
c
c
        subroutine h2soleva(rk,target,xys,sol,whts,npts,cpot)
        implicit real *8 (a-h,o-z)
        dimension target(2),xys(2,1),whts(1)
        complex *16 sol(npts),cpot,rk,ima,z,h0,h1
        data ima/(0.0d0,1.0d0)/
c
        cpot=0
        do 1200 i=1,npts
        dx=target(1)-xys(1,i)
        dy=target(2)-xys(2,i)
        r=sqrt(dx**2+dy**2)
c
        ifexpon=1
        z=rk*r
        call hank103(z,h0,h1,ifexpon)
c
        cpot=cpot+h0*sol(i)*whts(i)
 1200   continue
c        
        return
        end
c
c
c
c
c
        subroutine h2direva(rk,source,target,cpot)
        implicit real *8 (a-h,o-z)
        dimension source(2),target(2)
        complex *16 cpot,rk,ima,z,h0,h1
        data ima/(0.0d0,1.0d0)/
c
        dx=target(1)-source(1)
        dy=target(2)-source(2)
        r=sqrt(dx**2+dy**2)
c
        ifexpon=1
        z=rk*r
        call hank103(z,h0,h1,ifexpon)
c
        cpot=h0
c        
        return
        end
c
c
c
c
c
        subroutine chunkdiag(cmatr,npts,cd)
        implicit real *8 (a-h,o-z)
        complex *16 cmatr(npts,npts)
c
        do 1200 i=1,npts
        cmatr(i,i)=cmatr(i,i)+cd
 1200   continue
c
        return
        end
c
c
c
c
c
        subroutine h2dmulta(a,p1,p2,p3,p4,p5,p6,x,y,n)
        implicit real *8 (a-h,o-z)
        complex *16 a(n,n),x(n),y(n)
c
        do 1200 i=1,n
        y(i)=0
        do 1100 j=1,n
        y(i)=y(i)+a(i,j)*x(j)
 1100   continue
 1200   continue
c
        return
        end
c
c
c
c
c
        subroutine h2dmultb(a,p1,p2,p3,p4,p5,p6,x,y,n)
        implicit real *8 (a-h,o-z)
        complex *16 a(n,n),x(n),y(n)
c
        do 1200 i=1,n
        y(i)=0
        do 1100 j=1,n
        y(i)=y(i)+conjg(a(j,i))*x(j)
 1100   continue
 1200   continue
c
        return
        end
c
c
c
c
c
