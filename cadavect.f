cc Copyright (C) 2012: Zydrunas Gimbutas 
cc 
cc This software is being released under a modified FreeBSD license
cc (see COPYING in home directory). 
c
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c        this is the end of the debugging code and the 
c        start of the actual quadrature routines.
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
c
c
c
        subroutine cadavect(ier,a,b,fun,n,
     $      p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11,p12,
     $      m,eps,cints,maxrec,numint,ww,lenww)
        implicit real *8 (a-h,o-z)
        complex *16 cints(1),ww(1)
        dimension t(100),w(100),stack(400)
        external fun
c
c       this subroutine uses the adaptive Gaussian quadrature
c       to evaluate the integral of the user-supplied vector-valued
c       function fun: [a,b] \to R^n on the user-specified interval [a,b]
c
c                       input parameters:
c
c  a,b - the ends of the interval on which the integral is to 
c       be evaluated
c  fun - the user-supplied SUBROUTINE evaluating the function to 
c       be integrated. the calling sequence of fun must be 
c
c        callfun(x,n,fout,p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11,p12).                            (1)
c
c        in (1), x is a point on the interval [a,b] where
c        the function is to be evaluated, and par1, par2
c        are two parameters to be used by fun; they can be 
c        variables or arrays, real or integer, as desired. 
c        The output vector fout is assumed to be real *8 vector
c        of length n.
c  n - the dimension of the vector-values function fun to be 
c        integrated
c  par1, par2 - parameters to be used by the user-supplied 
c        subroutine fun (see above)
c  m - the order of the quadrature to me used on each subinterval
c  eps - the accuracy (absolute) to which the integral will be 
c       evaluated
c  lenww - the length of the user-supplied work array ww. While
c        n*400 is alway sufficient, the actual requirements are 
c        normally much less. If the array is too short, it will
c        limit the depth of recursion to which the subroutine 
c        can go. When it runs out of memory in ww, it sets ier 
c        to 4 and bombs.
c
c                       output parameters:
c
c  ier - error return code. 
c          ier=0 means normal conclusion
c          ier=4 means that the user-supplied work array ww is too 
c                short. this is a fatal error.
c          ier=8 means that at some point, one subinterval in the
c                subdivision was smaller than (b-a)/2**200. this 
c                is a fatal error.
c          ier=16 means that the total number of subintervals in the
c                adaptive subdivision of [a,b] turned out to be greater 
c                than 100000.  this is a fatal error.
c                
c  cints - the integral as evaluated
c  maxrec - the maximum depth to which the recursion went at its 
c         deepest point. can not be greater than 200, since at that
c         point ier is set to 8 and the execution of the subroutine
c         terminated.
c  numint - the total number of intervals in the subdivision. can not 
c         be greater than 100000,  since at that
c         point ier is set to 16 and the execution of the subroutine
c         terminated.
c
c                        work arrays:
c
c  ww - must be sufficiently long (see parameter lenww above).         
c
c          . . . check if the legendre quadrature has been
c                initialized at a preceeding call to this routine
c

        ifinit=1
        call legewhts(m,t,w,ifinit)

c
c        integrate the user-supplied function using the 
c        adaptive gaussian quadratures
c
        nnmax=100000
        maxdepth=200
c
        ivalue2=1
        lvalue2=n*2+4
c
        ivalue3=ivalue2+lvalue2
        lvalue3=n*2+4
c
        ifout=ivalue3+lvalue3
        lfout=n*2+4
c
        ival=ifout+lfout
c
        lleft=lenww-ival
c
        mmm=(lleft-100)/(n*2)
cccc        call prinf('mmm as calculated*',mmm,1)
        if(maxdepth .lt. mmm) maxdepth=mmm
c        
        call cadvecre(ier,stack,a,b,fun,n,
     $      p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11,p12,t,w,m,
     1      ww(ival),nnmax,eps,cints,maxdepth,maxrec,numint,
     2      ww(ivalue2),ww(ivalue3),ww(ifout) )
        return
        end
c
c
c
c
c
        subroutine cadvecre(ier,stack,a,b,fun,n,
     1      p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11,p12,
     $      t,w,m,vals,nnmax,eps,
     2      cints,maxdepth,maxrec,numint,value2,value3,fout)
        implicit real *8 (a-h,o-z)
        complex *16 cints(1),value2(1),value3(1),vals(n,1),fout(1)
        dimension stack(2,1),t(1),w(1)
        external fun
c
c       start the recursion
c
        stack(1,1)=a
        stack(2,1)=b
        call convecin(a,b,fun,n,
     $     p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11,p12,t,w,m,vals(1,1),fout)
c
c       recursively integrate the thing
c
        j=1
        do 1200 jjj=1,n
c        
        cints(jjj)=0
 1200 continue
c
        ier=0
        maxrec=0
        do 3000 i=1,nnmax
ccc        call prinf('i=*',i,1)
        numint=i
        if(j .gt. maxrec) maxrec=j
ccc        call prinf('j=*',j,1)
c
c       subdivide the current subinterval
c
        c=(stack(1,j)+stack(2,j))/2
        call convecin(stack(1,j),c,fun,n,
     1      p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11,p12,t,w,m,value2,fout)
c
        call convecin(c,stack(2,j),fun,n,
     1      p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11,p12,t,w,m,value3,fout)
c
c
        ifdone=1
c
c       ... estimate the maximum magnitude of integrals to be evaluated,
c       currently, this is done for each subdivision separately, which
c       may cause some difficulties if all functions are very small
c       inside a subregion. This should not be a big problem for
c       non-oscillatory kernels.
c
        ifrel=0
c
        dmax = 0
        do ij=1,n
        if( dmax .lt. abs(vals(ij,j)) ) dmax = abs(vals(ij,j))
        enddo
c
        do 1400 ij=1,n
c
c       ... check if absolute precision has been reached
c
        if( ifrel .eq. 0 ) then
        if( abs(value2(ij)+value3(ij)
     1      -vals(ij,j)) .gt. eps) ifdone=0
        endif
c
c       ... check if relative precision has been reached
c
        if( ifrel .eq. 1 ) then
        if( abs(value2(ij)+value3(ij)
     1      -vals(ij,j)) .gt. eps*dmax) ifdone=0
        endif
c
 1400 continue

c
c       if the function on this subinterval has been 
c       integrated with sufficient accuracy - add the 
c       value to that of the global integral and move up
c       in the stack
c
        if(ifdone  .eq. 0) goto 2000
c
        do 1600 jjj=1,n
c
        cints(jjj)=cints(jjj)+value2(jjj)+value3(jjj)
 1600 continue
        j=j-1
c
c        if the whole thing has been integrated - return
c
        if(j .eq. 0) return
        goto 3000
 2000 continue
c        
c       if the function on this subinterval has not been 
c       integrated with sufficient accuracy - move 
c       down the stack
c
        stack(1,j+1)=stack(1,j)
        stack(2,j+1)=(stack(1,j)+stack(2,j))/2
c
        do 2200 jjj=1,n
        vals(jjj,j+1)=value2(jjj)
 2200 continue
c
        stack(1,j)=(stack(1,j)+stack(2,j))/2
c
        do 2400 jjj=1,n
        vals(jjj,j)=value3(jjj)
 2400 continue
c
cccc        vals(j)=value3
c
        j=j+1
c     
c       if the depth of the recursion has become excessive - bomb
c
        if(j+1 .le. maxdepth) goto 3000
c
        ier=8
        if(maxdepth .lt. 200) ier=4
        return
 3000 continue
        ier=16
        return
        end
c
c
c
c
c
        subroutine convecin(a,b,fun,n,
     $      p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11,p12,
     $      t,w,m,cints,fout)
        implicit real *8 (a-h,o-z)
        dimension t(1),w(1)
        complex *16 cints(1),fout(1)
        external fun
c
c       integrate the function fun on the interval [a,b]
c 
        do 1100 i=1,n
c    
        cints(i)=0
 1100 continue
c
        u=(b-a)/2
        v=(b+a)/2
        do 1200 i=1,m
        tt=u*t(i)+v
c
        call fun(tt,n,fout,p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11,p12)
c
        do 1150 j=1,n
c
        cints(j)=cints(j)+fout(j)*w(i)
 1150 continue
c
 1200 continue
c
        do 1400 i=1,n
c
        cints(i)=cints(i)*u
 1400 continue
        return
        end
c
c
c
c
c
c
c
c
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c        this is the end of the debugging code and the 
c        start of the actual quadrature routines.
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
c
c
c
        subroutine cadachunk(ier,a,b,fun,n,
     $      p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11,p12,
     $      m,eps,cints,maxrec,numint,ww,lenww,
     $      ichunk,targinfo,info)
        implicit real *8 (a-h,o-z)
        complex *16 cints(1),ww(1)
        dimension t(100),w(100),stack(400),targinfo(6),info(2)
        external fun, p1
c
c       this subroutine uses the adaptive Gaussian quadrature
c       to evaluate the integral of the user-supplied vector-valued
c       function fun: [a,b] \to R^n on the user-specified interval [a,b]
c
c                       input parameters:
c
c  a,b - the ends of the interval on which the integral is to 
c       be evaluated
c  fun - the user-supplied SUBROUTINE evaluating the function to 
c       be integrated. the calling sequence of fun must be 
c
c        callfun(x,n,fout,p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11,p12).                            (1)
c
c        in (1), x is a point on the interval [a,b] where
c        the function is to be evaluated, and par1, par2
c        are two parameters to be used by fun; they can be 
c        variables or arrays, real or integer, as desired. 
c        The output vector fout is assumed to be real *8 vector
c        of length n.
c  n - the dimension of the vector-values function fun to be 
c        integrated
c  par1, par2 - parameters to be used by the user-supplied 
c        subroutine fun (see above)
c  m - the order of the quadrature to me used on each subinterval
c  eps - the accuracy (absolute) to which the integral will be 
c       evaluated
c  lenww - the length of the user-supplied work array ww. While
c        n*400 is alway sufficient, the actual requirements are 
c        normally much less. If the array is too short, it will
c        limit the depth of recursion to which the subroutine 
c        can go. When it runs out of memory in ww, it sets ier 
c        to 4 and bombs.
c
c                       output parameters:
c
c  ier - error return code. 
c          ier=0 means normal conclusion
c          ier=4 means that the user-supplied work array ww is too 
c                short. this is a fatal error.
c          ier=8 means that at some point, one subinterval in the
c                subdivision was smaller than (b-a)/2**200. this 
c                is a fatal error.
c          ier=16 means that the total number of subintervals in the
c                adaptive subdivision of [a,b] turned out to be greater 
c                than 100000.  this is a fatal error.
c                
c  cints - the integral as evaluated
c  maxrec - the maximum depth to which the recursion went at its 
c         deepest point. can not be greater than 200, since at that
c         point ier is set to 8 and the execution of the subroutine
c         terminated.
c  numint - the total number of intervals in the subdivision. can not 
c         be greater than 100000,  since at that
c         point ier is set to 16 and the execution of the subroutine
c         terminated.
c
c                        work arrays:
c
c  ww - must be sufficiently long (see parameter lenww above).         
c
c          . . . check if the legendre quadrature has been
c                initialized at a preceeding call to this routine
c

        ifinit=1
        call legewhts(m,t,w,ifinit)

c
c        integrate the user-supplied function using the 
c        adaptive gaussian quadratures
c
        nnmax=100000
        maxdepth=200
c
        ivalue2=1
        lvalue2=n*2+4
c
        ivalue3=ivalue2+lvalue2
        lvalue3=n*2+4
c
        ifout=ivalue3+lvalue3
        lfout=n*2+4
c
        ival=ifout+lfout
c
        lleft=lenww-ival
c
        mmm=(lleft-100)/(n*2)
cccc        call prinf('mmm as calculated*',mmm,1)
        if(maxdepth .lt. mmm) maxdepth=mmm
c        
        call cadpatchre(ier,stack,a,b,fun,n,
     $      p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11,p12,t,w,m,
     1      ww(ival),nnmax,eps,cints,maxdepth,maxrec,numint,
     2      ww(ivalue2),ww(ivalue3),ww(ifout),
     $      ichunk,targinfo,info)
        return
        end
c
c
c
c
c
        subroutine cadpatchre(ier,stack,a,b,fun,n,
     1      p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11,p12,
     $      t,w,m,vals,nnmax,eps,
     2      cints,maxdepth,maxrec,numint,value2,value3,fout,
     $      ichunk,targinfo,info)
        implicit real *8 (a-h,o-z)
        complex *16 cints(1),value2(1),value3(1),vals(n,1),fout(1)
        dimension stack(2,1),t(1),w(1),targinfo(6),info(2)
        dimension xy0(2),xy1(2),xy2(2)
        dimension dxydt0(2),dxydt1(2),dxydt2(2)
c
        dimension irec(207)
c
        external fun, p1
c
c       start the recursion
c
        info(1)=0
        info(2)=0
c
        stack(1,1)=a
        stack(2,1)=b
        call conpatchin(a,b,fun,n,
     $     p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11,p12,t,w,m,vals(1,1),fout)
c
c       recursively integrate the thing
c
        j=1
        do 1200 jjj=1,n
c        
        cints(jjj)=0
 1200 continue
c
        irec(j)=1
c
        ier=0
        maxrec=0
        do 3000 i=1,nnmax
ccc        call prinf('i=*',i,1)
        numint=i
        if(j .gt. maxrec) maxrec=j
c
c
c       first, check if the target point is in the far field 
c       and act accordingly
c
        iffar=0
c
        t1=stack(1,j)
        call p1(ichunk,t1,xy1,dxydt1,p2,p3,p4,p5)
        t2=stack(2,j)
        call p1(ichunk,t2,xy2,dxydt2,p2,p3,p4,p5)
        t0=(stack(1,j)+stack(2,j))/2
        call p1(ichunk,t0,xy0,dxydt0,p2,p3,p4,p5)
        dx=targinfo(1)-xy0(1)
        dy=targinfo(2)-xy0(2)
        r0=sqrt(dx*dx+dy*dy)
        dx=xy1(1)-xy0(1)
        dy=xy1(2)-xy0(2)
        r1=sqrt(dx*dx+dy*dy)
        dx=xy2(1)-xy0(1)
        dy=xy2(2)-xy0(2)
        r2=sqrt(dx*dx+dy*dy)
c
        rmax=r1
        if( r2 .gt. rmax ) rmax=r2
c
cc        dx=xy1(1)-xy2(1)
cc        dy=xy1(2)-xy2(2)
cc        rmax=sqrt(dx*dx+dy*dy)

        if( r0 .gt. 3*rmax) iffar=1 
c
ccc        write(*,*) r0, rmax, (r0 .gt. 3*rmax), iffar, maxrec
c
        if( iffar .eq. 1 ) then        
ccc        call prinf('iffar=*',iffar,1)
        endif

        if( irec(j) .ge. maxdepth ) then 
        ier=4 
        iffar=1
        endif

ccc        call prinf('j=*',j,1)

        if( iffar .eq. 0 ) then
c
c       subdivide the current subinterval
c
        c=(stack(1,j)+stack(2,j))/2
        call conpatchin(stack(1,j),c,fun,n,
     1      p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11,p12,t,w,m,value2,fout)
c
        call conpatchin(c,stack(2,j),fun,n,
     1      p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11,p12,t,w,m,value3,fout)
c
c
        ifdone=1
c
c       ... estimate the maximum magnitude of integrals to be evaluated,
c       currently, this is done for each subdivision separately, which
c       may cause some difficulties if all functions are very small
c       inside a subregion. This should not be a big problem for
c       non-oscillatory kernels.
c
        ifrel=1
c
        dmax = 0
        do ij=1,n
        if( dmax .lt. abs(vals(ij,j)) ) dmax = abs(vals(ij,j))
        enddo
c
        do 1400 ij=1,n
c
c       ... check if absolute precision has been reached
c
        if( ifrel .eq. 0 ) then
        if( abs(value2(ij)+value3(ij)
     1      -vals(ij,j)) .gt. eps) ifdone=0
        endif
c
c       ... check if relative precision has been reached
c
        if( ifrel .eq. 1 ) then
        if( abs(value2(ij)+value3(ij)
     1      -vals(ij,j)) .gt. eps*dmax) ifdone=0
        endif
c
 1400 continue

c
c       if the function on this subinterval has been 
c       integrated with sufficient accuracy - add the 
c       value to that of the global integral and move up
c       in the stack
c
        if(ifdone  .eq. 0) goto 2000
c
        do 1600 jjj=1,n
c
        cints(jjj)=cints(jjj)+value2(jjj)+value3(jjj)
 1600 continue
        j=j-1
c
        info(1)=info(1)+1
c
        endif
c

        if( iffar .eq. 1 ) then
c
c       ... target is in the far field, 
c       add the contribution due to the smooth quadrature
c
        do 1900 ij=1,n
        cints(ij)=cints(ij)+vals(ij,j)
 1900 continue
c
        j=j-1
        ifdone=1
c
        info(2)=info(2)+1
c
        endif
c
c        if the whole thing has been integrated - return
c
        if(j .eq. 0) return
        goto 3000
 2000 continue
c        
c       if the function on this subinterval has not been 
c       integrated with sufficient accuracy - move 
c       down the stack
c
        stack(1,j+1)=stack(1,j)
        stack(2,j+1)=(stack(1,j)+stack(2,j))/2
c
        do 2200 jjj=1,n
        vals(jjj,j+1)=value2(jjj)
 2200 continue
c
        stack(1,j)=(stack(1,j)+stack(2,j))/2
        stack(2,j)=stack(2,j)
c
        do 2400 jjj=1,n
        vals(jjj,j)=value3(jjj)
 2400 continue
c
cccc        vals(j)=value3
c
        irectmp=irec(j)+1
        irec(j)=irectmp
        irec(j+1)=irectmp
c
        j=j+1
c     
c       if the depth of the recursion has become excessive - bomb
c
        if(j+1 .le. maxdepth) goto 3000
c
        ier=8
        if(maxdepth .lt. 200) ier=4
        return
 3000 continue
        ier=16
        return
        end
c
c
c
c
c
        subroutine conpatchin(a,b,fun,n,
     $      p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11,p12,
     $      t,w,m,cints,fout)
        implicit real *8 (a-h,o-z)
        dimension t(1),w(1)
        complex *16 cints(1),fout(1)
        external fun, p1
c
c       integrate the function fun on the interval [a,b]
c 
        do 1100 i=1,n
c    
        cints(i)=0
 1100 continue
c
        u=(b-a)/2
        v=(b+a)/2
        do 1200 i=1,m
        tt=u*t(i)+v
c
        call fun(tt,n,fout,p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11,p12)
c
        do 1150 j=1,n
c
        cints(j)=cints(j)+fout(j)*w(i)
 1150 continue
c
 1200 continue
c
        do 1400 i=1,n
c
        cints(i)=cints(i)*u
 1400 continue
        return
        end
c
c
c
c
c
