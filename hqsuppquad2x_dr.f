        implicit real *8 (a-h,o-z)
        dimension w(90 000 000),xs0(1000),xs(1000),ws(1000),
     $     rints(1000),rints2(1000),errs(1000),
     $     xsupp(1000),wsupp(1000)
c
        external funuser,funuser2,funwht,fun1,fun2,fun3
c
        call prini(6,13)
C
C       SET ALL PARAMETERS
C
C
        PRINT *, 'ENTER npols, [1-6,8,10,12,16,20]'
        READ *,npols
        CALL PRINf('npols=*',npols,1 )
c
c        initialize the function computer
c
        ifquad = 0
c                      
        ax=0.0d0
        bx=+1.0d0

c
        nsupp = npols
c
        call prinf('nsupp=*',nsupp,1)
c
        ifwhts=1
        call legewhts(nsupp,xsupp,wsupp,ifwhts)
        do i=1,nsupp
        xsupp(i)=(xsupp(i)+1)/2
        wsupp(i)=wsupp(i)/2
        enddo
        call prin2('xsupp=*',xsupp,nsupp)
        call prin2('wsupp=*',wsupp,nsupp)
c
c
        nfuns0=2*npols
c
        ifhilbert=1
        ifhyper=1
c
c       
        do 5000 inode = 1,nsupp
        x07 = xsupp(inode)
        call prin2('x07=*',x07,1)
c
        call funuserini(nsupp,x07)
        call funuser2ini(nsupp,x07)
c
        call hqsuppquad2x(nsupp,inode,xs,ws,nnodes)

        m=20
        if( ifquad .eq. 0 ) eps=1d-15
        if( ifquad .eq. 1 ) eps=1d-15**2
        if( ifquad .eq. 2 ) eps=1d-18
c
        do 3000 i=1,nfuns0
        call adapgauss(ier,ax,bx,fun1,funuser,i,m,eps,
     1      rints(i),maxrec,numint)
 3000   continue
c
        do 3400 i=1,nfuns0
        rints2(i)=0
        do 3200 j=1,nnodes
        call funuser(xs(j),i,val)
        rints2(i)=rints2(i)+ws(j)*val
 3200   continue
        errs(i)=rints2(i)-rints(i)
 3400   continue
c
        call prin2('rints=*',rints,nfuns0)
        call prin2('rints2=*',rints2,nfuns0)
        call prin2('errors=*',errs,nfuns0)
c
        if( ifhilbert .eq. 1 ) then
c
        hint=0
        do 4100 i=1,nnodes
c
        x=xs(i)
        hint=hint+ws(i)/(x-x07)
 4100   continue
c
        hexact=log(bx-x07)-log(x07-ax)
c
        call prin2('via quadr, hilbert integral=*',hint,1)
        call prin2('exact,     hilbert integral=*',hexact,1)
        call prin2('error=*',hint-hexact,1)
        call prin2('rel error=*',(hint-hexact)/hexact,1)
c
        endif
c
c
        if( ifhyper .eq. 1 ) then
c
        qint=0
        do 4200 i=1,nnodes
c
        x=xs(i)
        qint=qint+ws(i)/(x-x07)**2
 4200   continue
c
        qexact=-1/(bx-x07)-1/(x07-ax)
c
        call prin2('via quadr, quadrupole integral=*',qint,1)
        call prin2('exact,     quadrupole integral=*',qexact,1)
        call prin2('error=*',qint-qexact,1)
        call prin2('rel error=*',(qint-qexact)/qexact,1)
c
        endif
c
 5000   continue
c
c        write(*,*) nnodes
c        do i=1,nnodes
c        write(*,*) 2*xs(i)-1,2*ws(i)
c        enddo
        
        stop
        end
c
c
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c      The user-specified functions    
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
c
        subroutine funwht(x,weight)
        implicit real *8 (a-h,o-z)
c
        save
c
        weight=1
c
        return
        end
c
c
c
c
c
        subroutine fun1(x,funuser,i,f)
        implicit real *8 (a-h,o-z)
        external funuser
        call funuser(x,i,f)
        f=f
        return
        end
c
c
c
c
c
c
        subroutine fun2(x,funuser,i,f)
        implicit real *8 (a-h,o-z)
        external funuser
        call funuser(x,i,f)
        f=f**2
        return
        end
c
c
c
c
c
c
        subroutine fun3(x,n,funuser,par2,f)
        implicit real *8 (a-h,o-z)
        dimension f(*)
        external funuser
        do 1200 i=1,n
        call funuser(x,i,f(i))
        f(i)=f(i)**2
 1200   continue
        return
        end
c
c
c
c
c
c
        subroutine funuser(x,i,f)
        implicit real *8 (a-h,o-z)
        save
c
        xx=2*x-1
c
        if(i .gt. npols) goto 1200
c
        call legepol(xx,i-1,f,der)
        f=f 
        return
c
 1200 continue
c
        if(i .gt. npols*2) goto 1400
c
        ii=i-npols
c
        call legepol(xx,ii-1,f,der)
        f=f*log(abs(x-x0)) 
        return
c
 1400 continue
c
        if(i .gt. npols*2+1) goto 1600
c
        h=1.0d-12
        scale=1d-6
c
        h=1.0d-14
        scale=1d-6
c
        f=(x-x0)/((x-x0)**2+h**2)
        f=f*scale
        return
c
 1600 continue
c
        if(i .gt. npols*2+2) goto 1800
c
        h=1.0d-8
        scale=1d-11
c
        f=((x-x0)**2-h**2)/((x-x0)**2+h**2)**2
        f=f*scale
        return
c
 1800 continue
c
c
        call prin2('disaster in funuser!!!*',i,0)
c
        stop
c
        return
c
c
c
c
        entry funuserini(npols7,x07)
c
        npols=npols7
        x0=x07
c
        return
        end
c
c
c
c
c
c
c
        subroutine funuser2(x,i,f)
        implicit real *8 (a-h,o-z)
        save
c
        xx=2*x-1
c
        if(i .gt. npols) goto 1200
c
        call legepol(xx,i-1,f,der)
        f=f 
        return
c
 1200 continue
c
        if(i .gt. npols*2) goto 1400
c
        ii=i-npols
c
        call legepol(xx,ii-1,f,der)
        f=f*log(abs(x-x0)) 
        return
c
 1400 continue
c
        if(i .gt. npols*2+1) goto 1600
c
        f=1/(x-x0)
        return
c
 1600 continue
c
        if(i .gt. npols*2+2) goto 1800
c
        f=1/(x-x0)**2
        return
c
 1800 continue
c
c
        call prin2('disaster in funuser2!!!*',i,0)
c
        stop
c
        return
c
c
c
c
        entry funuser2ini(npols7,x07)
c
        npols=npols7
        x0=x07
c
        return
        end

        
