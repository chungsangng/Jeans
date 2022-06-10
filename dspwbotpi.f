c----------------------------------------------------------------------2
      program dspwbotpi
c
c by C. S. Ng
c
c Version: 1/25/01
c
c calculate the root of the dispersion relation by using
c the plasma dispersion function  -- for Jeans instability
c with two maxwellians (bump on tail)
c
c
      integer mmax
      parameter (mmax=10000)
      complex pdf
      real lr,li,alpha,dlr,dli,c1,c2,c3,p1,p2,eps,xi,xi3,lid
      real lrmin,lrmax,limin,limax,errtol,evtol
      real lx(2),plx(3,2),yx(3),ftol,funk,ev(mmax,2)
      integer i,j,k,nlr,nli,iter,itmax,imode
      common/para/alpha,p1,p2,eps,xi,xi3
      namelist/in/lrmin,lrmax,limin,limax,nlr,nli,alpha,c1,c2,c3,
     .  ftol,errtol,evtol,itmax,eps,xi,xi3
c
      read(*,in)
      write(*,in)
c
      dlr = (lrmax-lrmin)/nlr
      dli = (limax-limin)/nli
      p1 = alpha/(1. + eps*xi*xi)
      p2 = p1*eps*xi**4
      read(*,50)imode
      do i = 1,4
        read(*,*)
      enddo
c      imode = 0
      do i = 1,imode
        read(*,10) ev(i,1),ev(i,2)
      enddo
      do j = 2,nli-1
        li = limin + j*dli
        lid = li-limin
        do k = 0,nlr-1
          lr = lrmin + k*dlr + c1*lid + c2*lid*lid
          lx(1) = lr - dlr*c3
          lx(2) = li - dli*c3
          plx(1,1) = lx(1)
          plx(1,2) = lx(2)
          yx(1) = funk(lx)
          lx(1) = lr - dlr*c3
          lx(2) = li + dli*c3
          plx(2,1) = lx(1)
          plx(2,2) = lx(2)
          yx(2) = funk(lx)
          lx(1) = lr + dlr*c3
          lx(2) = li
          plx(3,1) = lx(1)
          plx(3,2) = lx(2)
          yx(3) = funk(lx)
c        write(*,*)"step 1"
          call amoeba(plx,yx,3,2,2,ftol,funk,iter,itmax)
c        write(*,*)"step 3"
          lr = plx(1,1)
          li = plx(1,2)
          if ((lr.ge.-0.1e-5).and.(yx(1).lt.errtol)) then
            do i = 1,imode
              if ((ev(i,1)-lr)**2+(ev(i,2)-li)**2.lt.evtol
     .          *(lr**2+li**2)) go to 1
            enddo
            imode = imode + 1
            if (imode.gt.mmax) go to 20
            ev(imode,1) = lr
            ev(imode,2) = li
            write(*,10) lr,li
1         endif
        enddo
      enddo
10    format(e24.16,x,e24.16)
20    write(*,*)"imode = ",imode
      write(*,*)alpha
      write(*,*)eps
      write(*,*)xi
      write(*,*)xi3
      do i = 1,imode
        write(*,10) ev(i,1),ev(i,2)
      enddo
30    format(x,i4)
40    format(i4)
50    format(8x,i20)
      stop
      end
c----------------------------------------------------------------------2
      function funk(x)
c
c Version: 1/25/01
c
      complex k,pdf,z,k1,z1
      real alpha,zerr,p1,p2,eps,xi,xi3
      real x(2),funk
      integer i,nmax,nmin,icheck
      common/para/alpha,p1,p2,eps,xi,xi3
c
      k = cmplx(x(1),x(2))
      k1 = xi*((1.,0.)/k-(1.,0.)*xi3)
      z = pdf((1.,0.)/k)
      z1 = pdf(k1)
      zerr = zabs((1.,0.) + p1*((1.,0.)+z/k)
     .         + p2*((1.,0.) + k1*z1) )
      if (zerr.eq.0.) then
        funk = -120.
      else
        funk = log(zerr)
      endif
      return
      end
c----------------------------------------------------------------------2
      subroutine amoeba(p,y,mp,np,ndim,ftol,funk,iter,itmax)
c
c from Numerical Recipes
c
c Version: 11/4/98
c
c     Multidimensional minimization of the function funk(x) where x is
c  an ndim-dimensional vector, by the downhill simplex method of Neder
c  and Mead. Input is a matrix p whose ndim+1 rows are ndim-dimensional
c  vectors which are the vertices of the starting simplex. [Logical
c  dimensions of p are p(ndim+1,ndim); physical dimensions are input as
c  p(mp,np)]. Also input is the vector y of length ndim+1, whose compon-
c  ents must be pre-initialized to the values of funk evaluated at the 
c  ndim+1 vertices (rows) of p; and ftol the fractional convergence 
c  tolerance to be archieved in the function value (n.b.!). On output,
c  p and y will have been reset to ndim+1 new points all within ftol 
c  of a minimum function value, and iter gives the number of iterations 
c  taken.
c
c      parameter (nmax=20,alpha=1.0,beta=0.5,gamma=2.0,itmax=1000)
c
c     Expected maximum number of dimensions, three parameters which 
c  difine the expansions and contractions, and maximum allowed 
c  number of iterations.
c
      integer nmax,itmax,iter,ilo,ihi,inhi,i,j,mp,np,mpts,ndim
      real p(mp,np),y(mp),pr(20),prr(20),pbar(20)
      real alpha,beta,gamma,ftol,rtol,ypr,yprr,funk,diff
      real tmp1
      nmax = 20
      alpha = 1.0
      beta = 0.5
      gamma = 2.0
c      itmax = 1000
      mpts = ndim + 1
c
c     Note that mp is the physical dimension corresponding to the
c  logical dimension mpts, np to ndim.
c
      iter = 0
c
c     First we must determine which point is the highest (worst), 
c  next-highest, and lowest (best),
c
    1 ilo = 1
      if (y(1).gt.y(2)) then
        ihi = 1
        inhi = 2
      else
        ihi = 2
        inhi = 1
      endif
c        write(*,*)"step a"
   15 format(3(f8.5,x,f8.5,2x))
c
c     by looping over the points in the simplex.
c
      do i = 1,mpts
        if (y(i).lt.y(ilo)) ilo = i
        if (y(i).gt.y(ihi)) then
          inhi = ihi
          ihi = i
        else if (y(i).gt.y(inhi)) then
               if (i.ne.ihi) inhi = i
        endif
      enddo
c        write(*,*)"step b"

c
c     Compute the fractional range from highest to lowest and return
c  if satisfactory.
c
      diff = 0.
      do i = 1,np
        diff = diff + 2.*abs(p(ihi,i)-p(ilo,i))/(abs(p(ihi,i))+
     .         abs(p(ilo,i)))
      enddo
c        write(*,*)"step c"
      if (diff.lt.1.e-15) return
      tmp1 = abs(y(ihi))+abs(y(ilo))
      if (tmp1.ne.0.) then
        rtol = 2.*abs(y(ihi)-y(ilo))/(abs(y(ihi))+abs(y(ilo)))
      else
        rtol = abs(y(ihi)-y(ilo))
      endif
c        write(*,*)"step d"
      if (rtol.lt.ftol) return
      if (iter.eq.itmax) then
        write(*,*) 'Amoba exceeding maximum iterations.'
        return
      endif
c        write(*,*)"step e"
      iter = iter + 1
      do j = 1,ndim
        pbar(j) = 0.0
      enddo
c        write(*,*)"step f"
c
c     Begin a new iteration. Compute the vector average of all points
c  except the highest, i.e. the center of the "face" of the simplex 
c  across from the high point. We will subsequently explore along the
c  ray from the high point through that center.
c
      do i = 1,mpts
        if (i.ne.ihi) then
          do j = 1,ndim
            pbar(j) = pbar(j) + p(i,j)
          enddo
        endif
      enddo
c        write(*,*)"step g"
c
c     Extrapolate by a factor alpha through the face, i.e. reflect the
c  simplex across the ray from the high point.
c
      do j = 1,ndim
        pbar(j) = pbar(j)/ndim
        pr(j) = (1.0 + alpha)*pbar(j) - alpha*p(ihi,j)
      enddo
c        write(*,*)"step h"
c
c     Evaluate the function at the reflected point.
c
      ypr = funk(pr)
c
c     Gives a result better than the best point, so try an additional
c 
c        write(*,*)"step i"
    2 if (ypr.le.y(ilo)) then
        do j = 1,ndim
          prr(j) = gamma*pr(j) + (1.0 - gamma)*pbar(j)
        enddo
        yprr = funk(prr)
c        write(*,*)"step j"
    3   if (yprr.lt.y(ilo)) then
          do j = 1,ndim
            p(ihi,j) = prr(j)
          enddo
          y(ihi) = yprr
c        write(*,*)"step k"
    4   else
          do j = 1,ndim
            p(ihi,j) = pr(j)
          enddo
          y(ihi) = ypr
        endif
c        write(*,*)"step l"
  5   else if (ypr.ge.y(inhi)) then
  6    if (ypr.lt.y(ihi)) then
          do j = 1,ndim
            p(ihi,j) = pr(j)
          enddo
          y(ihi) = ypr
        endif
c        write(*,*)"step m"

       do j = 1,ndim
          prr(j) = beta*p(ihi,j) + (1.0 - beta)*pbar(j)
        enddo
        yprr = funk(prr)
c        write(*,*)"step n"
  7     if (yprr.lt.y(ihi)) then
          do j = 1,ndim
            p(ihi,j) = prr(j)
          enddo
          y(ihi) = yprr
  8     else
c        write(*,*)"step o"
          do i = 1,mpts
  9         if (i.ne.ilo) then
              do j = 1,ndim
                pr(j) = 0.5*(p(i,j) + p(ilo,j))
                p(i,j) = pr(j)
              enddo
              y(i) = funk(pr)
            endif
          enddo
        endif
c        write(*,*)"step p"
 10   else
        do j= 1,ndim
          p(ihi,j) = pr(j)
        enddo
        y(ihi) = ypr
      endif
c        write(*,*)"step q"
      go to 1
      end
c----------------------------------------------------------------------2
c        
C=============================================================
      complex function pdf(zeta)
C=============================================================
C
C  PLASMA DISPERSION FUNCTION
C
      complex z,zeta,a(25),b(25),acap(25),bcap(25),zn,zn1,z0,
     &rpim,zsq,p,tst
      real rd,d,rn,c
      real ax,x,y,ay,ab,sig,cx,cy,acc
      integer n,n1,n2
      data rpim/(0.,1.77245385090552)/
      acc = 1.e-14
      cx=1.
      cy=1.
      sig=-1.
      z0=zeta
      ab=cabs(zeta)
      x=zeta
      ax=abs(x)
      y=aimag(zeta)
      ay=abs(y)
      if(ab.gt.4.) go to 20
   10 if((ax.le.4.).and.((ay+ax*.2).lt.2.0)) go to 60
      go to 100
C.. ASYMPTOTIC SERIES
   20 zn=1.
      rn=-.5
      pdf=1.
      zsq=zeta*zeta
   30 rn=rn+1.
      zn1=rn/zsq
      c=cabs(zn1)
      if(c-1.)40,40,50
   40 zn=zn*zn1
      pdf=pdf+zn
      c=cabs(zn/pdf)
      if(c-acc)50,50,30
   50 pdf=-pdf/zeta
      go to 140
C... POWER SERIES
   60 z=cmplx(1.,0.)
      rd=.5
      p=z
      zsq=zeta*zeta
      tst=cmplx(1.,0.)
      if(ab-.7071)80,80,70
   70 tst=2.*zsq
   80 rd=rd+1.
      zn=-zsq/rd
      z=zn*z
      p=p+z
      cx=cabs(z)*cabs(tst)
      if(cx-acc) 90,90,80
   90 pdf=rpim*cexp(-zsq)-2.*zeta*p
      return
C.. CONTINUED FRACTION METHOD (Y SET POSITIVE)
  100 z=cmplx(x,ay)
      a(1)=z
      a(2)=-.5
      b(1)=-z*z+.5
      b(2)=b(1)+2.
      acap(1)=a(1)
      acap(2)=a(1)*b(2)
      bcap(1)=b(1)
      bcap(2)=a(2)+b(1)*b(2)
      zn=acap(2)/bcap(2)
      acc=acc*100.
      do 120 n=3,25
      n1=n-1
      n2=n-2
      a(n)=-n1*(n1-.5)
      b(n)=b(n1)+2.
      acap(n)=b(n)*acap(n1)+a(n)*acap(n2)
      bcap(n)=b(n)*bcap(n1)+a(n)*bcap(n2)
      pdf=acap(n)/bcap(n)
      if(cabs(pdf/zn-1.)-acc) 130,130,110
  110 zn=pdf
  120 continue
c      write(6,220) zeta
C.. REFLECTION FOR X AND/OR Y NEGATIVE
  130 if(y.gt.0.) go to 140
      pdf=conjg(pdf)
C.. DETERMINE SIGN OF SIGMA
  140 if(ay-ax) 150,170,170
  150 if(ay-.78539/ax) 160,170,170
  160 sig=1.
      go to 190
  170 if(y) 180,200,200
  180 sig=2.
  190 d=-zeta*zeta
      if(d.le.100.) go to 191
c      write(6,210)d,zeta
      pdf=pdf+rpim*sig*1.e-20
      go to 200
  191 if(100.+d)200,200,195
  195 pdf=pdf+rpim*sig*cexp(-zeta*zeta)
  200 continue
      return
  210 format("error in pdf, real (-zeta**2)=",e10.3," zeta=",2e10.3)
  220 format("error in pdf, cont. frac. did not converge. zeta=",2e10.3)
      end
c----------------------------------------------------------------------2
