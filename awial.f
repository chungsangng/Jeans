c----------------------------------------------------------------------2
      program awinu
c
c by C. S. Ng
c
c Version: 3/22/00
c
c calculate matching of the a_n/a_{n-1} ratio
c using the recurrion relation inversely and forwardly
c
      complex a0,a1,a2,a3,ci,sl,lambda,f0,f1
      real lr,li,nu,alpha,amax,eps,fal,dlr,dli
      real lr1,lr2,lr3,li1,li2,li3,wp,pi
      real lx(2),plx(3,2),yx(3),ftol,funk,emin
      real lrold,liold,alold,dlrda,dlida
      integer i,j,k,nmax,nlr,nli,nmin,iter,itmax,ialmax,iw
      common/para/ci,a0,nu,alpha,eps,nmin,nmax,amax
      namelist/in/lr,li,nu,fal,dlr,dli,alpha,nmax,nmin,amax,eps,
     .  ftol,itmax,ialmax,iw
c
      read(*,in)
      write(*,in)      
c
      if ((lr.eq.0.).and.(li.eq.0.)) then
        pi = 4.*atan(1.)
        lr = 0.
        li = (1. + 1./alpha)/sqrt(pi)
      endif 
      write(*,*)"wr = ",lr
      write(*,*)"wi = ",li
      ci = (0.,1.)
      a0 = (1.,0.)
      write(*,5)" alpha                   ",
     .          " wr                      ",
     .          " wi                      ",
     .          " nmax         ",
     .          " nmin         ",
     .          " funk                    "
      i = 0
      lx(1) = lr - dlr
      lx(2) = li - dli
      plx(1,1) = lx(1)
      plx(1,2) = lx(2)
      yx(1) = funk(lx)
      lx(1) = lr - dlr
      lx(2) = li + dli
      plx(2,1) = lx(1)
      plx(2,2) = lx(2)
      yx(2) = funk(lx)
      lx(1) = lr + dlr
      lx(2) = li
      plx(3,1) = lx(1)
      plx(3,2) = lx(2)
      yx(3) = funk(lx)
      call amoeba(plx,yx,3,2,2,ftol,funk,iter,itmax,emin)
      if (mod(i,iw).eq.0) then
        write(*,10)alpha,plx(1,1),plx(1,2),nmax,nmin,yx(1)
      endif
      lr = plx(1,1)
      li = plx(1,2)
      lrold = lr
      liold = li
      alold = alpha
      i = 1
      alpha = alpha*fal**0.01
      lx(1) = lr - dlr
      lx(2) = li - dli
      plx(1,1) = lx(1)
      plx(1,2) = lx(2)
      yx(1) = funk(lx)
      lx(1) = lr - dlr
      lx(2) = li + dli
      plx(2,1) = lx(1)
      plx(2,2) = lx(2)
      yx(2) = funk(lx)
      lx(1) = lr + dlr
      lx(2) = li
      plx(3,1) = lx(1)
      plx(3,2) = lx(2)
      yx(3) = funk(lx)
      call amoeba(plx,yx,3,2,2,ftol,funk,iter,itmax,emin)
      if (mod(i,iw).eq.0) then
        write(*,10)alpha,plx(1,1),plx(1,2),nmax,nmin,yx(1)
      endif
      lr = plx(1,1)
      li = plx(1,2)
      dlrda = (lr - lrold)/(alpha - alold)
      dlida = (li - liold)/(alpha - alold)
      lrold = lr
      liold = li
      do i = 2,ialmax
        alpha = alpha*fal
        lr = lr + dlrda*(alpha - alold)
        li = li + dlida*(alpha - alold)
        lx(1) = lr - dlr
        lx(2) = li - dli
        plx(1,1) = lx(1)
        plx(1,2) = lx(2)
        yx(1) = funk(lx)
        lx(1) = lr - dlr
        lx(2) = li + dli
        plx(2,1) = lx(1)
        plx(2,2) = lx(2)
        yx(2) = funk(lx)
        lx(1) = lr + dlr
        lx(2) = li
        plx(3,1) = lx(1)
        plx(3,2) = lx(2)
        yx(3) = funk(lx)
        call amoeba(plx,yx,3,2,2,ftol,funk,iter,itmax,emin)
        if (mod(i,iw).eq.0) then
          write(*,10)alpha,plx(1,1),plx(1,2),nmax,nmin,yx(1)
        endif
        lr = plx(1,1)
        li = plx(1,2)
        dlrda = (lr - lrold)/(alpha - alold)
        dlida = (li - liold)/(alpha - alold)
        lrold = lr
        liold = li
        alold = alpha
      enddo
      write(*,in)      
5     format(3a25,2a14,a16)
10    format(3(e24.16,x),2(i13,x),e16.8)
      stop
      end
c----------------------------------------------------------------------2
      function funk(x)
c
c Version: 11/5/98
c
      complex a0,a1,a2,a3,ci,sl,omeg,f0,f1
      real lr,li,nu,alpha,amax,eps,zerr
      real x(2),funk
      integer i,nmax,nmin,icheck
      common/para/ci,a0,nu,alpha,eps,nmin,nmax,amax
c
      omeg = cmplx(x(1),x(2))
      sl = a0*sqrt(2.)
1     a3 = a0*eps
      a2 = a3*sqrt(2.*nmax)*nu*ci
      icheck = 0
      nmin = 3
      do i = nmax,nmin,-1
        a1 = (sl*(omeg+ci*i*nu)*a2 - sqrt(1.*(i+1))*a3)
     .           /sqrt(1.*i)
        if (zabs(a1).gt.amax) then
          a2 = a2*eps/amax
          a1 = a1*eps/amax
          icheck = nmax - (i+1)
        endif
        a3 = a2
        a2 = a1
        if (zabs(a2).lt.zabs(a3)) then
          if ((icheck.ne.0).or.((icheck.eq.0).and.
     .        (zabs(a2)+zabs(a3).gt.1.))) then
            nmin = i
            go to 5
          endif 
        endif 
      enddo
5     if (icheck.ne.0) then
        nmax = icheck
      elseif (zabs(a2)+zabs(a3).lt.1.) then
        nmax = nmax*2
        go to 1
      endif
      f1 = -a3/a2
      a1 = sl*omeg*a0
      a2 = (sl*(omeg+ci*nu)*a1 - (1.+alpha)*a0)/sqrt(2.)
      a3 = (sl*(omeg+2.*ci*nu)*a2 - sqrt(2.)*a1)/sqrt(3.)
      a2 = a3/a2
      a1 = a2
      do i = 4,nmin
        a1 = (sl*(omeg+ci*(i-1)*nu) - sqrt(1.*(i-1))/a2)
     .         /sqrt(1.*i)
c        if (zabs(a1).gt.amax) go to 1
        a2 = a1
      enddo
      f0 = -a1
      zerr = zabs(f0-f1)
      if (zerr.lt.1.e-120) then
        funk = -120.
      else
        funk = log(2.*zerr/(zabs(f0)+zabs(f1)))
      endif
      return
      end
c----------------------------------------------------------------------2
      subroutine amoeba(p,y,mp,np,ndim,ftol,funk,iter,itmax)
c
c from Numerical Recipes
c
c Version: 12/15/98
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
      real p(mp,np),y(mp),pr(20),prr(20),pbar(20),emin
      real alpha,beta,gamma,ftol,rtol,ypr,yprr,funk,diff
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

c
c     Compute the fractional range from highest to lowest and return
c  if satisfactory.
c
      diff = 0.
      do i = 1,np
        diff = diff + 2.*abs(p(ihi,i)-p(ilo,i))/(abs(p(ihi,i))+
     .         abs(p(ilo,i)))
      enddo
      if (diff.lt.1.e-15) go to 100
      rtol = 2.*abs(y(ihi)-y(ilo))/(abs(y(ihi))+abs(y(ilo)))
      if (rtol.lt.ftol) go to 100
      if (iter.eq.itmax) then
c        write(*,*) 'Amoba exceeding maximum iterations.'
        return
      endif
      iter = iter + 1
      do j = 1,ndim
        pbar(j) = 0.0
      enddo
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
c
c     Extrapolate by a factor alpha through the face, i.e. reflect the
c  simplex across the ray from the high point.
c
      do j = 1,ndim
        pbar(j) = pbar(j)/ndim
        pr(j) = (1.0 + alpha)*pbar(j) - alpha*p(ihi,j)
      enddo
c
c     Evaluate the function at the reflected point.
c
      ypr = funk(pr)
c
c     Gives a result better than the best point, so try an additional
c 
    2 if (ypr.le.y(ilo)) then
        do j = 1,ndim
          prr(j) = gamma*pr(j) + (1.0 - gamma)*pbar(j)
        enddo
        yprr = funk(prr)
    3   if (yprr.lt.y(ilo)) then
          do j = 1,ndim
            p(ihi,j) = prr(j)
          enddo
          y(ihi) = yprr
    4   else
          do j = 1,ndim
            p(ihi,j) = pr(j)
          enddo
          y(ihi) = ypr
        endif
  5   else if (ypr.ge.y(inhi)) then
  6    if (ypr.lt.y(ihi)) then
          do j = 1,ndim
            p(ihi,j) = pr(j)
          enddo
          y(ihi) = ypr
        endif

       do j = 1,ndim
          prr(j) = beta*p(ihi,j) + (1.0 - beta)*pbar(j)
        enddo
        yprr = funk(prr)
  7     if (yprr.lt.y(ihi)) then
          do j = 1,ndim
            p(ihi,j) = prr(j)
          enddo
          y(ihi) = yprr
  8     else
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
 10   else
        do j= 1,ndim
          p(ihi,j) = pr(j)
        enddo
        y(ihi) = ypr
      endif
      go to 1
100   return 
      end
c----------------------------------------------------------------------2
c        
