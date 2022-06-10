c----------------------------------------------------------------------2
      program m2botp
c
c by C. S. Ng
c
c Version: 5/2/00
c
c plot the two maxwellians distribution function
c
      integer nsmax,nvmax
      parameter (nsmax=6553600,nvmax=5000)
c
      integer nmax,nv
      integer i,i1,j,jp1,jm1,k
      real b(nsmax)
      real h0(-nvmax:nvmax),h1(-nvmax:nvmax)
      real f(-nvmax:nvmax),df(-nvmax:nvmax),df1(-nvmax:nvmax)
      real f1(-nvmax:nvmax),f0,df0
      real dv,sq2,aos,v1
      real pis,piq,v,sj,sjp1,sjm1,evv(-nvmax:nvmax),evv1,v0
      real alpha,eps,xi,xi3,p1,p2,p3,p4,p5
c
      complex lambda,kk,z,z1
      namelist/in/nmax,nv,dv,v0,alpha,eps,xi,xi3

      read(*,in)
      write(*,in)
      write(*,*)nmax      
      write(*,*)nv      
      write(*,*)dv      
      write(*,*)v0      
      write(*,*)alpha
      write(*,*)eps
      write(*,*)xi
      write(*,*)xi3
c
      pis = sqrt(4.*atan(1.))
c***
c***        write(*,*)pis
c***
      piq = sqrt(pis)
c***
c***        write(*,*)piq
c***
      sq2 = sqrt(2.)
      aos = alpha/sq2/piq
      p1 = 1./(1. + eps*xi*xi)/pis
      p2 = eps*xi**3
      p3 = -alpha*p1
      p4 = xi*xi*p2
      p5 = piq*dv/3.
c***
c***        write(*,*)p1
c***        write(*,*)p2
c***        write(*,*)p3
c***        write(*,*)p4
c***        write(*,*)p5
c***
      i = -nv
      v = v0 + i*dv
      v1 = v - xi3
      evv(i) = exp(-v*v)
      evv1 = exp(-xi*xi*v1*v1)
c***
c***        write(*,*)evv1
c***
      h0(i) = 1./piq
      h1(i) = sq2*v*h0(i)
      f(i) = p1*(evv(i) + p2*evv1)
c***
c***        write(*,*)f(i)
c***
c***
c***        write(*,*)v*evv(i)
c***        write(*,*)v1**evv1
c***        write(*,*)v1**evv1
c***        write(*,*)p4*v1**evv1
c***
      df(i) = p3*(v*evv(i) + p4*v1*evv1)
c***
c***        write(*,*)df(i)
c***
      b(1) = f(i)*h1(i)
c***
c***        write(*,*)i
c***        write(*,*)evv(i)
c***        write(*,*)evv1
c***        write(*,*)f(i)
c***        write(*,*)df(i)
c***
      do i = -nv+1,nv-3,2
        v = v0 + i*dv
        v1 = v - xi3
        evv(i) = exp(-v*v)
        evv1 = exp(-xi*xi*v1*v1)
        h0(i) = 1./piq
        h1(i) = sq2*v*h0(i)
        f(i) = p1*(evv(i) + p2*evv1)
        df(i) = p3*(v*evv(i) + p4*v1*evv1)
        b(1) = b(1) + 4.*f(i)*h1(i) 
c***
c***        write(*,*)i
c***        write(*,*)evv(i)
c***        write(*,*)evv1
c***        write(*,*)f(i)
c***        write(*,*)df(i)
c***
        i1 = i + 1
        v = v + dv
        v1 = v - xi3
        evv(i1) = exp(-v*v)
        evv1 = exp(-xi*xi*v1*v1)
        h0(i1) = 1./piq
        h1(i1) = sq2*v*h0(i1)
        f(i1) = p1*(evv(i1) + p2*evv1)
        df(i1) = p3*(v*evv(i1) + p4*v1*evv1)
        b(1) = b(1) + 2.*f(i1)*h1(i1) 
c***
c***        write(*,*)i1
c***        write(*,*)evv(i1)
c***        write(*,*)evv1
c***        write(*,*)f(i1)
c***        write(*,*)df(i1)
c***
      enddo
      i = nv-1
      v = v0 + i*dv
      v1 = v - xi3
      evv(i) = exp(-v*v)
      evv1 = exp(-xi*xi*v1*v1)
      h0(i) = 1./piq
      h1(i) = sq2*v*h0(i)
      f(i) = p1*(evv(i) + p2*evv1)
      df(i) = p3*(v*evv(i) + p4*v1*evv1)
      b(1) = b(1) + 4.*f(i)*h1(i)
c***
c***        write(*,*)i
c***        write(*,*)evv(i)
c***        write(*,*)evv1
c***        write(*,*)f(i)
c***        write(*,*)df(i)
c***
      i = i + 1
      v = v + dv
      v1 = v - xi3
      evv(i) = exp(-v*v)
      evv1 = exp(-xi*xi*v1*v1)
      h0(i) = 1./piq
      h1(i) = sq2*v*h0(i)
      f(i) = p1*(evv(i) + p2*evv1)
      df(i) = p3*(v*evv(i) + p4*v1*evv1)
      b(1) = p5*(b(1) + f(i)*h1(i))
c***
c***        write(*,*)i
c***        write(*,*)evv(i)
c***        write(*,*)evv1
c***        write(*,*)f(i)
c***        write(*,*)df(i)
c***
      do i = -nv,nv
        f1(i) = (h0(i) + b(1)*h1(i))/piq*evv(i)
        df1(i) = -aos*h1(i)*evv(i)
      enddo
c
      do j = 2,nmax,2
        jp1 = j + 1
        jm1 = j - 1
        sj = sqrt(1.*j)
        sjp1 = sqrt(1.*(jp1))
        sjm1 = sqrt(1.*(jm1))
        i = -nv
        v = v0 + i*dv
        h0(i) = (sq2*v*h1(i) - sjm1*h0(i))/sj
        h1(i) = (sq2*v*h0(i) - sj*h1(i))/sjp1
        b(j) = f(i)*h0(i)
        b(jp1) = f(i)*h1(i)
        do i = -nv+1,nv-3,2
          v = v0 + i*dv
          h0(i) = (sq2*v*h1(i) - sjm1*h0(i))/sj
          h1(i) = (sq2*v*h0(i) - sj*h1(i))/sjp1
          b(j) = b(j) + 4.*f(i)*h0(i) 
          b(jp1) = b(jp1) + 4.*f(i)*h1(i) 
          i1 = i + 1
          v = v + dv
          h0(i1) = (sq2*v*h1(i1) - sjm1*h0(i1))/sj
          h1(i1) = (sq2*v*h0(i1) - sj*h1(i1))/sjp1
          b(j) = b(j) + 2.*f(i1)*h0(i1) 
          b(jp1) = b(jp1) + 2.*f(i1)*h1(i1) 
        enddo
        i = nv-1
        v = v0 + i*dv
        h0(i) = (sq2*v*h1(i) - sjm1*h0(i))/sj
        h1(i) = (sq2*v*h0(i) - sj*h1(i))/sjp1
        b(j) = b(j) + 4.*f(i)*h0(i) 
        b(jp1) = b(jp1) + 4.*f(i)*h1(i) 
        i = i + 1
        v = v + dv
        h0(i) = (sq2*v*h1(i) - sjm1*h0(i))/sj
        h1(i) = (sq2*v*h0(i) - sj*h1(i))/sjp1
        b(j) = p5*(b(j) + f(i)*h0(i)) 
        b(jp1) = p5*(b(jp1) + f(i)*h1(i)) 
        do i = -nv,nv
          v = v0 + i*dv
          f1(i) = f1(i) + (b(j)*h0(i) + b(jp1)*h1(i))/piq*evv(i)
          df1(i)=df1(i)-aos*(sj*b(jm1)*h0(i)+sjp1*b(j)*h1(i))*evv(i)
        enddo
      enddo
      do i = 1,nmax
        write(*,1) i,b(i)
      enddo
c
      write(*,5)" v                 ",
     .          " f0                ",
     .          " df0               ",
     .          " f                 ",
     .          " df                ",
     .          " f1                ",
     .          " df1               "
      do i = -nv,nv
        v = v0 + i*dv
        f0 = evv(i)/pis
        df0 = -alpha*v*f0
        write(*,10)v,f0,df0,f(i),df(i),f1(i),df1(i)
      enddo
1     format(i5,x,e24.16)     
5     format(7a17)
10    format(7(e16.8,x))
      stop
      end
c----------------------------------------------------------------------2
