      program main
      implicit double precision (a-h,o-z)
C
C String Method for Mueller potential
C by W. Ren
C
C 2d problem; 21 images along the string
      parameter (md=2,n=20)
      double precision x(md,0:n),x0(md,0:n)
      double precision s
      integer nx,ny,i
      character*20,fname0,fname1,PESdata
C
C time step
      dt=0.04d0
C
      fname0='data0'
      fname1='data1'
      PESdata='PES2.data'
C
C input initial path (From output of init.f)
      call input(fname0,md,n,x)
      call readPES(PESdata)
C
      ic=0
 100  continue
      ic=ic+1
      x0=x
C
C Evolve the images by steepest descent
      call mvper(dt,x)
C
C Reparameterization of the string for every 10 steps
      call mvtan(x)
C Calculating the error
C
C   By  Chen Xin
CCCCCCCCCCCCCCCCCCCCCCCCC
      s=0
      do 1 i=0,n
         s=(x(1,i)-x0(1,i))**2+(x(2,i)-x0(2,i))**2
 1    continue
      error=dsqrt(s)/(n-1)
C      print *,error

C
      if(mod(ic,10).eq.0) then
        write(*,*)'ic=',ic, '  error =',error
      endif



C
C output the string for every 1000 steps
      if(mod(ic,1000).eq.0) then
        call output(fname1,md,n,x)
      endif
C
C check convergence
      if(error.lt.1.e-15) then
        call output(fname1,md,n,x)
        stop
      endif
C
      goto 100
C
      end 
