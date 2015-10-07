      subroutine mvper(dt,x)
      implicit double precision (a-h,o-z)
C 
C Evlove each image by projected potential force
C by W. Ren
C
C input dt: time step
C       x: coordinates of images along the string
C Output x: new position of the string after one time step 
C        error: the L_\inf norm of the projected force
C
      parameter (md=2,n=20)
      double precision dt,x(md,0:n),t(md,n-1),err,err1
      double precision r(md),r2(md)
C
C calculate the unit tangent vector along the string
      call tangent(x,t)
C
      do 1 k=0,n
C
C calculate the potential force of the k-th image
C        print *,k,x(1,k)
        call force(md,x(1,k),r)
        do 3 i=1,md
          x(i,k)=x(i,k)+dt*r(i)
  3     continue
  1     continue
C 
C project the force onto the plane normal to the string
C      do 5 k=1,n-1
C        s=0.d0
C        call force(md,x(1,k),r)
C        do 4 i=1,md
C          s=s+r(i)*t(i,k)
C  4     continue
C        do 2 i=1,md
C          r2(i)=r(i)-s*t(i,k)
C  2     continue
C        tmp=0.d0
C        do 6 i=1,md
C          tmp=tmp+r2(i)**2.d0
C  6     continue
C        tmp=dsqrt(tmp)
C        error=dmax1(error,tmp)
C  5     continue
C Update the images
C
C
      end
